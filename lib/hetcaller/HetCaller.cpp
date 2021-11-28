// ***************************************************************************
// HetCaller.cpp (c) 2021 zhenhua yu <qasim0208@163.com>
// Health Informatics Lab, Ningxia University
// All rights reserved.

#include <iostream>
#include <cstdio>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <algorithm>
#include <string>
#include <limits>
#include <array>
#include <tuple>

#include "HetCaller.h"
#include "split.h"
#include "MyDefine.h"

using namespace std;

HetCaller::HetCaller() {
	pthread_mutex_init(&pm, NULL);
}

void HetCaller::loadData() {
	loadMutationData();
	loadRealData();
}

void HetCaller::loadMutationData() {
	string inputFile = config.getStringPara("input");
	ifstream ifs;
	ifs.open(inputFile.c_str());
	if(!ifs.is_open()) {
		cerr << "Error: cannot open file " << inputFile << endl;
		exit(-1);
	}
	
	FILE *fp;
	char buf[500];
	string cmd = "cat "+inputFile+" | wc -l";
	fp = popen(cmd.c_str(), "r");
	if(!fp) {
		cerr << "Error: cannot open file " << inputFile << endl;
		exit(-1);
	}
	char *p = fgets(buf, 500, fp);
	fclose(fp);
	int num_cell = atoi(buf);
	
	int i, j, k, indx;
	long line_num = 0;
	string line;
	int num_muta = 0;
	while(getline(ifs, line)) {
		line_num++;
		vector<string> fields = split(line, '\t');
		if(num_muta == 0) {
			num_muta = fields.size();
			obsData.resize(num_cell, num_muta, false);
		}
		if(num_muta != fields.size()) {
			cerr << "Error: malformed input file " << inputFile <<
                    ", inconsistent number of fields @line " << line_num << endl;
            cerr << buf << endl;
			exit(1);
		}
		for(i = 0; i < num_muta; i++) {
			indx = (line_num-1)*num_muta+i;
			k = atoi(fields[i].c_str());
			if(k != 0 && k != 1 && k != 3) {
				cerr << "unsupported genotype value " << k << " at line " << line_num << " in file " << inputFile << endl;
				exit(1);
			}
			obsData[indx] = k;
		}
	}
	ifs.close();
	assert(line_num == num_cell);
	
	cerr << "Total " << num_cell << " cells and " << num_muta << " mutations were loaded from file " << inputFile << endl;
}

void HetCaller::loadRealData() {
	string inputFile = config.getStringPara("rinput");
	if(inputFile.empty()) {
		return;
	}
	
	FILE *fp;
	char buf[500];
	string cmd = "cat "+inputFile+" | wc -l";
	fp = popen(cmd.c_str(), "r");
	if(!fp) {
		cerr << "Error: cannot open file " << inputFile << endl;
		exit(-1);
	}
	char *p = fgets(buf, 500, fp);
	fclose(fp);
	int num_cell = atoi(buf)-1;
	
	int i, j, k, indx;
	long line_num = 0;
	string line;
	int num_muta = 0;
	
	ifstream ifs;
	ifs.open(inputFile.c_str());
	if(!ifs.is_open()) {
		cerr << "Error: cannot open file " << inputFile << endl;
		exit(-1);
	}
	getline(ifs, line);
	vector<string> fields = split(line, ':');
	if(fields.size() > 1) {
		fields = split(fields[1], '\t');
		for(i = 0; i < fields.size(); i++) {
			doublet_indxs.push_back(atoi(fields[i].c_str())-1);
		}
	}
	
	while(getline(ifs, line)) {
		line_num++;
		vector<string> fields = split(line, '\t');
		if(num_muta == 0) {
			num_muta = fields.size();
			realData.resize(num_cell, num_muta, false);
		}
		if(num_muta != fields.size()) {
			cerr << "Error: malformed input file " << inputFile <<
                    ", inconsistent number of fields @line " << line_num << endl;
            cerr << buf << endl;
			exit(1);
		}
		for(i = 0; i < num_muta; i++) {
			indx = (line_num-1)*num_muta+i;
			k = atoi(fields[i].c_str());
			realData[indx] = k;
		}
	}
	ifs.close();
	assert(line_num == num_cell);
	
	cerr << "real mutation data were loaded from file " << inputFile << endl;
}

void HetCaller::preProcess() {
	int num_muta = obsData.getCOLS();
	int num_cell = obsData.getROWS();
	int i, j, k;
	
	Matrix<int> tmp = obsData;
	uniqueData = obsData;
	
	missing_rate = 0;
	for(i = 0; i < num_cell; i++) {
		for(j = 0; j < num_muta; j++) {
			k = obsData[i*num_muta+j];
			if(k == 3) {
				tmp[i*num_muta+j] = 0;
				missing_rate++;
			}
		}
	}
	missing_rate /= (num_cell*num_muta);
	
	int rindx = 0;
	for(i = 0; i < num_cell; i++) {
		int flag = 0;
		for(k = 0; k < rindx; k++) {
			for(j = 0; j < num_muta; j++) {
				if(uniqueData[k*num_muta+j] != tmp[i*num_muta+j]) {
					break;
				}
			}
			if(j == num_muta) {
				flag = 1;
				break;
			}
		}
		if(flag == 0) {
			for(j = 0; j < num_muta; j++) {
				uniqueData[rindx*num_muta+j] = tmp[i*num_muta+j];
			}
			rindx++;
		}
	}
	uniqueData.resize(rindx, num_muta, true);
	
}

void HetCaller::call() {
	int i, j, k, n;
	
	preProcess();
	
	int num_muta = obsData.getCOLS();
	int num_cell = obsData.getROWS();
	
	/*** find the optimal number of clusters based on BIC ***/
	int maxc = config.getIntPara("maxc");
	if(maxc <= 0) {
		maxc = num_cell/10;
	}
	else {
		maxc = min(maxc, num_cell);
	}
	
	int num_cluster = 1;
	int count = 0;
	int m_indx = -1;
	double pre_dist = numeric_limits<float>::max();
	while(num_cluster <= maxc) {
		/*** reasoning the mutation states of each cluster ***/
		num_cluster_c = num_cluster;
		predict(num_cluster);
		cerr << "--------------- screening report -----------------" << endl;
		printSolution(candi_s[num_cluster-1]);
		cerr << "--------------- screening report -----------------" << endl;
		
		if(candi_s[num_cluster-1].valid) {
			if(m_indx == -1) {
				m_indx = num_cluster-1;
			}
			else {
				
				if(candi_s[num_cluster-1].score-candi_s[m_indx].score >= 1e-5 && candi_s[num_cluster-1].ll > candi_s[m_indx].ll) {
					m_indx = num_cluster-1;
					count = 0;
				}
				else {
					count++;
					if(count >= 10) {
						break;
					}
				}
			}
		}
		num_cluster++;
	}
	
	for(i = candi_s.size()-1; i >= 1; i--) {
		if(!candi_s[i].valid)
			continue;
		double tmp = 0;
		for(j = i-1; j >= 0; j--) {
			tmp += 45;
			if((candi_s[i].ll-candi_s[j].ll) >= tmp) {
				candi_s[j].valid = 0;
			}
		}
	}
	
	
	best_s_indx = -1;
	for(i = 0; i < candi_s.size(); i++) {
		if(candi_s[i].valid) {
			if(best_s_indx == -1)	best_s_indx = i;
			
			else if(candi_s[i].score-candi_s[best_s_indx].score >= 1e-5) {
				best_s_indx = i;
			}
		}
	}
	if(best_s_indx == -1) {
		best_s_indx = 0;
		for(i = 1; i < candi_s.size(); i++) {
			if(candi_s[i].score-candi_s[best_s_indx].score >= 1e-5) {
				best_s_indx = i;
			}
		}
	}
	
	/*** best solution ***/
	cerr << "--------------- best solution -----------------" << endl;
	printSolution(candi_s[best_s_indx]);
	cerr << "-----------------------------------------------" << endl;
	
	/*** save results ***/
	saveResults();
}

void HetCaller::printSolution(Solution& s) {
	cerr << "#clusters: " << s.num_cluster << endl;
	//cerr << "Valid: " << s.valid << endl;
	cerr << "LL: " << s.ll << endl;
	cerr << "Score: " << s.score << endl;
	
	if(s.acc >= 0) {
		cerr << "Acc: " << s.acc << endl;
	}
	cerr << "alpha: " << s.alpha << endl;
	cerr << "beta: " << s.beta << endl;
}

void HetCaller::predict(int num_cluster) {
	int num_muta = obsData.getCOLS();
	int num_cell = obsData.getROWS();
	int num_cell_u = uniqueData.getROWS();
	int i, j, k, n;
	
	double alpha = config.getRealPara("alpha");
	double beta = config.getRealPara("beta");
	
	vector<double> betas;
	if(beta < 0) {
		beta = 0.01;
		while(beta <= 0.5) {
			betas.push_back(beta);
			beta += 0.05;
		}	
	}
	else {
		betas.push_back(beta);
	}
	
	int epoches = 50;
	epoches = min(epoches, num_cell_u);
	
	vector<int> indxs;
	for(i = 0; i < epoches; i++) {
		k = threadpool->randomInteger(0, num_cell_u);
		if(indxs.empty()) {
			indxs.push_back(k);
		}
		else {
			vector<int>::iterator it = find(indxs.begin(), indxs.end(), k);
			while(it != indxs.end()) {
				k = threadpool->randomInteger(0, num_cell_u);
				it = find(indxs.begin(), indxs.end(), k);
			}
			indxs.push_back(k);
		}
	}
	
	candi_s1.clear();
	vector<double*> paras_t;
	for(i = 0; i < epoches; i++) {
		k = indxs[i];
		for(j = 0; j < betas.size(); j++) {
			double* paras = new double[2];
			paras[0] = k; paras[1] = betas[j];
			threadpool->pool_add_work(&HetCaller::inferSolution, paras, i);
			paras_t.push_back(paras);
		}
	}
	threadpool->wait();
	for(i = 0; i < paras_t.size(); i++) {
		delete[] paras_t[i];
	}
	
	int best_indx = -1;
	for(i = 0; i < candi_s1.size(); i++) {
		Matrix<double> pie = candi_s1[i].pie;
		for(k = 0; k < num_cluster; k++) {
			if(pie[k]*num_cell < 5) {
				break;
			}
		}
		if(k < num_cluster)	continue;
		if(best_indx == -1)	best_indx = i;
		else if(candi_s1[i].ll > candi_s1[best_indx].ll) {
			best_indx = i;
		}
	}
	
	if(best_indx == -1) {
		best_indx = 0;
		for(i = 1; i < candi_s1.size(); i++) {
			if(candi_s1[i].ll > candi_s1[best_indx].ll) {
				best_indx = i;
			}
		}
	}
	
	evalAccuracy(candi_s1[best_indx]);
	evaluateSolution(candi_s1[best_indx]);
	candi_s1[best_indx].num_cluster = num_cluster;
	candi_s.push_back(candi_s1[best_indx]);
}

void* HetCaller::inferSolution(const void *arg) {
	double* tmp_p = (double*) arg;
	int c_i = tmp_p[0];
	double beta = tmp_p[1];
	int i, j, k, n;
	
	Matrix<int>& uniqueData = hetcaller.getUniqueData();
	int num_cluster = hetcaller.getNumOfClusters();
	int num_muta = uniqueData.getCOLS();
	int num_cell_u = uniqueData.getROWS();
	
	double alpha = config.getRealPara("alpha");	
	bool alpha_fixed, beta_fixed;
	if(alpha >= 0) {
		alpha_fixed = true;
	}
	else {
		alpha = 0.01;
		alpha_fixed = false;
	}
	if(config.getRealPara("beta") >= 0) {
		beta_fixed = true;
	}
	else {
		beta_fixed = false;
	}
	
	Matrix<int> states(num_cluster, num_muta);
	Matrix<double> pie(1, num_cluster);
	
	if(num_cluster == 1) {
		for(j = 0; j < num_muta; j++) {
			states[j] = uniqueData[c_i*num_muta+j];
		}
		pie.set(1.0);
	}
	else {
		Solution& s_pre = hetcaller.getSolution(num_cluster-2);
		double p = threadpool->randomDouble(0, 1);
		if(s_pre.valid || p > 0.3) {
			states.setRows(0, num_cluster-2, s_pre.states);
			for(j = 0; j < num_muta; j++) {
				states[(num_cluster-1)*num_muta+j] = uniqueData[c_i*num_muta+j];
			}
			pie.set(1.0/num_cluster);
		}
		else {
			vector<int> indxs;
			for(i = 0; i < num_cluster; i++) {
				k = threadpool->randomInteger(0, num_cell_u);
				if(indxs.empty()) {
					indxs.push_back(k);
				}
				else {
					vector<int>::iterator it = find(indxs.begin(), indxs.end(), k);
					while(it != indxs.end()) {
						k = threadpool->randomInteger(0, num_cell_u);
						it = find(indxs.begin(), indxs.end(), k);
					}
					indxs.push_back(k);
				}
			}
			
			uniqueData.Rows(indxs, states);
			pie.set(1.0/num_cluster);
		}
	}
	
	
	vector<int> para_updates(4, 1);
	para_updates[0] = !alpha_fixed;
	para_updates[1] = !beta_fixed;
	
	Solution s(num_cluster, alpha, beta, pie, states);
	hetcaller.inferParas(s, para_updates);
	//hetcaller.evaluateSolution(s);
	hetcaller.saveSolution(s);
	
	return NULL;
}

void HetCaller::inferParas(Solution& s, const vector<int>& para_updates) {
	int i, j, k, n;
	int num_cell = obsData.getROWS();
	int num_muta = obsData.getCOLS();
	
	int num_cluster = s.num_cluster;
	double& alpha = s.alpha;
	double& beta = s.beta;
	
	Matrix<int> states_b = s.states;
	Matrix<int>& states = s.states;
	Matrix<int> states_n(num_cluster, num_muta);
	Matrix<double> pie_b = s.pie;
	Matrix<double>& pie = s.pie;
	Matrix<double> pie_n(1, num_cluster);
	s.post_probs.resize(num_cell, num_cluster, false);
	Matrix<double>& post_probs = s.post_probs;
	Matrix<long double> probs(num_cell, num_cluster);
	
	double max_alpha = 0.1, max_beta = 0.99;
	double eps = numeric_limits<long double>::epsilon();
	
	double ll, pre_ll = numeric_limits<long>::min();
	int iter = 0, max_iter = config.getIntPara("max_iter");
	while(iter < max_iter) {
		
		obsProbs(states, probs, alpha, beta);
		// E-step: posterior probability calculation 
		ll = 0;
		for(i = 0; i < num_cell; i++) {
			long double sum = 0;
			for(k = 0; k < num_cluster; k++) {
				sum += probs[i*num_cluster+k]*pie[k];
			}
			for(k = 0; k < num_cluster; k++) {
				post_probs[i*num_cluster+k] = probs[i*num_cluster+k]*pie[k]/sum;
			}
			ll += logl(sum+eps);
		}
		if(isnan(ll)) {
			alpha = 0.01;
			beta = 0.01;
			ll = numeric_limits<int>::min();
			post_probs.set(1.0/num_cluster);
			states = states_b;
			pie = pie_b;
			break;
		}
		
		// M-step: parameter update
		
		// update alpha and beta
		
		double tmp1 = 0, tmp2 = 0, tmp3 = 0, tmp4 = 0, tmp5 = 0;
		for(i = 0; i < num_cell; i++) {
			for(k = 0; k < num_cluster; k++) {
				double gamma = post_probs[i*num_cluster+k];
				for(j = 0; j < num_muta; j++) {
					if(obsData[i*num_muta+j] == 3) {
						continue;
					}
					int s = obsData[i*num_muta+j];
					tmp1 += gamma*(1-states[k*num_muta+j])*s;
					tmp2 += gamma*(1-states[k*num_muta+j]);
					tmp3 += gamma*states[k*num_muta+j]*(1-s);
					tmp4 += gamma*states[k*num_muta+j];
				}
			}
		}
		double alpha_n = tmp1/(tmp2+eps);
		double beta_n = tmp3/(tmp4+eps);
		
		if(para_updates[0]) {
			alpha = isnan(alpha_n)? alpha:alpha_n;
			alpha = (alpha > max_alpha)? max_alpha:alpha;
		}
		if(para_updates[1]) {
			beta = isnan(beta_n)? beta:beta_n;
			beta = (beta > max_beta)? max_beta:beta;
		}
		
		// update states
		if(para_updates[2] == 1) {
			for(k = 0; k < num_cluster; k++) {
				for(j = 0; j < num_muta; j++) {
					double f1 = 0, f2 = 0;
					for(i = 0; i < num_cell; i++) {		
						if(obsData[i*num_muta+j] == 3) {
							continue;
						}
						int s = obsData[i*num_muta+j];
						f1 += post_probs[i*num_cluster+k]*(s*log(1-beta+eps)+(1-s)*log(beta+eps));
						f2 += post_probs[i*num_cluster+k]*(s*log(alpha+eps)+(1-s)*log(1-alpha+eps));			
					}
					states_n[k*num_muta+j] = (f1 > f2)? 1:0;
				}
			}
			states = states_n;
		}
		
		// update pie
		for(k = 0; k < num_cluster; k++) {
			double sum = 0;
			for(i = 0; i < num_cell; i++) {
				sum += post_probs[i*num_cluster+k];
			}
			pie_n[k] = sum/num_cell;
		}
		if(para_updates[3]) {
			pie = pie_n;
		}
		
		if(fabs(pre_ll-ll) < 1e-3) {
			break;
		}
		else {
			pre_ll = ll;
		}
		iter++;
	}
	s.ll = ll;
}

void HetCaller::evalAccuracy(Solution& s) {
	if(realData.getROWS() == 0) {
		s.acc = -1;
		return;
	}
	int i, j, k;
	int num_cell = obsData.getROWS();
	int num_muta = obsData.getCOLS();
	
	Matrix<int>& states = s.states;
	Matrix<double>& post_probs = s.post_probs;
	Matrix<int> s_indxs;
	post_probs.max(2, s_indxs);
	Matrix<double> results(num_cell, num_muta, 0.0);
	int num_cluster = s.num_cluster;
	for(i = 0; i < num_cell; i++) {
		k = s_indxs[i];
		for(j = 0; j < num_muta; j++) {
			results[i*num_muta+j] = states[k*num_muta+j];
		}
	}
	
	double acc = 0;
	int n = 0;
	for(i = 0; i < num_cell; i++) {
		for(j = 0; j < doublet_indxs.size(); j++) {
			if(doublet_indxs[j] == i) {
				break;
			}
		}
		if(j < doublet_indxs.size()) {
			continue;
		}
		n++;
		for(j = 0; j < num_muta; j++) {
			k = (results[i*num_muta+j] > 0.5)? 1:0;
			if(k == realData[i*num_muta+j])	acc++;
		}
	}
	s.acc = acc/(n*num_muta);
}

void HetCaller::evaluateSolution(Solution& s) {
	int i, j, k, n;
	int num_cluster = s.num_cluster;
	int num_cell = obsData.getROWS();
	int num_muta = obsData.getCOLS();
	
	
	Matrix<int>& states = s.states;
	s.valid = true;
	
	Matrix<int> tmp = states.sumRows();
	int unmutated_count = 0;
	for(i = 0; i < num_muta; i++) {
		if(tmp[i] == 0) {
			unmutated_count++;
		}
	}
	if(1.0*unmutated_count/num_muta >= 0.1) {
		s.valid = false;
	}
	
	calculateScore(s);
	
	double mean_score = 0, min_score = 1;
	if(num_cluster > 1) {
		for(i = 0, k = 0; i < num_cluster-1; i++) {
			for(j = i+1; j < num_cluster; j++) {
				if(s.scores[i*num_cluster+j] < 0) {
					continue;
				}
				mean_score += s.scores[i*num_cluster+j];
				k++;
				if(min_score > s.scores[i*num_cluster+j]) {
					min_score = s.scores[i*num_cluster+j];
				}
			}
		}
		mean_score /= k;
	}
	else {
		mean_score = 0;
	}
	
	s.score = mean_score;
}

void HetCaller::calculateScore(Solution& s) {
	int c1, c2;
	
	int i, j, k, n;
	
	int num_cluster = s.num_cluster;
	double beta = s.beta;
	double alpha = s.alpha;
	Matrix<int>& states = s.states;
	Matrix<double>& post_probs = s.post_probs;
	int num_cell = post_probs.getROWS();
	int num_muta = states.getCOLS();
	
	s.scores.resize(num_cluster, num_cluster, false);
	
	Matrix<int> indxs;
	post_probs.max(2, indxs);
	map<int, vector<int>> c_indxs_a;
	for(j = 0; j < num_cell; j++) {
		c_indxs_a[indxs[j]].push_back(j);
	}
	
	double p_00 = 2*alpha*(1-alpha)*pow(1-missing_rate,2)+2*missing_rate*(1-missing_rate);
	double p_01 = pow(1-missing_rate,2)*((1-alpha)*(1-beta)+alpha*beta)+2*missing_rate*(1-missing_rate);
	double p_10 = p_01;
	double p_11 = 2*beta*(1-beta)*pow(1-missing_rate,2)+2*missing_rate*(1-missing_rate);
	
	for(c1 = 0; c1 < num_cluster-1; c1++) {
		for(c2 = c1+1; c2 < num_cluster; c2++) {
			int counts[4] = {0};
			for(k = 0; k < num_muta; k++) {
				i = states[c1*num_muta+k]*2+states[c2*num_muta+k];
				counts[i]++;
			}
			double expected_dist = counts[0]*p_00+(counts[1]+counts[2])*p_01+counts[3]*p_11;
			double expected_dist0 = counts[0]*p_00;
			double expected_dist1 = (counts[1]+counts[2])*p_01;
			double expected_dist2 = counts[3]*p_11;
			
			vector<int>& c1_indxs = c_indxs_a[c1];
			vector<int>& c2_indxs = c_indxs_a[c2];
			
			if(c1_indxs.size()*c2_indxs.size() < 10) {
				s.scores[c1*num_cluster+c2] = -1;
				s.scores[c2*num_cluster+c1] = -1;
			}
			else {
				int num_c1 = c1_indxs.size();
				int num_c2 = c2_indxs.size();
				double ave_dist = 0;
				for(i = 0; i < num_c1; i++) {
					for(j = 0; j < num_c2; j++) {
						double dist = 0;
						for(k = 0; k < num_muta; k++) {
							if(obsData[c1_indxs[i]*num_muta+k] != obsData[c2_indxs[j]*num_muta+k]) {
								dist++;
							}
						}
						ave_dist += dist;
					}
				}				
				ave_dist /= (num_c1*num_c2);
				s.scores[c1*num_cluster+c2] = exp(-pow(expected_dist-ave_dist, 2));
				s.scores[c2*num_cluster+c1] = s.scores[c1*num_cluster+c2];
			}
		}
	}
}

void HetCaller::obsProbs(Matrix<int>& states, Matrix<long double>& probs, double alpha, double beta) {
	int i, j, k;
	int num_muta = obsData.getCOLS();
	int num_cell = obsData.getROWS();
	int num_cluster = states.getROWS();
	long double eps = numeric_limits<long double>::epsilon();
	
	probs.resize(num_cell, num_cluster, false);
	
	for(i = 0; i < num_cell; i++) {
		for(k = 0; k < num_cluster; k++) {
			long double prob = 1;
			for(j = 0; j < num_muta; j++) {
				if(states[k*num_muta+j] == 1) {
					if(obsData[i*num_muta+j] == 0) {
						prob *= beta;
					}
					else if(obsData[i*num_muta+j] == 1) {
						prob *= 1-beta;
					}
				}
				else {
					if(obsData[i*num_muta+j] == 0) {
						prob *= 1-alpha;
					}
					else if(obsData[i*num_muta+j] == 1) {
						prob *= alpha;
					}
				}
			}
			probs[i*num_cluster+k] = prob;
		}
	}
}

void HetCaller::saveSolution(Solution& s) {
	pthread_mutex_lock(&pm);
	candi_s1.push_back(s);
	pthread_mutex_unlock(&pm);
}

void HetCaller::saveResults() {
	int i, j, k;
	ofstream ofs;
	string fn, outputPrefix = config.getStringPara("output");
	
	int num_cell = obsData.getROWS();
	int num_muta = obsData.getCOLS();
	/*
	fn = outputPrefix+".solution.txt";
	ofs.open(fn.c_str());
	if(!ofs.is_open()) {
		cerr << "Error: cannot open file " << fn << endl;
		exit(-1);
	}
	ofs << "alpha = " << candi_s[best_s_indx].alpha << endl;
	ofs << "beta = " << candi_s[best_s_indx].beta << endl;
	//ofs << "K = " << candi_s[best_s_indx].num_cluster << endl;
	ofs.close();
	*/
	
	int num_cluster = candi_s[best_s_indx].num_cluster;
	Matrix<int>& states = candi_s[best_s_indx].states;
	fn = outputPrefix+".clone_genotypes";
	ofs.open(fn.c_str());
	if(!ofs.is_open()) {
		cerr << "Error: cannot open file " << fn << endl;
		exit(-1);
	}
	for(i = 0; i < num_cluster; i++) {
		for(j = 0; j < num_muta; j++) {
			if(j < num_muta-1) {
				ofs << states[i*num_muta+j] << '\t';
			}
			else {
				ofs << states[i*num_muta+j] << endl;
			}
		}
	}
	ofs.close();
	
	Matrix<double>& post_probs = candi_s[best_s_indx].post_probs;
	Matrix<int> s_indxs;
	post_probs.max(2, s_indxs);
	fn = outputPrefix+".cell_assignment";
	ofs.open(fn.c_str());
	if(!ofs.is_open()) {
		cerr << "Error: cannot open file " << fn << endl;
		exit(-1);
	}
	for(i = 0; i < num_cell-1; i++) {
		ofs << s_indxs[i] << " ";
	}
	ofs << s_indxs[num_cell-1] << endl;
	/*
	for(k = 0; k < num_cluster; k++) {
		for(j = 0; j < num_cell; j++) {
			if(s_indxs[j] != k) {
				continue;
			}
			ofs << j << " ";
		}
		ofs << endl;
	}
	*/
	ofs.close();
	
}


