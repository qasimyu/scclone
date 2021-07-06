// ***************************************************************************
// HetCaller.cpp (c) 2020 zhenhua yu <qasim0208@163.com>
// Health Informatics Lab, Ningxia University
// All rights reserved.

#include <iostream>
#include <cstdio>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <string>
#include <limits>
#include <array>
#include <tuple>

#include "HetCaller.h"
#include "split.h"
#include "MyDefine.h"
#include "PCA.h"
#include "dkm.hpp"

using namespace std;

HetCaller::HetCaller() {
	pthread_mutex_init(&pm, NULL);
}

void HetCaller::loadData() {
	loadMutationData();
	loadRealData();
	loadCellLabels();
	loadMutaLabels();
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
	homo_muta = false;
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
			if(k == 2) {
				homo_muta = true;
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

void HetCaller::loadCellLabels() {
	string inputFile = config.getStringPara("clabel");
	if(inputFile.empty()) {
		char buf[100];
		for(int i = 0; i < obsData.getROWS(); i++) {
			sprintf(buf, "%d", i+1);
			cLabels.push_back(buf);
		}
		return;
	}
	ifstream ifs;
	ifs.open(inputFile.c_str());
	if(!ifs.is_open()) {
		cerr << "Error: cannot open file " << inputFile << endl;
		exit(-1);
	}
	
	long line_num = 0;
	string line;
	while(getline(ifs, line)) {
		line_num++;
		if(line.empty()) {
			cerr << "Error: malformed label file " << inputFile <<
                    ", empty item is not allowed @line " << line_num << endl;
			exit(1);
		}
		cLabels.push_back(line);
	}
	ifs.close();
	
	if(cLabels.size() != obsData.getROWS()) {
		cerr << "Error: the number of cell labels is not consistent with the number of cells." << endl;
		exit(1);
	}
	
	cerr << "the labels of cells were loaded from file " << inputFile << endl;
}

void HetCaller::loadMutaLabels() {
	string inputFile = config.getStringPara("mlabel");
	if(inputFile.empty()) {
		char buf[100];
		for(int i = 1; i <= obsData.getCOLS(); i++) {
			sprintf(buf, "%d", i);
			mLabels.push_back(buf);
		}
		return;
	}
	ifstream ifs;
	ifs.open(inputFile.c_str());
	if(!ifs.is_open()) {
		cerr << "Error: cannot open file " << inputFile << endl;
		exit(-1);
	}
	
	long line_num = 0;
	string line;
	while(getline(ifs, line)) {
		line_num++;
		if(line.empty()) {
			cerr << "Error: malformed label file " << inputFile <<
                    ", empty item is not allowed @line " << line_num << endl;
			exit(1);
		}
		mLabels.push_back(line);
	}
	ifs.close();
	
	if(mLabels.size() != obsData.getCOLS()) {
		cerr << "Error: the number of mutation labels is not consistent with the number of mutations." << endl;
		exit(1);
	}
	
	cerr << "the labels of mutations were loaded from file " << inputFile << endl;
}

void HetCaller::fetchUniqueData() {
	int num_muta = obsData.getCOLS();
	int num_cell = obsData.getROWS();
	int i, j, k, rindx = 0;
	
	Matrix<int> tmp = obsData;
	uniqueData = obsData;
	
	missing_rate = 0;
	for(i = 0; i < num_cell; i++) {
		for(j = 0; j < num_muta; j++) {
			if(tmp[i*num_muta+j] == 3) {
				missing_rate++;
			}
			if(tmp[i*num_muta+j] == 1){
				positive_rate++;
			}
		}
	}
	positive_rate /= (num_cell*num_muta-missing_rate);
	missing_rate /= (num_cell*num_muta);
	
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
	replicate.resize(rindx, num_muta, false);
	for(i = 0; i < rindx; i++) {
		for(j = 0; j < num_muta; j++) {
			k = uniqueData[i*num_muta+j];
			if(homo_muta) {
				if(k == 3) {
					replicate[i*num_muta+j] = 0.5;
					//replicate[i*num_muta+j] = 0.01;
					uniqueData[i*num_muta+j] = 0;
				}
				else if(k == 2) {
					replicate[i*num_muta+j] = 0.99;
					uniqueData[i*num_muta+j] = 1;
				}
				else {
					replicate[i*num_muta+j] = k;
				}
			}
			else {
				if(k == 3) {
					replicate[i*num_muta+j] = 0.5;
					//replicate[i*num_muta+j] = 0.01;
					uniqueData[i*num_muta+j] = 0;
				}
				else {
					replicate[i*num_muta+j] = k;
				}
				
			}
		}
	}
	
}

void HetCaller::call() {
	int i, j, k, n;
	
	//derive unique data
	fetchUniqueData();
	
	Matrix<float> reduced_data;
	//if(uniqueData.getROWS() >= uniqueData.getCOLS()) {
		//perform PCA, keeping 90% energy
		PCA pca(replicate);
		pca.important_ONB(0.9);
		reduced_data = pca.getProj();
		//cerr << "Dimension reduces to " << reduced_data.getROWS() << "X" << reduced_data.getCOLS() << endl;
	/*}
	else {
		//reduced_data = replicate;
	}*/
	
	for(i = 0; i < reduced_data.getROWS(); i++) {
		vector<float> tmp;
		tmp.reserve(reduced_data.getCOLS());
		for(j = 0; j < reduced_data.getCOLS(); j++) {
			tmp.push_back(reduced_data[i*reduced_data.getCOLS()+j]);
		}
		data_for_cluster.push_back(tmp);
	}
	
	int num_cluster = obsData.getROWS()/10;
	int km_reps = config.getIntPara("km_reps");	
		
	min_dist = numeric_limits<float>::max();
	vector<int*> t_paras;
	for(i = 0; i < km_reps; i++) {
		int* paras = new int[2];
		paras[0] = num_cluster; paras[1] = i;
		t_paras.push_back(paras);
		threadpool->pool_add_work(&HetCaller::clusterCells, paras, i);
	}
	threadpool->wait();
	for(i = 0; i < km_reps; i++) {
		delete[] t_paras[i];
	}
	
	/*** estiamte cluster centers and parameters ***/
	int num_clone_max = config.getIntPara("maxc");
	if(num_clone_max == -1) {
		num_clone_max = min(50, obsData.getROWS()/20+2);
	}
	else {
		num_clone_max = min(num_clone_max, obsData.getROWS()/20+2);
	}
	for(i = 1; i <= num_clone_max; i++) {
		predict(i);
		cerr << "iteration " << i << " done." << endl;
	}
	
	for(i = 0; i < num_clone_max; i++) {
		if(candi_s[i].valid) {
			break;
		}
	}
	if(i == num_clone_max) {
		for(i = 0; i < num_clone_max; i++) {
			candi_s[i].valid = 1;
		}
	}
	
	for(i = num_clone_max-1; i >= 1; i--) {
		double sum = 0;
		for(j = i-1; j >= 0; j--) {
			sum += 0.03*fabs(candi_s[j].ll); // 0.005
			if((candi_s[i].ll-candi_s[j].ll) >= sum) {
				candi_s[j].valid = 0;
			}
			else if(fabs(candi_s[i].ll-candi_s[j].ll) < 1 && candi_s[j].beta-candi_s[i].beta >= 0.02) {
				candi_s[j].valid = 0;
			}
		}
	}
	/*
	for(i = 1; i < num_clone_max; i++) {
		if(candi_s[i-1].beta-candi_s[i].beta >= 0.02) {
			candi_s[i-1].valid = 0;
		}
	}
	*/
	
	//find solution with maximum score
	double min_ll = candi_s[0].ll, max_ll = candi_s[0].ll;
	for(i = 1; i < num_clone_max; i++) {
		if(candi_s[i].ll > max_ll) {
			max_ll = candi_s[i].ll;
		}
		if(candi_s[i].ll < min_ll) {
			min_ll = candi_s[i].ll;
		}
	}
	for(i = 0; i < num_clone_max; i++) {
		if(fabs(max_ll-min_ll) < 1.0e-3) {
			candi_s[i].score1 = 1;
		}
		else {
			candi_s[i].score1 = (candi_s[i].ll-min_ll)/(max_ll-min_ll);
		}
		candi_s[i].score = candi_s[i].score1+candi_s[i].score2;
		//cerr << i+1 << ": " << candi_s[i].score << ", " << candi_s[i].score1 << endl;
	}
	
	/*
	for(i = 1; i <= num_clone_max; i++) {
		cerr << "--------------- screening report -----------------" << endl;
		printSolution(candi_s[i-1]);
		cerr << "--------------- screening report -----------------" << endl;
		
	}
	*/
	
	int num_cell = obsData.getROWS();
	best_s_indx = -1;
	for(i = 0; i < num_clone_max; i++) {
		if(candi_s[i].valid) {
			if(best_s_indx == -1) {
				best_s_indx = i;
			}
			else if(candi_s[i].score > candi_s[best_s_indx].score) {
				best_s_indx = i;
			}
		}
	}
	if(best_s_indx == -1) {
		best_s_indx = 0;
		for(i = 1; i < num_clone_max; i++) {
			if(candi_s[i].score > candi_s[best_s_indx].score) {
				best_s_indx = i;
			}
		}
	}
	/*** best solution ***/
	cerr << "--------------- best solution -----------------" << endl;
	printSolution(candi_s[best_s_indx]);
	cerr << "-----------------------------------------------" << endl;
	
	int rindx = 0;
	int num_clone = candi_s[best_s_indx].num_clone;
	int num_muta = obsData.getCOLS();
	Matrix<int> uniq_clusters(num_clone, num_muta);
	Matrix<int>& clusters = candi_s[best_s_indx].clusters;
	for(i = 0; i < num_clone; i++) {
		int flag = 0;
		for(k = 0; k < rindx; k++) {
			for(j = 0; j < num_muta; j++) {
				if(uniq_clusters[k*num_muta+j] != clusters[i*num_muta+j]) {
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
				uniq_clusters[rindx*num_muta+j] = clusters[i*num_muta+j];
			}
			rindx++;
		}
	}
	uniq_clusters.resize(rindx, num_muta, true);
	
	if(rindx != num_clone) {
		candi_s[best_s_indx].num_clone = rindx;
		candi_s[best_s_indx].pie.resize(1, rindx, false);
		candi_s[best_s_indx].pie.set(1.0/rindx);
		candi_s[best_s_indx].clusters = uniq_clusters;
		vector<int> para_updates(4, 0);
		para_updates[3] = 1;
		
		inferParas(candi_s[best_s_indx], para_updates);
		evaluateSolution(candi_s[best_s_indx]);
		candi_s[best_s_indx].score = (candi_s[best_s_indx].ll-min_ll)/(max_ll-min_ll)+candi_s[best_s_indx].score2;
		
		/*** refined best solution ***/
		cerr << "--------------- refined best solution -----------------" << endl;
		printSolution(candi_s[best_s_indx]);
		cerr << "-----------------------------------------------------------------------------" << endl;
	}
	
	clonal_tree = buildTree(candi_s[best_s_indx]);
	
	saveResults();
}

void HetCaller::printSolution(Solution& s) {
	cerr << "clones: " << s.num_clone << endl;
	cerr << "LL: " << s.ll << endl;
	//cerr << "BIC: " << s.bic << endl;
	
	cerr << "Score: " << s.score << ", " << s.score1 << ", " << s.score2 << endl;
	//cerr << "Valid: " << s.valid << endl;
	if(s.acc >= 0) {
		cerr << "Acc: " << s.acc << endl;
	}
	cerr << "alpha: " << s.alpha << endl;
	cerr << "beta: " << s.beta << endl;
	cerr << "pie: ";
	s.pie.Print();
}

void HetCaller::predict(int num_clone) {
	int num_cell = obsData.getROWS();
	int num_muta = obsData.getCOLS();
	int num_cell_u = uniqueData.getROWS();
	int i, j, epoches;
	
	int num_cluster = num_cell/10;
	if(num_clone == num_cluster) {
		epoches = 50;
	}
	else {
		int N = max(1,num_cluster/num_clone);
		double tmp = 1.0;
		for(i = num_clone+1, j = 1; i <= num_cluster; i++, j++) {
			tmp *= 1.0*i/(i-num_clone);
			if(j <= num_clone) {
				tmp /= N;
			}
		}
		while(j <= num_clone) {
			tmp /= N;
			j++;
		}
		epoches = round(tmp);
		if(epoches > 200) {
			epoches = 200;
		}
		if(epoches < 100) {
			epoches = 100;
		}
	}
	
	candi_s1.clear();
	for(i = 0; i < epoches; i++) {
		threadpool->pool_add_work(&HetCaller::inferSolution, &num_clone, i);
	}
	threadpool->wait();
	
	int num_s = candi_s1.size();
	int num_thread = threadpool->getThreadNumber();
	int N = max(1, num_s/num_thread);
	i = 0;
	vector<int*> tparas;
	while(i < num_s) {
		int *paras = new int[2];
		paras[0] = i;
		if(i+N-1 < num_s) {
			paras[1] = i+N-1;
		}
		else {
			paras[1] = num_s-1;
		}
		threadpool->pool_add_work(&HetCaller::evaluateSolutions, paras, i);
		tparas.push_back(paras);
		i += N;
	}
	threadpool->wait();
	for(i = 0; i < tparas.size(); i++) {
		delete[] tparas[i];
	}
	
	int best_indx = -1;
	for(i = 0; i < num_s; i++) {
		if(candi_s1[i].valid) {
			if(best_indx == -1) {
				best_indx = i;
			}
			else if(candi_s1[i].score2 > candi_s1[best_indx].score2) {
				best_indx = i;
			}
		}
	}
	if(best_indx == -1) {
		best_indx = 0;
		for(i = 1; i < num_s; i++) {
			if(candi_s1[i].score2 > candi_s1[best_indx].score2) {
				best_indx = i;
			}
		}
	}
	candi_s1[best_indx].num_clone = num_clone;
	candi_s1[best_indx].bic = -2*candi_s1[best_indx].ll+0.1*log(num_cell)*(num_clone*(1+num_muta)+1);
	candi_s.push_back(candi_s1[best_indx]);
}

void* HetCaller::clusterCells(const void *arg) {
	int* tmp = (int*) arg;
	int num_clone = tmp[0];
	int seed = tmp[1];
	
	vector<vector<float>>& data_for_cluster = hetcaller.getDataForCluster();
	
	float dist, score;
	dkm::clustering_parameters<float> paras(num_clone);
	//paras.set_random_seed(threadpool->randomInteger(0, 10000));
	paras.set_random_seed(seed);
	paras.set_max_iteration(100);
	auto cluster_data = dkm::kmeans_lloyd(data_for_cluster, paras, dist, score);
	
	hetcaller.saveClusterResults(dist, score, get<1>(cluster_data));
	
	return NULL;
}

void HetCaller::saveClusterResults(float dist, float score, vector<uint32_t>& indices) {
	pthread_mutex_lock(&pm);
	
	if(min_dist > dist) {
		min_dist = dist;
		best_indices = indices;
	}
	pthread_mutex_unlock(&pm);
}

void* HetCaller::inferSolution(const void *arg) {
	int num_clone = *((int*) arg);
	int i, j, k, n;
	
	Matrix<int>& uniqueData = hetcaller.getUniqueData();
	vector<uint32_t>& indices = hetcaller.getBestIndices();
	vector<Matrix<int>>& clusters_m = hetcaller.getClusters_m();
	
	int num_muta = uniqueData.getCOLS();
	int num_cell_u = uniqueData.getROWS();
	
	int num_cluster_m = clusters_m.size();
	Matrix<int> clusters_p(num_clone, num_muta);
	
	map<int, vector<int>> cell_assignments;
	for(i = 0; i < indices.size(); i++) {
		cell_assignments[indices[i]].push_back(i);
	}
	
	int num_cluster = cell_assignments.size();
	vector<int> indxs;
	int *tmp = new int[num_cluster];
	for(i = 0; i < num_cluster; i++) {
		tmp[i] = i;
	}
	
	j = num_cluster;
	//cerr << "selected unique clusters:" << endl;
	for(i = 0; i < num_clone; i++, j--) {
		k = threadpool->randomInteger(0, j);
		//cerr << tmp[k] << ' ';
		indxs.push_back(tmp[k]);
		for(; k < num_cluster-1; k++) {
			tmp[k] = tmp[k+1];
		}
	}
	delete[] tmp;
	for(i = 0; i < num_clone; i++) {
		vector<int>& cells = cell_assignments[indxs[i]];
		if(cells.empty()) {
			k = threadpool->randomInteger(0, num_cell_u);
			for(n = 0; n < num_muta; n++) {
				clusters_p[i*num_muta+n] = uniqueData[k*num_muta+n];
			}
		}
		else if(1.0*cells.size()/num_cell_u >= 0.2) {
			for(n = 0; n < num_muta; n++) {
				Matrix<int> counts(1, 2, 0);
				for(j = 0; j < cells.size(); j++) {
					int c = cells[j];
					int state = uniqueData[c*num_muta+n];
					counts[state]++;
				}
				counts.max(k);
				clusters_p[i*num_muta+n] = k;
			}
		}
		else {
			k = threadpool->randomInteger(0, cells.size());
			for(n = 0; n < num_muta; n++) {
				clusters_p[i*num_muta+n] = uniqueData[cells[k]*num_muta+n];
			}
		}
	}
	
	// search for optimal initial values of beta
	double max_alpha = config.getRealPara("max_alpha");
	double max_beta = config.getRealPara("max_beta");
	double alpha = config.getRealPara("alpha");
	double beta = config.getRealPara("beta");
	
	vector<double> alphas, betas;
	bool alpha_fixed, beta_fixed;
	if(alpha >= 0) {
		alphas.push_back(alpha);
		alpha_fixed = true;
	}
	else {
		/*
		if(num_clone == 1) {
			int iter1 = 5;
			for(n = 0; n <= iter1; n++) {
				alpha = n*max_alpha/iter1+5e-6;
				alphas.push_back(alpha);
			}
		}
		else {
			Solution& s = hetcaller.getSolution(num_clone-2);
			alphas.push_back(s.alpha);
		}
		*/
		alphas.push_back(0.01);
		alpha_fixed = false;
	}
	if(beta >= 0) {
		betas.push_back(beta);
		beta_fixed = true;
	}
	else {
		/*
		if(num_clone == 1) {
			int iter1 = 5;
			for(n = 0; n <= iter1; n++) {
				beta = n*max_beta/iter1+0.01;
				betas.push_back(beta);
			}
		}
		else {
			Solution& s = hetcaller.getSolution(num_clone-2);
			betas.push_back(s.beta);
		}
		*/
		betas.push_back(0.01);
		beta_fixed = false;
	}
	
	vector<int> para_updates(4, 1);
	para_updates[0] = !alpha_fixed;
	para_updates[1] = !beta_fixed;
	//para_updates[2] = 0;
	
	//record best solution
	Solution best_s;
	best_s.ll = numeric_limits<long>::min();
	
	Matrix<double> pie(1, num_clone);
	pie.set(1.0/num_clone);
	
	for(i = 0; i < alphas.size(); i++) {
		for(j = 0; j < betas.size(); j++) {
			Solution s(num_clone, alphas[i], betas[j], pie, clusters_p);
			hetcaller.inferParas(s, para_updates);
			if(s.ll > best_s.ll) {
				best_s = s;
			}
		}
	}
	hetcaller.saveSolution(best_s);
	
	return NULL;
}

void HetCaller::inferParas(Solution& s, const vector<int>& para_updates) {
	int i, j, k, n;
	int num_cell = obsData.getROWS();
	int num_muta = obsData.getCOLS();
	int num_cell_u = uniqueData.getROWS();
	
	int num_clone = s.num_clone;
	double& alpha = s.alpha;
	double& beta = s.beta;
	
	Matrix<int> clusters_b = s.clusters;
	Matrix<int>& clusters = s.clusters;
	Matrix<int> clusters_n(num_clone, num_muta);
	Matrix<double>& pie = s.pie;
	Matrix<double> pie_n(1, num_clone);
	s.post_probs.resize(num_cell, num_clone, false);
	Matrix<double>& post_probs = s.post_probs;
	Matrix<long double> probs(num_cell, num_clone);
	
	double max_alpha = config.getRealPara("max_alpha");
	double max_beta = config.getRealPara("max_beta");
	double eps = numeric_limits<long double>::epsilon();
	
	double ll, pre_ll = numeric_limits<long>::min();
	int iter = 0, max_iter = config.getIntPara("max_iter");
	while(iter < max_iter) {
		obsProbs(clusters, probs, alpha, beta);
		/*** E-step: posterior probability calculation ***/
		ll = 0;
		for(i = 0; i < num_cell; i++) {
			long double sum = 0;
			for(k = 0; k < num_clone; k++) {
				sum += probs[i*num_clone+k]*pie[k];
			}
			//cerr << sum << endl;
			for(k = 0; k < num_clone; k++) {
				post_probs[i*num_clone+k] = probs[i*num_clone+k]*pie[k]/sum;
			}
			ll += logl(sum+eps);
		}
		if(isnan(ll)) {
			alpha = 0.01;
			beta = 0.01;
			ll = numeric_limits<int>::min();
			post_probs.set(1.0/num_clone);
			clusters = clusters_b;
			pie.set(1.0/num_clone);
			break;
		}
		/*** M-step: parameter update ***/
		
		// update alpha and beta
		double tmp1 = 0, tmp2 = 0, tmp3 = 0, tmp4 = 0, tmp5 = 0;
		for(i = 0; i < num_cell; i++) {
			for(k = 0; k < num_clone; k++) {
				double gamma = post_probs[i*num_clone+k];
				for(j = 0; j < num_muta; j++) {
					if(obsData[i*num_muta+j] == 3) {
						continue;
					}
					int s = obsData[i*num_muta+j];
					if(homo_muta) {
						tmp1 += gamma*(1-s)/2*(-s+clusters[k*num_muta+j]*(2-s));
						tmp2 += gamma*(1-clusters[k*num_muta+j])*(1-s)*(2-s)/2;
						tmp3 += gamma*clusters[k*num_muta+j]*s*(2-s);
						tmp4 += gamma*(1-clusters[k*num_muta+j])*s*(3-s);
						tmp5 += gamma*(1-clusters[k*num_muta+j]);
					}
					else {
						tmp1 += gamma*(1-clusters[k*num_muta+j])*s;
						tmp2 += gamma*(1-clusters[k*num_muta+j]);
						tmp3 += gamma*clusters[k*num_muta+j]*(1-s);
						tmp4 += gamma*clusters[k*num_muta+j];
					}
				}
			}
		}
		double alpha_n, beta_n;
		if(homo_muta) {
			if(alpha <= eps) {
				beta_n = tmp1/(tmp1+tmp3);
			}
			else {
				double a = alpha*(tmp1+tmp2+tmp3);
				double b = tmp1*(alpha-2)-tmp2*alpha+2*tmp3*(alpha-1);
				double c = 2*tmp1*(1-alpha);
				beta_n = 0.5*(-b-sqrt(b*b-4*a*c))/(a+eps);
			}
			alpha_n = tmp4/((2+beta_n)*tmp5);
		}
		else {
			alpha_n = tmp1/(tmp2+eps);
			beta_n = tmp3/(tmp4+eps);
		}
		
		// update clusters
		if(para_updates[2] == 1) {
			for(k = 0; k < num_clone; k++) {
				for(j = 0; j < num_muta; j++) {
					double f1 = 0, f2 = 0;
					for(i = 0; i < num_cell; i++) {
						if(obsData[i*num_muta+j] == 3) {
							continue;
						}
						int s = obsData[i*num_muta+j];
						if(homo_muta) {
							double tmp1 = pow(1-s, 2)*log(beta/2+eps)+s*(2-s)*log(1-beta+eps);
							f1 += post_probs[i*num_clone+k]*tmp1;
							double tmp2 = (1-s)*(2-s)/2*log(1-alpha-alpha*beta/2+eps)+s*(2-s)*log(alpha+eps)+s*(s-1)/2*log(alpha*beta/2+eps);
							f2 += post_probs[i*num_clone+k]*tmp2;
						}
						else {
							f1 += post_probs[i*num_clone+k]*(s*log(1-beta+eps)+(1-s)*log(beta+eps));
							f2 += post_probs[i*num_clone+k]*(s*log(alpha+eps)+(1-s)*log(1-alpha+eps));
						}
						
					}
					clusters_n[k*num_muta+j] = (f1 > f2)? 1:0;
				}
			}
			clusters = clusters_n;
		}
		
		if(para_updates[0]) {
			alpha = isnan(alpha_n)? alpha:alpha_n;
		}
		if(para_updates[1]) {
			beta = isnan(beta_n)? beta:beta_n;
		}
		
		
		if(alpha > max_alpha) {
			alpha = max_alpha;
		}
		if(beta > max_beta) {
			beta = max_beta;
		}
		
		// update pie
		for(k = 0; k < num_clone; k++) {
			double sum = 0;
			for(i = 0; i < num_cell; i++) {
				sum += post_probs[i*num_clone+k];
			}
			pie_n[k] = sum/num_cell;
		}
		if(para_updates[3]) {
			pie = pie_n;
		}
		
		if(fabs(pre_ll-ll) < 0.01) {
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
	
	Matrix<int>& clusters = s.clusters;
	Matrix<double>& post_probs = s.post_probs;
	Matrix<double> results(num_cell, num_muta, 0.0);
	int num_clone = s.num_clone;
	for(i = 0; i < num_cell; i++) {
		for(j = 0; j < num_muta; j++) {
			for(k = 0; k < num_clone; k++)
				results[i*num_muta+j] += post_probs[i*num_clone+k]*clusters[k*num_muta+j];
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
			k = (results[i*num_muta+j] >= 0.5)? 1:0;
			if(k == realData[i*num_muta+j])	acc++;
		}
	}
	s.acc = acc/(n*num_muta);
}

void* HetCaller::evaluateSolutions(const void *arg) {
	int *tmp = (int*) arg;
	int s = tmp[0], e = tmp[1];
	while(s <= e) {
		hetcaller.evaluateSolution(hetcaller.candi_s1[s]);
		s++;
	}
	return NULL;
}

void HetCaller::evaluateSolution(Solution& s) {
	int i, j, k, n;
	int num_clone = s.num_clone;
	Matrix<int>& clusters = s.clusters;
	int num_cell = obsData.getROWS();
	int num_muta = obsData.getCOLS();
	Matrix<int> tmp = clusters.sumRows();
	int unmutated_count = 0;
	for(i = 0; i < num_muta; i++) {
		if(tmp[i] == 0) {
			unmutated_count++;
		}
	}
	if(1.0*unmutated_count/num_muta >= 0.1) {
		s.valid = false;
	}
	else {
		s.valid = true;
	}
	
	calculateIntraScore(s);
	calculateInterScore(s);
	
	double mean_intra_score = 0;
	for(i = 0; i < num_clone; i++) {
		mean_intra_score += s.intra_scores[i];
	}
	mean_intra_score /= num_clone;
	double sigma = 0;
	for(i = 0; i < num_clone; i++) {
		sigma += pow(s.intra_scores[i]-mean_intra_score, 2);
	}
	sigma = sqrt(sigma/num_clone);
	
	double mean_inter_score = 0;
	Matrix<double> scores(1, (num_clone*(num_clone-1))/2);
	if(num_clone > 1) {
		for(i = 0, k = 0; i < num_clone-1; i++) {
			for(j = i+1; j < num_clone; j++) {
				mean_inter_score += s.inter_scores[i*num_clone+j];
				scores[k++] = s.inter_scores[i*num_clone+j];
			}
		}
		mean_inter_score *= 2.0/(num_clone*(num_clone-1));
	}
	else {
		mean_inter_score = mean_intra_score;
	}
	
	s.score2 = mean_inter_score;
	evalAccuracy(s);
}

void HetCaller::calculateIntraScore(Solution& s) {
	int i, j, k, n, m;
	
	int num_clone = s.num_clone;
	double beta = s.beta;
	double alpha = s.alpha;
	Matrix<int>& clusters = s.clusters;
	Matrix<double>& post_probs = s.post_probs;
	int num_cell = post_probs.getROWS();
	int num_muta = clusters.getCOLS();
	
	s.intra_scores.resize(1, num_clone, false);
	
	Matrix<int> indxs;
	post_probs.max(2, indxs);
	map<int, vector<int>> c_indxs_a;
	for(j = 0; j < num_cell; j++) {
		c_indxs_a[indxs[j]].push_back(j);
	}
	
	double p1 = 2*beta*(1-beta)*pow(1-missing_rate,2)+2*missing_rate*(1-missing_rate);
	double p0 = 2*alpha*(1-alpha)*pow(1-missing_rate,2)+2*missing_rate*(1-missing_rate);
	
	for(i = 0; i < num_clone; i++) {
		int count1 = 0;
		for(k = 0; k < num_muta; k++) {
			if(clusters[i*num_muta+k] == 1)	count1++;
		}
		double expected_dist = count1*p1+(num_muta-count1)*p0;
		
		vector<int>& c_indxs = c_indxs_a[i];
		double score = 0;
		if(c_indxs.size() > 1) {
			double ave_dist = 0;
			int num_c = c_indxs.size();
			Matrix<double> dists(1, num_c*(num_c-1)/2);
			for(j = 0, m = 0; j < num_c-1; j++) {
				for(k = j+1; k < num_c; k++) {
					double dist = 0;
					for(n = 0; n < num_muta; n++) {
						if(obsData[c_indxs[j]*num_muta+n] != obsData[c_indxs[k]*num_muta+n]) {
							dist++;
						}
					}
					dists[m++] = dist;
					ave_dist += dist;
				}
			}
			ave_dist *= 2.0/(c_indxs.size()*(c_indxs.size()-1));
			score = exp(-pow(expected_dist-ave_dist, 2));
			//score = 1/ave_dist;
			//double med_dist = median(dists.getEntrance(), num_c*(num_c-1)/2);
			//score = exp(-pow(expected_dist-med_dist, 2));
		}
		
		s.intra_scores[i] = score;
	}
	
}

void HetCaller::calculateInterScore(Solution& s) {
	int c1, c2;
	
	int i, j, k, n;
	
	int num_clone = s.num_clone;
	double beta = s.beta;
	double alpha = s.alpha;
	Matrix<int>& clusters = s.clusters;
	Matrix<double>& post_probs = s.post_probs;
	int num_cell = post_probs.getROWS();
	int num_muta = clusters.getCOLS();
	
	s.inter_scores.resize(num_clone, num_clone, false);
	
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
	
	for(c1 = 0; c1 < num_clone-1; c1++) {
		for(c2 = c1+1; c2 < num_clone; c2++) {
			int counts[4] = {0};
			for(k = 0; k < num_muta; k++) {
				i = clusters[c1*num_muta+k]*2+clusters[c2*num_muta+k];
				counts[i]++;
			}
			double expected_dist = counts[0]*p_00+(counts[1]+counts[2])*p_01+counts[3]*p_11;
			
			vector<int>& c1_indxs = c_indxs_a[c1];
			vector<int>& c2_indxs = c_indxs_a[c2];
			
			if(c1_indxs.empty() || c2_indxs.empty()) {
				s.inter_scores[c1*num_clone+c2] = 0;
				s.inter_scores[c2*num_clone+c1] = 0;
			}
			else {
				int num_c1 = c1_indxs.size();
				int num_c2 = c2_indxs.size();
				Matrix<double> dists(1, num_c1*num_c2);
				double ave_dist = 0;
				for(i = 0, n = 0; i < num_c1; i++) {
					for(j = 0; j < num_c2; j++) {
						double dist = 0;
						for(k = 0; k < num_muta; k++) {
							if(obsData[c1_indxs[i]*num_muta+k] != obsData[c2_indxs[j]*num_muta+k]) {
								dist++;
							}
						}
						dists[n++] = dist;
						ave_dist += dist;
					}
				}
				ave_dist /= (c1_indxs.size()*c2_indxs.size());
				//s.inter_scores[c1*num_clone+c2] = ave_dist;
				s.inter_scores[c1*num_clone+c2] = exp(-pow(expected_dist-ave_dist, 2));
				//double med_dist = median(dists.getEntrance(), num_c1*num_c2);
				//s.inter_scores[c1*num_clone+c2] = exp(-pow(expected_dist-med_dist, 2));
				s.inter_scores[c2*num_clone+c1] = s.inter_scores[c1*num_clone+c2];
			}
			
		}
	}
}

double HetCaller::getObsProb(int gt, int observed, double alpha, double beta) {
	double prob = 1, eps = numeric_limits<double>::epsilon();
	if(gt == 1) {
		if(observed == 0) {
			prob = (homo_muta)? beta/2+eps:beta+eps;
		}
		else if(observed == 1) {
			prob = 1-beta+eps;
		}
		else if(observed == 2) {
			prob = beta/2+eps;
		}
	}
	else {
		if(observed == 0) {
			prob = (homo_muta)? 1-alpha-alpha*beta/2+eps:1-alpha+eps;
		}
		else if(observed == 1) {
			prob = alpha+eps;
		}
		else if(observed == 2) {
			prob = alpha*beta/2+eps;
		}
	}
	return prob;
}

void HetCaller::obsProbs(Matrix<int>& clusters, Matrix<long double>& probs, double alpha, double beta) {
	int i, j, k;
	int num_muta = obsData.getCOLS();
	int num_cell = obsData.getROWS();
	int num_clone = clusters.getROWS();
	long double eps = numeric_limits<long double>::epsilon();
	
	probs.resize(num_cell, num_clone, false);
	
	for(i = 0; i < num_cell; i++) {
		for(k = 0; k < num_clone; k++) {
			double log_prob = 0;
			long double prob = 1;
			for(j = 0; j < num_muta; j++) {
				if(clusters[k*num_muta+j] == 1) {
					if(obsData[i*num_muta+j] == 0) {
						log_prob += (homo_muta)? log(beta/2+eps):log(beta+eps);
						prob *= (homo_muta)? (beta/2+eps):(beta+eps);
					}
					else if(obsData[i*num_muta+j] == 1) {
						log_prob += log(1-beta+eps);
						prob *= (1-beta+eps);
					}
					else if(obsData[i*num_muta+j] == 2) {
						log_prob += log(beta/2+eps);
						prob *= (beta/2+eps);
					}
				}
				else {
					if(obsData[i*num_muta+j] == 0) {
						log_prob += (homo_muta)? log(1-alpha-alpha*beta/2+eps):log(1-alpha+eps);
						prob *= (homo_muta)? (1-alpha-alpha*beta/2+eps):(1-alpha+eps);
					}
					else if(obsData[i*num_muta+j] == 1) {
						log_prob += log(alpha+eps);
						prob *= (alpha+eps);
					}
					else if(obsData[i*num_muta+j] == 2) {
						log_prob += log(alpha*beta/2+eps);
						prob *= (alpha*beta/2+eps);
					}
				}
			}
			//probs[i*num_clone+k] = exp(log_prob);
			probs[i*num_clone+k] = prob;
		}
	}
}

void HetCaller::saveSolution(Solution& s) {
	pthread_mutex_lock(&pm);
	candi_s1.push_back(s);
	pthread_mutex_unlock(&pm);
}

Matrix<int> HetCaller::buildTree(Solution& s) {
	int i, j, k, n;
	int num_clone = s.num_clone;
	Matrix<int>& clusters = s.clusters;
	int num_muta = obsData.getCOLS();
	int num_cell = obsData.getROWS();
	
	/*** inter-cluster distance ***/
	Matrix<int> distances(num_clone, num_clone, 0);
	for(i = 0; i < num_clone-1; i++) {
		for(j = i+1; j < num_clone; j++) {
			int dist = 0;
			for(k = 0; k < num_muta; k++) {
				if(clusters[i*num_muta+k] != clusters[j*num_muta+k]) {
					dist++;
				}
			}
			distances[i*num_clone+j] = dist;
			distances[j*num_clone+i] = dist;
		}
	}
	
	/*** reconstruct clonal tree as a minimum spanning tree ***/
	Matrix<int> clonal_tree(1, num_clone);
	clonal_tree.set(-1);
	int root = 0;
	int min_dist = num_muta+1;
	for(i = 0; i < num_clone; i++) {
		int dist = 0;
		for(k = 0; k < num_muta; k++) {
			if(clusters[i*num_muta+k] != 0) {
				dist++;
			}
		}
		if(dist < min_dist) {
			min_dist = dist;
			root = i;
		}
	}
	clonal_tree[root] = num_clone;
	
	Matrix<int> max_num_childs(1, num_clone, num_clone);
	if(min_dist == 0) {
		max_num_childs[root] = 1;
	}
	while(1) {
		int s, t;
		min_dist = num_muta+1;
		for(i = 0; i < num_clone; i++) {
			if(clonal_tree[i] < 0 || max_num_childs[i] == 0) {
				continue;
			}
			for(j = 0; j < num_clone; j++) {
				if(clonal_tree[j] >= 0) {
					continue;
				}
				if(min_dist > distances[i*num_clone+j]) {
					min_dist = distances[i*num_clone+j];
					s = i;
					t = j;
				}
			}
		}
		if(min_dist > num_muta) {
			break;
		}
		clonal_tree[t] = s;
		max_num_childs[s]--;
	}
	
	return clonal_tree;
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
	//ofs << "K = " << candi_s[best_s_indx].num_clone << endl;
	ofs.close();
	*/
	
	int num_clone = candi_s[best_s_indx].num_clone;
	Matrix<int>& clusters = candi_s[best_s_indx].clusters;
	fn = outputPrefix+".clone_genotypes";
	ofs.open(fn.c_str());
	if(!ofs.is_open()) {
		cerr << "Error: cannot open file " << fn << endl;
		exit(-1);
	}
	for(i = 0; i < num_clone; i++) {
		for(j = 0; j < num_muta; j++) {
			if(j < num_muta-1) {
				ofs << clusters[i*num_muta+j] << '\t';
			}
			else {
				ofs << clusters[i*num_muta+j] << endl;
			}
		}
	}
	ofs.close();
	
	Matrix<double>& post_probs = candi_s[best_s_indx].post_probs;
	Matrix<int> indxs;
	post_probs.max(2, indxs);
	map<int, string> clone_labels;
	map<int, int> counts;
	//int cells_per_line = 8;
	int cells_per_line = max(5, (int) sqrt(num_cell/num_clone));
	for(i = 0; i < num_cell; i++) {
		k = indxs[i];
		if(clone_labels[k].empty()) {
			clone_labels[k] = cLabels[i];
			counts[k] = 1;
		}
		else {
			counts[k]++;
			if(counts[k] == cells_per_line) {
				counts[k] = 0;
				clone_labels[k] = clone_labels[k]+" "+cLabels[i]+"\\n";
			}
			else {
				if(counts[k] == 1) {
					clone_labels[k] = clone_labels[k]+cLabels[i];
				}
				else {
					clone_labels[k] = clone_labels[k]+" "+cLabels[i];
				}
			}
		}
	}
	
	map<int, string> node_labels;
	for(i = 0; i < num_clone; i++) {
		char buf[100];
		sprintf(buf, "subclone%d", i+1);
		node_labels[i] = buf;
	}
	
	fn = outputPrefix+".dot";
	ofs.open(fn.c_str());
	if(!ofs.is_open()) {
		cerr << "Error: cannot open file " << fn << endl;
		exit(-1);
	}
	
	ofs << "digraph T {" << endl;
	ofs << "node [color=deeppink4, fontcolor=black, penwidth=2];" << endl;
	for(i = 0; i < num_clone; i++) {
		ofs << i << " [label=\"" << node_labels[i] << "\"];" << endl;
	}
	ofs << "node [color=lightgrey, fontcolor=black, penwidth=2.5];" << endl;
	for(i = 0; i < num_clone; i++) {
		if(!clone_labels[i].empty()) {
			ofs << i+num_clone << " [label=\"" << clone_labels[i] << "\"];" << endl;
		}
	}
	
	ofs << "edge [penwidth=1.5];" << endl;
	for(i = 0; i < num_clone; i++) {
		int p = clonal_tree[i];
		if(p != num_clone) {
			ofs << clonal_tree[i] << "->" << i << ";" << endl;
		}
	}
	for(i = 0; i < num_clone; i++) {
		if(!clone_labels[i].empty()) {
			ofs << i << "->" << i+num_clone << ";" << endl;
		}
	}
	
	ofs << "}" << endl;
	
	
	ofs.close();
	
}


