// ***************************************************************************
// HetCaller.h (c) 2020 Zhenhua Yu <qasim0208@163.com>
// Health Informatics Lab, Ningxia University
// All rights reserved.

#ifndef _HETCALLER_H
#define _HETCALLER_H

#include <vector>
#include <map>
#include <pthread.h>

#include "Matrix.h"

class Solution {
	public:
		Solution() {}
		Solution(int num_cluster, double alpha, double beta, Matrix<double>& pie, Matrix<int>& states) :
			num_cluster(num_cluster), alpha(alpha), beta(beta), pie(pie), states(states) {}
		
		int num_cluster; //number of clusters
		double ll; //log-likelihood
		double score;
		bool valid;
		double acc; //accuracy
		double alpha; //false positive rate
		double beta; //false negative rate
		Matrix<int> states; //mutation states of clusters
		Matrix<double> pie; //occurrence probability of each clone
		Matrix<double> post_probs; //posterior probability
		Matrix<double> scores; //inter-clone score
};

class HetCaller {
	private:
		Matrix<int> obsData;
		Matrix<int> uniqueData;
		Matrix<int> realData;
		
		double missing_rate;
		vector<int> doublet_indxs;
		
		int num_cluster_c;
		vector<Solution> candi_s;
		vector<Solution> candi_s1;
		int best_s_indx;
		
		pthread_mutex_t pm;
		
		void preProcess();
		
		void predict(int num_cluster);
		static void* inferSolution(const void *arg);
		void inferParas(Solution& s, const vector<int>& para_updates);
		
		void obsProbs(Matrix<int>& states, Matrix<long double>& probs, double alpha, double beta);
		void saveSolution(Solution& s);
		
		void printSolution(Solution& s);
		
		void loadMutationData();
		void loadRealData();
		
		void saveResults();
		void evalAccuracy(Solution& s);
		void evaluateSolution(Solution& s);
		
		void calculateScore(Solution& s);
		
	public:
		HetCaller();
	
		Matrix<int>& getObsData() {return obsData;}
		Matrix<int>& getUniqueData() {return uniqueData;}
		int getNumOfClusters() {return num_cluster_c;}
		
		Solution& getSolution(int i) {return candi_s[i];}
		
		void loadData();
		void call();
};

#endif
