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
		Solution(int num_clone, double alpha, double beta, Matrix<double>& pie, Matrix<int>& clusters) :
			num_clone(num_clone), alpha(alpha), beta(beta), pie(pie), clusters(clusters) {}
		
		int num_clone; //number of clones
		double ll; //log-likelihood
		double bic; //BIC value
		double score, score1, score2;
		bool valid;
		double acc; //accuracy
		double alpha; //false positive rate
		double beta; //false negative rate
		Matrix<double> pie; //occurrence probability of each clone
		Matrix<int> clusters; //mutation states of clones
		Matrix<double> post_probs; //posterior probability
		Matrix<double> intra_scores; //intra-clone score
		Matrix<double> inter_scores; //inter-clone score
};

class HetCaller {
	private:
		Matrix<int> obsData;
		Matrix<int> uniqueData;
		Matrix<float> replicate;
		Matrix<int> realData;
		vector<string> cLabels;
		vector<string> mLabels;
		bool homo_muta;
		double missing_rate;
		double positive_rate;
		
		vector<int> doublet_indxs;
		
		//k-means
		vector<vector<float>> data_for_cluster;
		float min_dist; 
		vector<uint32_t> best_indices;
		vector<Matrix<int>> clusters_m;
		
		Matrix<int> para_updates;
		vector<Solution> candi_s;
		vector<Solution> candi_s1;
		int best_s_indx;
		
		Matrix<int> clonal_tree;
		
		pthread_mutex_t pm;
		
		void fetchUniqueData();
		
		static void* clusterCells(const void *arg);
		void saveClusterResults(float dist, float score, vector<uint32_t>& indices);
		
		void predict(int num_clone);
		static void* inferSolution(const void *arg);
		void inferParas(Solution& s, const vector<int>& para_updates);
		
		void obsProbs(Matrix<int>& clusters, Matrix<long double>& probs, double alpha, double beta);
		double getObsProb(int gt, int observed, double alpha, double beta);
		void saveSolution(Solution& s);
		
		void printSolution(Solution& s);
		void calculateIntraScore(Solution& s);
		void calculateInterScore(Solution& s);
		
		void loadMutationData();
		void loadRealData();
		void loadCellLabels();
		void loadMutaLabels();
		
		Matrix<int> buildTree(Solution& s);
		void saveResults();
		
		static void* evaluateSolutions(const void *arg);
		void evaluateSolution(Solution& s);
		void evalAccuracy(Solution& s);
		
	public:
		HetCaller();
	
		Matrix<int>& getObsData() {return obsData;}
		Matrix<int>& getUniqueData() {return uniqueData;}
		vector<vector<float>>& getDataForCluster() {return data_for_cluster;}
		vector<uint32_t>& getBestIndices() {return best_indices;}
		
		Matrix<int>& getClusters() {return candi_s[best_s_indx].clusters;}
		Solution& getSolution(int i) {return candi_s[i];}
		Matrix<int>& getParaUpdates() {return para_updates;}
		vector<Matrix<int>>& getClusters_m() {return clusters_m;}
		
		bool getHomoMuta() {return homo_muta;}
		double getMissingRate() {return missing_rate;}
		double getPositiveRate() {return positive_rate;}
		
		void loadData();
		void call();
};

#endif
