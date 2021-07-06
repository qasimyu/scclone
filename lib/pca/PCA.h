#ifndef _PCA_H
#define _PCA_H

#include <cmath>
#include <eigen3/Eigen/SVD>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Eigen>

#include "Matrix.h"

using namespace std;
using Eigen::MatrixXf;
using Eigen::VectorXf;

class PCA {
	public:
		PCA(Matrix<float>& data);
		virtual ~PCA(){};
    
		void first_k_ONB(const int k);
		void important_ONB(const double energy_percent);
		void sole_k_ONB(const int k);
		
		VectorXf Proj2LowDim(VectorXf OriginalVector);
		MatrixXf Projection(MatrixXf data);
		Matrix<float> getProj();
		VectorXf Reconstruction();
		MatrixXf Whitening();
		
		VectorXf GetOriginalVector(const int k);
    
	private:
    
		MatrixXf OriginalData;
		MatrixXf ONB;
		
		VectorXf EigenVector;
		
		VectorXf Average;
		VectorXf proj;

		int n_Samples;
		int n_Dimensions;
		int n_Dimensions_p;
};

#endif
