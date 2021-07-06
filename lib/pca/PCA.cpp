#include "PCA.h"

PCA::PCA(Matrix<float>& data) {
    /*
     *  Description:
	 *  convert data to MatrixXf
     *  @param data: data matrix
     *
     */
	n_Samples = data.getROWS();
	n_Dimensions = data.getCOLS();
    OriginalData = Eigen::MatrixXf(n_Samples, n_Dimensions);
    for(int i = 0; i < n_Samples; i++) {
		for(int j = 0; j < n_Dimensions; j++) {
			OriginalData(i, j) = data[i*n_Dimensions+j];
		}
	}
    
    /* Make the variable mean to be zero */
    Average = OriginalData.row(0);
    for (int i = 1; i < n_Samples; i++)
    {
        Average += OriginalData.row(i);
    }
    Average /= n_Samples;
    
    VectorXf ones(n_Samples);
    ones.setOnes();
    OriginalData -= ones * Average.transpose();
}

void PCA::first_k_ONB(const int k) 
{
    /*
     *  Description:
     *  Find the first k largest eigen values and the corresponding
     *  eigen vector
     *
     *  @param k: first k eigen vector
     *
     */
    
    if (k > n_Dimensions || k > n_Samples)
    {
        cerr << "Too hign dimensions\n";
    }
    
    Eigen::JacobiSVD<MatrixXf> svd(OriginalData, Eigen::ComputeThinU | Eigen::ComputeThinV);
    
    EigenVector = svd.singularValues();
    cout << "The accuracy is " << EigenVector.head(k).sum() / EigenVector.sum() * 1.0 << endl;
    
    ONB =  svd.matrixV().block(0, 0, n_Dimensions, k);
}

void PCA::important_ONB(const double energy_percent) 
{
    /*
     *  Description:
     *  Find the first k largest eigen values and the corresponding
     *  eigen vector that account for a given percent of energy
     *
     *  @param energy_percent: percent of energy to keep
     *
     */
    
    if (energy_percent < 0)
    {
        cerr << "input error\n";
    }
    
    Eigen::JacobiSVD<MatrixXf> svd(OriginalData, Eigen::ComputeThinU | Eigen::ComputeThinV);
    
    EigenVector = svd.singularValues();
	
	int i;
	for(i = 1; i <= n_Dimensions && i <= n_Samples; i++) {
		//cerr << EigenVector(i-1) << endl;
		if(EigenVector.head(i).sum()/EigenVector.sum() > energy_percent) {
			//cerr << "The accuracy is " << EigenVector.head(i).sum()/EigenVector.sum() << endl;
			break;
		}
	}
	n_Dimensions_p = i;
	ONB =  svd.matrixV().block(0, 0, n_Dimensions, i);
}

void PCA::sole_k_ONB(const int k)
{
    /*
     *  Description:
     *  Find the k-th eigen vector
     *
     *  @param k: k-th eigen vector
     *
     */

    if (k > n_Dimensions || k > n_Samples)
    {
        cerr << "Too hign dimensions\n";
    }
    Eigen::JacobiSVD<MatrixXf> svd(OriginalData, Eigen::ComputeThinU | Eigen::ComputeThinV);
    
    EigenVector = svd.singularValues();
    cout << "The accuracy is " << EigenVector(k - 1) / EigenVector.sum() * 1.0 << endl;
    
    ONB =  svd.matrixV().block(0, k - 1, n_Dimensions, 1);
}

VectorXf PCA::Proj2LowDim(VectorXf OriginalVector)
{
    /*
     *  Description:
     *  Mapping the original vector to a lower space
     *
     *  @param OriginalVector: Could be one of the sample or test
     *                          set
     */
    
    proj = ONB.transpose() * OriginalVector;
    return proj;
}

MatrixXf PCA::Projection(MatrixXf data)
{
    /*
     *  Description:
     *  Mapping the original data to a lower space
     *
     *  @param data: a matrix containing original data
     */
    
    proj = data * ONB;
    return proj;
}

Matrix<float> PCA::getProj()
{
    /*
     *  Description:
     *  get the projected data
     *
     */
    Matrix<float> ret(n_Samples, n_Dimensions_p);
    MatrixXf tmp = OriginalData * ONB;
	for(int i = 0; i < n_Samples; i++) {
		for(int j = 0; j < n_Dimensions_p; j++) {
			ret[i*n_Dimensions_p+j] = tmp(i, j);
		}
	}
    return ret;
}

VectorXf PCA::Reconstruction()
{
    /*
     *  Description:
     *  Mapping the projected vector back to the original space
     *
     */
    return ONB * proj + Average;
}

MatrixXf PCA::Whitening()
{
    /*
     *  Description:
     *  Eliminate the correlation between the data
     *
     */
    
    first_k_ONB(n_Dimensions);

    MatrixXf whitenData(OriginalData.rows(), OriginalData.cols());
    for (int i = 0; i < n_Samples; i++)
    {
        for (int j = 0; j < n_Dimensions; j++)
        {
            whitenData(i, j) = (Proj2LowDim(OriginalData.row(i)))(j) / sqrt(EigenVector(j));
        }
    }
    
    return whitenData;
}

VectorXf PCA::GetOriginalVector(const int k)
{
    /*
     *  Description:
     *  Get the k-th original vector
     *
     */
    assert(k < n_Samples);
    return OriginalData.row(k);
}
