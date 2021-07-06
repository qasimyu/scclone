// ***************************************************************************
// MyDefine.h (c) 2020 Zhenhua Yu <qasim0208@163.com>
// Health Informatics Lab, Ningxia University
// All rights reserved.

#ifndef _MYDEFINE_H
#define _MYDEFINE_H

#include <iostream>
#include <cstring>
#include <cstdlib>

#include "Config.h"
#include "InputParser.h"
#include "MathFunc.h"
#include "HetCaller.h"
#include "Matrix.h"
#include "ThreadPool.h"

class MyDefine {
	public:
		MyDefine();
		~MyDefine();
};

/*** declaration of global vars ***/
extern string current_version;
extern Config config;
extern InputParser parser;
extern HetCaller hetcaller;
extern ThreadPool *threadpool;
/*** end of declaration ***/

/*** declaration of general functions ***/
void norm_trans(Matrix<double> &T, double thres);

double normpdf(double x, double mu, double sigma);
double stupdf(double x, double mu, double nu, double sigma);
double binopdf(double x, double n, double p);
double betabinopdf(double x, double n, double alpha, double beta);
double negabinopdf(double x, double lambda, double p);
double gammapdf(double x, double k, double theta);

//****** quick sort ******//
template <class numtype>
void Qsort(numtype *a, int low, int high){
    if(high <= low) {
		return;
	}
    int i = low;
    int j = high + 1;
    numtype key = a[low];
    while(true) {
        while(a[++i] < key) {
            if(i == high) {
                break;
            }
        }
        while(a[--j] > key) {
            if(j == low) {
                break;
            }
        }
        if(i >= j) {
			break;
		}
        numtype temp = a[i];
        a[i] = a[j];
        a[j] = temp;
    }
    numtype temp = a[low];
    a[low] = a[j];
    a[j] = temp;
    Qsort(a, low, j-1);
    Qsort(a, j+1, high);
}

template <class numtype>
void Qsort(numtype *a, int *indxs, int low, int high){
    if(high <= low) {
		return;
	}
    int i = low;
    int j = high+1;
    numtype pivot = a[low];
	int pivot_indx = indxs[low];
    while(true) {
        while(a[++i] < pivot) {
            if(i == high) {
                break;
            }
        }
        while(a[--j] > pivot) {
            if(j == low) {
                break;
            }
        }
        if(i >= j) {
			break;
		}
        numtype temp = a[i];
        a[i] = a[j];
        a[j] = temp;
		int tmp_indx = indxs[i];
		indxs[i] = indxs[j];
		indxs[j] = tmp_indx;
    }
    a[low] = a[j];
    a[j] = pivot;
	indxs[low] = indxs[j];
	indxs[j] = pivot_indx;
    Qsort(a, indxs, low, j-1);
    Qsort(a, indxs, j+1, high);
}

template <class numtype>
int* QuickSort(numtype *a, int n, bool retIndx){
	if(retIndx) {
		int *indxs = new int[n];
		for(int i = 0; i < n; i++) {
			indxs[i] = i;
		}
		Qsort(a, indxs, 0, n-1);
		return indxs;
	}
	else {
		Qsort(a, 0, n-1);
		return NULL;
	}
}

//****** calculate median value ******//
template <class numtype>
numtype median(numtype *p, int n) {
	QuickSort(p, n, false);
	if(n%2 != 0) {
		return p[n/2];
	}
	else {
		return (p[n/2]+p[n/2-1])/2;
	}
}

//****** calculate median value of a matrix ******//
template <class numtype>
numtype median(Matrix<numtype> &mat) {
	int n = mat.getROWS()*mat.getCOLS();
	Matrix<numtype> T = mat;
	numtype *p = T.getEntrance();
	return median(p, n);
}

//****** calculate median value of a matrix along a dimension ******//
template <class numtype>
Matrix<numtype> median(Matrix<numtype> &mat, int dim) {
	if(dim == 0) {
		Matrix<numtype> ret(1, mat.getCOLS());
		for(int i = 0; i < mat.getCOLS(); i++) {
			Matrix<numtype> T = mat.Col(i);
			numtype *p = T.getEntrance();
			ret.set(0, i, median(p, mat.getROWS()));
		}
		return ret;
	}
	else if(dim == 1) {
		Matrix<numtype> ret(mat.getROWS(), 1);
		for(int i = 0; i < mat.getROWS(); i++) {
			Matrix<numtype> T = mat.Row(i);
			numtype *p = T.getEntrance();
			ret.set(i, 0, median(p, mat.getCOLS()));
		}
		return ret;
	}
	else {
		cerr << "Error: unknown dimension \"" << dim << "\"" << endl;
		exit(1);
	}
}

//****** calculate median value of a vector ******//
template <class numtype>
numtype median(vector<numtype> &x) {
	int n = x.size();
	numtype *p = new numtype[n];
	for(int i = 0; i < n; i++) {
		p[i] = x[i];
	}
	numtype ret = median(p, n);
	delete[] p;
	return ret;
}

//****** calculate standard deviation ******//
template <class numtype>
double stdev(numtype *p, int n, int flag) {
	double average = 0;
	int i;
	for(i = 0; i < n; i++) {
		average += p[i];
	}
	average /= n;
	double sigma = 0;
	for(i = 0; i < n; i++) {
		sigma += pow(p[i]-average,2);
	}
	if(flag == 0) {
		sigma = sqrt(sigma/(n-1));
	}
	else {
		sigma = sqrt(sigma/n);
	}
	return sigma;
}

//****** calculate variance ******//
template <class numtype>
double var(vector<numtype> &x) {
	int n = x.size();
	numtype *p = new numtype[n];
	for(int i = 0; i < n; i++) {
		p[i] = x[i];
	}
	double s = stdev(p, n, 1);
	delete[] p;
	return s*s;
}

//****** minimum element of an array ******//
template <class numtype>
numtype min(numtype *p, int n,int &indx) {
	numtype minValue = p[0];
	indx = 0;
	for(int i = 1; i < n; i++) {
		if(p[i] < minValue) {
			minValue = p[i];
			indx = i;
		}
	}
	return minValue;
}

//****** maximum element of an array ******//
template <class numtype>
numtype max(numtype *p, int n, int &indx) {
	numtype maxValue = p[0];
	indx = 0;
	for(int i = 1; i < n; i++) {
		if(p[i] > maxValue) {
			maxValue = p[i];
			indx = i;
		}
	}
	return maxValue;
}

//****** percentile ******//
template <class numtype>
numtype prctile(numtype *p, int n, double perc) {
	int i, k;
	numtype result;
	if(perc < 0) {
		perc = 0;
	}
	if(perc > 100) {
		perc = 100;
	}
	numtype *tmp = new numtype[n+2];
	for(i = 0; i < n; i++) {
		tmp[i+1] = p[i];
	}
	QuickSort(tmp+1, n, false);
	tmp[0] = tmp[1];
	tmp[n+1] = tmp[n];
	
	vector<double> q;
	double pro = -0.5;
	q.push_back(0);
	while(true) {
		pro += 1;
		if(pro <= n-0.5) {
			q.push_back(pro/n*100);
		}
		else {
			break;
		}
		k = q.size();
		if(perc >= q[k-2] && perc <= q[k-1]) { //linear interpolation
			result = (perc-q[k-2])*(tmp[k-1]-tmp[k-2])/(q[k-1]-q[k-2])+tmp[k-2];
			break;
		}
	}
	if(pro > n-0.5) {
		
		result = tmp[n];
	}
	
	delete[] tmp;
	return result;
}

string trim(const string &str, const char *charlist = " \t\r\n");
string abbrOfChr(string chr);

/*** end of declaration ***/


#endif

