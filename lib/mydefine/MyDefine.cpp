// ***************************************************************************
// MyDefine.cpp (c) 2020 Zhenhua Yu <qasim0208@163.com>
// Health Informatics Lab, Ningxia University
// All rights reserved.

#include <cstdio>
#include <fstream>
#include <vector>
#include <cmath>
#include <limits>

#include "MyDefine.h"

using namespace std;

MyDefine::MyDefine() {
}

MyDefine::~MyDefine() {
}

string current_version = "1.0";

/*** definition of global vars ***/
Config config;
InputParser parser;
HetCaller hetcaller;
ThreadPool *threadpool;
/*** end of definition ***/

/*** definition of general functions ***/
void norm_trans(Matrix<double> &T, double thres) {
	double *p1 = T.getEntrance();
	double ZERO_FINAL = 2.2204e-16;
	int rows = T.getROWS();
	int cols = T.getCOLS();
	for(int i = 0; i < rows; i++) {
		Matrix<double> temp = T.Row(i);
		double tmp = temp.sum();
		if(p1[i*cols+i] < thres*tmp) {
			temp.set(0, i, 0);
			Matrix<double> a = temp*((1-thres)/temp.sum());
			a.set(0, i, thres);
			T.setRow(i, a);
		}
		else {
			Matrix<double> a = temp/(tmp+ZERO_FINAL);
			T.setRow(i, a);
		}
	}
}

//****** pdf of normal distribution ******//
double normpdf(double x, double mu, double sigma) {
	double PI = 3.1415926;
	return exp(-pow(x-mu,2)/(2*pow(sigma,2)))/(sqrt(2*PI)*sigma);
}

//****** pdf of student't distribution ******//
double stupdf(double x, double mu, double nu, double sigma) {
	double PI = 3.1415926;
	double ret = gammaln((nu+1)*0.5)-(nu+1)*0.5*log(1+pow(x-mu,2)/(nu*sigma*sigma))
				-gammaln(nu*0.5)-0.5*log(nu)-log(sigma)-0.5*log(PI);
	ret = exp(ret);
	return ret;
}

//****** pdf of binomial distribution ******//
double binopdf(double x, double n, double p) {
	double ret = gammaln(n+1)-gammaln(x+1)-gammaln(n-x+1)+x*log(p)+(n-x)*log(1-p);
	ret = exp(ret);
	return ret;
}

//****** pdf of beta-binomial distribution ******//
double betabinopdf(double x, double n, double alpha, double beta) {
	double ret = betaln(x+alpha, n-x+beta)-betaln(alpha, beta)+gammaln(n+1)-gammaln(x+1)-gammaln(n-x+1);
	ret = exp(ret);
	return ret;
}

//****** pdf of negative binomial distribution ******//
double negabinopdf(double x, double lambda, double p) {
	double ret = gammaln(x+lambda*(1-p)/p)+(lambda*(1-p)/p)*log(1-p)+x*log(p)
				-gammaln(x+1)-gammaln(lambda*(1-p)/p);
	ret = exp(ret);
	return ret;
}

//****** pdf of gamma distribution ******//
double gammapdf(double x, double k, double theta) {
	return exp((k-1)*log(x)-x/theta-gammaln(k)-k*log(theta));
}

double normrand(double mu, double sigma) {
	double eps = numeric_limits<double>::epsilon();
    double pi = 3.14159265358979323846;

    static double z0, z1;
    static bool flag = true;

    if(!flag) {
       return z1*sigma+mu;
	}
	flag = !flag;
	
    double u1, u2;
    do {
       u1 = rand()*(1.0/RAND_MAX);
       u2 = rand()*(1.0/RAND_MAX);
	}while(u1 <= eps);

    z0 = sqrt(-2.0*log(u1))*cos(2*pi*u2);
    z1 = sqrt(-2.0*log(u1))*sin(2*pi*u2);
	
    return z0*sigma+mu;
}

//****** trim string ******//
string trim(const string &str, const char *charlist) {
	string ret(str);
	size_t indx = ret.find_first_not_of(charlist);
	if(indx != string::npos) {
		ret.erase(0, indx);
		indx = ret.find_last_not_of(charlist);
		ret.erase(indx+1);
	}
	else {
		ret.erase();
	}
	return ret;
}

//****** abbreviation of chromosome name ******//
string abbrOfChr(string chr) {
	string abbr_chr = chr;
	size_t i = abbr_chr.find("chrom");
	if(i == string::npos) {
		i = abbr_chr.find("chr");
		if(i != string::npos) {
			abbr_chr = abbr_chr.substr(i+3,abbr_chr.size()-3);
		}
	}
	else {
		abbr_chr = abbr_chr.substr(i+5,abbr_chr.size()-5);
	}
	return abbr_chr;
}

