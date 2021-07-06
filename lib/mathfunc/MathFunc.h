// ***************************************************************************
// MathFunc.h (c) 2019 Zhenhua Yu <qasim0208@163.com>
// Health Informatics Lab, Ningxia University
// All rights reserved.

#ifndef _MATHFUNC_H
#define _MATHFUNC_H

/*** digamma function ***/
double digamma(double x);

/*** trigamma function ***/
double trigamma(double x);

/*** polygamma function ***/
double psi(int degree, double x);

/***  log gamma function ***/
double gammaln(double x);

/***  log beta function ***/
double betaln(double z, double w);

#endif

