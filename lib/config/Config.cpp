// ***************************************************************************
// Config.cpp (c) 2019 Zhenhua Yu <qasim0208@163.com>
// Health Informatics Lab, Ningxia University
// All rights reserved.

#include <iostream>
#include <fstream>
#include <cassert>
#include <cmath>
#include <limits>

#include "Config.h"

Config::Config() {
	string strParaNames[] = {"input", "rinput", "output", "binPath", "mlabel", "clabel"};
	
	/*---start default configuration---*/
	
	int n = sizeof(strParaNames)/sizeof(string);
	for(int i = 0; i < n; i++) {
		stringParas.insert(make_pair(strParaNames[i], ""));
	}
	
	intParas.insert(make_pair("maxc", -1));
	realParas.insert(make_pair("alpha", -1));
	realParas.insert(make_pair("beta", -1));
	
	intParas.insert(make_pair("max_iter", 30));
	intParas.insert(make_pair("threads", 100));
	intParas.insert(make_pair("verbose", 1));
	
	realParas.insert(make_pair("eps", numeric_limits<double>::epsilon()));
	
	/*---end default configuration---*/
}

string Config::getStringPara(string paraName) {
	if(stringParas.find(paraName) != stringParas.end()) {
		return stringParas[paraName];
	}
	else {
		cerr << "Error: unrecognized parameter name \"" << paraName << "\"" << endl;
		exit(1);
	}
	return "";
}

void Config::setStringPara(string paraName, string value) {
	if(stringParas.find(paraName) != stringParas.end()) {
		stringParas[paraName] = value;
	}
	else {
		cerr << "Warning: unrecognized parameter name \"" << paraName << "\"" << endl;
	}
}

long Config::getIntPara(string paraName) {
	if(intParas.find(paraName) != intParas.end()) {
		return intParas[paraName];
	}
	else {
		cerr << "Error: unrecognized parameter name \"" << paraName << "\"" << endl;
		exit(1);
	}
}

void Config::setIntPara(string paraName, long value) {
	if(intParas.find(paraName) != intParas.end()) {
		intParas[paraName] = value;
	}
	else {
		cerr << "Warning: unrecognized parameter name \"" << paraName << "\"" << endl;
	}
}

double Config::getRealPara(string paraName) {
	if(realParas.find(paraName) != realParas.end()) {
		return realParas[paraName];
	}
	else {
		cerr << "Error: unrecognized parameter name \"" << paraName << "\"" << endl;
		exit(1);
	}
}

void Config::setRealPara(string paraName, double value) {
	if(realParas.find(paraName) != realParas.end()) {
		realParas[paraName] = value;
	}
	else {
		cerr << "Warning: unrecognized parameter name \"" << paraName << "\"" << endl;
	}
}


