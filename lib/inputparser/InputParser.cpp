// ***************************************************************************
// InputParser.cpp (c) 2020 zhenhua yu <qasim0208@163.com>
// Health Informatics Lab, Ningxia University
// All rights reserved.

#include <iostream>
#include <string>
#include <unistd.h>
#include <getopt.h>
#include <sys/stat.h>
#include <limits.h>

#include "InputParser.h"
#include "MyDefine.h"

using namespace std;

void InputParser::parseArgs(int argc, char *argv[]) {
	string inputFile = "", outputPrefix = "";
	int maxc = -1;
	double alpha = -1, beta = -1;

	struct option long_options[] = {
		{"help", no_argument, 0, 'h'},
		{"version", no_argument, 0, 'v'},
		{"input", required_argument, 0, 'i'},
		{"output", required_argument, 0, 'o'},
		{"maxc", required_argument, 0, 'K'},
		{"threads", required_argument, 0, 't'},
		{"alpha", required_argument, 0, 'a'},
		{"beta", required_argument, 0, 'b'},
		{0, 0, 0, 0}
	};

	int c;
	//Parse command line parameters
	while((c = getopt_long(argc, argv, "hvi:o:K:a:b:", long_options, NULL)) != -1){
		switch(c){
			case 'h':
				usage(argv[0]);
				exit(0);
			case 'v':
				cerr << "SCClone version " << current_version << endl;
				exit(0);
			case 'i':
				inputFile = optarg;
				break;
			case 'o':
				outputPrefix = optarg;
				break;
			case 'K':
				maxc = atoi(optarg);
				break;
			case 'a':
				alpha = atof(optarg);
				break;
			case 'b':
				beta = atof(optarg);
				break;
			default :
				usage(argv[0]);
				exit(1);
		}
	}
	
	if(inputFile.empty()){
        cerr << "Use --input to specify the file containing mutation data." << endl;
		usage(argv[0]);
        exit(1);
    }
	
	if(outputPrefix.empty()){
		cerr << "Use --output to specify the prefix of result file names." << endl;
		usage(argv[0]);
		exit(1);
	}
	
	if(alpha > 1) {
		cerr << "Error: the value of parameter \"alpha\" should be in [0, 1]." << endl;
		usage(argv[0]);
		exit(1);
	}
	
	if(beta > 1) {
		cerr << "Error: the value of parameter \"beta\" should be in [0, 1]." << endl;
		usage(argv[0]);
		exit(1);
	}
	
	config.setStringPara("input", inputFile);
	config.setStringPara("output", outputPrefix);
	config.setIntPara("maxc", maxc);
	config.setRealPara("alpha", alpha);
	config.setRealPara("beta", beta);
	
	/*** check output directory ***/
	size_t i = outputPrefix.find_last_of('/');
	string outputDir = outputPrefix;
	if(i != string::npos) {
		outputDir = outputPrefix.substr(0, i);
	}
	bool is_exist = (access(outputDir.c_str(), F_OK) == 0);
	if(!is_exist) {
		cerr << "Error: the output directory " << outputDir << " does not exist!" << endl;
		exit(1);
	}
	//string cmd = "test ! -e "+outputDir+" && mkdir -m 755 -p "+outputDir;
	//system(cmd.c_str());
	
	/*** create thread pool ***/
	threadpool = new ThreadPool(config.getIntPara("threads"));
	threadpool->pool_init();
}

void InputParser::usage(const char* app) {
	cerr << "Usage: " << app << " [options]" << endl
		<< endl
		<< "Options:" << endl
		<< "    -h, --help                      give this information" << endl
		<< "    -v, --version                   print software version" << endl
		<< "    -i, --input <string>            input file containing mutation data" << endl
		<< "    -o, --output <string>           prefix of output file names" << endl
		<< "    -K, --maxc <int>                maximum number of clones to consider" << endl
		<< "    -a, --alpha <double>            false positive rate [default:inferred from data]" << endl
		<< "    -b, --beta <double>             false negative rate [default:inferred from data]" << endl
		<< endl
		<< "Example:" << endl
		<< app << " -i ./testdata/example.txt -o ./testdata/example" << endl
		<< endl
		<< "Author: Zhenhua Yu <qasim0208@163.com>\n" << endl;
}
