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
	string realInputFile = "";
	string clabelFile = "", mlabelFile = "";
	int threads = 1;
	int maxc = -1;
	double alpha = -1, beta = -1;
	double max_alpha = 0.05, max_beta = 0.5;

	struct option long_options[] = {
		{"help", no_argument, 0, 'h'},
		{"version", no_argument, 0, 'v'},
		{"input", required_argument, 0, 'i'},
		//{"rinput", required_argument, 0, 'r'},
		{"output", required_argument, 0, 'o'},
		{"clabel", required_argument, 0, 'c'},
		{"mlabel", required_argument, 0, 'm'},
		{"maxc", required_argument, 0, 'K'},
		{"threads", required_argument, 0, 't'},
		{"alpha", required_argument, 0, 'a'},
		{"beta", required_argument, 0, 'b'},
		{"max_alpha", required_argument, 0, 'A'},
		{"max_beta", required_argument, 0, 'B'},
		{0, 0, 0, 0}
	};

	int c;
	//Parse command line parameters
	while((c = getopt_long(argc, argv, "hvi:r:o:c:m:K:t:a:b:A:B:", long_options, NULL)) != -1){
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
			/*
			case 'r':
				realInputFile = optarg;
				break;
			*/
			case 'o':
				outputPrefix = optarg;
				break;
			case 'c':
				clabelFile = optarg;
				break;
			case 'm':
				mlabelFile = optarg;
				break;
			case 'K':
				maxc = atoi(optarg);
				break;
			case 't':
				threads = atoi(optarg);
				break;
			case 'a':
				alpha = atof(optarg);
				break;
			case 'b':
				beta = atof(optarg);
				break;
			case 'A':
				max_alpha = atof(optarg);
				break;
			case 'B':
				max_beta = atof(optarg);
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
	/*
	if(realInputFile.empty()){
        //cerr << "Warning: the file containing real mutation data was not specified." << endl;
		//usage(argv[0]);
    }
	*/
	
	if(outputPrefix.empty()){
		cerr << "Use --output to specify the prefix of result file names." << endl;
		usage(argv[0]);
		exit(1);
	}
	
	/*
	if(maxc < 1) {
		cerr << "Error: the value of parameter \"maxc\" should be a positive integer." << endl;
		usage(argv[0]);
		exit(1);
	}
	*/
	
	if(threads < 1) {
		cerr << "Error: the value of parameter \"threads\" should be a positive integer." << endl;
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
	
	if(max_alpha > 1) {
		cerr << "Error: the value of parameter \"max_alpha\" should be in [0, 1]." << endl;
		usage(argv[0]);
		exit(1);
	}
	
	if(max_beta > 1) {
		cerr << "Error: the value of parameter \"max_beta\" should be in [0, 1]." << endl;
		usage(argv[0]);
		exit(1);
	}
	
	config.setStringPara("input", inputFile);
	//config.setStringPara("rinput", realInputFile);
	config.setStringPara("output", outputPrefix);
	config.setStringPara("clabel", clabelFile);
	config.setStringPara("mlabel", mlabelFile);
	config.setIntPara("maxc", maxc);
	config.setIntPara("threads", threads);
	config.setRealPara("alpha", alpha);
	config.setRealPara("beta", beta);
	config.setRealPara("max_alpha", max_alpha);
	config.setRealPara("max_beta", max_beta);
	
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
	threadpool = new ThreadPool(threads);
	threadpool->pool_init();
	/*
	string binaryPath = argv[0];
	char abs_path_buff[PATH_MAX];
	realpath(binaryPath.c_str(), abs_path_buff);
	string binPath = abs_path_buff;
	int indx = binPath.rfind("bin");
	binPath = binPath.substr(0, indx-1);
	config.setStringPara("binPath", binPath);
	*/
}

void InputParser::usage(const char* app) {
	cerr << "Usage: " << app << " [options]" << endl
		<< endl
		<< "Options:" << endl
		<< "    -h, --help                      give this information" << endl
		<< "    -v, --version                   print software version" << endl
		<< "    -i, --input <string>            input file containing mutation data" << endl
		//<< "    -r, --rinput <string>           input file containing real mutation data" << endl
		<< "    -o, --output <string>           prefix of output file names" << endl
		<< "    -c, --clabel <string>           file defining labels of the cells" << endl
		<< "    -m, --mlabel <string>           file defining labels of the mutations" << endl
		<< "    -K, --maxc <int>                maximum number of clones to consider" << endl
		<< "    -t, --threads <int>             number of threads to use [default:1]" << endl
		<< "    -a, --alpha <double>            false positive rate [default:inferred from data]" << endl
		<< "    -b, --beta <double>             false negative rate [default:inferred from data]" << endl
		<< "    -A, --max_alpha <double>        maximum false positive rate [default:0.05]" << endl
		<< "    -B, --max_beta <double>         maximum false negative rate [default:0.5]" << endl
		<< endl
		<< "Example:" << endl
		<< app << " -i ./testdata/example.txt -K 10 -t 5 -o ./testdata/example" << endl
		<< endl
		<< "Author: Zhenhua Yu <qasim0208@163.com>\n" << endl;
}
