// ***************************************************************************
// InputParser.h (c) 2020 Zhenhua Yu <qasim0208@163.com>
// Health Informatics Lab, Ningxia University
// All rights reserved.

#ifndef _INPUTPARSER_H
#define _INPUTPARSER_H

class InputParser {
	private:
		void usage(const char* app);
	public:
		void parseArgs(int argc, char *argv[]);
};

#endif
