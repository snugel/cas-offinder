#include "config.h"

#include "cas-offinder.h"
#include "oclfunctions.h"
#include "read_fasta.h"
#include "read_twobit.h"

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#ifdef _WIN32
	#include <ctime>
#else
	#include <sys/time.h>
	#include <unistd.h>
	#include <sys/stat.h>
#endif

#include <cstdlib>
#include <cstdarg>

//#define DEBUG

using namespace std;

void error_exit(int nargs, ...) {
	va_list errmsgs;
	va_start(errmsgs, nargs);
	for (int i=1; i<=nargs; i++)
		cerr << (char*)(va_arg(errmsgs, char*));
	cerr << endl;
	va_end(errmsgs);
	exit(1);
}

int check_file(const char* filename) {
	struct stat file_stat;
	return (stat(filename, &file_stat) == 0);
}

void run_cas_offinder(Cas_OFFinder &s, const char* chromfilename, const char* outfilename, const char* summaryfilename, int *cnum) {
	string filepath = chromfilename;
	bool issummary = strlen(summaryfilename) > 0;

	cerr << "Reading " << filepath << "..." << endl;
	
	if (read_fasta(filepath, s.m_chrnames, s.m_chrdata, s.m_chrpos)) {
		if (read_twobit(filepath, s.m_chrnames, s.m_chrdata, s.m_chrpos)) {
			cerr << "Skipping non-acceptable file " << filepath << "..." << endl;
			return;
		}
	}
	cerr << "Sending data to devices..." << endl;
	s.setChrData();
	cerr << "Chunk load started." << endl;
	while (s.loadNextChunk()) {
		// Find patterns in the chunk
		cerr << "Finding pattern in chunk #" << ++(*cnum) << "..." << endl;
		s.findPattern();
		cerr << "Comparing patterns in chunk #" << *cnum << "..." << endl;
		s.compareAll(outfilename, issummary);
		s.releaseLociinfo();
	}
	if (issummary) {
		cerr << "Writing summary file..." << endl;
		s.writeSummaryTable(summaryfilename);
	}
}

int main(int argc, char *argv[]) {
#ifdef _WIN32
	clock_t start, end;
	start = clock();
#else
	struct timeval start, end;
	gettimeofday(&start, NULL);
#endif
	float seconds;
	string filepath, tmpstr, outfilename, summaryfilename = "";
	DIR* dir;
	dirent *ent;
	vector<string> positionalargs;

	Cas_OFFinder::initOpenCLPlatforms();

	for (size_t i=1; i<argc; i++) {
		tmpstr = string(argv[i]);
		if (tmpstr.size() > 2 && tmpstr[0] == '-' && tmpstr[1] == '-') {
			tmpstr = tmpstr.substr(2);
			if (tmpstr.compare("summary") == 0) {
				summaryfilename = argv[++i];
			} else {
				Cas_OFFinder::print_usage();
				exit(0);
			}
		} else {
			positionalargs.push_back(string(argv[i]));
		}
	}

	if (positionalargs.size() < 3) { // Not all option specified
		Cas_OFFinder::print_usage();
		exit(0);
	}

	cl_device_type devtype;
	switch (positionalargs[1][0]) {
	case 'C':
		devtype = CL_DEVICE_TYPE_CPU;
		break;
	case 'G':
		devtype = CL_DEVICE_TYPE_GPU;
		break;
	case 'A':
		devtype = CL_DEVICE_TYPE_ACCELERATOR;
		break;
	default:
		error_exit(2, "Unknown option: ", positionalargs[1]);
	}

	string devarg = positionalargs[1].substr(1);
	Cas_OFFinder s(devtype, devarg);

	cerr << "Loading input file..." << endl;
	s.readInputFile(positionalargs[0].c_str());

	outfilename = positionalargs[2];
	s.writeHeaders(outfilename.c_str());

	int cnum = 0;
	if ((dir = opendir(s.m_chrdir.c_str())) == NULL) {
		if (check_file(s.m_chrdir.c_str())) {
			run_cas_offinder(s, s.m_chrdir.c_str(), outfilename.c_str(), summaryfilename.c_str(), &cnum);
		} else {
			error_exit(2, "An error has occured while opening directory: ", s.m_chrdir.c_str());
			exit(1);
		}
	} else {
		while ((ent = readdir(dir)) != NULL) {
			if (ent->d_type == DT_REG || ent->d_type == DT_LNK) {
				filepath = s.m_chrdir + "/" + ent->d_name;
				run_cas_offinder(s, filepath.c_str(), outfilename.c_str(), summaryfilename.c_str(), &cnum);
			}
		}
	}
#ifdef _WIN32
	end = clock();
	seconds = (float)(end - start) / CLOCKS_PER_SEC;
#else
	gettimeofday(&end, NULL);
	seconds = (float)((end.tv_sec + (end.tv_usec / 1000000.0)) - (start.tv_sec + (start.tv_usec / 1000000.0)));
#endif
	cerr << seconds << " seconds elapsed." << endl;
	return 0;
}
