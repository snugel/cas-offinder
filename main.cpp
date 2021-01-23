#include "config.h"

#include "cas-offinder.h"
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

void run_cas_offinder(Cas_OFFinder &s, const char* chromfilename, const char* outfilename, int *cnum) {
	string filepath = chromfilename;

	cerr << "Reading " << filepath << "..." << endl;
	
	if (read_fasta(filepath, s.chrnames, s.chrdata, s.chrpos)) {
		if (read_twobit(filepath, s.chrnames, s.chrdata, s.chrpos)) {
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
		s.compareAll(outfilename);
		s.releaseLociinfo();
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
	string filepath, tmpstr, outfilename;
	DIR* dir;
	dirent *ent;

	int maxdevnum = MAX_DEVICE_NUM;

	Cas_OFFinder::init_platforms();

	if (argc < 4) { // Not all option specified
		Cas_OFFinder::print_usage();
		exit(0);
	}

	cl_device_type devtype;
	switch (argv[2][0]) {
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
		error_exit(2, "Unknown option: ", argv[2]);
	}

	string devarg = argv[2]+1;
	Cas_OFFinder s(devtype, devarg);

	cerr << "Loading input file..." << endl;
	s.readInputFile(argv[1]);

	outfilename = argv[3]; remove(argv[3]);

	int cnum = 0;
	if ((dir = opendir(s.chrdir.c_str())) == NULL) {
		if (check_file(s.chrdir.c_str())) {
			run_cas_offinder(s, s.chrdir.c_str(), outfilename.c_str(), &cnum);
		} else {
			error_exit(2, "An error has occured while opening directory: ", s.chrdir.c_str());
			exit(1);
		}
	} else {
		while ((ent = readdir(dir)) != NULL) {
			if (ent->d_type == DT_REG || ent->d_type == DT_LNK) {
				filepath = s.chrdir + "/" + ent->d_name;
				run_cas_offinder(s, filepath.c_str(), outfilename.c_str(), &cnum);
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
