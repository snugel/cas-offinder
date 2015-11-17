#include "config.h"

#include "cas-offinder.h"
#include "read_fasta.h"
#include "read_twobit.h"

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#ifdef _MSC_VER
#include <ctime>
#else
#include <sys/time.h>
#include <unistd.h>
#endif
#include <cstdlib>
#include <cstdarg>

//#define DEBUG

using namespace std;

void error_exit(int nargs, ...) {
	va_list errmsgs;
	va_start(errmsgs, nargs);
	for (int i=1; i<=nargs; i++)
		cout << (char*)(va_arg(errmsgs, char*));
	cout << endl;
	va_end(errmsgs);
	exit(1);
}

int main(int argc, char *argv[]) {
#ifdef _MSC_VER
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
	unsigned int cnt;

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

	Cas_OFFinder s(devtype);

	cout << "Loading input file..." << endl;
	s.readInputFile(argv[1]);

	outfilename = argv[3]; remove(argv[3]);

	int cnum = 0, pnum = 0;
	if ((dir = opendir(s.chrdir.c_str())) == NULL) {
		error_exit(2, "No such directory: ", s.chrdir.c_str());
		exit(1);
	} else {
		while ((ent = readdir(dir)) != NULL) {
			if (ent->d_type == DT_REG) {
				filepath = ent->d_name;
				filepath = s.chrdir + "/" + filepath;
				cout << "Reading " << filepath << "..." << endl;
				cnum = 0;
				cnt = 0;
				if (read_fasta(filepath, s.chrnames, s.chrdata, s.chrpos)) {
				    if (read_twobit(filepath, s.chrnames, s.chrdata, s.chrpos)) {
						cout << "Skipping non-acceptable file " << filepath << "..." << endl;
						continue;
					}
				}
				cout << "Sending data to devices..." << endl;
				s.setChrData();
				cout << "Chunk load started." << endl;
				while (s.loadNextChunk()) {
					// Find patterns in the chunk
					cout << "Finding pattern in chunk #" << ++cnum << "..." << endl;
					s.findPattern();
					cout << "Comparing patterns in chunk #" << cnum << "..." << endl;
					s.compareAll(outfilename.c_str());
					s.releaseLociinfo();
				}			
			}
		}
	}
#ifdef _MSC_VER
	end = clock();
	seconds = (float)(end - start) / CLOCKS_PER_SEC;
#else
	gettimeofday(&end, NULL);
	seconds = (float)((end.tv_sec + (end.tv_usec / 1000000.0)) - (start.tv_sec + (start.tv_usec / 1000000.0)));
#endif
	cout << seconds << " seconds elapsed." << endl;
	return 0;
}
