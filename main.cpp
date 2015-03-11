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

//#define DEBUG

using namespace std;

int main(int argc, char *argv[]) {
#ifdef _MSC_VER
	clock_t start, end;
#else
	struct timeval start, end;
#endif
	float seconds;
	string filepath, tmpstr, outfilename;
	DIR* dir;
	dirent *ent;
	unsigned int i, cnt;

	cl_platform_id platforms[MAX_PLATFORM_NUM];
	cl_uint platform_cnt;

	Cas_OFFinder::get_platforms(platforms, &platform_cnt);

	if (argc < 4) // Not all option specified
		Cas_OFFinder::print_usage(platforms, platform_cnt);

#ifdef _MSC_VER
	start = clock();
#else
	gettimeofday(&start, NULL);
#endif

	cl_device_type devtype;
	if (argv[2][0] == 'C') {
		devtype = CL_DEVICE_TYPE_CPU;
	}
	else if (argv[2][0] == 'G') {
		devtype = CL_DEVICE_TYPE_GPU;
	}
	else if (argv[2][0] == 'A') {
		devtype = CL_DEVICE_TYPE_ACCELERATOR;
	}
	else {
		cout << "Unknown option: " << argv[2] << endl;
		exit(-1);
	}

	Cas_OFFinder s(devtype, platforms, platform_cnt);

	string chrdir, pattern;
	vector <string> compares;
	vector <int> thresholds;

	cout << "Loading input file..." << endl;

	Cas_OFFinder::readInputFile(argv[1], chrdir, pattern, compares, thresholds);

	outfilename = argv[3];
	remove(argv[3]);

	int cnum = 0, pnum = 0;
	if ((dir = opendir(chrdir.c_str())) == NULL) {
		cout << "No such directory: " << chrdir << endl;
		exit(1);
	}
	else {
		while ((ent = readdir(dir)) != NULL) {
			if (ent->d_type == DT_REG) {
				filepath = ent->d_name;
				filepath = chrdir + "/" + filepath;
				cout << "Reading " << filepath << "..." << endl;
				cnum = 0;
				cnt = 0;
				if (read_fasta(filepath, s.chrnames, s.chrdata, s.chrpos)) {
				    if (read_twobit(filepath, &s.chrnames, &s.chrdata, &s.chrpos)) {
						cout << "Skipping non-acceptable file " << filepath << "..." << endl;
						continue;
					}
				}
				cout << "Sending data to devices..." << endl;
				s.setChrData();
				cout << "Setting pattern to devices..." << endl;
				s.setPattern(pattern.c_str());
				cout << "Chunk load started." << endl;
				while (s.loadNextChunk()) {
					// Find patterns in the chunk
					cout << "Finding pattern in chunk #" << ++cnum << "..." << endl;
					s.findPattern();
					cout << "Comparing patterns in chunk #" << cnum << "..." << endl;
					for (i = 0; i < thresholds.size(); i++)
						s.compareAll(compares[i].c_str(), thresholds[i], (outfilename).c_str());
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
