#include "oclkernels.h"

#include <map>
#include <vector>
#include <string>
#include <cstring>
#include <fstream>
#include <algorithm>
#include <utility>

#define CL_USE_DEPRECATED_OPENCL_1_2_APIS
#ifdef __APPLE__
#  include <OpenCL/cl.h>
#else
#  include <CL/cl.h>
#endif

#ifndef MIN
#define MIN(a,b) ( ((a)<(b))?(a):(b) )
#endif

using namespace std;

static cl_platform_id platforms[MAX_PLATFORM_NUM];
static cl_uint platform_cnt;

typedef pair<string, cl_ushort> compareinfo;
typedef pair<string, int> bulgeinfo;
typedef map<string, pair<compareinfo, vector<bulgeinfo> > > m_compare_t;

static inline bool isnumeric(string s) {
	return s.find_first_not_of("0123456789") == std::string::npos;
}

class Cas_OFFinder {
private:
	cl_device_type m_devtype;

	vector<cl_command_queue> m_queues;
	vector<cl_context> m_contexts;
	vector<cl_ulong> MAX_ALLOC_MEMORY; // on device, in bytes
	vector<size_t> MAX_LOCAL_SIZE;

	unsigned long long m_chrdatasize;
	m_compare_t m_compares;
	string m_pattern;
	unsigned int m_dnabulgesize;
	unsigned int m_rnabulgesize;

	cl_uint m_threshold;
	cl_uint m_patternlen;
	cl_uint m_entrycount;

	vector<cl_kernel> m_finderkernels;
	vector<cl_kernel> m_comparerkernels;

	vector<cl_mem> m_chrdatabufs;
	vector<cl_mem> m_patternbufs;
	vector<cl_mem> m_patternflagbufs;
	vector<cl_mem> m_flagbufs;
	vector<cl_mem> m_locibufs;
	vector<cl_mem> m_entrycountbufs;
	vector<cl_mem> m_comparebufs;
	vector<cl_mem> m_compareflagbufs;


	vector<cl_mem> m_mmlocibufs;
	vector<cl_mem> m_mmcountbufs;
	vector<cl_mem> m_directionbufs;

	vector <cl_uint> m_locicnts;
	vector <cl_uint *> m_locis;
	vector <cl_ushort *> m_mmcounts;
	vector <cl_char *> m_flags;
	vector <cl_char *> m_directions;
	vector <cl_uint *> m_mmlocis;

	vector<size_t> m_dicesizes;
	unsigned long long m_totalanalyzedsize;
	unsigned long long m_lasttotalanalyzedsize;
	vector<unsigned long long> m_worksizes;
	size_t m_devnum;
	unsigned int m_activedevnum;
	unsigned long long m_lastloci;

	unsigned long long m_linenum;
	unsigned long long m_filenum;

	void set_complementary_sequence(cl_char* seq, size_t seqlen);
	void set_seq_flags(int* seq_flags, const cl_char* seq, size_t seqlen);
	void initOpenCL(vector<unsigned int> dev_ids);
	void parseInput(istream& input);
	void indicate_mismatches(cl_char* seq, cl_char* comp);

public:
	string m_chrdir;
	string m_chrdata;
	vector<string> m_chrnames;
	vector<unsigned long long> m_chrpos;

	Cas_OFFinder(cl_device_type devtype, string devarg);
	~Cas_OFFinder();

	void writeHeaders(const char* outfilename);
	void setChrData();

	bool loadNextChunk();
	void findPattern();
	void releaseLociinfo();

	void compareAll(const char* outfilename);
	void readInputFile(const char* inputfile);

	static void print_usage();
	static void init_platforms();
};
