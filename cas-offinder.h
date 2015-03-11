#include <vector>
#include <string>
#include <cstring>
#include <fstream>
#include <algorithm>
#include <CL/cl.h>

#ifndef min
#define min(a,b) ( ((a)<(b))?(a):(b) )
#endif

using namespace std;

class Cas_OFFinder {
private:
	vector<cl_command_queue> m_queues;
	vector<cl_context> m_contexts;
	vector<cl_ulong> MAX_ALLOC_MEMORY; // on device, in bytes

	unsigned long long m_chrdatasize;
	vector<string> m_chrnames;
	vector<unsigned long long> m_chrpos;
	string m_chrdata;
	cl_char* m_pattern;

	cl_uint m_threshold;
	cl_uint m_patternlen;
	cl_uint m_entrycount;

	vector<cl_kernel> m_finderkernels;
	vector<cl_kernel> m_comparerkernels;

	vector<cl_mem> m_chrdatabufs;
	vector<cl_mem> m_patternbufs;
	vector<cl_mem> m_patternindexbufs;
	vector<cl_mem> m_flagbufs;
	vector<cl_mem> m_locibufs;
	vector<cl_mem> m_entrycountbufs;

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

	void set_complementary_sequence(cl_char* seq);
	void set_pattern_index(int* pattern_index, const cl_char* pattern);
	void initOpenCL(cl_device_type devtype, cl_platform_id* platforms, cl_uint platform_cnt);

public:
	vector<string> chrnames;
	string chrdata;
	vector<unsigned long long> chrpos;

	Cas_OFFinder(cl_device_type devtype, cl_platform_id* platforms, cl_uint platform_cnt);
	~Cas_OFFinder();

	void setChrData();
	void setPattern(const char* pattern);

	bool loadNextChunk();
	void findPattern();
	void releaseLociinfo();

	void indicate_mismatches(cl_char* seq, cl_char* comp);

	void compareAll(const char *arg_compare, unsigned short threshold, const char* outfilename);

	static void get_platforms(cl_platform_id *platforms, cl_uint *platform_cnt);
	static void print_usage(cl_platform_id *platforms, cl_uint platform_cnt);
	static void readInputFile(const char* inputfile, string &chrdir, string &pattern, vector<string> &compares, vector<int> &thresholds);
};
