#include "config.h"

#include <iostream>
#include <fstream>
#include <streambuf>
#include <string>
#include <vector>
#include <cassert>
#include <ctime>
#include <climits>
#include <algorithm>
#include "kseq.h"

#define myread(a,b,c) fread(b, 1, c, a)
KSEQ_INIT(FILE*, myread)

#ifndef min
#define min(a,b) ( ((a)<(b))?(a):(b) )
#endif

using namespace std;

class Genos {

private:
	vector<cl::CommandQueue> m_queues;
	vector<cl::Context> m_contexts;
	vector<cl_ulong> MAX_ALLOC_MEMORY; // on device, in bytes

	size_t m_chrdatasize;
	cl_char* m_chrdata;
	cl_char* m_chrname;
	cl_char* m_pattern;

	cl_uint m_threshold;
	cl_uint m_patternlen;
	cl_uint m_entrycount;

	vector<cl::Kernel> m_finderkernels;
	vector<cl::Kernel> m_comparerkernels;

	vector<cl::Buffer> m_chrdatabufs;
	vector<cl::Buffer> m_patternbufs;
	vector<cl::Buffer> m_patternindexbufs;
	vector<cl::Buffer> m_flagbufs;
	vector<cl::Buffer> m_locibufs;
	vector<cl::Buffer> m_entrycountbufs;

	vector<cl::Buffer> m_mmlocibufs;
	vector<cl::Buffer> m_mmcountbufs;
	vector<cl::Buffer> m_directionbufs;


	vector<cl::Event> m_finderevts;
	vector<cl::Event> m_comparerevts;

	vector <cl_uint> m_locicnts;
	vector <cl_uint *> m_locis;
	vector <cl_ushort *> m_mmcounts;
	vector <cl_char *> m_flags;
	vector <cl_char *> m_directions;
	vector <cl_uint *> m_mmlocis;

	vector<unsigned int> m_dicesizes;
	unsigned int m_totalanalyzedsize;
	unsigned int m_lasttotalanalyzedsize;
	vector<unsigned int> m_worksizes;
	unsigned int m_devnum;
	unsigned int m_activedevnum;
	unsigned int m_lastloci;

	unsigned int m_linenum;
	unsigned int m_filenum;

	void set_complementary_sequence(cl_char* seq)
	{
		size_t i, l = 0;
		cl_char tmp;

		for (i = 0; seq[i] != 0; i++) {
			if (seq[i] == 'A') seq[i] = 'T';
			else if (seq[i] == 'T') seq[i] = 'A';
			else if (seq[i] == 'G') seq[i] = 'C';
			else if (seq[i] == 'C') seq[i] = 'G';
			else if (seq[i] == 'R') seq[i] = 'Y';
			else if (seq[i] == 'Y') seq[i] = 'R';
			else if (seq[i] == 'M') seq[i] = 'K';
			else if (seq[i] == 'K') seq[i] = 'M';
			else if (seq[i] == 'H') seq[i] = 'D';
			else if (seq[i] == 'D') seq[i] = 'H';
			else if (seq[i] == 'B') seq[i] = 'V';
			else if (seq[i] == 'V') seq[i] = 'B';
			l++;
		}
		for (i = 0; i < l / 2; i++) {
			tmp = seq[i];
			seq[i] = seq[l - i - 1];
			seq[l - i - 1] = tmp;
		}
	}

	void set_pattern_index(int* pattern_index, const cl_char* pattern) {
		int i, n = 0;
		for (i = 0; pattern[i] != 0; i++) {
			if (pattern[i] != 'N') {
				pattern_index[n] = i;
				n++;
			}
		}
		if (i != n)
			pattern_index[n] = -1;
	}

	void initOpenCL(cl_device_type devtype) {
		int i, j;

		vector<cl::Platform> platforms;
		vector<cl::Device> devices;

		cl::Platform::get(&platforms);
		int number_of_platforms = platforms.size();

		vector<cl::Device> devices_per_platform;
		for (i = 0; i < number_of_platforms; i++) {
			devices_per_platform.clear();
			platforms[i].getDevices(devtype, &devices_per_platform);
			for (j = 0; j < devices_per_platform.size(); j++) {
				devices.push_back(devices_per_platform[j]);
			}
		}
		m_devnum = devices.size();
		if (m_devnum == 0) {
			cout << "No OpenCL devices found." << endl;
			exit(1);
		}
		const char* program_src =
			"#pragma OPENCL EXTENSION cl_khr_global_int32_base_atomics : enable\n"
			"#pragma OPENCL EXTENSION cl_khr_local_int32_base_atomics : enable\n"
			"#pragma OPENCL EXTENSION cl_khr_global_int32_extended_atomics : enable\n"
			"#pragma OPENCL EXTENSION cl_khr_local_int32_extended_atomics : enable\n"
			"__kernel void finder(__global char* chr,"
			"                     __global char* pat, __global int* pat_index, unsigned int patternlen,"
			"                     __global char* flag, __global unsigned int* entrycount, __global unsigned int* loci)"
			"{"
			"	unsigned int i = get_global_id(0);"
			"	unsigned int j;"
			"	unsigned int old;"
			"	int k;"
			"	char localflag = 0;"
			"	for (j=0; j<patternlen; j++) {"
			"		k = pat_index[j];"
			"		if (k == -1)"
			"			break;"
			"		if ( (pat[k] == 'R' && (chr[i+k] == 'C' || chr[i+k] == 'T')) ||"
			"		     (pat[k] == 'Y' && (chr[i+k] == 'A' || chr[i+k] == 'G')) ||"
			"		     (pat[k] == 'K' && (chr[i+k] == 'A' || chr[i+k] == 'C')) ||"
			"		     (pat[k] == 'M' && (chr[i+k] == 'G' || chr[i+k] == 'T')) ||"
			"		     (pat[k] == 'W' && (chr[i+k] == 'C' || chr[i+k] == 'G')) ||"
			"		     (pat[k] == 'S' && (chr[i+k] == 'A' || chr[i+k] == 'T')) ||"
			"		     (pat[k] == 'H' && (chr[i+k] == 'G')) ||"
			"		     (pat[k] == 'B' && (chr[i+k] == 'A')) ||"
			"		     (pat[k] == 'V' && (chr[i+k] == 'T')) ||"
			"		     (pat[k] == 'D' && (chr[i+k] == 'C')) ||"
			"		     (pat[k] == 'A' && (chr[i+k] != 'A')) ||"
			"		     (pat[k] == 'G' && (chr[i+k] != 'G')) ||"
			"		     (pat[k] == 'C' && (chr[i+k] != 'C')) ||"
			"		     (pat[k] == 'T' && (chr[i+k] != 'T')))"
			"			localflag |= 2;"
			"		k = pat_index[patternlen + j];"
			"		if ( (pat[k + patternlen] == 'R' && (chr[i+k] == 'C' || chr[i+k] == 'T')) ||"
			"		     (pat[k + patternlen] == 'Y' && (chr[i+k] == 'A' || chr[i+k] == 'G')) ||"
			"		     (pat[k + patternlen] == 'K' && (chr[i+k] == 'A' || chr[i+k] == 'C')) ||"
			"		     (pat[k + patternlen] == 'M' && (chr[i+k] == 'G' || chr[i+k] == 'T')) ||"
			"		     (pat[k + patternlen] == 'W' && (chr[i+k] == 'C' || chr[i+k] == 'G')) ||"
			"		     (pat[k + patternlen] == 'S' && (chr[i+k] == 'A' || chr[i+k] == 'T')) ||"
			"		     (pat[k + patternlen] == 'H' && (chr[i+k] == 'G')) ||"
			"		     (pat[k + patternlen] == 'B' && (chr[i+k] == 'A')) ||"
			"		     (pat[k + patternlen] == 'V' && (chr[i+k] == 'T')) ||"
			"		     (pat[k + patternlen] == 'D' && (chr[i+k] == 'C')) ||"
			"		     (pat[k + patternlen] == 'A' && (chr[i+k] != 'A')) ||"
			"		     (pat[k + patternlen] == 'G' && (chr[i+k] != 'G')) ||"
			"		     (pat[k + patternlen] == 'C' && (chr[i+k] != 'C')) ||"
			"		     (pat[k + patternlen] == 'T' && (chr[i+k] != 'T')))"
			"			localflag |= 1;"
			"		if (localflag == 3)"
			"			break;"
			"	}"
			"	if (localflag != 3) {"
			"		old = atomic_inc(entrycount);"
			"		loci[old] = i;"
			"		flag[old] = localflag;"
			"	}"
			"}"
			"__kernel void comparer(__global char* chr, __global unsigned int* loci, __global unsigned int* mm_loci,"
			"                       __global char* comp, __global int* comp_index, unsigned int patternlen, unsigned short threshold,"
			"                       __global char* flag, __global unsigned short* mm_count, __global char* direction, __global unsigned int* entrycount)"
			"{"
			"	unsigned int i = get_global_id(0);"
			"	unsigned int j, loopcnt, lmm_count, old;"
			"	int k;"
			"	if (flag[i] == 0 || flag[i] == 1) {"
			"		lmm_count = 0;"
			"		for (j=0; j<patternlen; j++) {"
			"			k = comp_index[j]; "
			"			if (k == -1) break;"
			"			if ( (comp[k] == 'R' && (chr[loci[i]+k] == 'C' || chr[loci[i]+k] == 'T')) ||"
			"			     (comp[k] == 'Y' && (chr[loci[i]+k] == 'A' || chr[loci[i]+k] == 'G')) ||"
			"			     (comp[k] == 'K' && (chr[loci[i]+k] == 'A' || chr[loci[i]+k] == 'C')) ||"
			"			     (comp[k] == 'M' && (chr[loci[i]+k] == 'G' || chr[loci[i]+k] == 'T')) ||"
			"			     (comp[k] == 'W' && (chr[loci[i]+k] == 'C' || chr[loci[i]+k] == 'G')) ||"
			"			     (comp[k] == 'S' && (chr[loci[i]+k] == 'A' || chr[loci[i]+k] == 'T')) ||"
			"			     (comp[k] == 'H' && (chr[loci[i]+k] == 'G')) ||"
			"			     (comp[k] == 'B' && (chr[loci[i]+k] == 'A')) ||"
			"			     (comp[k] == 'V' && (chr[loci[i]+k] == 'T')) ||"
			"			     (comp[k] == 'D' && (chr[loci[i]+k] == 'C')) ||"
			"				 (comp[k] == 'A' && (chr[loci[i]+k] != 'A')) ||"
			"			     (comp[k] == 'G' && (chr[loci[i]+k] != 'G')) ||"
			"			     (comp[k] == 'C' && (chr[loci[i]+k] != 'C')) ||"
			"			     (comp[k] == 'T' && (chr[loci[i]+k] != 'T'))) {"
			"				lmm_count++;"
			"				if (lmm_count > threshold) break;"
			"			}"
			"		}"
			"		if (lmm_count <= threshold) {"
			"			old = atomic_inc(entrycount);"
			"			mm_count[old] = lmm_count;"
			"			direction[old] = '+';"
			"			mm_loci[old] = loci[i];"
			"		}"
			"	}"
			"	if (flag[i] == 0 || flag[i] == 2) {"
			"		lmm_count = 0;"
			"		for (j=0; j<patternlen; j++) {"
			"			k = comp_index[patternlen + j]; "
			"			if (k == -1) break;"
			"			if ( (comp[k+patternlen] == 'R' && (chr[loci[i]+k] == 'C' || chr[loci[i]+k] == 'T')) ||"
			"			     (comp[k+patternlen] == 'Y' && (chr[loci[i]+k] == 'A' || chr[loci[i]+k] == 'G')) ||"
			"			     (comp[k+patternlen] == 'K' && (chr[loci[i]+k] == 'A' || chr[loci[i]+k] == 'C')) ||"
			"			     (comp[k+patternlen] == 'M' && (chr[loci[i]+k] == 'G' || chr[loci[i]+k] == 'T')) ||"
			"			     (comp[k+patternlen] == 'W' && (chr[loci[i]+k] == 'C' || chr[loci[i]+k] == 'G')) ||"
			"			     (comp[k+patternlen] == 'S' && (chr[loci[i]+k] == 'A' || chr[loci[i]+k] == 'T')) ||"
			"			     (comp[k+patternlen] == 'H' && (chr[loci[i]+k] == 'G')) ||"
			"			     (comp[k+patternlen] == 'B' && (chr[loci[i]+k] == 'A')) ||"
			"			     (comp[k+patternlen] == 'V' && (chr[loci[i]+k] == 'T')) ||"
			"			     (comp[k+patternlen] == 'D' && (chr[loci[i]+k] == 'C')) ||"
			"			     (comp[k+patternlen] == 'A' && (chr[loci[i]+k] != 'A')) ||"
			"			     (comp[k+patternlen] == 'G' && (chr[loci[i]+k] != 'G')) ||"
			"			     (comp[k+patternlen] == 'C' && (chr[loci[i]+k] != 'C')) ||"
			"				 (comp[k+patternlen] == 'T' && (chr[loci[i]+k] != 'T'))) {"
			"				lmm_count++;"
			"				if (lmm_count > threshold) break;"
			"			}"
			"		}"
			"		if (lmm_count <= threshold) {"
			"			old = atomic_inc(entrycount);"
			"			mm_count[old] = lmm_count;"
			"			direction[old] = '-';"
			"			mm_loci[old] = loci[i];"
			"		}"
			"	}"
			"}";

		cl::Program::Sources source(1, std::make_pair(program_src, strlen(program_src)));;
		vector<cl::Device> devvec;
		for (i = 0; i < devices.size(); i++) {
			// Create completely separate contexts per device to avoid unknown errors
			devvec.clear(); devvec.push_back(devices[i]);
			m_contexts.push_back(cl::Context(devvec));
			cl::Program program(m_contexts[i], source);
			program.build(devvec);
			m_finderkernels.push_back(cl::Kernel(program, "finder"));
			m_comparerkernels.push_back(cl::Kernel(program, "comparer"));
			m_queues.push_back(cl::CommandQueue(m_contexts[i], devices[i]));
			MAX_ALLOC_MEMORY.push_back(0);
			devices[i].getInfo(CL_DEVICE_MAX_MEM_ALLOC_SIZE, &MAX_ALLOC_MEMORY[i]);
		}
		cout << "Total " << m_devnum << " device(s) found." << endl;
	}

public:
	Genos(cl_device_type devtype) {
		initOpenCL(devtype);
	}

	void setChrData(char *chrname, char *chrdata, size_t chrdatasize) {
		m_chrname = (cl_char *)chrname;
		m_chrdata = (cl_char *)chrdata;
		m_chrdatasize = chrdatasize;
		m_totalanalyzedsize = 0;
		m_lasttotalanalyzedsize = 0;
		m_lastloci = 0;
	}

	void setPattern(const char* pattern) {
		int dev_index;
		m_pattern = (cl_char*)pattern;
		m_patternlen = strlen(pattern);

		m_dicesizes.clear();
		m_chrdatabufs.clear();
		m_patternbufs.clear();
		m_patternindexbufs.clear();
		m_flagbufs.clear();
		m_locibufs.clear();
		m_entrycountbufs.clear();

		for (dev_index = 0; dev_index < m_devnum; dev_index++) {
			m_dicesizes.push_back(
				min(
				(MAX_ALLOC_MEMORY[dev_index] - sizeof(cl_char)* (3 * m_patternlen - 1) - sizeof(cl_uint)* (2 * m_patternlen + 3) - sizeof(cl_ushort)) / (4 * sizeof(cl_char)+3 * sizeof(cl_uint)+2 * sizeof(cl_ushort)),
				((m_chrdatasize / m_devnum) + ((m_chrdatasize%m_devnum == 0) ? 0 : 1))
				)
				); // No more than maximum allocation per device
			// cout << "Dicesize: " << m_dicesizes[dev_index] << endl;
			m_chrdatabufs.push_back(cl::Buffer(m_contexts[dev_index], CL_MEM_READ_ONLY, sizeof(cl_char)* (m_dicesizes[dev_index] + m_patternlen - 1)));
			m_patternbufs.push_back(cl::Buffer(m_contexts[dev_index], CL_MEM_READ_ONLY, sizeof(cl_char)* m_patternlen * 2));
			m_patternindexbufs.push_back(cl::Buffer(m_contexts[dev_index], CL_MEM_READ_ONLY, sizeof(cl_int)* m_patternlen * 2));
			m_flagbufs.push_back(cl::Buffer(m_contexts[dev_index], CL_MEM_WRITE_ONLY, sizeof(cl_char)* m_dicesizes[dev_index]));
			m_entrycountbufs.push_back(cl::Buffer(m_contexts[dev_index], CL_MEM_READ_WRITE, sizeof(cl_uint)));
			m_locibufs.push_back(cl::Buffer(m_contexts[dev_index], CL_MEM_WRITE_ONLY, sizeof(cl_uint)* m_dicesizes[dev_index]));
		}
	}

	bool loadNextChunk() {
		if (m_totalanalyzedsize == m_chrdatasize)
			return false;

		int dev_index;
		int tailsize;

		cl_char *c_pattern = (cl_char *)malloc(sizeof(cl_char)* (m_patternlen + 1)); c_pattern[m_patternlen] = 0; memcpy(c_pattern, m_pattern, m_patternlen); set_complementary_sequence(c_pattern);
		int *c_pattern_index = (cl_int *)malloc(sizeof(cl_int)* m_patternlen); set_pattern_index(c_pattern_index, c_pattern);
		int *pattern_index = (cl_int *)malloc(sizeof(cl_int)* m_patternlen); set_pattern_index(pattern_index, m_pattern);

		m_activedevnum = 0;
		m_worksizes.clear();
		m_lasttotalanalyzedsize = m_totalanalyzedsize;

		for (dev_index = 0; dev_index < m_devnum; dev_index++) {
			tailsize = m_chrdatasize - m_totalanalyzedsize;
			m_activedevnum++;
			if (tailsize <= m_dicesizes[dev_index]) {
				m_queues[dev_index].enqueueWriteBuffer(m_chrdatabufs[dev_index], CL_FALSE, 0, sizeof(cl_char)* (tailsize + m_patternlen - 1), m_chrdata + m_totalanalyzedsize, NULL);
				m_totalanalyzedsize += tailsize;
				m_worksizes.push_back(tailsize);
				// cout << "Worksize: " <<  m_worksizes[dev_index] << ", Tailsize: " << tailsize << endl;
				break;
			}
			else {
				m_queues[dev_index].enqueueWriteBuffer(m_chrdatabufs[dev_index], CL_FALSE, 0, sizeof(cl_char)* (m_dicesizes[dev_index] + m_patternlen - 1), m_chrdata + m_totalanalyzedsize, NULL);
				m_totalanalyzedsize += m_dicesizes[dev_index];
				m_worksizes.push_back(m_dicesizes[dev_index]);
				// cout << "Worksize: " <<  m_worksizes[dev_index] << ", Tailsize: " << tailsize << endl;
			}
		}
		cout << m_activedevnum << " devices selected to analyze..." << endl;

		for (dev_index = 0; dev_index < m_activedevnum; dev_index++) {
			cl_uint zero = 0;
			m_queues[dev_index].enqueueWriteBuffer(m_patternbufs[dev_index], CL_FALSE, 0, sizeof(cl_char)* m_patternlen, m_pattern, NULL);
			m_queues[dev_index].enqueueWriteBuffer(m_patternbufs[dev_index], CL_FALSE, sizeof(cl_char)* m_patternlen, sizeof(cl_char)* m_patternlen, c_pattern, NULL);
			m_queues[dev_index].enqueueWriteBuffer(m_patternindexbufs[dev_index], CL_FALSE, 0, sizeof(cl_int)* m_patternlen, pattern_index, NULL);
			m_queues[dev_index].enqueueWriteBuffer(m_patternindexbufs[dev_index], CL_FALSE, sizeof(cl_int)* m_patternlen, sizeof(cl_int)* m_patternlen, c_pattern_index, NULL);
			m_queues[dev_index].enqueueWriteBuffer(m_entrycountbufs[dev_index], CL_FALSE, 0, sizeof(cl_uint), &zero, NULL);
			m_queues[dev_index].finish();

			m_finderkernels[dev_index].setArg(0, m_chrdatabufs[dev_index]);
			m_finderkernels[dev_index].setArg(1, m_patternbufs[dev_index]);
			m_finderkernels[dev_index].setArg(2, m_patternindexbufs[dev_index]);
			m_finderkernels[dev_index].setArg(3, sizeof(cl_uint), &m_patternlen);
			m_finderkernels[dev_index].setArg(4, m_flagbufs[dev_index]);
			m_finderkernels[dev_index].setArg(5, m_entrycountbufs[dev_index]);
			m_finderkernels[dev_index].setArg(6, m_locibufs[dev_index]);
		}

		free((void*)pattern_index);
		free((void*)c_pattern);
		free((void*)c_pattern_index);

		return true;
	}

	void findPattern() {
		int dev_index, i;

		m_finderevts.clear();
		for (dev_index = 0; dev_index < m_activedevnum; dev_index++) {
			m_finderevts.push_back(cl::Event());
			m_queues[dev_index].enqueueNDRangeKernel(m_finderkernels[dev_index], cl::NullRange, cl::NDRange(m_worksizes[dev_index]), cl::NullRange, NULL, &m_finderevts[dev_index]);
		}

		for (dev_index = 0; dev_index < m_activedevnum; dev_index++) {
			m_finderevts[dev_index].wait();
			m_locicnts.push_back(0);

			m_queues[dev_index].enqueueReadBuffer(m_entrycountbufs[dev_index], CL_TRUE, 0, sizeof(cl_uint), &m_locicnts[dev_index], NULL);
			m_flags.push_back((cl_char *)malloc(sizeof(cl_char)* m_locicnts[dev_index]));
			m_queues[dev_index].enqueueReadBuffer(m_flagbufs[dev_index], CL_TRUE, 0, sizeof(cl_char)*m_locicnts[dev_index], m_flags[dev_index], NULL);

			m_mmcounts.push_back((cl_ushort *)malloc(sizeof(cl_ushort)* m_locicnts[dev_index] * 2)); // Maximum numbers of mismatch counts
			m_directions.push_back((cl_char *)malloc(sizeof(cl_char)* m_locicnts[dev_index] * 2));
			m_mmlocis.push_back((cl_uint *)malloc(sizeof(cl_uint)* m_locicnts[dev_index] * 2));

			m_mmlocibufs.push_back(cl::Buffer(m_contexts[dev_index], CL_MEM_WRITE_ONLY, sizeof(cl_uint)* m_locicnts[dev_index] * 2));
			m_mmcountbufs.push_back(cl::Buffer(m_contexts[dev_index], CL_MEM_WRITE_ONLY, sizeof(cl_ushort)* m_locicnts[dev_index] * 2));
			m_directionbufs.push_back(cl::Buffer(m_contexts[dev_index], CL_MEM_WRITE_ONLY, sizeof(cl_char)* m_locicnts[dev_index] * 2));

			m_comparerkernels[dev_index].setArg(1, m_locibufs[dev_index]);
			m_comparerkernels[dev_index].setArg(2, m_mmlocibufs[dev_index]);
			m_comparerkernels[dev_index].setArg(7, m_flagbufs[dev_index]);
		}
	}

	void releaseLociinfo() {
		int dev_index;

		for (dev_index = 0; dev_index < m_activedevnum; dev_index++) {
			free((void *)m_mmcounts[dev_index]);
			free((void *)m_flags[dev_index]);
			free((void *)m_directions[dev_index]);
			free((void *)m_mmlocis[dev_index]);
		}
		m_directions.clear();
		m_mmlocis.clear();
		m_mmcounts.clear();
		m_locicnts.clear();
		m_mmlocibufs.clear();
		m_flags.clear();
		m_mmcountbufs.clear();
		m_directionbufs.clear();
	}

	void indicate_mismatches(cl_char* seq, cl_char* comp) {
		int k;
		for (k = 0; k < m_patternlen; k++)
		if ((comp[k] == 'R' && (seq[k] == 'C' || seq[k] == 'T')) ||
			(comp[k] == 'Y' && (seq[k] == 'A' || seq[k] == 'G')) ||
			(comp[k] == 'K' && (seq[k] == 'A' || seq[k] == 'C')) ||
			(comp[k] == 'M' && (seq[k] == 'G' || seq[k] == 'T')) ||
			(comp[k] == 'W' && (seq[k] == 'C' || seq[k] == 'G')) ||
			(comp[k] == 'S' && (seq[k] == 'A' || seq[k] == 'T')) ||
			(comp[k] == 'H' && (seq[k] == 'G')) ||
			(comp[k] == 'B' && (seq[k] == 'A')) ||
			(comp[k] == 'V' && (seq[k] == 'T')) ||
			(comp[k] == 'D' && (seq[k] == 'C')) ||
			(comp[k] == 'A' && (seq[k] != 'A')) ||
			(comp[k] == 'G' && (seq[k] != 'G')) ||
			(comp[k] == 'C' && (seq[k] != 'C')) ||
			(comp[k] == 'T' && (seq[k] != 'T')))
			seq[k] += 32;
	}

	void compareAll(const char *arg_compare, unsigned short threshold, const char* outfilename) {
		int i, j, k, dev_index;

		vector <cl::Buffer> comparebufs;
		vector <cl::Buffer> compareindexbufs;

		cl_char *compare = (cl_char *)arg_compare;
		cl_char *c_compare = (cl_char *)malloc(sizeof(cl_char)* (m_patternlen + 1)); c_compare[m_patternlen] = 0; memcpy(c_compare, compare, m_patternlen); set_complementary_sequence(c_compare);
		int *c_compare_index = (cl_int *)malloc(sizeof(cl_int)* m_patternlen); set_pattern_index(c_compare_index, c_compare);
		int *compare_index = (cl_int *)malloc(sizeof(cl_int)* m_patternlen); set_pattern_index(compare_index, compare);

		m_comparerevts.clear();
		for (dev_index = 0; dev_index<m_activedevnum; dev_index++) {
			if (m_locicnts[dev_index] > 0) {
				const cl_char *compare = (const cl_char*)arg_compare;

				comparebufs.push_back(cl::Buffer(m_contexts[dev_index], CL_MEM_READ_ONLY, sizeof(cl_char)* m_patternlen * 2));
				compareindexbufs.push_back(cl::Buffer(m_contexts[dev_index], CL_MEM_READ_ONLY, sizeof(cl_uint)* m_patternlen * 2));

				cl_uint zero = 0;
				m_queues[dev_index].enqueueWriteBuffer(comparebufs[dev_index], CL_FALSE, 0, sizeof(cl_char)* m_patternlen, compare, NULL);
				m_queues[dev_index].enqueueWriteBuffer(comparebufs[dev_index], CL_FALSE, sizeof(cl_char)* m_patternlen, sizeof(cl_char)* m_patternlen, c_compare, NULL);
				m_queues[dev_index].enqueueWriteBuffer(compareindexbufs[dev_index], CL_FALSE, 0, sizeof(cl_int)* m_patternlen, compare_index, NULL);
				m_queues[dev_index].enqueueWriteBuffer(compareindexbufs[dev_index], CL_FALSE, sizeof(cl_int)* m_patternlen, sizeof(cl_int)* m_patternlen, c_compare_index, NULL);
				m_queues[dev_index].enqueueWriteBuffer(m_entrycountbufs[dev_index], CL_FALSE, 0, sizeof(cl_uint), &zero, NULL);
				m_queues[dev_index].finish();

				cl_ushort cl_threshold = threshold;
				m_comparerkernels[dev_index].setArg(0, m_chrdatabufs[dev_index]);
				//m_comparerkernels[dev_index].setArg(1, m_locibufs[dev_index]);
				//m_comparerkernels[dev_index].setArg(2, m_mmlocibufs[dev_index]);
				m_comparerkernels[dev_index].setArg(3, comparebufs[dev_index]);
				m_comparerkernels[dev_index].setArg(4, compareindexbufs[dev_index]);
				m_comparerkernels[dev_index].setArg(5, sizeof(cl_uint), &m_patternlen);
				m_comparerkernels[dev_index].setArg(6, sizeof(cl_ushort), &cl_threshold);
				//m_comparerkernels[dev_index].setArg(7, m_flagbufs[dev_index]);
				m_comparerkernels[dev_index].setArg(8, m_mmcountbufs[dev_index]);
				m_comparerkernels[dev_index].setArg(9, m_directionbufs[dev_index]);
				m_comparerkernels[dev_index].setArg(10, m_entrycountbufs[dev_index]);

				m_comparerevts.push_back(cl::Event());
				m_queues[dev_index].enqueueNDRangeKernel(m_comparerkernels[dev_index], cl::NullRange, cl::NDRange(m_locicnts[dev_index]), cl::NullRange, NULL, &m_comparerevts[dev_index]);
			}
		}

		unsigned int loci;

		char comp_symbol[2] = { '+', '-' };
		char *strbuf = (char *)malloc(sizeof(char)* (m_patternlen + 1)); strbuf[m_patternlen] = 0;

		ofstream fo(outfilename, ios::out | ios::app);
		unsigned int localanalyzedsize = 0;
		unsigned int cnt = 0;
		for (dev_index = 0; dev_index<m_activedevnum; dev_index++) {
			if (m_locicnts[dev_index] > 0) {
				m_comparerevts[dev_index].wait();
				m_queues[dev_index].enqueueReadBuffer(m_entrycountbufs[dev_index], CL_TRUE, 0, sizeof(cl_uint), &cnt, NULL);
				m_queues[dev_index].enqueueReadBuffer(m_mmcountbufs[dev_index], CL_TRUE, 0, sizeof(cl_ushort)* cnt, m_mmcounts[dev_index], NULL);
				m_queues[dev_index].enqueueReadBuffer(m_directionbufs[dev_index], CL_TRUE, 0, sizeof(cl_char)* cnt, m_directions[dev_index], NULL);
				m_queues[dev_index].enqueueReadBuffer(m_mmlocibufs[dev_index], CL_TRUE, 0, sizeof(cl_uint)* cnt, m_mmlocis[dev_index], NULL);

				for (i = 0; i < cnt; i++) {
					loci = m_mmlocis[dev_index][i] + m_lasttotalanalyzedsize + localanalyzedsize;
					if (m_mmcounts[dev_index][i] <= threshold) {
						strncpy(strbuf, (char *)(m_chrdata + loci), m_patternlen);
						if (m_directions[dev_index][i] == '-') set_complementary_sequence((cl_char *)strbuf);
						indicate_mismatches((cl_char*)strbuf, compare);
						fo << compare << "\t" << m_chrname << "\t" << loci << "\t" << strbuf << "\t" << m_directions[dev_index][i] << "\t" << m_mmcounts[dev_index][i] << endl;
					}
				}
			}
			localanalyzedsize += m_worksizes[dev_index];
		}
		fo.close();

		free((void *)strbuf);
		free((void *)c_compare);
		free((void *)c_compare_index);
		free((void *)compare_index);
	}
};

int StrUpr(char *str)
{
	int loop = 0;
	while (str[loop] != '\0')
	{
		str[loop] = (char)toupper(str[loop]);
		loop++;
	}
	return loop;
}

int StrUpr(string &str)
{
	int loop = 0;
	for (loop = 0; loop < str.length(); loop++)
	{
		str[loop] = (int)toupper((int)str[loop]);
		loop++;
	}
	return loop;
}

int main(int argc, char *argv[]) {
	clock_t start, end;
	float seconds;
	string fastapath, chrname, pattern, tmpstr, chrdir, outfilename;
	DIR* dir;
	dirent *ent;
	int i, j, l, tmpint;

	FILE *f_fasta;

	kseq_t *seq_f;

	vector<cl::Platform> platforms;

	cl::Platform::get(&platforms);
	int number_of_platforms = platforms.size();
	if (number_of_platforms == 0) {
		cout << "No OpenCL platforms found. Check OpenCL installation!" << endl;
		exit(1);
	}

	if (argc < 4) { // Not all option specified
		cout << "Cas-OFFinder v2.1 (2014-03-23)" << endl <<
			endl <<
			"Copyright (c) 2013 Jeongbin Park and Sangsu Bae" << endl <<
			"Website: http://github.com/snugel/cas-offinder" << endl <<
			endl <<
			"Usage: cas-offinder {input_file} {C|G|A} {output_file}" << endl <<
			"(C: using CPUs, G: using GPUs, A: using accelerators)" << endl <<
			endl <<
			"Example input file:" << endl <<
			"/var/chromosomes/human_hg19" << endl <<
			"NNNNNNNNNNNNNNNNNNNNNRG" << endl <<
			"GGCCGACCTGTCGCTGACGCNNN 5" << endl <<
			"CGCCAGCGTCAGCGACAGGTNNN 5" << endl <<
			"ACGGCGCCAGCGTCAGCGACNNN 5" << endl <<
			"GTCGCTGACGCTGGCGCCGTNNN 5" << endl <<
			endl <<
			"Available device list:" << endl;
		vector<cl::Device> devices_per_platform;
		for (i = 0; i < number_of_platforms; i++) {
			devices_per_platform.clear();
			platforms[i].getDevices(CL_DEVICE_TYPE_CPU, &devices_per_platform);
			for (j = 0; j < devices_per_platform.size(); j++) {
				cout << "Type: CPU, '" << devices_per_platform[j].getInfo<CL_DEVICE_NAME>() << "'" << endl;
			}
			devices_per_platform.clear();
			platforms[i].getDevices(CL_DEVICE_TYPE_GPU, &devices_per_platform);
			for (j = 0; j < devices_per_platform.size(); j++) {
				cout << "Type: GPU, '" << devices_per_platform[j].getInfo<CL_DEVICE_NAME>() << "'" << endl;
			}
			platforms[i].getDevices(CL_DEVICE_TYPE_ACCELERATOR, &devices_per_platform);
			for (j = 0; j < devices_per_platform.size(); j++) {
				cout << "Type: ACCELERATOR, '" << devices_per_platform[j].getInfo<CL_DEVICE_NAME>() << "'" << endl;
			}
		}
		exit(0);
	}

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
	Genos s(devtype);

	vector <string> compare;
	vector <int> threshold;

	start = clock();
	cout << "Loading input file..." << endl;

	ifstream fi(argv[1], ios::in);
	if (!fi.good()) {
		cout << "No such file: " << argv[1] << endl;
		exit(0);
	}
	if (!fi.eof())
		fi >> chrdir;
	if (!fi.eof())
		fi >> pattern;
	while (true) {
		if (fi.eof()) break;
		fi >> tmpstr;
		StrUpr(tmpstr);
		compare.push_back(tmpstr);
		if (fi.eof()) break;
		fi >> tmpint;
		threshold.push_back(tmpint);
	}
	fi.close();

	StrUpr(pattern);

	if (argc > 3) {
		outfilename = argv[3];
	}
	else {
		outfilename = pattern;
	}
	remove(outfilename.c_str());

	int cnum = 0, pnum = 0;
	if ((dir = opendir(chrdir.c_str())) == NULL) {
		cout << "No such directory: " << chrdir << endl;
		exit(1);
	}
	else {
		while ((ent = readdir(dir)) != NULL) {
			if (ent->d_type == DT_REG) {
				fastapath = ent->d_name;
				fastapath = chrdir + "/" + fastapath;
				cout << "Opening " << fastapath << "..." << endl;
				cnum = 0;

				f_fasta = fopen(fastapath.c_str(), "r");
				seq_f = kseq_init(f_fasta);
				while ((l = kseq_read(seq_f)) >= 0) {
					cout << "Analyzing " << seq_f->name.s << "..." << endl;
					StrUpr(seq_f->seq.s);
					cout << "Sending data to devices..." << endl;
					s.setChrData(seq_f->name.s, seq_f->seq.s, seq_f->seq.l);
					cout << "Setting pattern to devices..." << endl;
					s.setPattern(pattern.c_str());
					cout << "Chunk load started." << endl;
					while (s.loadNextChunk()) {
						// Find patterns in the chunk
						cout << "Finding pattern in chunk #" << ++cnum << "..." << endl;
						s.findPattern();
						pnum = 0;
						for (i = 0; i < threshold.size(); i++) {
							cout << "Comparing pattern #" << ++pnum << " in chunk #" << cnum << "..." << endl;
							s.compareAll(compare[i].c_str(), threshold[i], (outfilename).c_str());
						}
						s.releaseLociinfo();
					}
				}
				kseq_destroy(seq_f);
				fclose(f_fasta);
			}
		}
	}
	end = clock(); seconds = (float)(end - start) / CLOCKS_PER_SEC; cout << seconds << " seconds elapsed." << endl;
	return 0;
}
