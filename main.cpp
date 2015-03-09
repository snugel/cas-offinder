#include "config.h"
#include "main.h"

#include <iostream>
#include <fstream>
#include <streambuf>
#include <string>
#include <cstring>
#include <vector>
#include <cassert>
#include <ctime>
#include <climits>
#include <algorithm>

#ifndef min
#define min(a,b) ( ((a)<(b))?(a):(b) )
#endif

//#define DEBUG

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
		unsigned int i, j;

		cl_platform_id platforms[10]; // Maximum 10 platforms
		cl_device_id devices_per_platform[10]; // Maximum 10 devices

		cl_uint platform_cnt;
		cl_uint device_cnt;

		vector<cl_device_id> devices;

		oclGetPlatformIDs(10, platforms, &platform_cnt);

		for (i = 0; i < platform_cnt; i++) {
			oclGetDeviceIDs(platforms[i], devtype, 10, devices_per_platform, &device_cnt);
			for (j = 0; j < device_cnt; j++) {
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
			"__kernel void finder(__global char* chr,\n"
			"                     __constant char* pat, __constant int* pat_index, unsigned int patternlen,\n"
			"                     __global char* flag, __global unsigned int* entrycount, __global unsigned int* loci)\n"
			"{\n"
			"	unsigned int i = get_global_id(0);\n"
			"	unsigned int j;\n"
			"	unsigned int old;\n"
			"	int k;\n"
			"	char localflag = 0;\n"
			"	for (j=0; j<patternlen; j++) {\n"
			"		k = pat_index[j];\n"
			"		if (k == -1)\n"
			"			break;\n"
			"		if ( (pat[k] == 'R' && (chr[i+k] == 'C' || chr[i+k] == 'T')) ||\n"
			"		     (pat[k] == 'Y' && (chr[i+k] == 'A' || chr[i+k] == 'G')) ||\n"
			"		     (pat[k] == 'K' && (chr[i+k] == 'A' || chr[i+k] == 'C')) ||\n"
			"		     (pat[k] == 'M' && (chr[i+k] == 'G' || chr[i+k] == 'T')) ||\n"
			"		     (pat[k] == 'W' && (chr[i+k] == 'C' || chr[i+k] == 'G')) ||\n"
			"		     (pat[k] == 'S' && (chr[i+k] == 'A' || chr[i+k] == 'T')) ||\n"
			"		     (pat[k] == 'H' && (chr[i+k] == 'G')) ||\n"
			"		     (pat[k] == 'B' && (chr[i+k] == 'A')) ||\n"
			"		     (pat[k] == 'V' && (chr[i+k] == 'T')) ||\n"
			"		     (pat[k] == 'D' && (chr[i+k] == 'C')) ||\n"
			"		     (pat[k] == 'A' && (chr[i+k] != 'A')) ||\n"
			"		     (pat[k] == 'G' && (chr[i+k] != 'G')) ||\n"
			"		     (pat[k] == 'C' && (chr[i+k] != 'C')) ||\n"
			"		     (pat[k] == 'T' && (chr[i+k] != 'T')) )\n"
			"			localflag |= 2;\n"
			"		k = pat_index[patternlen + j];\n"
			"		if ( (pat[k + patternlen] == 'R' && (chr[i+k] == 'C' || chr[i+k] == 'T')) ||\n"
			"		     (pat[k + patternlen] == 'Y' && (chr[i+k] == 'A' || chr[i+k] == 'G')) ||\n"
			"		     (pat[k + patternlen] == 'K' && (chr[i+k] == 'A' || chr[i+k] == 'C')) ||\n"
			"		     (pat[k + patternlen] == 'M' && (chr[i+k] == 'G' || chr[i+k] == 'T')) ||\n"
			"		     (pat[k + patternlen] == 'W' && (chr[i+k] == 'C' || chr[i+k] == 'G')) ||\n"
			"		     (pat[k + patternlen] == 'S' && (chr[i+k] == 'A' || chr[i+k] == 'T')) ||\n"
			"		     (pat[k + patternlen] == 'H' && (chr[i+k] == 'G')) ||\n"
			"		     (pat[k + patternlen] == 'B' && (chr[i+k] == 'A')) ||\n"
			"		     (pat[k + patternlen] == 'V' && (chr[i+k] == 'T')) ||\n"
			"		     (pat[k + patternlen] == 'D' && (chr[i+k] == 'C')) ||\n"
			"		     (pat[k + patternlen] == 'A' && (chr[i+k] != 'A')) ||\n"
			"		     (pat[k + patternlen] == 'G' && (chr[i+k] != 'G')) ||\n"
			"		     (pat[k + patternlen] == 'C' && (chr[i+k] != 'C')) ||\n"
			"		     (pat[k + patternlen] == 'T' && (chr[i+k] != 'T')) )\n"
			"			localflag |= 1;\n"
			"		if (localflag == 3)\n"
			"			break;\n"
			"	}\n"
			"	if (localflag != 3) {\n"
			"		printf(\"%d\\n\");\n"
			"		for (j=0; j<patternlen; j++)\n"
			"			if (chr[i+j] == ';') return;\n"
			"		old = atomic_inc(entrycount);\n"
			"		loci[old] = i;\n"
			"		flag[old] = localflag;\n"
			"	}\n"
			"}\n"
			"__kernel void comparer(__global char* chr, __global unsigned int* loci, __global unsigned int* mm_loci,\n"
			"                       __constant char* comp, __constant int* comp_index, unsigned int patternlen, unsigned short threshold,\n"
			"                       __global char* flag, __global unsigned short* mm_count, __global char* direction, __global unsigned int* entrycount)\n"
			"{\n"
			"	unsigned int i = get_global_id(0);\n"
			"	unsigned int j, lmm_count, old;\n"
			"	int k;\n"
			"	if (flag[i] == 0 || flag[i] == 1) {\n"
			"		lmm_count = 0;\n"
			"		for (j=0; j<patternlen; j++) {\n"
			"			k = comp_index[j];\n"
			"			if (k == -1) break;\n"
			"			if ( (comp[k] == 'R' && (chr[loci[i]+k] == 'C' || chr[loci[i]+k] == 'T')) ||\n"
			"			     (comp[k] == 'Y' && (chr[loci[i]+k] == 'A' || chr[loci[i]+k] == 'G')) ||\n"
			"			     (comp[k] == 'K' && (chr[loci[i]+k] == 'A' || chr[loci[i]+k] == 'C')) ||\n"
			"			     (comp[k] == 'M' && (chr[loci[i]+k] == 'G' || chr[loci[i]+k] == 'T')) ||\n"
			"			     (comp[k] == 'W' && (chr[loci[i]+k] == 'C' || chr[loci[i]+k] == 'G')) ||\n"
			"			     (comp[k] == 'S' && (chr[loci[i]+k] == 'A' || chr[loci[i]+k] == 'T')) ||\n"
			"			     (comp[k] == 'H' && (chr[loci[i]+k] == 'G')) ||\n"
			"			     (comp[k] == 'B' && (chr[loci[i]+k] == 'A')) ||\n"
			"			     (comp[k] == 'V' && (chr[loci[i]+k] == 'T')) ||\n"
			"			     (comp[k] == 'D' && (chr[loci[i]+k] == 'C')) ||\n"
			"				 (comp[k] == 'A' && (chr[loci[i]+k] != 'A')) ||\n"
			"			     (comp[k] == 'G' && (chr[loci[i]+k] != 'G')) ||\n"
			"			     (comp[k] == 'C' && (chr[loci[i]+k] != 'C')) ||\n"
			"			     (comp[k] == 'T' && (chr[loci[i]+k] != 'T'))) {\n"
			"				lmm_count++;\n"
			"				if (lmm_count > threshold) break;\n"
			"			}\n"
			"		}\n"
			"		if (lmm_count <= threshold) {\n"
			"			old = atomic_inc(entrycount);\n"
			"			mm_count[old] = lmm_count;\n"
			"			direction[old] = '+';\n"
			"			mm_loci[old] = loci[i];\n"
			"		}\n"
			"	}\n"
			"	if (flag[i] == 0 || flag[i] == 2) {\n"
			"		lmm_count = 0;\n"
			"		for (j=0; j<patternlen; j++) {\n"
			"			k = comp_index[patternlen + j];\n"
			"			if (k == -1) break;\n"
			"			if ( (comp[k+patternlen] == 'R' && (chr[loci[i]+k] == 'C' || chr[loci[i]+k] == 'T')) ||\n"
			"			     (comp[k+patternlen] == 'Y' && (chr[loci[i]+k] == 'A' || chr[loci[i]+k] == 'G')) ||\n"
			"			     (comp[k+patternlen] == 'K' && (chr[loci[i]+k] == 'A' || chr[loci[i]+k] == 'C')) ||\n"
			"			     (comp[k+patternlen] == 'M' && (chr[loci[i]+k] == 'G' || chr[loci[i]+k] == 'T')) ||\n"
			"			     (comp[k+patternlen] == 'W' && (chr[loci[i]+k] == 'C' || chr[loci[i]+k] == 'G')) ||\n"
			"			     (comp[k+patternlen] == 'S' && (chr[loci[i]+k] == 'A' || chr[loci[i]+k] == 'T')) ||\n"
			"			     (comp[k+patternlen] == 'H' && (chr[loci[i]+k] == 'G')) ||\n"
			"			     (comp[k+patternlen] == 'B' && (chr[loci[i]+k] == 'A')) ||\n"
			"			     (comp[k+patternlen] == 'V' && (chr[loci[i]+k] == 'T')) ||\n"
			"			     (comp[k+patternlen] == 'D' && (chr[loci[i]+k] == 'C')) ||\n"
			"			     (comp[k+patternlen] == 'A' && (chr[loci[i]+k] != 'A')) ||\n"
			"			     (comp[k+patternlen] == 'G' && (chr[loci[i]+k] != 'G')) ||\n"
			"			     (comp[k+patternlen] == 'C' && (chr[loci[i]+k] != 'C')) ||\n"
			"				 (comp[k+patternlen] == 'T' && (chr[loci[i]+k] != 'T'))) {\n"
			"				lmm_count++;\n"
			"				if (lmm_count > threshold) break;\n"
			"			}\n"
			"		}\n"
			"		if (lmm_count <= threshold) {\n"
			"			old = atomic_inc(entrycount);\n"
			"			mm_count[old] = lmm_count;\n"
			"			direction[old] = '-';\n"
			"			mm_loci[old] = loci[i];\n"
			"		}\n"
			"	}\n"
			"}";

		cl_context context;
		cl_program program;

		const size_t src_len = strlen(program_src);
		for (i = 0; i < devices.size(); i++) {
			// Create completely separate contexts per device to avoid unknown errors
			context = oclCreateContext(0, 1, &devices[i], 0, 0);
			m_contexts.push_back(context);
			program = oclCreateProgramWithSource(context, 1, &program_src, &src_len);
			oclBuildProgram(program, 1, &devices[i], "", 0, 0);
			m_finderkernels.push_back(oclCreateKernel(program, "finder"));
			m_comparerkernels.push_back(oclCreateKernel(program, "comparer"));
			m_queues.push_back(oclCreateCommandQueue(m_contexts[i], devices[i], 0));
			MAX_ALLOC_MEMORY.push_back(0);
			oclGetDeviceInfo(devices[i], CL_DEVICE_MAX_MEM_ALLOC_SIZE, sizeof(cl_ulong), &MAX_ALLOC_MEMORY[i], 0);
		}
		cout << "Total " << m_devnum << " device(s) found." << endl;
	}

public:
	Cas_OFFinder(cl_device_type devtype) {
		initOpenCL(devtype);
	}

	~Cas_OFFinder() {
		unsigned int i;
		for (i = 0; i < m_finderkernels.size(); i++)
			oclReleaseKernel(m_finderkernels[i]);
		for (i = 0; i < m_comparerkernels.size(); i++)
			oclReleaseKernel(m_comparerkernels[i]);
		for (i = 0; i < m_devnum; i++) {
			oclReleaseCommandQueue(m_queues[i]);
			oclReleaseContext(m_contexts[i]);
		}
	}

	void setChrData(vector<string> chrnames, string chrdata, vector<unsigned long long> chrpos) {
		m_chrnames = chrnames;
		m_chrpos = chrpos;
		m_chrdata = chrdata;
		m_chrdatasize = chrdata.size();
		m_totalanalyzedsize = 0;
		m_lasttotalanalyzedsize = 0;
		m_lastloci = 0;
	}

	void setPattern(const char* pattern) {
		unsigned int dev_index;
		m_pattern = (cl_char*)pattern;
		m_patternlen = (cl_uint)strlen(pattern);

		m_dicesizes.clear();
		clearbufvec(&m_chrdatabufs);
		clearbufvec(&m_patternbufs);
		clearbufvec(&m_patternindexbufs);
		clearbufvec(&m_flagbufs);
		clearbufvec(&m_locibufs);
		clearbufvec(&m_entrycountbufs);

		for (dev_index = 0; dev_index < m_devnum; dev_index++) {
			m_dicesizes.push_back(
				(size_t)min(
				(MAX_ALLOC_MEMORY[dev_index] - sizeof(cl_char)* (3 * m_patternlen - 1) - sizeof(cl_uint)* (2 * m_patternlen + 3) - sizeof(cl_ushort)) / (4 * sizeof(cl_char)+3 * sizeof(cl_uint)+2 * sizeof(cl_ushort)),
				((m_chrdatasize / m_devnum) + ((m_chrdatasize%m_devnum == 0) ? 0 : 1))
				)
				); // No more than maximum allocation per device
			// cout << "Dicesize: " << m_dicesizes[dev_index] << endl;
			m_chrdatabufs.push_back(oclCreateBuffer(m_contexts[dev_index], CL_MEM_READ_ONLY, sizeof(cl_char)* (m_dicesizes[dev_index] + m_patternlen - 1), 0));
			m_patternbufs.push_back(oclCreateBuffer(m_contexts[dev_index], CL_MEM_READ_ONLY, sizeof(cl_char)* m_patternlen * 2, 0));
			m_patternindexbufs.push_back(oclCreateBuffer(m_contexts[dev_index], CL_MEM_READ_ONLY, sizeof(cl_int)* m_patternlen * 2, 0));
			m_flagbufs.push_back(oclCreateBuffer(m_contexts[dev_index], CL_MEM_WRITE_ONLY, sizeof(cl_char)* m_dicesizes[dev_index], 0));
			m_entrycountbufs.push_back(oclCreateBuffer(m_contexts[dev_index], CL_MEM_READ_WRITE, sizeof(cl_uint), 0));
			m_locibufs.push_back(oclCreateBuffer(m_contexts[dev_index], CL_MEM_WRITE_ONLY, sizeof(cl_uint)* m_dicesizes[dev_index], 0));
		}
	}

	bool loadNextChunk() {
		if (m_totalanalyzedsize == m_chrdatasize)
			return false;

		unsigned int dev_index;
		unsigned long long tailsize;

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
				oclEnqueueWriteBuffer(m_queues[dev_index], m_chrdatabufs[dev_index], CL_FALSE, 0, (size_t)(sizeof(cl_char)* (tailsize + m_patternlen - 1)), (cl_char *)m_chrdata.c_str() + m_totalanalyzedsize, 0, 0, 0);
				m_totalanalyzedsize += tailsize;
				m_worksizes.push_back(tailsize);
				#ifdef DEBUG
				cout << "Worksize: " <<  m_worksizes[dev_index] << ", Tailsize: " << tailsize << endl;
				#endif
				break;
			}
			else {
				oclEnqueueWriteBuffer(m_queues[dev_index], m_chrdatabufs[dev_index], CL_FALSE, 0, sizeof(cl_char)* (m_dicesizes[dev_index] + m_patternlen - 1), (cl_char *)m_chrdata.c_str() + m_totalanalyzedsize, 0, 0, 0);
				m_totalanalyzedsize += m_dicesizes[dev_index];
				m_worksizes.push_back(m_dicesizes[dev_index]);
				#ifdef DEBUG
				cout << "Worksize: " <<  m_worksizes[dev_index] << ", Tailsize: " << tailsize << endl;
				#endif
			}
		}
		cout << m_activedevnum << " devices selected to analyze..." << endl;

		for (dev_index = 0; dev_index < m_activedevnum; dev_index++) {
			cl_uint zero = 0;
			#ifdef DEBUG
			cout << "Writing buffer for find..." << endl;
			#endif
			oclEnqueueWriteBuffer(m_queues[dev_index], m_patternbufs[dev_index], CL_FALSE, 0, sizeof(cl_char)* m_patternlen, m_pattern, 0, 0, 0);
			oclEnqueueWriteBuffer(m_queues[dev_index], m_patternbufs[dev_index], CL_FALSE, sizeof(cl_char)* m_patternlen, sizeof(cl_char)* m_patternlen, c_pattern, 0, 0, 0);
			oclEnqueueWriteBuffer(m_queues[dev_index], m_patternindexbufs[dev_index], CL_FALSE, 0, sizeof(cl_int)* m_patternlen, pattern_index, 0, 0, 0);
			oclEnqueueWriteBuffer(m_queues[dev_index], m_patternindexbufs[dev_index], CL_FALSE, sizeof(cl_int)* m_patternlen, sizeof(cl_int)* m_patternlen, c_pattern_index, 0, 0, 0);
			oclEnqueueWriteBuffer(m_queues[dev_index], m_entrycountbufs[dev_index], CL_FALSE, 0, sizeof(cl_uint), &zero, 0, 0, 0);
			oclFinish(m_queues[dev_index]);
			#ifdef DEBUG
			cout << "Done." << endl;
			#endif
			oclSetKernelArg(m_finderkernels[dev_index], 0, sizeof(cl_mem), &m_chrdatabufs[dev_index]);
			oclSetKernelArg(m_finderkernels[dev_index], 1, sizeof(cl_mem), &m_patternbufs[dev_index]);
			oclSetKernelArg(m_finderkernels[dev_index], 2, sizeof(cl_mem), &m_patternindexbufs[dev_index]);
			oclSetKernelArg(m_finderkernels[dev_index], 3, sizeof(cl_uint), &m_patternlen);
			oclSetKernelArg(m_finderkernels[dev_index], 4, sizeof(cl_mem), &m_flagbufs[dev_index]);
			oclSetKernelArg(m_finderkernels[dev_index], 5, sizeof(cl_mem), &m_entrycountbufs[dev_index]);
			oclSetKernelArg(m_finderkernels[dev_index], 6, sizeof(cl_mem), &m_locibufs[dev_index]);
		}

		free((void*)pattern_index);
		free((void*)c_pattern);
		free((void*)c_pattern_index);

		return true;
	}

	void findPattern() {
		unsigned int dev_index;
		for (dev_index = 0; dev_index < m_activedevnum; dev_index++) {
			const size_t worksize = (size_t)m_worksizes[dev_index];
			oclEnqueueNDRangeKernel(m_queues[dev_index], m_finderkernels[dev_index], 1, 0, &worksize, 0, 0, 0, 0);
		}

		for (dev_index = 0; dev_index < m_activedevnum; dev_index++) {
			oclFinish(m_queues[dev_index]);
			m_locicnts.push_back(0);
			oclEnqueueReadBuffer(m_queues[dev_index], m_entrycountbufs[dev_index], CL_TRUE, 0, sizeof(cl_uint), &m_locicnts[dev_index], 0, 0, 0);
			if (m_locicnts[dev_index] > 0) {
				m_flags.push_back((cl_char *)malloc(sizeof(cl_char)* m_locicnts[dev_index]));
				oclEnqueueReadBuffer(m_queues[dev_index], m_flagbufs[dev_index], CL_TRUE, 0, sizeof(cl_char)*m_locicnts[dev_index], m_flags[dev_index], 0, 0, 0);

				m_mmcounts.push_back((cl_ushort *)malloc(sizeof(cl_ushort)* m_locicnts[dev_index] * 2)); // Maximum numbers of mismatch counts
				m_directions.push_back((cl_char *)malloc(sizeof(cl_char)* m_locicnts[dev_index] * 2));
				m_mmlocis.push_back((cl_uint *)malloc(sizeof(cl_uint)* m_locicnts[dev_index] * 2));

				m_mmlocibufs.push_back(oclCreateBuffer(m_contexts[dev_index], CL_MEM_WRITE_ONLY, sizeof(cl_uint)* m_locicnts[dev_index] * 2, 0));
				m_mmcountbufs.push_back(oclCreateBuffer(m_contexts[dev_index], CL_MEM_WRITE_ONLY, sizeof(cl_ushort)* m_locicnts[dev_index] * 2, 0));
				m_directionbufs.push_back(oclCreateBuffer(m_contexts[dev_index], CL_MEM_WRITE_ONLY, sizeof(cl_char)* m_locicnts[dev_index] * 2, 0));

				oclSetKernelArg(m_comparerkernels[dev_index], 1, sizeof(cl_mem), &m_locibufs[dev_index]);
				oclSetKernelArg(m_comparerkernels[dev_index], 2, sizeof(cl_mem), &m_mmlocibufs[dev_index]);
				oclSetKernelArg(m_comparerkernels[dev_index], 7, sizeof(cl_mem), &m_flagbufs[dev_index]);
			} else {
				m_flags.push_back(0);
				m_mmcounts.push_back(0);
				m_directions.push_back(0);
				m_mmlocis.push_back(0);
				m_mmlocibufs.push_back(0);
				m_mmcountbufs.push_back(0);
				m_directionbufs.push_back(0);
			}
		}
	}

	void releaseLociinfo() {
		unsigned int dev_index;

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
		clearbufvec(&m_mmlocibufs);
		m_flags.clear();
		clearbufvec(&m_mmcountbufs);
		clearbufvec(&m_directionbufs);
	}

	void indicate_mismatches(cl_char* seq, cl_char* comp) {
		unsigned int k;
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
		unsigned int i, j, dev_index;
		cl_int err;
		cl_uint zero = 0;

		vector <cl_mem> comparebufs;
		vector <cl_mem> compareindexbufs;

		cl_char *compare = (cl_char *)arg_compare;
		cl_char *c_compare = (cl_char *)malloc(sizeof(cl_char)* (m_patternlen + 1)); c_compare[m_patternlen] = 0; memcpy(c_compare, compare, m_patternlen); set_complementary_sequence(c_compare);
		int *c_compare_index = (cl_int *)malloc(sizeof(cl_int)* m_patternlen); set_pattern_index(c_compare_index, c_compare);
		int *compare_index = (cl_int *)malloc(sizeof(cl_int)* m_patternlen); set_pattern_index(compare_index, compare);

		for (dev_index = 0; dev_index<m_activedevnum; dev_index++) {
			if (m_locicnts[dev_index] > 0) {
				comparebufs.push_back(clCreateBuffer(m_contexts[dev_index], CL_MEM_READ_ONLY, sizeof(cl_char)* m_patternlen * 2, 0, &err));
				compareindexbufs.push_back(clCreateBuffer(m_contexts[dev_index], CL_MEM_READ_ONLY, sizeof(cl_uint)* m_patternlen * 2, 0, &err));

				const cl_char *compare = (const cl_char*)arg_compare;

				#ifdef DEBUG
				cout << "Writing compare buffer (frontside)..." << endl;
				#endif
				oclEnqueueWriteBuffer(m_queues[dev_index], comparebufs[dev_index], CL_FALSE, 0, sizeof(cl_char)* m_patternlen, compare, 0, 0, 0);
				#ifdef DEBUG
				cout << "Writing compare buffer (backside)..." << endl;
				#endif
				oclEnqueueWriteBuffer(m_queues[dev_index], comparebufs[dev_index], CL_FALSE, sizeof(cl_char)* m_patternlen, sizeof(cl_char)* m_patternlen, c_compare, 0, 0, 0);
				#ifdef DEBUG
				cout << "Writing index buffer (frontside)..." << endl;
				#endif
				oclEnqueueWriteBuffer(m_queues[dev_index], compareindexbufs[dev_index], CL_FALSE, 0, sizeof(cl_int)* m_patternlen, compare_index, 0, 0, 0);
				#ifdef DEBUG
				cout << "Writing index buffer (backside)..." << endl;
				#endif
				oclEnqueueWriteBuffer(m_queues[dev_index], compareindexbufs[dev_index], CL_FALSE, sizeof(cl_int)* m_patternlen, sizeof(cl_int)* m_patternlen, c_compare_index, 0, 0, 0);
				#ifdef DEBUG
				cout << "Writing entry count buffer..." << endl;
				#endif
				oclEnqueueWriteBuffer(m_queues[dev_index], m_entrycountbufs[dev_index], CL_FALSE, 0, sizeof(cl_uint), &zero, 0, 0, 0);
				oclFinish(m_queues[dev_index]);
				#ifdef DEBUG
				cout << "Done." << endl;
				#endif
				cl_ushort cl_threshold = threshold;
				oclSetKernelArg(m_comparerkernels[dev_index], 0, sizeof(cl_mem), &m_chrdatabufs[dev_index]);
				//m_comparerkernels[dev_index].setArg(1, m_locibufs[dev_index]);
				//m_comparerkernels[dev_index].setArg(2, m_mmlocibufs[dev_index]);
				oclSetKernelArg(m_comparerkernels[dev_index], 3, sizeof(cl_mem), &comparebufs[dev_index]);
				oclSetKernelArg(m_comparerkernels[dev_index], 4, sizeof(cl_mem), &compareindexbufs[dev_index]);
				oclSetKernelArg(m_comparerkernels[dev_index], 5, sizeof(cl_uint), &m_patternlen);
				oclSetKernelArg(m_comparerkernels[dev_index], 6, sizeof(cl_ushort), &cl_threshold);
				//m_comparerkernels[dev_index].setArg(7, m_flagbufs[dev_index]);
				oclSetKernelArg(m_comparerkernels[dev_index], 8, sizeof(cl_mem), &m_mmcountbufs[dev_index]);
				oclSetKernelArg(m_comparerkernels[dev_index], 9, sizeof(cl_mem), &m_directionbufs[dev_index]);
				oclSetKernelArg(m_comparerkernels[dev_index], 10, sizeof(cl_mem), &m_entrycountbufs[dev_index]);
				const size_t locicnts = m_locicnts[dev_index];
				oclEnqueueNDRangeKernel(m_queues[dev_index], m_comparerkernels[dev_index], 1, 0, &locicnts, 0, 0, 0, 0);
			} else {
				comparebufs.push_back(0);
				compareindexbufs.push_back(0);
			}
		}

		unsigned long long loci;

		char comp_symbol[2] = { '+', '-' };
		char *strbuf = (char *)malloc(sizeof(char)* (m_patternlen + 1)); strbuf[m_patternlen] = 0;

		ofstream fo(outfilename, ios::out | ios::app);
		unsigned long long localanalyzedsize = 0;
		unsigned int cnt = 0;
		unsigned int idx;
		for (dev_index = 0; dev_index<m_activedevnum; dev_index++) {
			 if (m_locicnts[dev_index] > 0) {
				oclFinish(m_queues[dev_index]);
				oclEnqueueReadBuffer(m_queues[dev_index], m_entrycountbufs[dev_index], CL_TRUE, 0, sizeof(cl_uint), &cnt, 0, 0, 0);
				if (cnt > 0) {
					oclEnqueueReadBuffer(m_queues[dev_index], m_mmcountbufs[dev_index], CL_TRUE, 0, sizeof(cl_ushort)* cnt, m_mmcounts[dev_index], 0, 0, 0);
					oclEnqueueReadBuffer(m_queues[dev_index], m_directionbufs[dev_index], CL_TRUE, 0, sizeof(cl_char)* cnt, m_directions[dev_index], 0, 0, 0);
					oclEnqueueReadBuffer(m_queues[dev_index], m_mmlocibufs[dev_index], CL_TRUE, 0, sizeof(cl_uint)* cnt, m_mmlocis[dev_index], 0, 0, 0);
					for (i = 0; i < cnt; i++) {
						loci = m_mmlocis[dev_index][i] + m_lasttotalanalyzedsize + localanalyzedsize;
						if (m_mmcounts[dev_index][i] <= threshold) {
							strncpy(strbuf, (char *)(m_chrdata.c_str() + loci), m_patternlen);
							if (m_directions[dev_index][i] == '-') set_complementary_sequence((cl_char *)strbuf);
							indicate_mismatches((cl_char*)strbuf, compare);
							for (j = 0; ((j < m_chrpos.size()) && (loci > m_chrpos[j])); j++) idx = j;
							fo << compare << "\t" << m_chrnames[idx] << "\t" << loci - m_chrpos[idx] << "\t" << strbuf << "\t" << m_directions[dev_index][i] << "\t" << m_mmcounts[dev_index][i] << endl;
						}
					}
				}
			}
			localanalyzedsize += m_worksizes[dev_index];
		}
		fo.close();
		clearbufvec(&comparebufs);
		clearbufvec(&compareindexbufs);
		free((void *)strbuf);
		free((void *)c_compare);
		free((void *)c_compare_index);
		free((void *)compare_index);
	}
};

int main(int argc, char *argv[]) {
	clock_t start, end;
	float seconds;
	string fastapath, chrname, pattern, tmpstr, chrdir, outfilename;
	DIR* dir;
	dirent *ent;
	unsigned int i, j, cnt, tmpint;

	std::ifstream input;
	std::string line, name, content;

	cl_platform_id platforms[10];
	cl_uint platform_cnt;

	oclGetPlatformIDs(10, platforms, &platform_cnt);

	if (platform_cnt == 0) {
		cout << "No OpenCL platforms found. Check OpenCL installation!" << endl;
		exit(1);
	}

	if (argc < 4) { // Not all option specified
		cout << "Cas-OFFinder v2.3 (2015-03-08)" << endl <<
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

		cl_device_id devices_per_platform[10];
		cl_uint device_cnt;
		cl_char devname[255] = {0, };
		for (i = 0; i < platform_cnt; i++) {
			oclGetDeviceIDs(platforms[i], CL_DEVICE_TYPE_CPU, 10, devices_per_platform, &device_cnt);
			for (j = 0; j < device_cnt; j++) {
				oclGetDeviceInfo(devices_per_platform[j], CL_DEVICE_NAME, 255, &devname, 0);
				cout << "Type: CPU, '" << devname << "'" << endl;
			}
			oclGetDeviceIDs(platforms[i], CL_DEVICE_TYPE_GPU, 10, devices_per_platform, &device_cnt);
			for (j = 0; j < device_cnt; j++) {
				oclGetDeviceInfo(devices_per_platform[j], CL_DEVICE_NAME, 255, &devname, 0);
				cout << "Type: GPU, '" << devname << "'" << endl;
			}
			oclGetDeviceIDs(platforms[i], CL_DEVICE_TYPE_ACCELERATOR, 10, devices_per_platform, &device_cnt);
			for (j = 0; j < device_cnt; j++) {
				oclGetDeviceInfo(devices_per_platform[j], CL_DEVICE_NAME, 255, &devname, 0);
				cout << "Type: ACCELERATOR, '" << devname << "'" << endl;
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
	Cas_OFFinder s(devtype);

	vector <string> compare;
	vector <string> chrnames;
	vector <unsigned long long> chrpos;
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
		transform(tmpstr.begin(), tmpstr.end(), tmpstr.begin(), ::toupper);
		compare.push_back(tmpstr);
		if (fi.eof()) break;
		fi >> tmpint;
		threshold.push_back(tmpint);
	}
	fi.close();

	transform(pattern.begin(), pattern.end(), pattern.begin(), ::toupper);

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
				cout << "Reading " << fastapath << "..." << endl;
				cnum = 0;
				input.open(fastapath.c_str());
				cnt = 0;
				chrnames.clear();
				content.clear();
				chrpos.clear();
				while (getline(input, line).good()){
					if (!line.empty() && line[0] == '>') { // Identifier marker
						name = line.substr(1);
						//if (((cnt++) % 10000) == 0) cout << "Reading " << name << endl;
						chrnames.push_back(name);
						if (chrpos.size() != 0) content += ";"; // seperator
						chrpos.push_back(content.size());
					}
					else if (!name.empty()){
						transform(line.begin(), line.end(), line.begin(), ::toupper);
						try {
							content += line;
						} catch (const std::bad_alloc&) {
							cout << "Not enough memory!" << endl << "Split fasta file, or try 64bit version of Cas-OFFinder." << endl;
							exit(1);
						}
					}
				}
				input.close();
				cout << "Sending data to devices..." << endl;
				s.setChrData(chrnames, content, chrpos);
				cout << "Setting pattern to devices..." << endl;
				s.setPattern(pattern.c_str());
				cout << "Chunk load started." << endl;
				while (s.loadNextChunk()) {
					// Find patterns in the chunk
					cout << "Finding pattern in chunk #" << ++cnum << "..." << endl;
					s.findPattern();
					cout << "Comparing patterns in chunk #" << cnum << "..." << endl;
					for (i = 0; i < threshold.size(); i++) {
						s.compareAll(compare[i].c_str(), threshold[i], (outfilename).c_str());
					}
					s.releaseLociinfo();
				}			
			}
		}
	}
	end = clock(); seconds = (float)(end - start) / CLOCKS_PER_SEC; cout << seconds << " seconds elapsed." << endl;
	return 0;
}
