#pragma OPENCL EXTENSION cl_khr_global_int32_base_atomics : enable
#pragma OPENCL EXTENSION cl_khr_local_int32_base_atomics : enable
#pragma OPENCL EXTENSION cl_khr_global_int32_extended_atomics : enable
#pragma OPENCL EXTENSION cl_khr_local_int32_extended_atomics : enable

__kernel void finder(__global char* chr,
                     __constant char* pat, __constant int* pat_index, unsigned int patternlen,
                     __global char* flag, __global unsigned int* entrycount, __global unsigned int* loci)
{
	unsigned int i = get_global_id(0);
	unsigned int j;
	unsigned int old;
	int k;

	char localflag = 0;
	for (j=0; j<patternlen; j++) {
		k = pat_index[j];
		if (k == -1)
			break;
		if ( (pat[k] == 'R' && (chr[i+k] == 'C' || chr[i+k] == 'T')) ||
		     (pat[k] == 'Y' && (chr[i+k] == 'A' || chr[i+k] == 'G')) ||
		     (pat[k] == 'K' && (chr[i+k] == 'A' || chr[i+k] == 'C')) ||
		     (pat[k] == 'M' && (chr[i+k] == 'G' || chr[i+k] == 'T')) ||
		     (pat[k] == 'W' && (chr[i+k] == 'C' || chr[i+k] == 'G')) ||
		     (pat[k] == 'S' && (chr[i+k] == 'A' || chr[i+k] == 'T')) ||
		     (pat[k] == 'H' && (chr[i+k] == 'G')) ||
		     (pat[k] == 'B' && (chr[i+k] == 'A')) ||
		     (pat[k] == 'V' && (chr[i+k] == 'T')) ||
		     (pat[k] == 'D' && (chr[i+k] == 'C')) ||
		     (pat[k] == 'A' && (chr[i+k] != 'A')) ||
		     (pat[k] == 'G' && (chr[i+k] != 'G')) ||
		     (pat[k] == 'C' && (chr[i+k] != 'C')) ||
		     (pat[k] == 'T' && (chr[i+k] != 'T')) )
			localflag |= 2;
		k = pat_index[patternlen + j];
		if ( (pat[k + patternlen] == 'R' && (chr[i+k] == 'C' || chr[i+k] == 'T')) ||
		     (pat[k + patternlen] == 'Y' && (chr[i+k] == 'A' || chr[i+k] == 'G')) ||
		     (pat[k + patternlen] == 'K' && (chr[i+k] == 'A' || chr[i+k] == 'C')) ||
		     (pat[k + patternlen] == 'M' && (chr[i+k] == 'G' || chr[i+k] == 'T')) ||
		     (pat[k + patternlen] == 'W' && (chr[i+k] == 'C' || chr[i+k] == 'G')) ||
		     (pat[k + patternlen] == 'S' && (chr[i+k] == 'A' || chr[i+k] == 'T')) ||
		     (pat[k + patternlen] == 'H' && (chr[i+k] == 'G')) ||
		     (pat[k + patternlen] == 'B' && (chr[i+k] == 'A')) ||
		     (pat[k + patternlen] == 'V' && (chr[i+k] == 'T')) ||
		     (pat[k + patternlen] == 'D' && (chr[i+k] == 'C')) ||
		     (pat[k + patternlen] == 'A' && (chr[i+k] != 'A')) ||
		     (pat[k + patternlen] == 'G' && (chr[i+k] != 'G')) ||
		     (pat[k + patternlen] == 'C' && (chr[i+k] != 'C')) ||
		     (pat[k + patternlen] == 'T' && (chr[i+k] != 'T')) )
			localflag |= 1;
		if (localflag == 3)
			break;
	}
	if (localflag != 3) {
		for (j=0; j<patternlen; j++)
			if (chr[i+j] == ';') return;
		old = atomic_inc(entrycount);
		loci[old] = i;
		flag[old] = localflag;
	}
}
__kernel void comparer(__global char* chr, __global unsigned int* loci, __global unsigned int* mm_loci,
                       __constant char* comp, __constant int* comp_index, unsigned int patternlen, unsigned short threshold,
                       __global char* flag, __global unsigned short* mm_count, __global char* direction, __global unsigned int* entrycount,
					   __local char* l_comp, __local int* l_comp_index)
{
	unsigned int i = get_global_id(0);
	unsigned int j, lmm_count, old;
	int k;

	unsigned int li = i - get_group_id(0)*get_local_size(0);
	if (li < patternlen*2) {
		l_comp[li] = comp[li];
		l_comp_index[li] = comp_index[li];
	}
	barrier(CLK_LOCAL_MEM_FENCE);
	if (flag[i] == 0 || flag[i] == 1) {
		lmm_count = 0;
		for (j=0; j<patternlen; j++) {
			k = l_comp_index[j];
			if (k == -1) break;
			if ( (l_comp[k] == 'R' && (chr[loci[i]+k] == 'C' || chr[loci[i]+k] == 'T')) ||
			     (l_comp[k] == 'Y' && (chr[loci[i]+k] == 'A' || chr[loci[i]+k] == 'G')) ||
			     (l_comp[k] == 'K' && (chr[loci[i]+k] == 'A' || chr[loci[i]+k] == 'C')) ||
			     (l_comp[k] == 'M' && (chr[loci[i]+k] == 'G' || chr[loci[i]+k] == 'T')) ||
			     (l_comp[k] == 'W' && (chr[loci[i]+k] == 'C' || chr[loci[i]+k] == 'G')) ||
			     (l_comp[k] == 'S' && (chr[loci[i]+k] == 'A' || chr[loci[i]+k] == 'T')) ||
			     (l_comp[k] == 'H' && (chr[loci[i]+k] == 'G')) ||
			     (l_comp[k] == 'B' && (chr[loci[i]+k] == 'A')) ||
			     (l_comp[k] == 'V' && (chr[loci[i]+k] == 'T')) ||
			     (l_comp[k] == 'D' && (chr[loci[i]+k] == 'C')) ||
				 (l_comp[k] == 'A' && (chr[loci[i]+k] != 'A')) ||
			     (l_comp[k] == 'G' && (chr[loci[i]+k] != 'G')) ||
			     (l_comp[k] == 'C' && (chr[loci[i]+k] != 'C')) ||
			     (l_comp[k] == 'T' && (chr[loci[i]+k] != 'T'))) {
				lmm_count++;
				if (lmm_count > threshold) break;
			}
		}
		if (lmm_count <= threshold) {
			old = atomic_inc(entrycount);
			mm_count[old] = lmm_count;
			direction[old] = '+';
			mm_loci[old] = loci[i];
		}
	}
	if (flag[i] == 0 || flag[i] == 2) {
	lmm_count = 0;
		for (j=0; j<patternlen; j++) {
			k = l_comp_index[patternlen + j];
			if (k == -1) break;
			if ( (l_comp[k+patternlen] == 'R' && (chr[loci[i]+k] == 'C' || chr[loci[i]+k] == 'T')) ||
			     (l_comp[k+patternlen] == 'Y' && (chr[loci[i]+k] == 'A' || chr[loci[i]+k] == 'G')) ||
			     (l_comp[k+patternlen] == 'K' && (chr[loci[i]+k] == 'A' || chr[loci[i]+k] == 'C')) ||
			     (l_comp[k+patternlen] == 'M' && (chr[loci[i]+k] == 'G' || chr[loci[i]+k] == 'T')) ||
			     (l_comp[k+patternlen] == 'W' && (chr[loci[i]+k] == 'C' || chr[loci[i]+k] == 'G')) ||
			     (l_comp[k+patternlen] == 'S' && (chr[loci[i]+k] == 'A' || chr[loci[i]+k] == 'T')) ||
			     (l_comp[k+patternlen] == 'H' && (chr[loci[i]+k] == 'G')) ||
			     (l_comp[k+patternlen] == 'B' && (chr[loci[i]+k] == 'A')) ||
			     (l_comp[k+patternlen] == 'V' && (chr[loci[i]+k] == 'T')) ||
			     (l_comp[k+patternlen] == 'D' && (chr[loci[i]+k] == 'C')) ||
			     (l_comp[k+patternlen] == 'A' && (chr[loci[i]+k] != 'A')) ||
			     (l_comp[k+patternlen] == 'G' && (chr[loci[i]+k] != 'G')) ||
			     (l_comp[k+patternlen] == 'C' && (chr[loci[i]+k] != 'C')) ||
				 (l_comp[k+patternlen] == 'T' && (chr[loci[i]+k] != 'T'))) {
				lmm_count++;
				if (lmm_count > threshold) break;
			}
		}
		if (lmm_count <= threshold) {
			old = atomic_inc(entrycount);
			mm_count[old] = lmm_count;
			direction[old] = '-';
			mm_loci[old] = loci[i];
		}
	}
}