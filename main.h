#include <CL/cl.h>
#include <iostream>
#include <vector>

using namespace std;

cl_mem oclCreateBuffer(cl_context context, cl_mem_flags flags, size_t size, void *host_ptr) {
	cl_int err;
	cl_mem mem;
	mem = clCreateBuffer(context, flags, size, host_ptr, &err);
	if (err != CL_SUCCESS) {
		cout << "clCreateBuffer Failed: " << err << endl;
		exit(1);
	}
	return mem;
}

void clearbufvec(vector<cl_mem> *bufvec) {
	unsigned int i;
	for (i = 0; i < bufvec->size(); i++) {
		if ((*bufvec)[i] != 0) clReleaseMemObject((*bufvec)[i]);
	}
	bufvec->clear();
}