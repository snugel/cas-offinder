#ifdef __APPLE__
#  include <OpenCL/cl.h>
#  define CL_CALLBACK
#else
#  include <CL/cl.h>
#endif
#include <iostream>
#include <vector>
#pragma warning (disable : 4996)
using namespace std;

cl_mem oclCreateBuffer(cl_context context,
	cl_mem_flags flags,
	size_t size,
	void *host_ptr)
{
	cl_int err;
	cl_mem mem = clCreateBuffer(context, flags, size, host_ptr, &err);
	if (err != CL_SUCCESS) {
		cerr << "clCreateBuffer Failed: " << err << endl;
		exit(1);
	}
	return mem;
}

void oclGetPlatformIDs(cl_uint num_entries,
	cl_platform_id *platforms,
	cl_uint *num_platforms)
{
	cl_int err = clGetPlatformIDs(num_entries, platforms, num_platforms);
	if (err != CL_SUCCESS) {
		cerr << "clGetPlatformIDs Failed: " << err << endl;
		exit(1);
	}
}

void oclGetDeviceIDs(cl_platform_id platform,
	cl_device_type device_type, cl_uint num_entries,
	cl_device_id *devices,
	cl_uint *num_devices)
{
	cl_int err = clGetDeviceIDs(platform, device_type, num_entries, devices, num_devices);
	if (err == CL_DEVICE_NOT_FOUND) (*num_devices) = 0;
	if (err != CL_SUCCESS && err != CL_DEVICE_NOT_FOUND) {
		cerr << "clGetDeviceIDs Failed: " << err << endl;
		exit(1);
	}
}

cl_context oclCreateContext(cl_context_properties *properties,
	cl_uint num_devices,
	const cl_device_id *devices,
	void (CL_CALLBACK *pfn_notify)(
	const char *errinfo,
	const void *private_info,
	size_t cb,
	void *user_data
	),
	void *user_data)
{
	cl_int err;
	cl_context context = clCreateContext(properties, num_devices, devices, pfn_notify, user_data, &err);
	if (err != CL_SUCCESS) {
		cerr << "clCreateContext Failed: " << err << endl;
		exit(1);
	}
	return context;
}

cl_program oclCreateProgramWithSource(cl_context context,
	cl_uint count,
	const char **strings,
	const size_t *lengths)
{
	cl_int err;
	cl_program program = clCreateProgramWithSource(context, count, strings, lengths, &err);
	if (err != CL_SUCCESS) {
		cerr << "clCreateProgramWithSource Failed: " << err << endl;
		exit(1);
	}
	return program;
}

void oclBuildProgram(cl_program program,
	cl_uint num_devices,
	const cl_device_id *device_list,
	const char *options,
	void(CL_CALLBACK *pfn_notify)(cl_program, void *user_data),
	void *user_data) {
	cl_int err = clBuildProgram(program, num_devices, device_list, options, pfn_notify, user_data);
	if (err != CL_SUCCESS) {
		cerr << "clBuildProgram Failed: " << err << endl;
        if (err == CL_BUILD_PROGRAM_FAILURE) {
            for (int i = 0; i < num_devices; i++) {
                // Determine the size of the log
                size_t log_size;
                clGetProgramBuildInfo(program, device_list[i], CL_PROGRAM_BUILD_LOG, 0, NULL, &log_size);

                // Allocate memory for the log
                char *log = (char *)malloc(log_size);

                // Get the log
                clGetProgramBuildInfo(program, device_list[i], CL_PROGRAM_BUILD_LOG, log_size, log, NULL);

                // Print the log
                printf("Log from device %d:\n", i);
                printf("%s\n", log);
                free(log);
            }
        }
        exit(1);
	}
}

cl_kernel oclCreateKernel(cl_program  program,
	const char *kernel_name)
{
	cl_int err;
	cl_kernel kernel = clCreateKernel(program, kernel_name, &err);
	if (err != CL_SUCCESS) {
		cerr << "clCreateKernel Failed: " << err << endl;
		exit(1);
	}
	return kernel;
}

cl_command_queue oclCreateCommandQueue(cl_context context,
	cl_device_id device,
	cl_command_queue_properties properties)
{
	cl_int err;
	cl_command_queue command_queue = clCreateCommandQueue(context, device, properties, &err);
	if (err != CL_SUCCESS) {
		cerr << "clCreateCommandQueue Failed: " << err << endl;
		exit(1);
	}
	return command_queue;
}

void oclGetPlatformInfo(cl_platform_id platform,
    cl_platform_info param_name,
    size_t param_value_size,
    void *param_value,
    size_t *param_value_size_ret)
{
    cl_int err = clGetPlatformInfo(platform, param_name, param_value_size, param_value, param_value_size_ret);
    if (err != CL_SUCCESS) {
        cerr << "clGetDeviceInfo Failed: " << err << endl;
        exit(1);
    }
}

void oclGetDeviceInfo(cl_device_id device,
	cl_device_info param_name,
	size_t param_value_size,
	void *param_value,
	size_t *param_value_size_ret)
{
	cl_int err = clGetDeviceInfo(device, param_name, param_value_size, param_value, param_value_size_ret);
	if (err != CL_SUCCESS) {
		cerr << "clGetDeviceInfo Failed: " << err << endl;
		exit(1);
	}
}

void oclReleaseKernel(cl_kernel kernel)
{
	cl_int err;
	if (kernel != 0) {
		err = clReleaseKernel(kernel);
		if (err != CL_SUCCESS) {
			cerr << "clReleaseKernel Failed: " << err << endl;
			exit(1);
		}
	}
}

void oclReleaseCommandQueue(cl_command_queue command_queue)
{
	cl_int err;
	if (command_queue != 0) {
		err = clReleaseCommandQueue(command_queue);
		if (err != CL_SUCCESS) {
			cerr << "clReleaseCommandQueue Failed: " << err << endl;
			exit(1);
		}
	}
}

void oclReleaseContext(cl_context context)
{
	cl_int err;
	if (context != 0) {
		err = clReleaseContext(context);
		if (err != CL_SUCCESS) {
			cerr << "clReleaseContext Failed: " << err << endl;
			exit(1);
		}
	}
}

void oclEnqueueWriteBuffer(cl_command_queue command_queue,
	cl_mem buffer,
	cl_bool blocking_write,
	size_t offset,
	size_t cb,
	const void *ptr,
	cl_uint num_events_in_wait_list,
	const cl_event *event_wait_list,
	cl_event *event)
{
	cl_int err = clEnqueueWriteBuffer(command_queue, buffer, blocking_write, offset, cb, ptr, num_events_in_wait_list, event_wait_list, event);
	if (err != CL_SUCCESS) {
		cerr << "clEnqueueWriteBuffer Failed: " << err << endl;
		exit(1);
	}
}

void oclFinish(cl_command_queue command_queue)
{
	cl_int err = clFinish(command_queue);
	if (err != CL_SUCCESS) {
		cerr << "clFinish Failed: " << err << endl;
		exit(1);
	}
}

void oclSetKernelArg(cl_kernel kernel,
	cl_uint arg_index,
	size_t arg_size,
	const void *arg_value)
{
	cl_int err = clSetKernelArg(kernel, arg_index, arg_size, arg_value);
	if (err != CL_SUCCESS) {
		cerr << "clSetKernelArg Failed: " << err << endl;
		exit(1);
	}
}

void oclEnqueueNDRangeKernel(cl_command_queue command_queue,
	cl_kernel kernel,
	cl_uint work_dim,
	const size_t *global_work_offset,
	const size_t *global_work_size,
	const size_t *local_work_size,
	cl_uint num_events_in_wait_list,
	const cl_event *event_wait_list,
	cl_event *event)
{
	cl_int err = clEnqueueNDRangeKernel(command_queue, kernel, work_dim, global_work_offset, global_work_size, local_work_size, num_events_in_wait_list, event_wait_list, event);
	if (err != CL_SUCCESS) {
		cerr << "clEnqueueNDRangeKernel Failed: " << err << endl;
		exit(1);
	}
}

void oclEnqueueReadBuffer(cl_command_queue command_queue,
	cl_mem buffer,
	cl_bool blocking_read,
	size_t offset,
	size_t cb,
	void *ptr,
	cl_uint num_events_in_wait_list,
	const cl_event *event_wait_list,
	cl_event *event)
{
	cl_int err = clEnqueueReadBuffer(command_queue, buffer, blocking_read, offset, cb, ptr, num_events_in_wait_list, event_wait_list, event);
	if (err != CL_SUCCESS) {
		cerr << "clEnqueueReadBuffer Failed: " << err << endl;
		exit(1);
	}
}

void clearbufvec(vector<cl_mem> *bufvec)
{
	unsigned int i;
	cl_int err;
	for (i = 0; i < bufvec->size(); i++) {
		if ((*bufvec)[i] != 0) {
			err = clReleaseMemObject((*bufvec)[i]);
			if (err != CL_SUCCESS) {
				cerr << "clReleaseMemObject Failed: " << err << endl;
				exit(1);
			}
		}
	}
	bufvec->clear();
}
