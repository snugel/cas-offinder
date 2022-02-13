#pragma once
#define CL_TARGET_OPENCL_VERSION 120
#define CL_USE_DEPRECATED_OPENCL_1_2_APIS

#include <CL/cl.h>

cl_mem oclCreateBuffer(cl_context context,
                       cl_mem_flags flags,
                       size_t size,
                       void* host_ptr);

void oclGetPlatformIDs(cl_uint num_entries,
                       cl_platform_id* platforms,
                       cl_uint* num_platforms);

void oclGetDeviceIDs(cl_platform_id platform,
                     cl_device_type device_type,
                     cl_uint num_entries,
                     cl_device_id* devices,
                     cl_uint* num_devices);

cl_context oclCreateContext(
  const cl_context_properties* properties,
  cl_uint num_devices,
  const cl_device_id* devices,
  void(CL_CALLBACK* pfn_notify)(const char* errinfo,
                                const void* private_info,
                                size_t cb,
                                void* user_data),
  void* user_data);

cl_program oclCreateProgramWithSource(cl_context context,
                                      cl_uint count,
                                      const char** strings,
                                      const size_t* lengths);

void oclBuildProgram(cl_program program,
                     cl_uint num_devices,
                     const cl_device_id* device_list,
                     const char* options,
                     void(CL_CALLBACK* pfn_notify)(cl_program, void* user_data),
                     void* user_data);

cl_kernel oclCreateKernel(cl_program program, const char* kernel_name);

cl_command_queue oclCreateCommandQueue(cl_context context,
                                       cl_device_id device,
                                       cl_command_queue_properties properties);

void oclGetPlatformInfo(cl_platform_id platform,
                        cl_platform_info param_name,
                        size_t param_value_size,
                        void* param_value,
                        size_t* param_value_size_ret);

void oclGetDeviceInfo(cl_device_id device,
                      cl_device_info param_name,
                      size_t param_value_size,
                      void* param_value,
                      size_t* param_value_size_ret);

void oclReleaseKernel(cl_kernel kernel);

void oclReleaseProgram(cl_program program);

void oclReleaseCommandQueue(cl_command_queue command_queue);

void oclReleaseContext(cl_context context);

void oclEnqueueWriteBuffer(cl_command_queue command_queue,
                           cl_mem buffer,
                           cl_bool blocking_write,
                           size_t offset,
                           size_t cb,
                           const void* ptr,
                           cl_uint num_events_in_wait_list,
                           const cl_event* event_wait_list,
                           cl_event* event);

void oclFinish(cl_command_queue command_queue);

void oclSetKernelArg(cl_kernel kernel,
                     cl_uint arg_index,
                     size_t arg_size,
                     const void* arg_value);

void oclEnqueueNDRangeKernel(cl_command_queue command_queue,
                             cl_kernel kernel,
                             cl_uint work_dim,
                             const size_t* global_work_offset,
                             const size_t* global_work_size,
                             const size_t* local_work_size,
                             cl_uint num_events_in_wait_list,
                             const cl_event* event_wait_list,
                             cl_event* event);

void oclEnqueueReadBuffer(cl_command_queue command_queue,
                          cl_mem buffer,
                          cl_bool blocking_read,
                          size_t offset,
                          size_t cb,
                          void* ptr,
                          cl_uint num_events_in_wait_list,
                          const cl_event* event_wait_list,
                          cl_event* event);

void oclReleaseMemObject(cl_mem buf);
