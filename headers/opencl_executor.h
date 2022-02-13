#include <cassert>
#include <string>
#include <memory>
#include <vector>
#include "oclfunctions.h"

#define CL_TARGET_OPENCL_VERSION 120

//#include <cl.hpp>
#include <CL/cl.h>
//#include <CL/cl_platform.h>

class CLKernelArg
{
  public:
    CLKernelArg(size_t in_size, char* in_ptr)
      : size(in_size)
      , ptr(in_ptr)
    {}
    template<class baseval>
    CLKernelArg(baseval val):
        size(sizeof(val)),
        ptr(new char[sizeof(val)])
    {
        *reinterpret_cast<baseval*>(ptr.get()) = val;
    }
    size_t size;
    std::shared_ptr<char[]> ptr;
};
template<typename item_ty>
class CLBuffer
{
  protected:
    size_t bufsize;
    cl_mem buf;
    cl_context mycontext;
    cl_command_queue myqueue;

  public:
    CLBuffer(cl_context context, cl_command_queue queue, size_t size)
    {
        bufsize = size;
        mycontext = context;
        myqueue = queue;

        buf = oclCreateBuffer(context, CL_MEM_READ_WRITE, bytes(), nullptr);
    }
    CLBuffer(){}

    void write_buffer(std::vector<item_ty>& data)
    {
        assert(data.size() == bufsize);
        write_buffer(&data[0], data.size());
    }
    void write_buffer(item_ty* ptr, size_t size)
    {
        assert(size <= bufsize);
        oclEnqueueWriteBuffer(myqueue,
                                        buf,
                                        CL_TRUE,
                                        0,
                                        size * sizeof(item_ty),
                                        ptr,
                                        0,
                                        nullptr,
                                        nullptr);
        oclFinish(myqueue);
    }
    void read_buffer(std::vector<item_ty>& read_into)
    {
        assert(read_into.size() == bufsize);
        read_buffer(&read_into[0], read_into.size());
    }
    void read_buffer(item_ty* ptr, size_t size)
    {
        assert(size <= bufsize);
        clFinish(myqueue);
        clEnqueueReadBuffer(myqueue,
                                       buf,
                                       CL_TRUE,
                                       0,
                                       size * sizeof(item_ty),
                                       ptr,
                                       0,
                                       nullptr,
                                       nullptr);
        clFinish(myqueue);
    }
    CLKernelArg k_arg() { return CLKernelArg(buf); }
    size_t bytes() { return bufsize * sizeof(item_ty); }
    void copy_buffer(CLBuffer<item_ty> src_buf)
    {
        assert(src_buf.bufsize == this->bufsize);
        assert(src_buf.myqueue == this->myqueue);
        assert(src_buf.mycontext == this->mycontext);
        clFinish(myqueue);
        clEnqueueCopyBuffer(
          myqueue, src_buf.buf, this->buf, 0, 0, bytes(), 0, nullptr, nullptr);
    }
    void clear_buffer()
    {
        int zero = 0;
        clEnqueueFillBuffer(myqueue,
                                       buf,
                                       &zero,
                                       1,
                                       0,
                                       bufsize * sizeof(item_ty),
                                       0,
                                       nullptr,
                                       nullptr);
        clFinish(myqueue);
    }
};
class CL_NDRange
{
  public:
    size_t x;
    size_t y;
    size_t z;
    cl_uint ndim;
    CL_NDRange(size_t in_x, size_t in_y, size_t in_z)
    {
        x = in_x;
        y = in_y;
        z = in_z;
        ndim = 3;
         assert(in_x != 0);
         assert(in_y != 0);
         assert(in_z != 0);
    }
    CL_NDRange(size_t in_x, size_t in_y)
    {
        x = in_x;
        y = in_y;
        z = -1;
        ndim = 2;
         assert(in_x != 0);
         assert(in_y != 0);
    }
    CL_NDRange(size_t in_x)
    {
        x = in_x;
        y = -1;
        z = -1;
        ndim = 1;
         assert(in_x != 0);
    }
    CL_NDRange()
    {
        x = -1;
        y = -1;
        z = -1;
        ndim = 0;
    }
    size_t* array_view()
    {
        return ndim == 0 ? nullptr : reinterpret_cast<size_t*>(this);
    }
    cl_uint dim() { return ndim; }
};
inline CL_NDRange div_nd(CL_NDRange range, CL_NDRange divisor)
{
    size_t* arr = range.array_view();
    size_t* div = divisor.array_view();
    for (int i = 0; i < 3; i++) {
        arr[i] /= div[i];
    }
    return range;
}
class CLKernel
{
  protected:
    cl_command_queue myqueue;
    cl_program program;
    cl_kernel kern;

  public:
    CLKernel(cl_program in_prog,
             cl_command_queue in_queue,
             const char* kern_name,
             std::vector<CLKernelArg> args)
    {
        myqueue = in_queue;
        program = in_prog;

        // http://www.khronos.org/registry/cl/sdk/1.1/docs/man/xhtml/clCreateKernel.html
        kern = oclCreateKernel(program, kern_name);
        set_args(args);
    }
    void set_args(std::vector<CLKernelArg> args){
        int idx = 0;
        for (CLKernelArg& b_info : args) {
            oclSetKernelArg(kern, idx, b_info.size, b_info.ptr.get());
            idx++;
        }
    }
    CLKernel(){}
    void run(CL_NDRange run_range,
             CL_NDRange group_range,
             CL_NDRange exec_range)
    {
        assert(run_range.dim() > 0 &&
               "run_range needs to have at least 1 dimention specified");

        oclFinish(myqueue);
        if (exec_range.dim() == 0) {
            exec_range = CL_NDRange(1, 1, 1);
        }
        CL_NDRange glob_range = div_nd(run_range, exec_range);
        for (size_t x = 0; x < exec_range.x; x++) {
            for (size_t y = 0; y < exec_range.y; y++) {
                for (size_t z = 0; z < exec_range.z; z++) {
                    CL_NDRange global_offset(
                      x * glob_range.x, y * glob_range.y, z * glob_range.z);

                      oclEnqueueNDRangeKernel(myqueue,
                                             kern,
                                             run_range.dim(),
                                             global_offset.array_view(),
                                             glob_range.array_view(),
                                             group_range.array_view(),
                                             0,
                                             nullptr,
                                             nullptr);
                }
            }
        }
        oclFinish(myqueue);
    }
};

class OpenCLExecutor
{
  protected:
    std::string source;
    cl_platform_id platform;
    cl_device_id device;
    cl_context context;
    cl_program program;
    cl_command_queue queue;
    std::string arguments;

  public:
    OpenCLExecutor(std::string in_source,
                   cl_platform_id platid,
                   cl_device_id devid,
                   std::string in_arguments = "")
    {
        device = devid;
        platform = platid;
        source = in_source;
        arguments = in_arguments;
        build_program();
    }
    OpenCLExecutor(){}
    ~OpenCLExecutor()
    {
        oclReleaseProgram(program);
        oclReleaseContext(context);
        oclReleaseCommandQueue(queue);
    }
    template<typename item_ty>
    CLBuffer<item_ty> new_clbuffer(size_t size)
    {
        return CLBuffer<item_ty>(context, queue, size);
    }
    CLKernel new_clkernel(const char* kern_name,
                          std::vector<CLKernelArg> buflist={})
    {
        return CLKernel(program, queue, kern_name, buflist);
    }

  protected:
    void create_context()
    {
        // http://www.khronos.org/registry/cl/sdk/1.1/docs/man/xhtml/clCreateContext.html
        const cl_context_properties contextProperties[] = {
            CL_CONTEXT_PLATFORM,
            reinterpret_cast<cl_context_properties>(this->platform),
            0,
            0
        };

        this->context = oclCreateContext(
          contextProperties, 1, &device, nullptr, nullptr);
    }
    void create_queue()
    {
        // http://www.khronos.org/registry/cl/sdk/1.1/docs/man/xhtml/clCreateCommandQueue.html

        this->queue = oclCreateCommandQueue(context, this->device, 0);
    }
    void build_program()
    {
        create_context();
        create_queue();
        CreateProgram();
    }

    void CreateProgram()
    {
        // http://www.khronos.org/registry/cl/sdk/1.1/docs/man/xhtml/clCreateProgramWithSource.html
        size_t lengths[1] = { source.size() };
        const char* sources[1] = { source.data() };

        this->program =
          oclCreateProgramWithSource(this->context, 1, sources, lengths);

        oclBuildProgram(
          program, 1, &this->device, arguments.c_str(), nullptr, nullptr);
    }
};

struct ExecutorTemplate{
    cl_platform_id plat;
    cl_device_id device;
};

inline std::vector<ExecutorTemplate> get_executor_templates(cl_device_type device_ty=CL_DEVICE_TYPE_DEFAULT){
    std::vector<ExecutorTemplate> templates;

    constexpr size_t MAX_PLATFORM_NUM = 50;
    constexpr size_t MAX_DEVICE_NUM = 1000;

    cl_platform_id platforms[MAX_PLATFORM_NUM];
    cl_uint platform_cnt = 0;
    oclGetPlatformIDs(MAX_PLATFORM_NUM, platforms, &platform_cnt);

    cl_device_id devices[MAX_DEVICE_NUM];
    cl_uint device_cnt;
    for (size_t i = 0; i < platform_cnt; i++) {
        oclGetDeviceIDs(platforms[i], device_ty, MAX_DEVICE_NUM, devices, &device_cnt);
        for(size_t j = 0; j < device_cnt; j++){
            templates.push_back(ExecutorTemplate{
                                    .plat=platforms[i],
                                    .device=devices[j],
                                });
        }
    }
    return templates;
}
inline std::string get_template_info(ExecutorTemplate temp){
    std::string full_result;

    size_t size = 0;
    oclGetPlatformInfo(temp.plat, CL_PLATFORM_NAME, 0, nullptr, &size);
    std::string result;
    result.resize(size);
    oclGetPlatformInfo(temp.plat,
                      CL_PLATFORM_NAME,
                      size,
                      const_cast<char*>(result.data()),
                      nullptr);
    full_result = "Platform: '" + result + "',  ";

    clGetDeviceInfo(temp.device, CL_DEVICE_NAME, 0, nullptr, &size);

    result.resize(size);
    clGetDeviceInfo(
      temp.device, CL_DEVICE_NAME, size, const_cast<char*>(result.data()), nullptr);

    full_result += "Device: " + result;
    return full_result;
}

