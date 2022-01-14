#include <chrono>

template<typename fn_ty>
double time_spent(fn_ty fn){
    auto start = std::chrono::system_clock::now();
    fn();
    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> diff = end-start;
    return diff.count();
}