/*
Desynchronizes reader and writer, allowing arbitrary amounts of 
data to be written to the stream at once without the reader throttling 
the writer's capacity.
*/
#include <thread>
#include <iostream>
#include <vector>
#include <mutex>

std::mutex buffer_lock;
std::string buffer;

void read_thread(){
    std::string line;
    while(std::getline(std::cin, line)){
        line.push_back('\n'); //add back endline
        buffer_lock.lock();
        buffer += line;
        buffer_lock.unlock();
    }
}
void write_data(){
    buffer_lock.lock();
    std::string local_copy = buffer;
    buffer_lock.unlock();
    std::cout << local_copy;
}

int main(){
    std::thread thread(read_thread);
    
    while(!thread.joinable()){       
        std::this_thread::sleep_for(std::chrono::milliseconds(20));
        write_data();
    }    
    thread.join();
    write_data();
    return 0;
}
