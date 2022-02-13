#pragma once

#include <deque>
#include <memory>
#include <mutex>
#include <vector>
#include <condition_variable>

template<class ItemTy>
class Channel{
    /*
        A channel in a data processing pipeline that sends a bit 
        of information from one worker to another.
    */
protected:
    std::deque<ItemTy> queue;
    int count_until_done;
    size_t max_size;
    std::mutex lock;
    std::condition_variable can_rec;
    std::condition_variable can_send;
public:
    static constexpr size_t NO_MAX = ~size_t(0);
    Channel(size_t in_max_size=NO_MAX, int num_terminates=1){
        max_size = in_max_size;
        count_until_done=num_terminates;
    }
    void send(ItemTy item){
        std::unique_lock<std::mutex> lck(lock);
        while(queue.size() >= max_size){
            can_send.wait(lck);
        }
        queue.push_back(item);
        can_rec.notify_one();
    }
    void terminate(){
        std::unique_lock<std::mutex> lck(lock);
        count_until_done--;
        can_rec.notify_all();
    }
    bool receive(ItemTy & item){
        std::unique_lock<std::mutex> lck(lock);
        while(count_until_done && !queue.size()){
            can_rec.wait(lck);
        }
        if(!queue.size()){
            return false;
        }
        else{
            item = queue.front();
            queue.pop_front();
            if(max_size != NO_MAX){
                can_send.notify_one();
            }
            return true;
        }
    }
};
template<class ItemTy>
inline std::vector<ItemTy> get_stream(Channel<ItemTy> & channel){
    std::vector<ItemTy> result;
    ItemTy item;
    while(channel.receive(item)){
        result.push_back(item);
    }
    return result;
}
