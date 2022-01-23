#include <deque>
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
    bool is_done = false;
    std::mutex lock;
    std::condition_variable cv;
public:
    void send(ItemTy item){
        std::unique_lock<std::mutex> lck(lock);
        queue.push_back(item);
        cv.notify_one();
    }
    void terminate(){
        std::unique_lock<std::mutex> lck(lock);
        is_done = true;
        cv.notify_all();
    }
    bool receive(ItemTy & item){
        std::unique_lock<std::mutex> lck(lock);
        cv.wait(lck);
        
        bool ret_val = is_done;
        if(!ret_val){
            item = queue.pop_front();
        }
        
        return ret_val;
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