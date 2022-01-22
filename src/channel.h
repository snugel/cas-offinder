#include <deque>
#include <mutex>

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
        lock.lock();
        queue.push_back(item);
        lock.unlock();
        cv.notify_one();
    }
    void terminate(){
        lock.lock();
        is_done = true;
        lock.unlock();
        cv.notify_all();
    }
    bool receive(ItemTy & item){
        cv.wait();
        
        lock.lock();
        bool ret_val = is_done;
        if(!ret_val){
            item = queue.pop_front();
        }
        lock.unlock();
        
        return ret_val;
    }
};