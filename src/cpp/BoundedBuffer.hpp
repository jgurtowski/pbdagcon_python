
///
/// Thread-safe buffer container, lifted from boost.
///
template <class T>
class BoundedBuffer {
public:
    typedef boost::circular_buffer<T> container_type;
    typedef typename container_type::size_type size_type;
    typedef typename container_type::value_type value_type;
    typedef typename boost::call_traits<value_type>::param_type param_type;

    explicit BoundedBuffer(size_type capacity) : unread_(0), container_(capacity) {}

    void push(param_type item) {
        boost::mutex::scoped_lock lock(mutex_);
        not_full_.wait(lock, boost::bind(&BoundedBuffer<value_type>::is_not_full, this));
        container_.push_front(item);
        ++unread_;
        lock.unlock();
        not_empty_.notify_one();
    }

    void pop(value_type* pItem) {
        boost::mutex::scoped_lock lock(mutex_);
        not_empty_.wait(lock, boost::bind(&BoundedBuffer<value_type>::is_not_empty, this));
        *pItem = container_[--unread_];
        lock.unlock();
        not_full_.notify_one();
    }

private:
    BoundedBuffer(const BoundedBuffer&);              // Disabled copy constructor
    BoundedBuffer& operator = (const BoundedBuffer&); // Disabled assign operator

    bool is_not_empty() const { return unread_ > 0; }
    bool is_not_full() const { return unread_ < container_.capacity(); }

    size_type unread_;
    container_type container_;
    boost::mutex mutex_;
    boost::condition not_empty_;
    boost::condition not_full_;
};
