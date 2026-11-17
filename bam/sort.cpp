#include "sort.hpp"
#include "bgzf_writer.hpp"
#include <algorithm>
#include <atomic>
#include <deque>
#include <linux/falloc.h>
#include <numeric>

namespace nustar::star2711b::bam {
namespace {
using Clock = std::chrono::steady_clock;
double seconds(Clock::time_point start) { return std::chrono::duration<double>(Clock::now()-start).count(); }
void check(const std::atomic<bool>& stopped) {
    if (stopped.load(std::memory_order_relaxed)) throw std::runtime_error("BAM sorting cancelled");
}
void addMetrics(SortMetrics& a, const SortMetrics& b) {
    a.temporary_reads.add(b.temporary_reads); a.temporary_writes.add(b.temporary_writes);
    a.read_seconds += b.read_seconds; a.parse_seconds += b.parse_seconds; a.write_seconds += b.write_seconds;
    a.radix_seconds += b.radix_seconds; a.pack_seconds += b.pack_seconds;
}
struct Key {
    std::uint64_t coordinate, read, order;
    auto operator<=>(const Key&) const = default;
    unsigned digit(unsigned byte) const {
        const auto word = byte < 8 ? order : byte < 16 ? read : coordinate;
        return (word >> (8*(byte%8))) & 255;
    }
};
Key key(const std::byte* p, std::uint64_t read, std::uint64_t order) {
    return {(std::uint64_t(read32(p+4))<<32)|read32(p+8),read,order};
}
std::size_t recordSize(std::span<const std::byte> bytes) {
    if (bytes.size()<4) throw std::runtime_error("truncated BAM sort record length");
    const std::uint64_t n=std::uint64_t(read32(bytes.data()))+4;
    if (n<36 || n>bytes.size()) throw std::runtime_error("invalid BAM sort record length");
    return n;
}
struct Entry {
    Key key;
    const std::byte* data;
    std::size_t length;
    Entry() noexcept {}
    Entry(Key k,const std::byte* p,std::size_t n) : key(k),data(p),length(n) {}
};
struct Index {
    std::uint64_t active, index;
    Index() noexcept {}
    Index(std::uint64_t a,std::uint64_t i) : active(a),index(i) {}
};
static_assert(sizeof(Index)==16);
// Metadata + two index arrays + output range bookkeeping, rounded upward.
constexpr std::uint64_t entryCharge=sizeof(Entry)+2*sizeof(Index)+8;
struct Block {
    std::uint64_t sequence=0;
    Buffer data;
    std::vector<Entry> entries;
    void add(std::span<const std::byte> record,std::uint64_t read,std::uint64_t order) {
        entries.emplace_back(key(record.data(),read,order),reinterpret_cast<const std::byte*>(data.size),record.size());
        data.append(record);
    }
    void seal() {
        for (auto& e:entries) e.data=data.data.get()+reinterpret_cast<std::uintptr_t>(e.data);
    }
    std::uint64_t cost() const { return data.capacity+entryCharge*entries.capacity(); }
};
struct Lease {
    std::shared_ptr<std::atomic<std::uint64_t>> used;
    std::uint64_t cost;
    ~Lease() { used->fetch_sub(cost,std::memory_order_relaxed); }
};
struct OutputChunk { std::size_t record,offset,bytes; };
struct Batch {
    // Released last, after the arenas and metadata have been destroyed.
    std::shared_ptr<Lease> lease;
    std::vector<Block> blocks;
    std::vector<Entry> entries;
    std::vector<Index> order;
    std::vector<OutputChunk> chunks;
    SortMetrics metrics;
};
void keySort(std::span<Index> values,std::span<Index> scratch,const std::atomic<bool>& stopped,WorkQueue* pool=nullptr) {
    if (values.size()<2) return;
    if (values.size()<=32) {
        for (std::size_t i=1;i<values.size();++i) {
            auto value=values[i]; auto j=i;
            while (j && values[j-1].active>value.active) { values[j]=values[j-1]; --j; }
            values[j]=value;
        }
        return;
    }
    if (values.size()<=4096) {
        std::stable_sort(values.begin(),values.end(),[](const Index& a,const Index& b){return a.active<b.active;});
        return;
    }
    const unsigned workers=pool ? pool->concurrency() : 1;
    const auto run=[&](auto task) { if (pool) pool->parallel(workers,task); else task(0); };
    using Histogram=std::array<std::size_t,256>;
    std::vector<Histogram> counts(workers),offsets(workers);
    std::vector<std::uint64_t> differences(workers);
    const auto first=values.front().active;
    run([&](unsigned w) {
        for (auto i=values.size()*w/workers;i<values.size()*(w+1)/workers;++i) differences[w]|=values[i].active^first;
    });
    std::uint64_t varying=0;
    for (auto d:differences) varying|=d;
    auto source=values,destination=scratch;
    for (unsigned shift=0;shift<64;shift+=8) {
        if (!((varying>>shift)&255)) continue;
        check(stopped);
        run([&](unsigned w) {
            counts[w].fill(0);
            for (auto i=source.size()*w/workers;i<source.size()*(w+1)/workers;++i) ++counts[w][(source[i].active>>shift)&255];
        });
        std::size_t next=0;
        for (unsigned digit=0;digit<256;++digit) for (unsigned w=0;w<workers;++w) {
            offsets[w][digit]=next; next+=counts[w][digit];
        }
        run([&](unsigned w) {
            auto positions=offsets[w];
            for (auto i=source.size()*w/workers;i<source.size()*(w+1)/workers;++i)
                destination[positions[(source[i].active>>shift)&255]++]=source[i];
        });
        std::swap(source,destination);
    }
    if (source.data()!=values.data()) run([&](unsigned w) {
        const auto begin=source.size()*w/workers,end=source.size()*(w+1)/workers;
        std::copy(source.begin()+begin,source.begin()+end,values.begin()+begin);
    });
}
void sortBatch(Batch& batch,const std::atomic<bool>& stopped,WorkQueue* pool=nullptr) {
    auto start=Clock::now();
    auto& entries=batch.entries;
    batch.order.resize(entries.size());
    std::vector<Index> scratch(entries.size());
    const unsigned workers=pool ? pool->concurrency() : 1;
    const auto run=[&](auto task) { if(pool) pool->parallel(workers,task); else task(0); };
    run([&](unsigned w) {
        for (auto i=entries.size()*w/workers;i<entries.size()*(w+1)/workers;++i) batch.order[i]={entries[i].key.coordinate,i};
    });
    keySort(batch.order,scratch,stopped,pool);
    // Coordinate groups are disjoint. Small groups are bundled into tasks;
    // only the coordinator invokes team-wide passes for a large group.
    const auto refine=[&](std::size_t first,std::size_t last,WorkQueue* team) {
        for(auto begin=first;begin<last;) {
            check(stopped);
            auto end=begin+1;
            while(end<last && batch.order[end].active==batch.order[begin].active) ++end;
            if(end-begin>1) {
                auto group=std::span(batch.order).subspan(begin,end-begin);
                auto work=std::span(scratch).subspan(begin,end-begin);
                for(auto& index:group) index.active=entries[index.index].key.read;
                keySort(group,work,stopped,team && group.size()>=65536 ? team : nullptr);
                for(std::size_t b=0;b<group.size();) {
                    auto e=b+1;
                    while(e<group.size() && group[e].active==group[b].active) ++e;
                    if(e-b>1) {
                        for(auto i=b;i<e;++i) group[i].active=entries[group[i].index].key.order;
                        keySort(group.subspan(b,e-b),work.subspan(b,e-b),stopped,team && e-b>=65536 ? team : nullptr);
                    }
                    b=e;
                }
            }
            begin=end;
        }
    };
    if(!pool) refine(0,batch.order.size(),nullptr);
    else {
        std::vector<std::future<void>> pending; pending.reserve(2*workers);
        const auto drain=[&] {
            std::exception_ptr failure;
            for(auto& f:pending) try { f.get(); } catch(...) { if(!failure) failure=std::current_exception(); }
            pending.clear(); if(failure) std::rethrow_exception(failure);
        };
        const auto dispatch=[&](std::size_t begin,std::size_t end) {
            if(begin==end) return;
            if(pending.size()==2*workers) drain();
            pending.push_back(pool->submit([&,begin,end]{refine(begin,end,nullptr);}));
        };
        try {
            std::size_t chunk=0;
            for(std::size_t begin=0;begin<batch.order.size();) {
                auto end=begin+1;
                while(end<batch.order.size() && batch.order[end].active==batch.order[begin].active) ++end;
                if(end-begin>=65536) {
                    dispatch(chunk,begin); drain(); refine(begin,end,pool); chunk=end;
                } else if(end-chunk>=65536) { dispatch(chunk,end); chunk=end; }
                begin=end;
            }
            dispatch(chunk,batch.order.size()); drain();
        } catch(...) { for(auto& f:pending) f.wait(); throw; }
    }
    batch.metrics.radix_seconds+=seconds(start); start=Clock::now();
    // Only small range descriptors are materialized, not an ordered payload.
    std::size_t begin=0,offset=0,bytes=0;
    for (std::size_t i=0;i<batch.order.size();++i) {
        if (!(i&4095)) check(stopped);
        std::size_t at=0,n=entries[batch.order[i].index].length;
        while (at<n) {
            const auto take=std::min<std::size_t>((1ULL<<20)-bytes,n-at);
            bytes+=take; at+=take;
            if (bytes==(1ULL<<20)) {
                batch.chunks.push_back({begin,offset,bytes}); bytes=0;
                begin=at==n ? i+1 : i; offset=at==n ? 0 : at;
            }
        }
    }
    if (bytes) batch.chunks.push_back({begin,offset,bytes});
    batch.metrics.pack_seconds+=seconds(start);
}
struct Extent { unsigned file; std::uint64_t offset,bytes,records,sequence; };
struct Bucket {
    std::mutex mutex;
    std::vector<Extent> extents;
    std::uint64_t bytes=0,records=0;
    Key low{},high{};
    std::uint64_t cost() const { return bytes+entryCharge*records; }
};
using Buckets=std::vector<std::unique_ptr<Bucket>>;
Buckets buckets(std::size_t n) {
    Buckets result;
    for(std::size_t i=0;i<n;++i) result.push_back(std::make_unique<Bucket>());
    return result;
}
class Spool {
    struct Packet { Buffer data; Bucket* bucket; std::uint64_t records,sequence; Key low,high; };
    struct File { int fd=-1; std::uint64_t offset=0,queued=0; std::vector<Packet> packets; };
    struct Write { int fd=-1; unsigned file=0; std::uint64_t offset=0,bytes=0; std::vector<Packet> packets; };
    std::array<File,64> files_;
    std::filesystem::path directory_;
    const std::uint64_t cap_,target_;
    std::uint64_t held_=0;
    std::mutex mutex_;
    std::condition_variable available_;
    std::exception_ptr failure_;
    int open(unsigned id) {
        auto& f=files_[id];
        if(f.fd<0) {
            f.fd=::open((directory_/("bucket-"+std::to_string(id))).c_str(),O_RDWR|O_CREAT|O_EXCL|O_CLOEXEC,0600);
            if(f.fd<0) throw std::runtime_error("cannot create BAM sort bucket");
        }
        return f.fd;
    }
    Write take(unsigned id) {
        auto& f=files_[id];
        if(f.queued>std::uint64_t(std::numeric_limits<off_t>::max())-f.offset)
            throw std::runtime_error("BAM sort bucket offset overflow");
        Write w{f.fd,id,f.offset,f.queued,{}};
        f.offset+=f.queued; f.queued=0; w.packets.swap(f.packets); return w;
    }
    void execute(Write w,SortMetrics& stats,bool charged) {
        const auto start=Clock::now();
        try {
            std::vector<iovec> vectors;
            for(auto& p:w.packets) vectors.push_back({p.data.data.get(),p.data.size});
            vectorIo(w.fd,std::move(vectors),w.offset,true,true,stats.temporary_writes);
            auto offset=w.offset;
            for(auto& p:w.packets) {
                auto& b=*p.bucket;
                std::lock_guard lock(b.mutex);
                b.extents.push_back({w.file,offset,p.data.size,p.records,p.sequence});
                if(!b.records) { b.low=p.low; b.high=p.high; }
                else { b.low=std::min(b.low,p.low); b.high=std::max(b.high,p.high); }
                b.bytes+=p.data.size; b.records+=p.records; offset+=p.data.size;
            }
        } catch (...) {
            { std::lock_guard lock(mutex_); if(!failure_) failure_=std::current_exception();
              if(charged) held_-=w.bytes; }
            available_.notify_all(); throw;
        }
        stats.write_seconds+=seconds(start);
        if(charged) { std::lock_guard lock(mutex_); held_-=w.bytes; available_.notify_all(); }
    }
  public:
    explicit Spool(std::filesystem::path directory,std::uint64_t budget)
        : directory_(std::move(directory)),cap_(std::min<std::uint64_t>(64ULL<<20,budget/16)),
          target_(std::min<std::uint64_t>(1ULL<<20,cap_)) {}
    ~Spool() { for(auto& f:files_) if(f.fd>=0) ::close(f.fd); }
    void write(unsigned id,Buffer packet,Bucket& bucket,std::uint64_t records,std::uint64_t sequence,Key low,Key high,SortMetrics& stats) {
        const auto n=packet.size;
        std::unique_lock lock(mutex_);
        if(failure_) std::rethrow_exception(failure_);
        auto& f=files_[id]; open(id);
        if(n>=target_) {
            if(n>std::uint64_t(std::numeric_limits<off_t>::max())-f.offset) throw std::runtime_error("BAM sort bucket overflow");
            Write w{f.fd,id,f.offset,n,{}}; f.offset+=n;
            w.packets.push_back({std::move(packet),&bucket,records,sequence,low,high});
            lock.unlock(); execute(std::move(w),stats,false); return;
        }
        while(held_+n>cap_) {
            auto largest=std::max_element(files_.begin(),files_.end(),[](const File& a,const File& b){return a.queued<b.queued;});
            if(largest->queued) {
                auto w=take(largest-files_.begin()); lock.unlock(); execute(std::move(w),stats,true); lock.lock();
            } else available_.wait(lock,[&]{return failure_ || held_+n<=cap_;});
            if(failure_) std::rethrow_exception(failure_);
        }
        f.packets.push_back({std::move(packet),&bucket,records,sequence,low,high}); f.queued+=n; held_+=n;
        if(f.queued>=target_) { auto w=take(id); lock.unlock(); execute(std::move(w),stats,true); }
    }
    void flush(SortMetrics& stats) {
        for(unsigned id=0;id<files_.size();++id) {
            std::unique_lock lock(mutex_);
            if(failure_) std::rethrow_exception(failure_);
            if(files_[id].queued) { auto w=take(id); lock.unlock(); execute(std::move(w),stats,true); }
        }
        std::unique_lock lock(mutex_);
        available_.wait(lock,[&]{return failure_ || !held_;});
        if(failure_) std::rethrow_exception(failure_);
    }
    void read(unsigned file,std::uint64_t offset,std::vector<iovec> vectors,SortMetrics& stats) {
        const auto start=Clock::now();
        vectorIo(files_[file].fd,std::move(vectors),offset,false,true,stats.temporary_reads);
        stats.read_seconds+=seconds(start);
    }
    void discard(const Bucket& bucket) {
        // Called only after every child extent is written and parent readers
        // have joined. Exact byte ranges preserve adjacent live bucket data.
        auto extents=bucket.extents;
        std::sort(extents.begin(),extents.end(),[](const Extent& a,const Extent& b) {
            return std::tie(a.file,a.offset)<std::tie(b.file,b.offset);
        });
        for(std::size_t i=0;i<extents.size();) {
            const auto& first=extents[i++]; auto bytes=first.bytes;
            while(i<extents.size() && extents[i].file==first.file && extents[i].offset==first.offset+bytes)
                bytes+=extents[i++].bytes;
            int result;
            do { result=::fallocate(files_[first.file].fd,FALLOC_FL_PUNCH_HOLE|FALLOC_FL_KEEP_SIZE,first.offset,bytes); }
            while(result<0 && errno==EINTR);
            // Reclamation is optional on filesystems without hole punching.
            // A later failed write retains its actual errno for diagnosis.
            if(result<0) return;
        }
    }
    void prefetch(const Bucket& b,std::uint64_t budget) {
        auto left=std::min<std::uint64_t>(32ULL<<20,budget/32);
        // Hint at most one range per physical file; hints are not correctness operations.
        std::array<bool,64> seen{};
        for(const auto& e:b.extents) if(left && !seen[e.file]) {
            seen[e.file]=true;
            const auto n=std::min<std::uint64_t>(left,files_[e.file].offset-e.offset);
            (void)::posix_fadvise(files_[e.file].fd,e.offset,n,POSIX_FADV_WILLNEED); left-=n;
        }
    }
};
void unpack(std::span<const std::byte> bytes,std::span<Entry> entries) {
    std::size_t p=0;
    for(auto& entry:entries) {
        if(bytes.size()-p<52) throw std::runtime_error("truncated BAM sort record");
        std::uint64_t read,order;
        std::memcpy(&read,bytes.data()+p,8); std::memcpy(&order,bytes.data()+p+8,8); p+=16;
        const std::uint64_t n=std::uint64_t(read32(bytes.data()+p))+4;
        if(n<36 || n>bytes.size()-p) throw std::runtime_error("invalid BAM sort record length");
        entry={key(bytes.data()+p,read,order),bytes.data()+p,std::size_t(n)}; p+=n;
    }
    if(p!=bytes.size()) throw std::runtime_error("BAM bucket count mismatch");
}
}

struct CoordinateSorter::Impl {
    const std::uint64_t budget,block_target;
    const unsigned threads;
    std::mutex mutex;
    std::atomic<bool> stopped{false};
    bool finished=false,spilling=false;
    std::exception_ptr failure;
    std::uint64_t retained_cost=0,records=0,partitions=0;
    struct Sample { Key key; std::uint64_t weight; };
    std::vector<Sample> sample;
    std::uint64_t sample_stride=1;
    std::deque<Block> retained,flush_queue;
    std::vector<Key> boundaries;
    Buckets initial;
    Spool spool;
    SortMetrics totals;
    Impl(std::filesystem::path directory,std::uint64_t memory,unsigned n)
        : budget(std::max<std::uint64_t>(memory,8ULL<<20)),
          block_target(std::max<std::uint64_t>(4096,std::min<std::uint64_t>(8ULL<<20,budget/(16ULL*std::max(1U,n))))),
          threads(std::max(1U,n)),spool(std::move(directory),budget) {}
    void fail() {
        std::lock_guard lock(mutex); if(!failure) failure=std::current_exception();
        stopped.store(true,std::memory_order_relaxed);
    }
    void mergeMetrics(const SortMetrics& stats) { std::lock_guard lock(mutex); addMetrics(totals,stats); }
    void sampleBlock(const Block& block,std::uint64_t first) {
        const auto maximum=std::max<std::uint64_t>(1,std::min<std::uint64_t>(1ULL<<20,budget/32/sizeof(Sample)));
        auto at=(first+sample_stride-1)/sample_stride*sample_stride;
        const auto end=first+block.entries.size();
        while(at<end) {
            if(sample.size()==maximum) {
                std::size_t out=0;
                for(std::size_t i=0;i<sample.size();i+=2) sample[out++]=sample[i];
                sample.resize(out); sample_stride*=2;
                at=(at+sample_stride-1)/sample_stride*sample_stride;
                if(at>=end) break;
            }
            const auto& e=block.entries[at-first];
            sample.push_back({e.key,e.length+16+entryCharge}); at+=sample_stride;
        }
    }
    void chooseBoundaries() {
        std::sort(sample.begin(),sample.end(),[](const Sample& a,const Sample& b){return a.key<b.key;});
        const auto n=std::clamp(8ULL*threads,16ULL,256ULL);
        std::uint64_t total=0,accumulated=0,next=1;
        for(const auto& s:sample) total+=s.weight;
        for(const auto& s:sample) {
            accumulated+=s.weight;
            if(next<n && accumulated>=total*next/n) {
                if(boundaries.empty() || boundaries.back()<s.key) boundaries.push_back(s.key);
                ++next;
            }
        }
        initial=buckets(boundaries.size()+1);
        std::vector<Sample>().swap(sample);
    }
    template<class Route> void distribute(const Block& block,Buckets& destination,Route route,SortMetrics& stats) {
        check(stopped);
        const auto n=destination.size(),count=block.entries.size();
        std::vector<std::size_t> counts(n),offsets(n+1),cursor;
        auto routes=std::make_unique_for_overwrite<unsigned[]>(count);
        auto indices=std::make_unique_for_overwrite<std::size_t[]>(count);
        for(std::size_t i=0;i<count;++i) { routes[i]=route(block.entries[i].key); ++counts[routes[i]]; }
        std::partial_sum(counts.begin(),counts.end(),offsets.begin()+1); cursor=offsets;
        for(std::size_t i=0;i<count;++i) indices[cursor[routes[i]]++]=i;
        for(unsigned b=0;b<n;++b) {
            if(!counts[b]) continue;
            check(stopped);
            Key low=block.entries[indices[offsets[b]]].key,high=low;
            std::size_t bytes=0;
            for(auto j=offsets[b];j<offsets[b+1];++j) bytes+=16+block.entries[indices[j]].length;
            Buffer packet; packet.resize(bytes);
            std::size_t p=0;
            for(auto j=offsets[b];j<offsets[b+1];++j) {
                const auto& e=block.entries[indices[j]];
                std::memcpy(packet.data.get()+p,&e.key.read,8);
                std::memcpy(packet.data.get()+p+8,&e.key.order,8);
                std::memcpy(packet.data.get()+p+16,e.data,e.length); p+=16+e.length;
                low=std::min(low,e.key); high=std::max(high,e.key);
            }
            spool.write(b/((n+63)/64),std::move(packet),*destination[b],counts[b],block.sequence,low,high,stats);
        }
    }
    void route(const Block& block,SortMetrics& stats) {
        distribute(block,initial,[&](const Key& k){return std::upper_bound(boundaries.begin(),boundaries.end(),k)-boundaries.begin();},stats);
    }
    void accept(Block block) {
        if(block.entries.empty()) return;
        block.seal(); Block old;
        bool stored=false;
        {
            std::lock_guard lock(mutex); check(stopped);
            if(finished) throw std::runtime_error("BAM submit after finish");
            const auto first=records; block.sequence=first; records+=block.entries.size();
            if(!spilling) {
                sampleBlock(block,first);
                retained_cost+=block.cost(); retained.push_back(std::move(block)); stored=true;
                if(retained_cost<=budget/2) return;
                chooseBoundaries(); spilling=true; flush_queue.swap(retained);
            }
            if(!flush_queue.empty()) { old=std::move(flush_queue.front()); flush_queue.pop_front(); }
        }
        SortMetrics stats;
        if(!old.entries.empty()) route(old,stats);
        if(!stored) route(block,stats);
        mergeMetrics(stats);
    }
    void finish() {
        { std::lock_guard lock(mutex);
          if(failure) std::rethrow_exception(failure);
          check(stopped); if(finished) return; finished=true; }
        if(spilling) {
            // Mapping workers have joined. Their unprocessed retained slabs can
            // now use the output pool without oversubscribing mapping threads.
            WorkQueue pool(threads);
            std::vector<SortMetrics> stats(pool.concurrency());
            pool.parallel(pool.concurrency(),[&](unsigned w) {
                for(;;) {
                    Block block;
                    { std::lock_guard lock(mutex);
                      if(flush_queue.empty()) break;
                      block=std::move(flush_queue.front()); flush_queue.pop_front(); }
                    route(block,stats[w]);
                }
            });
            for(const auto& s:stats) addMetrics(totals,s);
            spool.flush(totals);
        }
    }
    Block load(const Bucket& bucket,SortMetrics& stats) {
        Block block; block.data.resize(bucket.bytes); block.entries.resize(bucket.records);
        struct Request { Extent extent; std::size_t bytes,records; };
        std::vector<Request> requests; requests.reserve(bucket.extents.size());
        std::size_t byte_offset=0,record_offset=0;
        for(const auto& e:bucket.extents) {
            if(e.bytes>bucket.bytes-byte_offset || e.records>bucket.records-record_offset)
                throw std::runtime_error("BAM extent exceeds bucket");
            requests.push_back({e,byte_offset,record_offset}); byte_offset+=e.bytes; record_offset+=e.records;
        }
        if(byte_offset!=bucket.bytes || record_offset!=bucket.records) throw std::runtime_error("BAM bucket extent totals disagree");
        std::sort(requests.begin(),requests.end(),[](const Request& a,const Request& b){
            return std::tie(a.extent.file,a.extent.offset)<std::tie(b.extent.file,b.extent.offset);
        });
        for(std::size_t i=0;i<requests.size();) {
            check(stopped);
            const auto& first=requests[i]; auto end=i; std::uint64_t bytes=0;
            std::vector<iovec> vectors;
            do {
                const auto& r=requests[end]; vectors.push_back({block.data.data.get()+r.bytes,std::size_t(r.extent.bytes)});
                bytes+=r.extent.bytes; ++end;
            } while(end<requests.size() && bytes<(4ULL<<20) && requests[end].extent.file==first.extent.file && requests[end].extent.offset==first.extent.offset+bytes);
            spool.read(first.extent.file,first.extent.offset,std::move(vectors),stats);
            const auto start=Clock::now();
            for(auto j=i;j<end;++j) {
                const auto& r=requests[j];
                unpack(block.data.bytes().subspan(r.bytes,r.extent.bytes),std::span(block.entries).subspan(r.records,r.extent.records));
            }
            stats.parse_seconds+=seconds(start); i=end;
        }
        return block;
    }
    template<class Visit> void chunks(const Bucket& bucket,WorkQueue& pool,Visit visit) {
        std::vector<SortMetrics> stats(pool.concurrency());
        pool.parallel(pool.concurrency(),[&](unsigned w) {
            for(auto i=bucket.extents.size()*w/pool.concurrency();i<bucket.extents.size()*(w+1)/pool.concurrency();++i) {
                check(stopped);
                const auto& e=bucket.extents[i];
                Block block; block.sequence=e.sequence; block.data.resize(e.bytes); block.entries.resize(e.records);
                spool.read(e.file,e.offset,{{block.data.data.get(),std::size_t(e.bytes)}},stats[w]);
                const auto start=Clock::now(); unpack(block.data.bytes(),block.entries); stats[w].parse_seconds+=seconds(start);
                visit(block,stats[w]);
            }
        });
        for(const auto& s:stats) addMetrics(totals,s);
    }
    void split(Bucket& bucket,std::uint64_t target,WorkQueue& pool,std::vector<std::shared_ptr<Bucket>>& leaves) {
        check(stopped);
        if(!bucket.records) return;
        std::stable_sort(bucket.extents.begin(),bucket.extents.end(),[](const Extent& a,const Extent& b){return a.sequence<b.sequence;});
        if(bucket.cost()<=target) {
            auto leaf=std::make_shared<Bucket>();
            leaf->extents=std::move(bucket.extents); leaf->bytes=bucket.bytes; leaf->records=bucket.records;
            leaf->low=bucket.low; leaf->high=bucket.high; leaves.push_back(std::move(leaf)); return;
        }
        if(bucket.low==bucket.high) {
            // Preserve the original logical extent order even though physical
            // reads may be reordered later. A single long record is indivisible.
            for(const auto& e:bucket.extents) {
                auto leaf=std::make_shared<Bucket>(); leaf->extents.push_back(e);
                leaf->bytes=e.bytes; leaf->records=e.records; leaf->low=bucket.low; leaf->high=bucket.high;
                leaves.push_back(std::move(leaf));
            }
            return;
        }
        unsigned digit=23;
        while(digit && bucket.low.digit(digit)==bucket.high.digit(digit)) --digit;
        auto children=buckets(256); ++partitions;
        chunks(bucket,pool,[&](const Block& block,SortMetrics& stats) {
            distribute(block,children,[&](const Key& k){return k.digit(digit);},stats);
        });
        spool.flush(totals);
        spool.discard(bucket);
        for(auto& child:children) split(*child,target,pool,leaves);
    }
    template<class Emit,class Progress> void run(WorkQueue& pool,Emit emit,Progress progress) {
        finish();
        if(!spilling) {
            auto batch=std::make_shared<Batch>(); batch->entries.resize(records);
            std::vector<std::size_t> offsets{0};
            while(!retained.empty()) {
                offsets.push_back(offsets.back()+retained.front().entries.size());
                batch->blocks.push_back(std::move(retained.front())); retained.pop_front();
            }
            pool.parallel(pool.concurrency(),[&](unsigned w) {
                for(auto i=batch->blocks.size()*w/pool.concurrency();i<batch->blocks.size()*(w+1)/pool.concurrency();++i) {
                    auto& block=batch->blocks[i];
                    std::copy(block.entries.begin(),block.entries.end(),batch->entries.begin()+offsets[i]);
                    std::vector<Entry>().swap(block.entries);
                }
            });
            sortBatch(*batch,stopped,&pool); addMetrics(totals,batch->metrics); emit(std::move(batch)); return;
        }
        std::vector<std::shared_ptr<Bucket>> leaves;
        const auto target=std::max<std::uint64_t>(block_target,budget/(std::max(4U,2*threads)));
        for(auto& b:initial) split(*b,target,pool,leaves);
        const auto limit=budget*3/4;
        auto used=std::make_shared<std::atomic<std::uint64_t>>(0);
        std::deque<std::future<std::shared_ptr<Batch>>> pending;
        std::size_t next=0;
        try {
            while(next<leaves.size() || !pending.empty()) {
                while(next<leaves.size() && pending.size()<2*pool.concurrency()) {
                    auto bucket=leaves[next]; const auto cost=bucket->cost();
                    const auto current=used->load(std::memory_order_relaxed);
                    if(current && (cost>limit || current>limit-cost)) break;
                    auto lease=std::make_shared<Lease>(); lease->used=used; lease->cost=cost;
                    used->fetch_add(cost,std::memory_order_relaxed);
                    spool.prefetch(*bucket,budget);
                    auto future=pool.submit([this,bucket,lease]() mutable {
                        check(stopped);
                        auto batch=std::make_shared<Batch>(); batch->lease=std::move(lease);
                        auto block=load(*bucket,batch->metrics); batch->entries=std::move(block.entries);
                        batch->blocks.push_back(std::move(block)); sortBatch(*batch,stopped); return batch;
                    });
                    try { pending.push_back(std::move(future)); }
                    catch(...) { if(future.valid()) future.wait(); throw; }
                    ++next;
                }
                if(pending.empty()) {
                    // Compression futures can retain the last admitted arena.
                    // Drain them rather than waiting for credits while holding
                    // completed output, which would deadlock a small budget.
                    if(!progress()) throw std::runtime_error("BAM memory credits cannot make progress");
                    continue;
                }
                const auto start=Clock::now();
                auto future=std::move(pending.front()); pending.pop_front();
                auto batch=future.get();
                totals.wait_seconds+=seconds(start); addMetrics(totals,batch->metrics); emit(std::move(batch));
            }
        } catch(...) {
            stopped.store(true,std::memory_order_relaxed);
            for(auto& future:pending) future.wait();
            throw;
        }
    }
};

CoordinateSorter::CoordinateSorter(std::filesystem::path directory,std::uint64_t budget,unsigned threads)
    : impl_(std::make_unique<Impl>(std::move(directory),budget,threads)) {}
CoordinateSorter::~CoordinateSorter()=default;
void CoordinateSorter::submit(const OutputBatch& batch) {
    auto& s=*impl_;
    try {
        check(s.stopped); Block block; std::uint64_t count=0;
        for(const auto& fragment:batch.fragments) {
            std::size_t p=0; std::uint64_t order=0;
            while(p<fragment.sam.size()) {
                if(!(count++&4095)) check(s.stopped);
                const auto record=std::span(fragment.sam).subspan(p); const auto n=recordSize(record);
                block.add(record.first(n),fragment.fragment.ordinal,order++); p+=n;
                if(block.data.size+entryCharge*block.entries.size()>=s.block_target) {
                    s.accept(std::move(block)); block={};
                }
            }
        }
        s.accept(std::move(block));
    } catch(...) { s.fail(); throw; }
}
void CoordinateSorter::finish() {
    try { impl_->finish(); } catch(...) { impl_->fail(); throw; }
}
void CoordinateSorter::emit(const std::function<void(std::uint64_t,std::span<const std::byte>)>& output,
                            const std::function<void(std::span<const std::byte>)>& blocks) {
    auto& s=*impl_; finish(); WorkQueue pool(s.threads);
    s.run(pool,[&](std::shared_ptr<Batch> batch) {
        const auto start=Clock::now(); std::size_t count=0;
        for(const auto& i:batch->order) {
            if(!(count++&4095)) check(s.stopped);
            const auto& e=batch->entries[i.index];
            if(blocks) blocks({e.data,e.length}); else output(e.key.read,{e.data,e.length});
        }
        s.totals.output_seconds+=seconds(start);
    },[]{return false;});
}
void CoordinateSorter::writeBam(const std::string& path,std::span<const std::byte> header,int compression,
    const std::function<std::pair<std::string,std::string>(std::uint64_t)>& correction) {
    auto& s=*impl_; finish(); WorkQueue pool(s.threads);
    BgzfWriter writer(path,compression,s.budget,pool,s.stopped);
    try {
        writer.append(header);
        s.run(pool,[&](std::shared_ptr<Batch> batch) {
            if(!correction) {
                for(const auto& c:batch->chunks) writer.appendChunk(batch,c.record,c.offset,c.bytes,
                    [](const Batch& b,std::size_t i) {
                        const auto& e=b.entries[b.order[i].index]; return std::span<const std::byte>(e.data,e.length);
                    });
            } else {
                Bytes corrected; Buffer pending;
                for(const auto& i:batch->order) {
                    const auto& e=batch->entries[i.index]; check(s.stopped);
                    corrected.assign(e.data,e.data+e.length);
                    auto [cb,umi]=correction(e.key.read);
                    appendStringTag(corrected,"CB",cb); appendStringTag(corrected,"UB",umi);
                    const auto n=corrected.size()-4;
                    if(n>UINT32_MAX) throw std::runtime_error("corrected BAM record too large");
                    for(unsigned j=0;j<4;++j) corrected[j]=std::byte(n>>(8*j));
                    if(pending.size+corrected.size()>(1ULL<<20)) { writer.append(pending.bytes()); pending.size=0; }
                    if(corrected.size()>=(1ULL<<20)) writer.append(corrected);
                    else pending.append(corrected);
                }
                writer.append(pending.bytes());
            }
        },[&]{return writer.drainOne();});
        writer.finish();
        const auto& m=writer.metrics();
        s.totals.output_writes=m.output; s.totals.pack_seconds+=m.prepare_seconds;
        s.totals.compression_seconds=m.compression_seconds; s.totals.wait_seconds+=m.wait_seconds;
        s.totals.sync_seconds=m.sync_seconds; s.totals.output_seconds=m.write_seconds;
    } catch(...) { s.stopped.store(true,std::memory_order_relaxed); throw; }
}
void CoordinateSorter::cancel() noexcept { impl_->stopped.store(true,std::memory_order_relaxed); }
SortMetrics CoordinateSorter::metrics() const {
    auto result=impl_->totals;
    result.records=impl_->records; result.partitions=impl_->partitions;
    result.spilled_bytes=result.temporary_writes.bytes; result.read_bytes=result.temporary_reads.bytes;
    return result;
}
}
