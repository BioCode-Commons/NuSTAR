#include "../../bam/bgzf_writer.hpp"
#include <fstream>
#include <filesystem>
#include <iostream>
#include <random>
#include <zlib.h>
using namespace nustar::star2711b::bam;
void require(bool b,const char* message) { if(!b) throw std::runtime_error(message); }
ssize_t shortWrite(int fd,const iovec* v,int) {
    static thread_local unsigned count=0;
    if(!(count++%5)) { errno=EINTR; return -1; }
    iovec first=*v; first.iov_len=std::min<std::size_t>(first.iov_len,19);
    return ::writev(fd,&first,1);
}
ssize_t shortPwrite(int fd,const iovec* v,int,off_t offset) {
    static thread_local unsigned count=0;
    if(!(count++%5)) { errno=EINTR; return -1; }
    iovec first=*v; first.iov_len=std::min<std::size_t>(first.iov_len,19);
    return ::pwritev(fd,&first,1,offset);
}
ssize_t shortRead(int fd,const iovec* v,int,off_t offset) {
    static thread_local unsigned count=0;
    if(!(count++%5)) { errno=EINTR; return -1; }
    iovec first=*v; first.iov_len=std::min<std::size_t>(first.iov_len,17);
    return ::preadv(fd,&first,1,offset);
}
ssize_t badWrite(int,const iovec*,int) { errno=ENOSPC; return -1; }
int badSync(int) { errno=EIO; return -1; }
void verify(const std::filesystem::path& path,std::span<const std::byte> expected) {
    std::ifstream in(path,std::ios::binary);
    std::vector<unsigned char> file((std::istreambuf_iterator<char>(in)),{});
    const unsigned char eof[28]={31,139,8,4,0,0,0,0,0,255,6,0,66,67,2,0,27,0,3,0,0,0,0,0,0,0,0,0};
    require(file.size()>=28 && std::equal(std::begin(eof),std::end(eof),file.end()-28),"missing EOF");
    std::size_t pos=0,total=0,eofs=0;
    while(pos<file.size()) {
        require(file.size()-pos>=28,"truncated block");
        auto* b=file.data()+pos;
        require(b[0]==31 && b[1]==139 && b[3]==4 && b[12]==66 && b[13]==67,"bad BGZF header");
        const std::size_t n=1+unsigned(b[16])+(unsigned(b[17])<<8);
        require(n>=28 && n<=65536 && n<=file.size()-pos,"bad block length");
        std::array<unsigned char,65536> decoded;
        z_stream stream{};
        require(inflateInit2(&stream,31)==Z_OK,"inflate init");
        stream.next_in=b; stream.avail_in=n; stream.next_out=decoded.data(); stream.avail_out=decoded.size();
        const auto code=inflate(&stream,Z_FINISH); const auto got=stream.total_out;
        const auto consumed=stream.total_in; inflateEnd(&stream);
        require(code==Z_STREAM_END && consumed==n,"gzip CRC or size failure");
        require(total+got<=expected.size(),"extra bytes");
        require(!std::memcmp(decoded.data(),expected.data()+total,got),"output differs");
        if(!got) ++eofs;
        total+=got; pos+=n;
    }
    require(total==expected.size() && eofs==1,"lost bytes or duplicate EOF");
}
int main(int argc,char** argv) {
    try {
        require(argc==2,"scratch directory required");
        const std::filesystem::path root=argv[1]; std::filesystem::create_directories(root);
        std::vector<std::byte> input((2ULL<<20)+73);
        std::mt19937 rng(109);
        for(auto& c:input) c=std::byte(rng());
        for(unsigned threads:{1U,4U,32U}) for(int level:{-1,0,1,9}) {
            for(std::size_t n:{std::size_t(0),std::size_t(1),std::size_t(0xff00),std::size_t(0xff01),input.size()}) {
                const auto path=root/"blocks.bgzf";
                WorkQueue pool(threads); std::atomic<bool> stopped{false};
                BgzfWriter writer(path.string(),level,8ULL<<20,pool,stopped);
                writer.append(std::span(input).first(n)); writer.finish(); verify(path,std::span(input).first(n));
            }
        }
        {
            WorkQueue pool(4); std::atomic<bool> stopped{false};
            IoCalls calls; calls.write=shortWrite;
            BgzfWriter writer((root/"short.bgzf").string(),1,8ULL<<20,pool,stopped,calls);
            writer.append(std::span(input).first(10000)); writer.finish(); verify(root/"short.bgzf",std::span(input).first(10000));
        }
        {
            const int fd=::open((root/"vectors").c_str(),O_RDWR|O_CREAT|O_TRUNC,0600);
            require(fd>=0,"open vector scratch");
            IoCalls calls; calls.read=shortRead; calls.pwrite=shortPwrite; IoCounters stats;
            vectorIo(fd,{{input.data(),31},{input.data()+31,6000},{input.data()+6031,96}},0,true,true,stats,calls);
            std::vector<std::byte> actual(6127);
            vectorIo(fd,{{actual.data(),5},{actual.data()+5,4000},{actual.data()+4005,2122}},0,false,true,stats,calls);
            require(std::equal(actual.begin(),actual.end(),input.begin()),"short vector transfer differs");
            bool failed=false;
            try { vectorIo(fd,{{actual.data(),actual.size()}},6127,false,true,stats,calls); }
            catch(const std::exception&) { failed=true; }
            ::close(fd); require(failed,"truncated read accepted");
        }
        for(bool sync:{false,true}) {
            WorkQueue pool(4); std::atomic<bool> stopped{false}; IoCalls calls;
            if(sync) calls.sync=badSync; else calls.write=badWrite;
            bool failed=false;
            try {
                BgzfWriter writer((root/"failure.bgzf").string(),1,8ULL<<20,pool,stopped,calls);
                writer.append(input); writer.finish();
            } catch(const std::system_error& e) {
                failed=true; require(!sync && e.code().value()==ENOSPC,"write errno lost");
            } catch(const std::exception&) { failed=true; require(sync,"write errno missing"); }
            require(failed,"output failure accepted");
        }
        {
            WorkQueue pool(4); std::atomic<bool> stopped{false};
            BgzfWriter writer((root/"ordered.bgzf").string(),1,64ULL<<20,pool,stopped);
            auto owner=std::make_shared<std::vector<std::byte>>(input);
            for(std::size_t offset=0;offset<input.size();offset+=100000) {
                const auto n=std::min<std::size_t>(100000,input.size()-offset);
                writer.appendChunk(owner,0,offset,n,[offset](const auto& bytes,std::size_t) {
                    if(!offset) std::this_thread::sleep_for(std::chrono::milliseconds(20));
                    return std::span<const std::byte>(bytes);
                });
            }
            writer.finish(); verify(root/"ordered.bgzf",input);
            stopped=true; bool failed=false;
            try { writer.append(input); } catch(const std::exception&) { failed=true; }
            require(failed,"cancel accepted");
        }
        for(bool cancel:{false,true}) {
            WorkQueue pool(4); std::atomic<bool> stopped{false};
            std::promise<void> release; auto gate=release.get_future().share();
            auto owner=std::make_shared<std::vector<std::byte>>(input);
            bool failed=false;
            try {
                BgzfWriter writer((root/"worker-failure.bgzf").string(),1,64ULL<<20,pool,stopped);
                writer.appendChunk(owner,0,0,100000,[gate,cancel](const auto& bytes,std::size_t) {
                    gate.wait();
                    if(!cancel) throw std::runtime_error("injected compression task failure");
                    return std::span<const std::byte>(bytes);
                });
                writer.appendChunk(owner,0,100000,100000,[](const auto& bytes,std::size_t) {
                    return std::span<const std::byte>(bytes);
                });
                stopped=cancel; release.set_value(); writer.finish();
            } catch(const std::exception&) { failed=true; }
            require(failed,"worker failure accepted");
        }
        std::filesystem::remove_all(root);
        std::cout<<"PASS BGZF levels/threads/long input/CRC/EOF, ordered completion, short I/O, EINTR, failures\n";
    } catch(const std::exception& e) { std::cerr<<e.what()<<'\n';return 1; }
}
