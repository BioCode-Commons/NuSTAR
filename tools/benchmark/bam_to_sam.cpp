// Benchmark adapter for STAR's textual SAM reader. NuSTAR does not invoke this
// program: its product reader consumes BAM directly through HTSlib.
#include "../../third_party/htslib/htslib/sam.h"
#include <cstdio>
int main(int argc, char** argv) {
    if (argc != 2) return 2;
    auto* input=sam_open(argv[1],"r");
    if(!input)return 1;
    auto* header=sam_hdr_read(input);
    auto* output=sam_open("-","w");
    auto* record=bam_init1();
    int result=-2;
    bool ok=header && output && record && sam_hdr_write(output,header)>=0;
    if(ok) {
        while((result=sam_read1(input,header,record))>=0)
            if(sam_write1(output,header,record)<0){ok=false;break;}
    }
    bam_destroy1(record);
    bam_hdr_destroy(header);
    if(output && sam_close(output)!=0)ok=false;
    if(sam_close(input)!=0)ok=false;
    return ok && result==-1 ? 0:1;
}
