"""Check input mate reordering and STAR's fragment filter flag source.
Run bam.py first; arguments are NuSTAR, STAR, fixtures, bam-test directory.
"""
from pathlib import Path
import collections,gzip,struct,subprocess,sys,zlib
from bam import decode,sam_text
binary,star,fixtures,work=(Path(x).resolve() for x in sys.argv[1:5])
# Generate BGZF independently, so the fixture does not use the product encoder.
def bgzf(path,data):
    with path.open('wb') as out:
        for i in range(0,len(data),60000):
            block=data[i:i+60000];compressed=zlib.compress(block,wbits=-15)
            out.write(bytes.fromhex('1f8b08040000000000ff060042430200')+
                struct.pack('<H',len(compressed)+25)+compressed+
                struct.pack('<II',zlib.crc32(block),len(block)))
        out.write(bytes.fromhex('1f8b08040000000000ff0600424302001b0003000000000000000000'))
data=gzip.open(work/'native-pe-1.Aligned.out.bam','rb').read()
n,=struct.unpack_from('<i',data,4);offset=8+n
n,=struct.unpack_from('<i',data,offset);offset+=4
for _ in range(n):
    size,=struct.unpack_from('<i',data,offset);offset+=4+size+4
header=data[:offset];records=[]
while offset<len(data):
    size,=struct.unpack_from('<i',data,offset);records.append(bytearray(data[offset:offset+size+4]));offset+=size+4
for i in range(0,len(records),2):
    for record in records[i:i+2]:
        word,=struct.unpack_from('<I',record,16);flag=word>>16
        if bool(flag&128)==bool((i//2)%2):flag|=2048
        else:flag&=~2048
        struct.pack_into('<I',record,16,(flag<<16)|(word&65535))
    records[i],records[i+1]=records[i+1],records[i]
bam=work/'reordered-filter.bam';bgzf(bam,header+b''.join(records))
sam=work/'reordered-filter.sam';sam.write_text(sam_text(bam))
for exe,label,source in ((binary,'native',bam),(star,'star',sam)):
    with (work/f'{label}-pair-flags.log').open('w') as log:
        result=subprocess.run([str(exe),'--genomeDir',str(fixtures/(label+'-index')),
            '--readFilesIn',str(source),'--readFilesType','SAM','PE','--outSAMtype','BAM','Unsorted',
            '--outSAMunmapped','Within','--runThreadN','4',
            '--outFileNamePrefix',str(work/f'{label}-pair-flags')+'.'],stdout=log,stderr=subprocess.STDOUT)
    assert result.returncode==0,(label,(work/f'{label}-pair-flags.log').read_text())
a=decode(work/'native-pair-flags.Aligned.out.bam')[2]
b=decode(work/'star-pair-flags.Aligned.out.bam')[2]
assert collections.Counter(a)==collections.Counter(b)
assert any(r[1]&512 for r in a) and any(not r[1]&512 for r in a)
(work/'pair-flags-validation.txt').write_text(f'{len(a)} records: exact decoded equality with STAR\n')
print('PASS reordered BAM mates and fragment filter flags',len(a))
