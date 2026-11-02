"""Compare native BAM IO, sorting and Solo tags with STAR 2.7.11b.
Usage: bam.py NUSTAR STAR existing-solo-fixtures output-directory
The decoder is independent of the product's HTSlib and codec.
"""
from pathlib import Path
import collections, gzip, json, struct, subprocess, sys


def decode(path):
    data=gzip.open(path,'rb').read()
    assert data[:4]==b'BAM\1'
    n,=struct.unpack_from('<i',data,4);header=data[8:8+n].decode();offset=8+n
    n,=struct.unpack_from('<i',data,offset);offset+=4;refs=[]
    for _ in range(n):
        size,=struct.unpack_from('<i',data,offset);offset+=4
        name=data[offset:offset+size-1].decode();offset+=size
        length,=struct.unpack_from('<i',data,offset);offset+=4;refs.append((name,length))
    records=[]
    while offset<len(data):
        size,=struct.unpack_from('<i',data,offset);offset+=4
        b=data[offset:offset+size];offset+=size
        tid,pos,bmn,fn,lseq,mtid,mpos,tlen=struct.unpack_from('<iiIIiiii',b)
        qlen=bmn&255;mapq=(bmn>>8)&255;flag=fn>>16;nc=fn&65535
        qname=b[32:32+qlen-1].decode();p=32+qlen
        cigar=''.join(str(c>>4)+'MIDNSHP=XB'[c&15] for c in struct.unpack_from('<'+'I'*nc,b,p)) or '*';p+=4*nc
        seq=''.join('=ACMGRSVTWYHKDBN'[(b[p+i//2]>>(4*(1-i%2)))&15] for i in range(lseq)) or '*';p+=(lseq+1)//2
        qual='*' if not lseq or b[p]==255 else ''.join(chr(x+33) for x in b[p:p+lseq]);p+=lseq
        tags=[];fmts={'c':'b','C':'B','s':'h','S':'H','i':'i','I':'I','f':'f'}
        while p<len(b):
            tag=b[p:p+2].decode();typ=chr(b[p+2]);p+=3
            if typ in fmts:
                value,=struct.unpack_from('<'+fmts[typ],b,p);p+=struct.calcsize(fmts[typ]);typ='f' if typ=='f' else 'i';value=str(value)
            elif typ in 'ZH':
                end=b.index(0,p);value=b[p:end].decode();p=end+1
            elif typ=='A':value=chr(b[p]);p+=1
            elif typ=='B':
                sub=chr(b[p]);n,=struct.unpack_from('<i',b,p+1);p+=5
                values=struct.unpack_from('<'+fmts[sub]*n,b,p);p+=struct.calcsize(fmts[sub])*n
                value=sub+','+','.join(map(str,values))
            else:raise AssertionError(('aux',typ))
            tags.append(f'{tag}:{typ}:{value}')
        records.append((qname,flag,tid,pos,mapq,cigar,mtid,mpos,tlen,seq,qual,tuple(tags)))
    assert offset==len(data)
    return header,refs,records


def sam_text(bam):
    header,refs,records=decode(bam)
    lines=[header]
    for q,f,t,p,m,c,mt,mp,tl,s,qual,tags in records:
        lines.append('\t'.join(map(str,[q,f,refs[t][0] if t>=0 else '*',p+1,m,c,
            '=' if mt==t and t>=0 else refs[mt][0] if mt>=0 else '*',mp+1,tl,s,qual,*tags]))+'\n')
    return ''.join(lines)


def main():
    binary,star,fixtures,work=(Path(x).resolve() for x in sys.argv[1:5]);work.mkdir(parents=True,exist_ok=True)
    def run(exe,args,name,success=True):
        with (work/(name+'.log')).open('w') as log:
            p=subprocess.run([str(exe),*map(str,args)],cwd=work,stdout=log,stderr=subprocess.STDOUT)
        assert (p.returncode==0)==success,(name,p.returncode,(work/(name+'.log')).read_text()[-4000:])
    def same_solo(left,right):
        a={p.relative_to(left):p.read_bytes() for p in left.rglob('*') if p.is_file()}
        b={p.relative_to(right):p.read_bytes() for p in right.rglob('*') if p.is_file()}
        assert a.keys()==b.keys()
        for key in a:assert a[key]==b[key],key
    base=['--soloType','CB_UMI_Simple','--soloCBstart',1,'--soloCBlen',8,'--soloUMIstart',9,
          '--soloUMIlen',6,'--soloCBwhitelist',fixtures/'whitelist.txt','--soloCellFilter','None']
    attrs=['NH','HI','AS','nM','CR','CY','UR','UY','CB','UB','GX','GN','sM','sS','sQ','sF']
    cases=[('se-tags',False,False,['--outSAMattributes','NH','HI','AS','nM','NM','MD','jM','jI','XS']),
           ('compression0',False,False,['--outBAMcompression',0]),
           ('se',False,False,[]),('pe',True,False,[]),('solo',False,True,[]),('solo-pe',True,True,[]),
           ('solo-full',False,True,['--soloFeatures','GeneFull']),
           ('solo-cr',False,True,['--soloUMIdedup','1MM_CR','--soloUMIfiltering','MultiGeneUMI_CR']),
           ('two-pass',False,True,['--twopassMode','Basic']),
           ('profile',False,True,[]),
           ('solo-cb-only',False,True,[]),('solo-ub-only',False,True,[]),
           ('tag-only',False,True,['--soloType','CB_samTagOut'])]
    passed=[]
    for name,paired,solo,extra in cases:
      for threads in (1,4):
        case=f'{name}-{threads}'
        for exe,label in ((binary,'native'),(star,'star')):
            reads=[fixtures/'reads.fq']+([fixtures/'mates.fq'] if paired else [])+([fixtures/'barcodes.fq'] if solo else [])
            args=['--genomeDir',fixtures/(label+'-index'),'--readFilesIn',*reads,
                '--outFileNamePrefix',str(work/(label+'-'+case))+'.','--runThreadN',threads,
                '--outSAMtype','BAM','Unsorted','SortedByCoordinate','--outSAMunmapped','Within',
                '--limitBAMsortRAM',4096 if label=='native' else 100000000,*([] if name=='tag-only' else extra)]
            if solo:
                caseattrs=attrs
                if name=='solo-cb-only':caseattrs=['NH','CB']
                if name=='solo-ub-only':caseattrs=['NH','UB']
                if name=='tag-only':caseattrs=['NH','CR','CY','UR','UY','CB','sM']
                args+=base+['--outSAMattributes',*caseattrs]
                if name=='tag-only':
                    args[args.index('CB_UMI_Simple')]='CB_samTagOut'
                    args+=['--soloCBmatchWLtype','1MM']
            if name=='profile' and label=='native':args+=['--profile-stages']
            run(exe,args,label+'-'+case)
        for suffix in ('Aligned.out.bam','Aligned.sortedByCoord.out.bam'):
            a=decode(work/f'native-{case}.{suffix}');b=decode(work/f'star-{case}.{suffix}')
            assert a[1]==b[1],(case,'references')
            # BAM integer widths are storage choices. Everything in each decoded record, including tag order, must match.
            ca,cb=collections.Counter(a[2]),collections.Counter(b[2])
            assert ca==cb,(case,suffix,'records',list((ca-cb).items())[:1],list((cb-ca).items())[:1])
            if 'sorted' in suffix:
                keys=[(r[2]&0xffffffff,r[3]&0xffffffff) for r in a[2]]
                assert keys==sorted(keys),(case,'coordinate order')
                assert 'SO:coordinate' in a[0]
                assert a[2]==b[2],(case,'equal-coordinate tie order')
        if solo:same_solo(work/f'native-{case}.Solo.out',work/f'star-{case}.Solo.out')
        assert not (work/f'native-{case}._NuSTARbamtmp').exists()
        passed.append(case);print('PASS',case,flush=True)
    # Use exactly the same records for native BAM input and STAR SAM input. This
    # exercises reverse strands, secondary records, unmapped reads and aux tags.
    for label,solo,paired in [('se',False,False),('solo',True,False),('pe',False,True),('solo-pe',True,True)]:
        source=work/f'native-{label}-1.Aligned.out.bam'
        sam=work/f'{label}-input.sam';sam.write_text(sam_text(source))
        for exe,engine,inputpath in ((binary,'native',source),(star,'star',sam)):
            args=['--genomeDir',fixtures/(engine+'-index'),'--readFilesIn',inputpath,
                '--readFilesType','SAM','PE' if paired else 'SE','--outFileNamePrefix',str(work/f'{engine}-input-{label}')+'.',
                '--runThreadN',4,'--outSAMtype','BAM','SortedByCoordinate','--outSAMunmapped','Within','--limitBAMsortRAM',100000000]
            if solo:args+=base+['--soloInputSAMattrBarcodeSeq','CR','UR','--soloInputSAMattrBarcodeQual','CY','UY','--outSAMattributes',*attrs]
            run(exe,args,engine+'-input-'+label)
        a=decode(work/f'native-input-{label}.Aligned.sortedByCoord.out.bam')[2]
        b=decode(work/f'star-input-{label}.Aligned.sortedByCoord.out.bam')[2]
        assert collections.Counter(a)==collections.Counter(b),('input',label,list((collections.Counter(a)-collections.Counter(b)).items())[:1])
        if solo:same_solo(work/f'native-input-{label}.Solo.out',work/f'star-input-{label}.Solo.out')
        passed.append('input-'+label);print('PASS','input-'+label,flush=True)
        # Also exercise the text SAM parser, including packed CIGAR/array tags.
        args[args.index('--genomeDir')+1]=fixtures/'native-index'
        args[args.index('--outFileNamePrefix')+1]=str(work/f'native-text-input-{label}')+'.'
        run(binary,args,'native-text-input-'+label)
        text_records=decode(work/f'native-text-input-{label}.Aligned.sortedByCoord.out.bam')[2]
        assert collections.Counter(text_records)==collections.Counter(b),('text-input',label)
        if solo:same_solo(work/f'native-text-input-{label}.Solo.out',work/f'star-input-{label}.Solo.out')
        passed.append('text-input-'+label)
    # Attribute retention is a BAM-specific STAR input option.
    for keep in (['None'],['CR','UR']):
        name='keep-'+'-'.join(keep)
        for exe,label,source in ((binary,'native',work/'native-solo-1.Aligned.out.bam'),(star,'star',work/'solo-input.sam')):
            run(exe,['--genomeDir',fixtures/(label+'-index'),'--readFilesIn',source,
                '--readFilesType','SAM','SE','--readFilesSAMattrKeep',*keep,
                '--outSAMtype','BAM','Unsorted','--outSAMunmapped','Within',
                '--outFileNamePrefix',str(work/f'{label}-{name}')+'.'],label+'-'+name)
        a=decode(work/f'native-{name}.Aligned.out.bam')[2];b=decode(work/f'star-{name}.Aligned.out.bam')[2]
        assert collections.Counter(a)==collections.Counter(b),name
        passed.append(name)
    bad=work/'bad-pair.sam'
    bad.write_text('@HD\tVN:1.4\n@SQ\tSN:chr1\tLN:30000\n'+
        'a\t65\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII\n'+
        'b\t129\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII\n')
    for exe,label in ((binary,'native'),(star,'star')):
        run(exe,['--genomeDir',fixtures/(label+'-index'),'--readFilesIn',bad,
            '--readFilesType','SAM','PE','--outSAMtype','None',
            '--outFileNamePrefix',str(work/f'{label}-bad-pair')+'.'],label+'-bad-pair',False)
    passed.append('reject-mismatched-mates')
    (work/'validation.json').write_text(json.dumps({'passed':passed},indent=2))
if __name__=='__main__':main()
