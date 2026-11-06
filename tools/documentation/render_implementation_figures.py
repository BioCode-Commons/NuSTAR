"""Generate the two README figures and their local walkthrough, using stdlib only.

Run from the repository root:
    python3 tools/documentation/render_implementation_figures.py

SVG is the editable figure source. PNG exports are browser captures of these
same SVGs; export_implementation_animation.py adds a panel highlight for GIF.
The interactive walkthrough exposes renderAt(milliseconds) for deterministic frames.
"""
from html import escape
from pathlib import Path
import json

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'docs/assets'
OUT.mkdir(parents=True, exist_ok=True)
STAR = 'https://github.com/alexdobin/STAR/blob/b1edc1208d91a53bf40ebae8669f71d50b994851/source/'
NUSTAR = 'https://github.com/BioCode-Commons/NuSTAR/blob/bdf6bfa7d465c965edac091c6f13d292a81be1ef/'
STYLE = '''text{font-family:Arial,Helvetica,sans-serif;font-size:19px;fill:#e6edf3}
.title{font-size:31px;font-weight:600}.heading{font-size:22px;font-weight:600}
.mono{font-family:Menlo,Consolas,monospace;font-size:17px}.small{font-size:16px}
.muted{fill:#9da7b3}.green{fill:#7ee787}.blue{fill:#79c0ff}.amber{fill:#e3b341}
.eyebrow{font-size:15px;letter-spacing:1.6px;fill:#9da7b3}.nodeLabel{font-size:19px;font-weight:600}
'''

def txt(x, y, label, cls='', anchor=None):
    attr = f' text-anchor="{anchor}"' if anchor else ''
    return f'<text x="{x}" y="{y}" class="{cls}"{attr}>{escape(label)}</text>'

def rect(x, y, w, h, fill, stroke=None, rx=8, extra=''):
    return f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="{rx}" fill="{fill}"'+(f' stroke="{stroke}"' if stroke else '')+f' {extra}/>'

def start(name, description, height):
    return [f'<svg xmlns="http://www.w3.org/2000/svg" width="1120" height="{height}" viewBox="0 0 1120 {height}" role="img" aria-labelledby="{name}-title {name}-desc">',f'<title id="{name}-title">{escape(description)}</title>',f'<desc id="{name}-desc">{escape(description)}. Source revisions: STAR 2.7.11b and NuSTAR bdf6bfa.</desc>',f'<style>{STYLE}</style>',rect(0,0,1120,height,'#0d1117',rx=14),'<defs><marker id="arrow" markerWidth="7" markerHeight="7" refX="6" refY="3.5" orient="auto"><path d="M0 0L7 3.5L0 7Z" fill="#6e7681"/></marker></defs>']

def node(x,y,title,lines,w=408,h=70,native=False):
    parts=[rect(x,y,w,h,'#12251f' if native else '#1c212b','#3a4c48' if native else '#303b4c')]
    parts.append(txt(x+w/2,y+25,title,'nodeLabel', 'middle'))
    parts += [txt(x+w/2,y+49+i*23,line,'small muted','middle') for i,line in enumerate(lines)]
    return parts

def path(d,color='#6e7681',arrow=True):
    return f'<path d="{d}" fill="none" stroke="{color}" stroke-width="1.6"'+(' marker-end="url(#arrow)"' if arrow else '')+'/>'

flow=start('flow','STAR and NuSTAR: end-to-end flow and processing boundaries',1138)
flow += [txt(32,37,'01 / END-TO-END FLOW','eyebrow'),txt(32,80,'The same alignment strategy, with explicit stage boundaries.','title'),txt(32,111,'Standard spliced alignment; the NuSTAR lane expands its normal single-end pipeline.','small muted')]
flow += [txt(66,153,'STAR 2.7.11b','heading'),txt(622,153,'NuSTAR','heading green')]
# Index generation is a setup step, outside the read/batch loops.
for x,native,title,lines in [(94,False,'Build or load the reference index',['Prefix chunks → qsort → packed suffix array','Genome / SA / SAindex']), (650,True,'Build or load the native index',['libsais64 → STAR rank normalization → pack','reference.bin / index.sa'])]:
    flow += node(x,173,title,lines,h=96,native=native)
    flow.append(path(f'M{x+204} 269V319'))
flow += [txt(52,301,'WORKER / CHUNK LOOP','eyebrow'),txt(608,301,'WORKER / NUMBERED BATCH LOOP','eyebrow')]
OLD=[('Read next chunk',['Input mutex protects shared streams']),('Prepare the next read',['ReadAlign::oneRead()']),('Find mappable seeds',['mapOneRead() updates ReadAlign arrays']),('Build windows and stitch',['stitchPieces() → stitchWindowAligns()']),('Select, filter, and emit',['multMapSelect → mappedFilter → output']),('Finish chunk output',['Chunk buffers / per-worker statistics'])]
NEW=[('Read next numbered batch',['Input mutex; reuse worker workspace']),('Prepare the next fragment',['prepareSingleFastq()']),('Seed search → SeedBatch',['seed.find(read, workspace, seeds)']),('WindowBatch → CandidateSet',['window.build() → stitchSingle()']),('Selection → FragmentEffects',['filter / project / Solo record* / encode']),('Commit batches in input order',['Drain by batch ID; reduce junctions / stats'])]
for x,rows,native in [(94,OLD,False),(650,NEW,True)]:
    for i,(title,lines) in enumerate(rows):
        y=320+i*88
        flow+=node(x,y,title,lines,native=native)
        if i<5:flow.append(path(f'M{x+204} {y+70}V{y+88}'))
    # Per-read and per-chunk loops are present in both versions.
    flow.append(path(f'M{x} 707H{x-22}V443H{x}'))
    flow.append(f'<text x="{x-32}" y="575" class="small muted" text-anchor="middle" transform="rotate(-90 {x-32} 575)">'+('next read' if not native else 'next fragment')+'</text>')
    flow.append(path(f'M{x} 795H{x-54}V355H{x}'))
    flow.append(f'<text x="{x-66}" y="575" class="small muted" text-anchor="middle" transform="rotate(-90 {x-66} 575)">'+('next chunk' if not native else 'next batch')+'</text>')
    flow.append(path(f'M{x+204} 830V882'))
    flow+=node(x,882,'After mapping',['Solo reduction*, then sorted BAM*'],native=native)
flow += [rect(32,978,1056,68,'#152032'),txt(50,1006,'Two-pass mode: discover junctions → extend index → reread input for pass two',''),txt(50,1031,'NuSTAR now merges inserted suffixes in parallel output-rank chunks.','small muted')]
flow += [txt(32,1080,'* Optional. Profiling paths, error exits, and unsupported STAR features are omitted.','small muted'),txt(32,1107,'Both implementations already process reads in worker-owned chunks; NuSTAR makes stage data explicit.','small muted'),'</svg>']
flow='\n'.join(flow)
(OUT/'implementation-flow.svg').write_text(flow)

DIFFS=[
    dict(title='Suffix-array construction',kind='sorting algorithm + compatibility repair',
         old=['for each prefix-range chunk:', '  qsort(..., funCompareSuffixes)', '  write SA_i to disk', 'reload chunks → pack SA'],
         new=['libsais64(...) / libsais64_omp(...)', 'filter non-base suffix starts', 'normalize STAR-equal suffix groups', 'pack + release → build prefix ranges'],
         note='STAR-equal groups still use std::sort; Linux source pages are released during packing.',
         star='Genome_genomeGenerate.cpp#L213-L325',native='index/suffix_array_index.cpp#L271-L291'),
    dict(title='Per-read window reset',kind='full-array clearing → touched-bin log',
         old=['memset(winBin[0], 255, full_size)', 'memset(winBin[1], 255, full_size)', 'then scan seed hits / build windows'],
         new=['setBin(): log bins on first mutation', 'for (bin : touched_bins[strand])', '  bin_to_window[strand][bin] = no_window', 'clear the touched-bin log'],
         note='The bin table stays allocated. Reset work follows touched bins instead of the full table.',
         star='ReadAlign_stitchPieces.cpp#L12-L16',native='window/star_window.cpp#L39-L57'),
    dict(title='Read loop and intermediate state',kind='shared object fields → explicit stage results',
         old=['oneRead() → mapOneRead()', 'ReadAlign fields: PC / WC / WA / trAll', 'stitchPieces() builds + stitches windows', 'multMapSelect() → outputAlignments()'],
         new=['seed.find() → SeedBatch', 'window.build() → WindowBatch', 'stitchSingle() → CandidateSet', 'select / filter / project → effects'],
         note='The recursive stitching kernel is retained. The change is its interface and ownership.',
         star='ReadAlign_mapOneRead.cpp#L34-L118',native='pipeline/single_end/mapping_composition.cpp#L72-L98'),
    dict(title='Coordinate-sorted BAM',kind='coordinate bins → external merge runs',
         old=['write per-worker coordinate-bin files', 'load one bin; qsort(startPos)', 'write each sorted bin', 'bam_cat() concatenates bins in order'],
         new=['buffer records; std::sort by key', 'spill sorted runs to disk', 'priority_queue → k-way merge', 'merge fan-in ≤ 64; repeat as needed'],
         note='Small outputs stay in memory. The sort budget does not cap total process memory.',
         star='BAMbinSortByCoordinate.cpp#L12-L52',native='bam/sink.cpp#L86-L121'),
]
DIFFS.extend([
    dict(title='Two-pass suffix-array merge',kind='serial merge loop → independent output-rank chunks',
         old=['walk the old suffix array in order', 'insert ranked junction suffixes', 'writePacked() into the extended SA'],
         new=['partition output ranks into chunks', 'binary-search first insertion per chunk', 'OpenMP: merge with a RangeWriter', 'discard consumed old-SA pages'],
         note='Chunks are aligned so packed writes share no bytes. No speedup is asserted here.',
         star='sjdbBuildIndex.cpp',native='index/junction_insert.cpp'),
    dict(title='Overlapping paired-end reads',kind='STAR overlap strategy → separate typed mapper workspace',
         old=['peOverlapMergeMap() finds overlap', 'merge mates → mapOneRead()', 'convert alignment back to mates', 'recalculate paired alignment score'],
         new=['map paired candidates first', 'overlap helper → single-read mapper', 'project merged blocks back to mates', 'rescore → paired selection / output'],
         note='NuSTAR defaults to overlap disabled. It keeps original mate sequences for output.',
         star='ReadAlign_peOverlapMergeMap.cpp',native='pipeline/paired_end/star_overlap.cpp'),
    dict(title='Annotation strand and winner selection',kind='STAR preparation rules → explicit proposal passes',
         old=['sjdbPrepare(): order by strand / locus', 'deduplicate normalized coordinates', 'resolve motif / strand / priority', 'build annotated junction sequences'],
         new=['stable_sort proposals by strand / locus', 'first-pass normalized deduplication', 'second-pass priority / strand rules', 'emit selected junction entries'],
         note='The updated rules have synthetic regression fixtures; real-data reruns remain pending.',
         star='sjdbPrepare.cpp',native='annotation/junction_database.cpp'),
])
DIFF_HEIGHT = 216 + len(DIFFS) * 236
diff=start('differences','Source-level implementation differences between STAR and NuSTAR',DIFF_HEIGHT)
diff += [txt(32,37,'02 / IMPLEMENTATION DIFFERENCES','eyebrow'),txt(32,80,'What actually changed inside the pipeline.','title'),txt(32,111,'Abbreviated code and control flow; source links point to the compared revisions.','small muted')]
for i,d in enumerate(DIFFS):
    y=140+i*236
    diff += [f'<g id="difference-{i}" data-difference="{i}">',rect(24,y,1072,220,'#161b22','#30363d'),txt(44,y+30,f'@@ {i+1:02}  {d["title"]} @@','heading'),txt(44,y+55,d['kind'],'small muted'),rect(25,y+67,535,112,'#231b23',rx=0),rect(560,y+67,535,112,'#12251f',rx=0)]
    for row,line in enumerate(d['old']):diff.append(txt(45,y+91+row*25,'− '+line,'mono'))
    for row,line in enumerate(d['new']):diff.append(txt(578,y+91+row*25,'+ '+line,'mono green'))
    diff += [txt(44,y+203,d['note'],'small muted'),f'<a href="{STAR+d["star"]}" target="_blank">',txt(882,y+29,'STAR source','small blue'),'</a>',f'<a href="{NUSTAR+d["native"]}" target="_blank">',txt(984,y+29,'NuSTAR','small green'),'</a>','</g>']
diff += ['<g id="difference-focus" opacity="0">'+rect(25,140,1070,220,'none','#58a6ff',extra='stroke-width="2"')+rect(25,140,4,220,'#58a6ff',rx=0)+'</g>','<circle id="difference-travel" cx="543" cy="263" r="4" fill="#79c0ff" opacity="0"/>',txt(32,DIFF_HEIGHT-38,'Source changes explain mechanisms. They do not measure the speedup of an individual stage.','small muted'),'</svg>']
diff='\n'.join(diff)
(OUT/'implementation-differences.svg').write_text(diff)

html='''<!doctype html>
<html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1"><title>NuSTAR — flows and implementation differences</title>
<style>:root{color-scheme:dark}*{box-sizing:border-box}body{margin:0;background:#0d1117;color:#e6edf3;font:15px/1.5 -apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif}main{max-width:1160px;margin:24px auto;padding:0 20px}svg{display:block;width:100%;height:auto}nav{display:flex;gap:12px;align-items:center;flex-wrap:wrap;padding:12px 8px}button,select{font:inherit;color:inherit;background:#21262d;border:1px solid #3d444d;border-radius:6px;padding:8px 12px;cursor:pointer}button:focus-visible,select:focus-visible{outline:2px solid #58a6ff;outline-offset:3px}button[aria-pressed=true]{border-color:#58a6ff}#status{color:#9da7b3}p{padding:0 8px;color:#9da7b3}a{color:#79c0ff}[hidden]{display:none}body.capture main{max-width:none;margin:0;padding:0}body.capture nav,body.capture p{display:none}body.capture svg{width:1120px}</style>
<main><nav aria-label="Figures"><button type="button" id="show-flow" aria-pressed="true">End-to-end flows</button><button type="button" id="show-diff" aria-pressed="false">Implementation differences</button></nav><section id="flow-panel">FLOW</section><section id="diff-panel" hidden>DIFFERENCES</section><nav id="walkthrough" aria-label="Difference walkthrough" hidden><button id="play" type="button">Play walkthrough</button><label for="selection">Focus</label><select id="selection"><option value="-1">All differences</option>OPTIONS</select><span id="status" role="status" aria-live="polite">All differences visible</span></nav><p>Source-backed descriptions and comparison scope: <a href="../IMPLEMENTATION_COMPARISON.md">implementation comparison</a>. Animation order has no timing scale.</p></main>
<script>
const titles=TITLES;
const flowPanel=document.getElementById('flow-panel');
const diffPanel=document.getElementById('diff-panel');
const flowButton=document.getElementById('show-flow');
const diffButton=document.getElementById('show-diff');
const controls=document.getElementById('walkthrough');
const focus=document.getElementById('difference-focus');
const travel=document.getElementById('difference-travel');
const selection=document.getElementById('selection');
const status=document.getElementById('status');
const play=document.getElementById('play');
let animation;
let playing=false;
function stop(){cancelAnimationFrame(animation);playing=false;play.textContent='Play walkthrough'}
function choose(index,progress=0.5){
 const valid=index>=0&&index<titles.length;
 focus.setAttribute('opacity',Number(valid));
 focus.setAttribute('transform',`translate(0 ${Math.max(0,index)*236})`);
 travel.setAttribute('opacity',Number(valid));
 travel.setAttribute('cx',543+34*progress);
 travel.setAttribute('cy',263+Math.max(0,index)*236);
}
function show(mode){stop();const isFlow=mode==='flow';flowPanel.hidden=!isFlow;diffPanel.hidden=isFlow;controls.hidden=isFlow;flowButton.setAttribute('aria-pressed',String(isFlow));diffButton.setAttribute('aria-pressed',String(!isFlow))}
function renderAt(time){
 const index=Math.floor(time/2400),p=(time%2400)/2400;
 choose(index,p);
 const f=Math.min(1,p/0.15),smooth=f*f*(3-2*f),row=Math.max(0,index-1+smooth);
 focus.setAttribute('transform',`translate(0 ${row*236})`);
 travel.setAttribute('cy',263+row*236);
 return index;
}
window.renderAt=renderAt;
window.showFigure=show;
function run(){
 stop();playing=true;play.textContent='Stop walkthrough';selection.value='-1';let start;
 function frame(now){start??=now;const i=renderAt(now-start);if(i<titles.length){status.textContent=titles[i];animation=requestAnimationFrame(frame)}else{stop();choose(-1);status.textContent='All differences visible'}}
 animation=requestAnimationFrame(frame);
}
flowButton.addEventListener('click',()=>show('flow'));
diffButton.addEventListener('click',()=>show('diff'));
play.addEventListener('click',()=>playing?stop():run());
selection.addEventListener('change',()=>{stop();const i=Number(selection.value);choose(i);status.textContent=i<0?'All differences visible':titles[i]});
const params=new URLSearchParams(location.search);
show(params.get('figure')||'flow');
if(params.has('capture'))document.body.classList.add('capture');
window.captureReady=true;
</script></html>
'''
html=html.replace('FLOW',flow,1).replace('DIFFERENCES',diff,1).replace('OPTIONS',''.join(f'<option value="{i}">{escape(d["title"])}</option>' for i,d in enumerate(DIFFS)),1).replace('TITLES',json.dumps([d['title'] for d in DIFFS]),1)
(OUT/'implementation-comparison.html').write_text(html)
print('Generated implementation-flow.svg, implementation-differences.svg, and implementation-comparison.html')
