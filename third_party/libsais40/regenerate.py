#!/usr/bin/env python3
"""Regenerate the SA-only five-byte adaptation from the pinned vendored source."""
from pathlib import Path
import hashlib
import re
root=Path(__file__).resolve().parents[2]
upstream=root/'third_party/libsais/libsais64.c'
assert hashlib.sha256(upstream.read_bytes()).hexdigest()=='088019c5fc5648284291f96657d0cfd2c3c9bf53a4d12514773fdab7f3e452bf', 'Upstream changed: audit the port before regenerating'
s=upstream.read_text()
s=s[:s.index('static sa_sint_t libsais64_main_long(')]
# Disable layouts tied to native 64-bit entries. Software prefetch remains enabled.
s=s.replace('#include "libsais.h"\n#include "libsais64.h"','#include "slot.hpp"\n#include "native_recursion.hpp"\n#define LIBSAIS_CPU_FEATURE 0')
# All declarations remain within a private namespace, after system headers.
idx=s.index('typedef int64_t')
s=s[:idx]+'namespace nustar::star2711b::detail::libsais40_impl {\n'+s[idx:]
s=re.sub(r'#define SAINT_BIT\s+\(64\)', '#define SAINT_BIT (40)',s)
s=re.sub(r'#define SAINT_MAX\s+INT64_MAX', '#define SAINT_MAX ((INT64_C(1) << 39) - 1)',s)
s=re.sub(r'#define SAINT_MIN\s+INT64_MIN', '#define SAINT_MIN (-(INT64_C(1) << 39))',s)
# Do not reinterpret five-byte storage as the upstream 8-byte/4-byte hybrid.
a=s.index('    if (n <= INT32_MAX && k <= INT32_MAX)',s.index('static sa_sint_t libsais64_main_32s_recursion'))
b=s.index('    if (k > 0 &&',a)
s=s[:a]+'#if !defined(NUSTAR_LIBSAIS40_FORCE_COMPACT_RECURSION)\n    if (n >= 64 && n <= INT32_MAX && k <= INT32_MAX)\n    {\n        return nativeRecursion(T, SA, n, k, fs, static_cast<int>(threads));\n    }\n#endif\n\n'+s[b:]
a=s.index('static void libsais64_convert_32u_to_64u(')
b=s.index('static sa_sint_t libsais64_main_32s_recursion(',a)
s=s[:a]+s[b:]
a=s.index('{',s.index('static sa_sint_t libsais64_induce_final_order_8u_omp('))
b=s.index('\nstatic void libsais64_induce_final_order_32s_6k',a)
s=s[:a]+'''{
    (void)flags; (void)r; (void)I;
    libsais64_final_sorting_scan_left_to_right_8u_omp(T, SA, n, k, &buckets[6 * ALPHABET_SIZE], threads, thread_state);
    if (threads > 1 && n >= 65536) { libsais64_clear_lms_suffixes_omp(SA, n, ALPHABET_SIZE, &buckets[6 * ALPHABET_SIZE], &buckets[7 * ALPHABET_SIZE], threads); }
    libsais64_final_sorting_scan_right_to_left_8u_omp(T, SA, 0, n, k, &buckets[7 * ALPHABET_SIZE], threads, thread_state);
    return 0;
}
''' + s[b:]
# Array storage, including recursive bucket arrays, is five bytes. Scalar locals
# and fixed-size thread-cache records intentionally remain int64_t.
s=re.sub(r'\bsa_sint_t\s*\*', 'Slot *', s)
s=re.sub(r'\bsa_uint_t\s*\*', 'Slot *', s)
s=re.sub(r'\bsa_sint_t(\s+\w+\s*\[)', r'Slot\1',s)
s=s.replace('sizeof(sa_sint_t)','sizeof(Slot)')
s=re.sub(r'(typedef struct LIBSAIS_NONUNIQUE_SUFFIX\s*\{)(.*?)(\} LIBSAIS_NONUNIQUE_SUFFIX;)',lambda m:m[1]+m[2].replace('sa_sint_t','Slot')+m[3],s,flags=re.S)
# Slot alignment is one byte; do not align a pointer away from its slot grid.
s=re.sub(r'sa_sint_t alignment = [^;]+;', 'sa_sint_t alignment = 0;',s)
s=re.sub(r'\(Slot \*\)libsais64_align_up\((&SA\[[^\]]+\]), \(size_t\)alignment \* sizeof\(Slot\)\)',r'\1',s)
s=s.replace('Slot * threads, LIBSAIS_THREAD_STATE', 'sa_sint_t * threads, LIBSAIS_THREAD_STATE')
# Preserve native mask objects instead of aliasing compressed bucket storage.
s=s.replace('LIBSAIS_THREAD_CACHE *          cache;', 'LIBSAIS_THREAD_CACHE *          cache;\n        uint64_t *                     masks;')
s=s.replace('    if (thread_state != NULL && thread_buckets != NULL && thread_cache != NULL)',
'''    uint64_t * thread_masks = (uint64_t *)libsais64_alloc_aligned((size_t)threads * (LIBSAIS_PER_THREAD_CACHE_SIZE / 64) * sizeof(uint64_t), 4096);
    if (thread_state != NULL && thread_buckets != NULL && thread_cache != NULL && thread_masks != NULL)''')
s=s.replace('thread_state[t].state.cache     = thread_cache;     thread_cache    += LIBSAIS_PER_THREAD_CACHE_SIZE;',
'''thread_state[t].state.cache     = thread_cache;     thread_cache    += LIBSAIS_PER_THREAD_CACHE_SIZE;
            thread_state[t].state.masks = thread_masks; thread_masks += LIBSAIS_PER_THREAD_CACHE_SIZE / 64;''')
s=s.replace('    libsais64_free_aligned(thread_cache);', '    libsais64_free_aligned(thread_masks);\n    libsais64_free_aligned(thread_cache);')
s=s.replace('        libsais64_free_aligned(thread_state[0].state.cache);', '        libsais64_free_aligned(thread_state[0].state.masks);\n        libsais64_free_aligned(thread_state[0].state.cache);')
s=s.replace('(uint64_t *)(void *)thread_state[0].state.buckets','thread_state[0].state.masks')
# Scratch record fields alias ordinary slots through an explicit view.
s=re.sub(r'typedef struct LIBSAIS_NONUNIQUE_SUFFIX\s*\{.*?\} LIBSAIS_NONUNIQUE_SUFFIX;', '',s,flags=re.S)
s=s.replace('LIBSAIS_NONUNIQUE_SUFFIX * RESTRICT suffixes','RecordView suffixes')
s=s.replace('LIBSAIS_NONUNIQUE_SUFFIX s =','Record s =')
s=s.replace('(Slot *)(void *)&suffixes[group_size]', '(suffixes.slots + 3 * group_size)')
s=re.sub(r'&suffixes\[([^\]]+)\]',r'(suffixes + (\1))',s)
s=s.replace('(LIBSAIS_NONUNIQUE_SUFFIX *)(void *)scratch', 'RecordView{scratch}')
s=s.replace('(LIBSAIS_NONUNIQUE_SUFFIX *)(void *)&scratch[omp_thread_num * omp_capacity]', 'RecordView{&scratch[omp_thread_num * omp_capacity]}')
# Casts in the original implement signed/unsigned 64-bit arithmetic; preserve
# those semantics at logical width 40 even though register arithmetic is 64-bit.
def skip(p):
 while p<len(s) and s[p].isspace():p+=1
 return p
def group(p):
 end={ '(':')','[':']','{':'}'}[s[p]];p+=1
 while p<len(s):
  if s[p]==end:return p+1
  if s[p] in '([{':p=group(p)
  else:p+=1
 raise ValueError('unclosed')
def unary(p):
 p=skip(p)
 if s.startswith(('++','--'),p):return unary(p+2)
 if s[p] in '+-!~*&':return unary(p+1)
 if s[p]=='(':
  q=group(p)
  if re.fullmatch(r'(?:const\s+)?\w+(?:\s*\*)*',s[p+1:q-1].strip()) and s[p+1:q-1].strip() in ['sa_sint_t','sa_uint_t','fast_sint_t','fast_uint_t','uint64_t','int64_t']:
   return unary(q)
  p=q
 else:
  m=re.match(r'\w+',s[p:]); assert m,s[p:p+80];p+=len(m[0])
 while True:
  q=skip(p)
  if q<len(s) and s[q] in '[(':p=group(q)
  elif s.startswith('->',q) or s[q:q+1]=='.':
   q=skip(q+(2 if s.startswith('->',q) else 1));m=re.match(r'\w+',s[q:]);p=q+len(m[0])
  elif s.startswith(('++','--'),q):p=q+2
  else:return p
# Replace innermost/rightmost casts first so nested expressions remain valid.
for token,fn in [('sa_uint_t','unsigned40'),('sa_sint_t','signed40')]:
 matches=list(re.finditer(r'\('+token+r'\)',s))
 for m in reversed(matches):
  a=skip(m.end());b=unary(a);s=s[:m.start()]+fn+'('+s[a:b]+')'+s[b:]
s+='\n} // namespace nustar::star2711b::detail::libsais40_impl\n'
s=s.replace('libsais64_', 'libsais40_')
s=s.replace('This file is a part of libsais,','This file is adapted from libsais v2.12.1 for NuSTAR five-byte storage,')
# Retain only functions reachable from the suffix-array entry point. Preserve
# preprocessor alternatives for a retained function; unused BWT kernels disappear.
functions=[]
for match in re.finditer(r'^static [^\n]*\b(libsais40_\w+)\([^\n]*\)\n\{',s,re.M):
 end=s.index('\n}',match.end())+2
 functions.append((match[1],match.start(),end,s[match.start():end]))
reachable={'libsais40_main'}
while True:
 deps=set().union(*(set(re.findall(r'\b(libsais40_\w+)\s*\(',body)) for name,_,_,body in functions if name in reachable))
 if deps<=reachable:break
 reachable|=deps
for name,start,end,_ in reversed(functions):
 if name not in reachable:s=s[:start]+s[end:]
s=re.sub(r'\n{4,}', '\n\n\n',s)
(root/'third_party/libsais40/libsais40.inc').write_text(s)
print('wrote',len(s.splitlines()),'lines')
