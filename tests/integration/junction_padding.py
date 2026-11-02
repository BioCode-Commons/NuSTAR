"""Compare chromosome-edge junction sequences, indexes and two-pass mapping with STAR.

Usage: junction_padding.py NUSTAR STAR COMPARE_INDEX WORK_DIRECTORY
The introns stay inside chromosomes; only their 150-base flanks enter padding.
"""
from pathlib import Path
import json
import random
import subprocess
import sys

from compare_genome_generate import nustar_genome

binary, star, comparator, work = (Path(x).resolve() for x in sys.argv[1:5])
work.mkdir(parents=True, exist_ok=True)
programs = {"star": star, "nustar": binary}


def run(command, name):
    with (work / (name + ".log")).open("w") as log:
        result = subprocess.run(list(map(str, command)), cwd=work, stdout=log,
                                stderr=subprocess.STDOUT, timeout=180)
    assert result.returncode == 0, (name, (work / (name + ".log")).read_text()[-4000:])


def records(path):
    return sorted(line for line in path.read_text().splitlines() if not line.startswith("@"))


def stats(path):
    return [line for line in path.read_text().splitlines() if "|" in line
            and not any(word in line.lower() for word in ("time", "speed", "started", "finished"))]


def compare_index(star_index, native_index, name):
    assert nustar_genome(native_index / "reference.bin") == (star_index / "Genome").read_bytes()
    for filename in ("sjdbInfo.txt", "sjdbList.out.tab"):
        assert (star_index / filename).read_bytes() == (native_index / filename).read_bytes(), filename
    run([comparator, star_index, native_index], name)


rng = random.Random(527)
sequences = {name: list("".join(rng.choice("ACGT") for _ in range(length)))
             for name, length in (("guard", 2000), ("XII", 17000), ("chrM", 16569))}
# The left-edge cases use coordinates found in the real failing datasets.
events = [("XII", 75, 5709), ("XII", 75, 11227), ("chrM", 144, 15989),
          ("chrM", 400, 16540)]
for chromosome, begin, end in events:
    sequences[chromosome][begin:begin + 2] = "GT"
    sequences[chromosome][end - 2:end] = "AG"
sequences = {name: "".join(seq) for name, seq in sequences.items()}
fasta, gtf, fastq = (work / name for name in ("genome.fa", "genes.gtf", "reads.fq"))
fasta.write_text("".join(f">{name}\n{seq}\n" for name, seq in sequences.items()))
annotation, reads = [], []
for number, (chromosome, begin, end) in enumerate(events):
    sequence = sequences[chromosome]
    for start, stop in ((max(0, begin - 100), begin), (end, min(len(sequence), end + 100))):
        annotation.append(f'{chromosome}\ttest\texon\t{start + 1}\t{stop}\t.\t+\t.\t'
                          f'gene_id "g{number}"; transcript_id "t{number}";\n')
    for anchor in (20, 35, 50, 70, 75, 80):
        if anchor > begin or end + 100 - anchor > len(sequence):
            continue
        read = sequence[begin - anchor:begin] + sequence[end:end + 100 - anchor]
        for repeat in range(6):
            reads.append(read if repeat % 2 == 0 else read.translate(str.maketrans("ACGT", "TGCA"))[::-1])
gtf.write_text("".join(annotation))
fastq.write_text("".join(f"@read{i}\n{read}\n+\n{'I' * len(read)}\n" for i, read in enumerate(reads)))

for annotated in (False, True):
    for label, exe in programs.items():
        run([exe, "--runMode", "genomeGenerate", "--genomeDir", work / f"{label}-{annotated}",
             "--genomeFastaFiles", fasta, "--genomeChrBinNbits", 18, "--genomeSAindexNbases", 3,
             "--runThreadN", 2,
             *(["--sjdbGTFfile", gtf, "--sjdbOverhang", 150] if annotated else [])], f"build-{label}-{annotated}")
compare_index(work / "star-True", work / "nustar-True", "annotated-index")
# Confirm the oracle fixture really contains padding inside the junction flanks.
genome = (work / "star-True/Genome").read_bytes()
real_end = int((work / "star-True/chrStart.txt").read_text().splitlines()[-1])
blocks = [genome[i:i + 301] for i in range(real_end, len(genome), 301)]
assert len(blocks) == len(events)
assert any(5 in block[:150] for block in blocks)
assert any(5 in block[150:300] for block in blocks)

checks = []
for annotated in (False, True):
    for threads in (1, 4):
        suffix = f"{annotated}-{threads}"
        for label, exe in programs.items():
            prefix = str(work / f"{label}-{suffix}.")
            run([exe, "--genomeDir", work / f"{label}-{annotated}", "--readFilesIn", fastq,
                 "--runThreadN", threads, "--twopassMode", "Basic", "--sjdbOverhang", 150,
                 "--sjdbInsertSave", "All", "--outSAMunmapped", "Within",
                 "--outFileNamePrefix", prefix], f"map-{label}-{suffix}")
            run([exe, "--genomeDir", prefix + "_STARgenome", "--readFilesIn", fastq,
                 "--runThreadN", threads, "--outSAMunmapped", "Within",
                 "--outFileNamePrefix", prefix + "reload."], f"reload-{label}-{suffix}")
        a, b = (str(work / f"{label}-{suffix}.") for label in ("star", "nustar"))
        for name in ("Aligned.out.sam", "SJ.out.tab", "Log.final.out", "_STARpass1/SJ.out.tab",
                     "_STARpass1/Log.final.out", "reload.Aligned.out.sam"):
            reader = records if name.endswith(".sam") else stats if name.endswith("Log.final.out") else Path.read_bytes
            assert reader(Path(a + name)) == reader(Path(b + name)), (suffix, name)
        assert records(Path(b + "Aligned.out.sam")) == records(Path(b + "reload.Aligned.out.sam"))
        found = {tuple(line.split()[:3]) for line in Path(a + "_STARpass1/SJ.out.tab").read_text().splitlines()}
        assert all((chromosome, str(begin + 1), str(end)) in found for chromosome, begin, end in events), found
        compare_index(Path(a + "_STARgenome"), Path(b + "_STARgenome"), "inserted-index-" + suffix)
        checks.append(suffix)
(work / "validation.json").write_text(json.dumps({"passed": checks, "reads": len(reads),
    "boundary_junctions": len(events), "checks": "exact padded sequences, SJ tables, SA ranks, prefixes, SAM and reload"}, indent=2) + "\n")
print("Chromosome padding: sequence/index equality and", len(checks), "two-pass/reload cases passed")
