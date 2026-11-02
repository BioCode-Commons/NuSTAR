"""Exercise STAR-style arguments through real index building and SE/PE mapping."""
from pathlib import Path
import random
import subprocess
import sys
import tempfile


binary = Path(sys.argv[1]).resolve()
star = Path(sys.argv[2]).resolve() if len(sys.argv) > 2 else None


def run(*args, ok=True, cwd=None):
    result = subprocess.run([str(binary), *map(str, args)], cwd=cwd,
                            capture_output=True, text=True, timeout=60)
    assert (result.returncode == 0) == ok, (args, result.stdout, result.stderr)
    return result


def records(path):
    return [line for line in path.read_text().splitlines() if not line.startswith("@")]


def statistics(path):
    return [line.split("|", 1) for line in path.read_text().splitlines()
            if "|" in line and not any(key in line for key in ("time", "speed", "Started", "Finished"))]


with tempfile.TemporaryDirectory(prefix="nustar-star-cli-") as temp:
    root = Path(temp)
    rng = random.Random(23)
    genome = "".join(rng.choice("ACGT") for _ in range(4000))
    fasta = root / "genome.fa"
    fasta.write_text(">chr1\n" + genome + "\n")
    mates = [root / "mate1.fastq", root / "mate2.fastq"]
    revcomp = lambda s: s.translate(str.maketrans("ACGT", "TGCA"))[::-1]
    for mate, seq in zip(mates, (genome[100:200], revcomp(genome[250:350]))):
        mate.write_text("@mapped\n" + seq + "\n+\n" + "I" * 100 +
                        "\n@unmapped\n" + "N" * 100 + "\n+\n" + "I" * 100 + "\n")
    index = root / "GenomeDir"
    # runMode need not be the first option; omitted genomeDir uses STAR's default.
    run("--genomeFastaFiles", fasta, "--runMode", "genomeGenerate",
        "--genomeSAindexNbases", 3, "--genomeChrBinNbits", 8, "--runThreadN", 2,
        "--sjdbGTFtagExonParentGeneType", "gene_type", "gene_biotype",
        "--sjdbGTFfile", "-", "--sjdbOverhang", 100, "--sjdbGTFchrPrefix", "-", cwd=root)
    if star:
        subprocess.run([str(star), "--runMode", "genomeGenerate", "--genomeDir", str(root / "star-index"),
                        "--genomeFastaFiles", str(fasta), "--genomeSAindexNbases", "3",
                        "--genomeChrBinNbits", "8", "--runThreadN", "2"],
                       cwd=root, check=True, capture_output=True, timeout=60)
    for paired in (False, True):
        reads = mates if paired else mates[:1]
        baseline_prefix = str(root / ("paired." if paired else "single."))
        baseline = Path(baseline_prefix + "Aligned.out.sam")
        run("--runMode", "alignReads", "--genomeDir", index, "--readFilesIn", *reads,
            "--outFileNamePrefix", baseline_prefix, "--outSAMunmapped", "Within",
            "--runThreadN", 1)
        for threads in (1, 2):
            prefix = str(root / "nested" / f"{paired}-{threads}.")
            args = ["--readFilesIn", *reads, "--genomeDir", index,
                    "--outFileNamePrefix", prefix, "--runThreadN", threads,
                    "--outSAMtype", "SAM", "--outSAMunmapped", "Within",
                    "--outSAMorder", "PairedKeepInputOrder"]
            run(*args, *([] if paired else ["--profile-stages"]))
            assert records(Path(prefix + "Aligned.out.sam")) == records(baseline)
            assert Path(prefix + "SJ.out.tab").read_bytes() == Path(baseline_prefix + "SJ.out.tab").read_bytes()
            assert statistics(Path(prefix + "Log.final.out")) == statistics(Path(baseline_prefix + "Log.final.out"))
        # Default prefix './', default mode alignReads and unmapped mode None.
        out = root / ("default-pe" if paired else "default-se")
        out.mkdir()
        run("--genomeDir", index, "--readFilesIn", *reads, cwd=out)
        expected = [line for line in records(baseline) if not int(line.split("\t")[1]) & 4]
        assert expected and records(out / "Aligned.out.sam") == expected
        assert statistics(out / "Log.final.out") == statistics(Path(baseline_prefix + "Log.final.out"))
        assert (out / "SJ.out.tab").exists()
        if star:
            for mode in ("None", "Within"):
                prefix = str(root / f"star-{paired}-{mode}.")
                subprocess.run([str(star), "--genomeDir", str(root / "star-index"),
                                "--readFilesIn", *map(str, reads), "--runThreadN", "2",
                                "--outSAMunmapped", mode, "--outSAMorder", "PairedKeepInputOrder",
                                "--outFileNamePrefix", prefix],
                               cwd=root, check=True, capture_output=True, timeout=60)
                actual = baseline if mode == "Within" else out / "Aligned.out.sam"
                assert records(Path(prefix + "Aligned.out.sam")) == records(actual)
                assert statistics(Path(prefix + "Log.final.out")) == statistics(out / "Log.final.out")
                assert Path(prefix + "SJ.out.tab").read_bytes() == (out / "SJ.out.tab").read_bytes()
    for args in (["--genomeDir", index, "--readFilesIn", mates[0], "--threads", "1"],
                 ["--runMode", "genomeGenerate", "--genomeFastaFiles", fasta,
                  "--genomeDir", root / "old-index", "--genomeSAindexNbases", "3", "--threads", "1"], ["alignSingleSam", index, mates[0], root / "old.sam"],
                 ["alignPairedSam", index, *mates, root / "old.sam"],
                 ["genomeGenerate"], ["build-index"], ["--runThreadN", "0"], ["--runThreadN", "x"], ["--genomeDir"],
                 ["--readFilesIn"], ["--readFilesIn", "a", "b", "c"],
                 ["--readFilesIn", "a,b"], ["--outSAMtype", "BAM"],
                 ["--outSAMunmapped", "Within", "KeepPairs"],
                 ["--runMode", "unknown"], ["--unknown", "x"],
                 ["--readFilesIn", *mates, "--profile-stages"]):
        run(*args, ok=False, cwd=root)
print("STAR CLI: SE/PE, profiling, threads, prefixes, defaults and errors passed")
