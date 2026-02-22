#!/usr/bin/env python3
import argparse
import csv
import gzip
import hashlib
import shutil
import subprocess
from collections import defaultdict
from pathlib import Path
from typing import Dict, Optional


# ----------------- utilities -----------------

def die(msg: str) -> None:
    raise SystemExit(f"[ERROR] {msg}")


def which_or_die(exe: str) -> None:
    if shutil.which(exe) is None:
        die(f"Required executable not found on PATH: {exe}")


def open_maybe_gz(path: str):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "r")


def normalize_read_id(s: str) -> str:
    s = s.strip()
    if s.startswith("@"):
        s = s[1:]
    return s.split()[0]


def rank_prefix(rank: str) -> Optional[str]:
    return {
        "D": "d__",
        "K": "k__",
        "P": "p__",
        "C": "c__",
        "O": "o__",
        "F": "f__",
        "G": "g__",
        "S": "s__",
    }.get(rank)


def run(cmd, *, cwd=None) -> None:
    print("[CMD]", " ".join(cmd))
    subprocess.run(cmd, cwd=cwd, check=True)


def shlex_quote(s: str) -> str:
    return "'" + s.replace("'", "'\"'\"'") + "'"


# ----------------- lineage summariser -----------------

def parse_report_build_label_map(report_path: str) -> Dict[str, str]:
    taxid_to_lineage: Dict[str, str] = {}
    lineage_by_rank: Dict[str, str] = {}

    with open_maybe_gz(report_path) as f:
        for line in f:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 6:
                continue

            rank = parts[3]
            taxid = parts[4]
            name = parts[5].strip()

            pref = rank_prefix(rank)
            if pref:
                lineage_by_rank[rank] = f"{pref}{name}"

            if rank in ("G", "S"):
                ordered = ["D", "K", "P", "C", "O", "F", "G"]
                if rank == "S":
                    ordered.append("S")
                lineage_str = ";".join([lineage_by_rank[r] for r in ordered if r in lineage_by_rank])
                taxid_to_lineage[taxid] = lineage_str

    return taxid_to_lineage


def load_kraken_read_map(kraken_path: str, classified_only: bool) -> Dict[str, str]:
    read_to_taxid: Dict[str, str] = {}
    with open_maybe_gz(kraken_path) as f:
        for line in f:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            status = parts[0]
            if classified_only and status != "C":
                continue
            read_id = normalize_read_id(parts[1])
            taxid = parts[2]
            read_to_taxid[read_id] = taxid
    return read_to_taxid


def gc_fraction(seq: str) -> float:
    seq = seq.strip().upper()
    if not seq:
        return 0.0
    gc = seq.count("G") + seq.count("C")
    return gc / len(seq)


def process_fastq(fastq_path, read_to_taxid, taxid_to_lineage, file_key, counts, gc_sum, gc_n):
    with open_maybe_gz(fastq_path) as f:
        while True:
            h = f.readline()
            if not h:
                break
            seq = f.readline()
            f.readline()
            qual = f.readline()
            if not qual:
                break

            read_id = normalize_read_id(h)
            taxid = read_to_taxid.get(read_id)
            lineage = taxid_to_lineage.get(taxid)
            if not lineage:
                continue

            counts[lineage]["total"] += 1
            counts[lineage][file_key] += 1
            gc_sum[lineage] += gc_fraction(seq)
            gc_n[lineage] += 1


def write_attrition_csv(report, kraken, mapped_fastq, unmapped_fastq, out_csv, classified_only, min_reads):
    taxid_to_lineage = parse_report_build_label_map(report)
    read_to_taxid = load_kraken_read_map(kraken, classified_only)

    counts = defaultdict(lambda: {"total": 0, "mapped": 0, "unmapped": 0})
    gc_sum = defaultdict(float)
    gc_n = defaultdict(int)

    process_fastq(mapped_fastq, read_to_taxid, taxid_to_lineage, "mapped", counts, gc_sum, gc_n)
    process_fastq(unmapped_fastq, read_to_taxid, taxid_to_lineage, "unmapped", counts, gc_sum, gc_n)

    with open(out_csv, "w", newline="") as out:
        w = csv.writer(out)
        w.writerow(["Lineage", "Total Count", "Mapped Read Proportion", "Unmapped Read Proportion", "Average GC Content"])

        for lineage, c in sorted(counts.items(), key=lambda kv: kv[1]["total"], reverse=True):
            if c["total"] < min_reads:
                continue
            avg_gc = gc_sum[lineage] / gc_n[lineage] if gc_n[lineage] else 0
            w.writerow([lineage, c["total"], c["mapped"]/c["total"], c["unmapped"]/c["total"], avg_gc])


# ----------------- bin reference caching -----------------

def bins_cache_key(bins_dir: Path) -> str:
    h = hashlib.sha256()
    for p in sorted(bins_dir.glob("*")):
        if p.is_file():
            st = p.stat()
            h.update(p.name.encode())
            h.update(str(st.st_size).encode())
            h.update(str(int(st.st_mtime)).encode())
    return h.hexdigest()[:16]


def build_or_get_binned_reference(bins_dir: Path, cache_dir: Path) -> Path:
    cache_dir.mkdir(exist_ok=True)
    ref = cache_dir / f"binned_{bins_cache_key(bins_dir)}.fasta"

    if ref.exists():
        return ref

    with open(ref, "w") as out:
        for bf in sorted(bins_dir.glob("*.fa*")):
            with open(bf) as f:
                shutil.copyfileobj(f, out)

    return ref


# ----------------- mapping + separation (gzipped FASTQs) -----------------

def map_and_separate(reference, r1, r2, outdir, threads, keep):
    outdir.mkdir(parents=True, exist_ok=True)

    if not Path(str(reference) + ".bwt").exists():
        run(["bwa", "index", str(reference)])

    bam = outdir / "aln.bam"
    sorted_bam = outdir / "aln.sorted.bam"

    run([
        "bash", "-lc",
        f"bwa mem -t {threads} {shlex_quote(str(reference))} {shlex_quote(str(r1))} {shlex_quote(str(r2))} "
        f"| samtools view -@ {threads} -bS - > {shlex_quote(str(bam))}"
    ])

    run(["samtools", "sort", "-@", str(threads), str(bam), "-o", str(sorted_bam)])

    mapped = outdir / "mapped.bam"
    unmapped = outdir / "unmapped.bam"

    run(["samtools", "view", "-b", "-F", "12", str(sorted_bam), "-o", str(mapped)])
    run(["samtools", "view", "-b", "-f", "12", str(sorted_bam), "-o", str(unmapped)])

    mapped_name = outdir / "mapped.name.bam"
    unmapped_name = outdir / "unmapped.name.bam"

    run(["samtools", "sort", "-@", str(threads), "-n", str(mapped), "-o", str(mapped_name)])
    run(["samtools", "sort", "-@", str(threads), "-n", str(unmapped), "-o", str(unmapped_name)])

    run([
        "bash", "-lc",
        f"samtools fastq -@ {threads} -n "
        f"-1 {shlex_quote(str(outdir/'mapped_R1.fq.gz'))} "
        f"-2 {shlex_quote(str(outdir/'mapped_R2.fq.gz'))} "
        f"-0 /dev/null -s /dev/null "
        f"{shlex_quote(str(mapped_name))}"
    ])

    run([
        "bash", "-lc",
        f"samtools fastq -@ {threads} -n "
        f"-1 {shlex_quote(str(outdir/'unmapped_R1.fq.gz'))} "
        f"-2 {shlex_quote(str(outdir/'unmapped_R2.fq.gz'))} "
        f"-0 /dev/null -s /dev/null "
        f"{shlex_quote(str(unmapped_name))}"
    ])

    if not keep:
        for f in [bam, sorted_bam, mapped, unmapped, mapped_name, unmapped_name]:
            if f.exists():
                f.unlink()

    return outdir/'mapped_R1.fq.gz', outdir/'unmapped_R1.fq.gz'


# ----------------- main -----------------

def main():
    which_or_die("bwa")
    which_or_die("samtools")

    ap = argparse.ArgumentParser()
    ap.add_argument("--sheet", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--classified-only", action="store_true")
    ap.add_argument("--min-reads", type=int, default=5)
    ap.add_argument("--cache-dir", default=".cache_binned_refs")
    ap.add_argument("--keep-mapping", action="store_true")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(exist_ok=True)

    with open(args.sheet) as f:
        reader = csv.DictReader(f)
        for row in reader:
            sample = row["sample_id"]
            print(f"\n[INFO] Processing {sample}")

            ref = build_or_get_binned_reference(Path(row["bins_dir"]), Path(args.cache_dir))

            mapped_fq, unmapped_fq = map_and_separate(
                ref,
                Path(row["r1_fastq"]),
                Path(row["r2_fastq"]),
                outdir/sample/"mapping",
                args.threads,
                args.keep_mapping
            )

            write_attrition_csv(
                row["kraken2_report_tsv"],
                row["kraken2_tsv"],
                mapped_fq,
                unmapped_fq,
                outdir/sample/f"{sample}_attrition_lineage.csv",
                args.classified_only,
                args.min_reads
            )


if __name__ == "__main__":
    main()