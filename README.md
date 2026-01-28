# UltraQ30

**UltraQ30** is a high-performance Swift command-line tool for analyzing FASTQ files.  
It computes Q30 statistics, GC content, read length metrics, and exact duplication rates using **streaming I/O**, **SIMD**, and **multi-core processing**.

Designed for speed and low memory overhead on large sequencing datasets.

---

## Features

- FASTQ (`.fastq`, `.fq`) support — single file or entire directory
- Metrics per file:
  - Total reads and bases
  - Q30 bases and Q30 percentage (Phred+33)
  - GC bases and GC percentage *(optional)*
  - Min / max / mean read length
  - Exact duplicate rate using 128-bit hashing *(optional)*
- Optimized for performance:
  - Streaming file parsing (no full file in memory)
  - Large batch processing (default: 32k reads)
  - Multi-threaded processing (DispatchQueue)
  - SIMD-accelerated Q30 and GC counting
  - Unaligned-safe SIMD loads
- TSV output per file + global summary

---

## How It Works

UltraQ30 processes FASTQ files in **two passes**:

1. **Mean read length estimation (seed pass)**  
   Reads a limited number of records to estimate average read length (informational only).

2. **Main analysis pass**
   - FASTQ records are parsed in streaming mode
   - Reads are packed into large batches
   - Batches are distributed across worker threads
   - Each worker computes local statistics
   - Results are merged into final per-file statistics

Exact duplicates are detected using a **128-bit hash** (two independent FNV-1a 64-bit hashes) of the full sequence.

---

## Installation

Requires **Swift 5.9+** and macOS or Linux.


git clone https://github.com/<your-org>/ultraQ30.git
cd ultraQ30
swift build -c release

---


## UltraQ30 basic usage

ultraQ30 <folder|fastq>
  [--out-dir DIR] [--threads N] [--batch READS] [--chunk BYTES]
  [--seed-reads N] [--no-gc] [--no-dup] [--no-validate] [--extra]

---


## Performance Notes

Uses SIMD (SIMD16<UInt8>) for Q30 and GC counting
Uses unaligned memory loads for safe streaming buffers
Avoids unnecessary sequence storage when GC / duplication checks are disabled
Designed for large FASTQ files (tens to hundreds of GB)

---

## Limitations

Exact duplicate detection uses hashing (not cryptographic)
Duplicate sets may consume significant memory on extremely large datasets
Assumes standard FASTQ format with Phred+33 encoding


