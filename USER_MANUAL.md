# UltraQ30 — User Manual

**Version 1.0**  
© 2026 Bluelynx

---

## 1. Overview

**UltraQ30** is a high-performance FASTQ quality analysis engine written in Swift.  
It is designed for speed, scalability, and deterministic results, using:

- SIMD acceleration
- Multithreaded processing
- Streaming file IO
- Minimal memory overhead

UltraQ30 can be used as:
- A command-line tool
- The core engine of an iOS/macOS application
- A Swift Package embedded into other bioinformatics workflows

The engine is optimized for modern Apple silicon CPUs.

---

## 2. Supported Input

### 2.1 File formats
- `.fastq`
- `.fq`

Only **plain (uncompressed)** FASTQ files are supported.

---

### 2.2 FASTQ format requirements

UltraQ30 expects standard **4-line FASTQ records**:

