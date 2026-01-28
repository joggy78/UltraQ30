/*
 * Copyright (c) 2026 Bluelynx
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

import Foundation
import simd

// =====================================================
// MARK: - Options (AGRESSIF)
// =====================================================

struct UltraQ30Options {
    // agressif par défaut
    var chunkSize: Int = 8_388_608 // 8MB
    var threads: Int = max(1, ProcessInfo.processInfo.activeProcessorCount)
    var batchReads: Int = 32_768   // gros batch

    // seed pour estimer la longueur moyenne (pour info uniquement)
    var seedReadsForMean: Int = 20_000

    // toggles
    var validateSeq: Bool = true
    var doGC: Bool = true
    var doDup: Bool = true

    // output
    var extraColumns: Bool = false
}

// =====================================================
// MARK: - Hash128 (pour duplications exactes)
// =====================================================

struct Hash128: Hashable {
    let a: UInt64
    let b: UInt64
}

// FNV-1a 64-bit paramétré (seed custom)
@inline(__always)
func fnv1a64(_ ptr: UnsafeBufferPointer<UInt8>, seed: UInt64) -> UInt64 {
    var h: UInt64 = seed
    for i in 0..<ptr.count {
        h ^= UInt64(ptr[i])
        h &*= 1099511628211
    }
    return h
}

// Hash "128-bit" = 2 hashes 64-bit indépendants
@inline(__always)
func hashFullRead128(_ seq: UnsafeBufferPointer<UInt8>) -> Hash128 {
    // 2 seeds différentes (offset basis standard + une variante)
    let h1 = fnv1a64(seq, seed: 14695981039346656037)
    let h2 = fnv1a64(seq, seed: 1099511628211 ^ 0x9E3779B97F4A7C15)
    return Hash128(a: h1, b: h2)
}

// =====================================================
// MARK: - Stats
// =====================================================

struct FastqStats {
    var reads = 0
    var bases = 0
    var q30Bases = 0
    var gcBases = 0

    var minReadLen = Int.max
    var maxReadLen = 0
    var sumReadLen: Int64 = 0

    // duplications exactes: set des uniques (par worker puis union)
    var uniques: Set<Hash128>? = nil

    mutating func enableDupSet() {
        self.uniques = Set<Hash128>()
    }

    mutating func addRead(len: Int) {
        reads += 1
        bases += len
        sumReadLen += Int64(len)

        if len < minReadLen { minReadLen = len }
        if len > maxReadLen { maxReadLen = len }
    }

    mutating func merge(_ other: FastqStats) {
        reads += other.reads
        bases += other.bases
        q30Bases += other.q30Bases
        gcBases += other.gcBases

        minReadLen = min(minReadLen, other.minReadLen)
        maxReadLen = max(maxReadLen, other.maxReadLen)
        sumReadLen += other.sumReadLen

        if var mine = uniques, let oth = other.uniques {
            // union exact
            mine.formUnion(oth)
            uniques = mine
        } else if uniques == nil, let oth = other.uniques {
            uniques = oth
        }
    }

    var meanReadLength: Double {
        reads > 0 ? Double(sumReadLen) / Double(reads) : 0
    }

    var q30Percent: Double {
        bases > 0 ? Double(q30Bases) / Double(bases) * 100.0 : 0
    }

    var gcPercent: Double {
        bases > 0 ? Double(gcBases) / Double(bases) * 100.0 : 0
    }

    var uniqueExact: Int? {
        uniques?.count
    }

    var dupPercentExact: Double? {
        guard let u = uniqueExact else { return nil }
        guard reads > 0 else { return 0 }
        let uniq = min(reads, max(0, u))
        return max(0.0, (Double(reads - uniq) / Double(reads)) * 100.0)
    }
}

// =====================================================
// MARK: - Safe fast helpers
// =====================================================

@inline(__always)
func isValidSequence(_ ptr: UnsafeBufferPointer<UInt8>) -> Bool {
    for b in ptr {
        if !(b == 65 || b == 67 || b == 71 || b == 84 || b == 78) {
            return false // A C G T N
        }
    }
    return true
}

// Unaligned-safe SIMD for Q30 (ASCII-33: Q30 => char >= 63)
@inline(__always)
func countQ30UnalignedSIMD(_ ptr: UnsafeBufferPointer<UInt8>) -> Int {
    let threshold = SIMD16<UInt8>(repeating: 63)
    var count = 0
    var i = 0

    while i + 16 <= ptr.count {
        let v: SIMD16<UInt8> = UnsafeRawPointer(ptr.baseAddress!.advanced(by: i))
            .loadUnaligned(as: SIMD16<UInt8>.self)
        let mask = (v .>= threshold)
        for lane in 0..<16 where mask[lane] {
            count += 1
        }
        i += 16
    }

    while i < ptr.count {
        if ptr[i] >= 63 { count += 1 }
        i += 1
    }
    return count
}

// Unaligned-safe SIMD for GC (sur toute la séquence)
@inline(__always)
func countGCUnalignedSIMD(_ ptr: UnsafeBufferPointer<UInt8>) -> Int {
    let g = SIMD16<UInt8>(repeating: 71) // 'G'
    let c = SIMD16<UInt8>(repeating: 67) // 'C'
    var count = 0
    var i = 0

    while i + 16 <= ptr.count {
        let v: SIMD16<UInt8> = UnsafeRawPointer(ptr.baseAddress!.advanced(by: i))
            .loadUnaligned(as: SIMD16<UInt8>.self)
        let mg = (v .== g)
        let mc = (v .== c)
        for lane in 0..<16 where (mg[lane] || mc[lane]) {
            count += 1
        }
        i += 16
    }

    while i < ptr.count {
        let b = ptr[i]
        if b == 71 || b == 67 { count += 1 }
        i += 1
    }
    return count
}

// =====================================================
// MARK: - Streaming line reader
// =====================================================

@inline(__always)
func readLine(from buffer: inout [UInt8], startIndex: inout Int) -> ArraySlice<UInt8>? {
    if startIndex >= buffer.count { return nil }
    var i = startIndex
    while i < buffer.count {
        if buffer[i] == 10 { // '\n'
            let line = buffer[startIndex..<i]
            startIndex = i + 1
            return line
        }
        i += 1
    }
    return nil
}

@inline(__always)
func compactBuffer(_ buffer: inout [UInt8], consumed: Int) {
    guard consumed > 0 else { return }
    if consumed >= buffer.count {
        buffer.removeAll(keepingCapacity: true)
        return
    }
    buffer.removeFirst(consumed)
}

// =====================================================
// MARK: - Packed batch
// NOTE: On packe toujours la QUAL.
// On packe la SEQ seulement si besoin (GC / DUP / VALIDATE).
// =====================================================

struct PackedBatch {
    var data: [UInt8] = [] // qual always

    var qualOff: [Int] = []
    var qualLen: [Int] = []

    // seq optional
    var hasSeq: Bool = true
    var seqOff: [Int] = []
    var seqLen: [Int] = []

    mutating func reserve(reads: Int, bytes: Int, hasSeq: Bool) {
        self.hasSeq = hasSeq
        data.reserveCapacity(bytes)
        qualOff.reserveCapacity(reads)
        qualLen.reserveCapacity(reads)
        if hasSeq {
            seqOff.reserveCapacity(reads)
            seqLen.reserveCapacity(reads)
        }
    }

    var count: Int { qualOff.count }
}

// =====================================================
// MARK: - Pass 1: estimate mean read length
// (needs SEQ lines, so only depends on file not flags)
// =====================================================

func estimateMeanReadLength(file: URL, opt: UltraQ30Options) throws -> Int {
    let h = try FileHandle(forReadingFrom: file)
    defer { try? h.close() }

    var buffer: [UInt8] = []
    buffer.reserveCapacity(opt.chunkSize / 2)

    var startIndex = 0
    var state = 0
    var seqLine: ArraySlice<UInt8>? = nil

    var seen = 0
    var sum = 0

    while seen < opt.seedReadsForMean {
        let chunk = try h.read(upToCount: min(opt.chunkSize, 1_048_576)) ?? Data()
        if chunk.isEmpty { break }

        buffer.append(contentsOf: chunk)

        while seen < opt.seedReadsForMean {
            guard let line = readLine(from: &buffer, startIndex: &startIndex) else { break }

            switch state {
            case 0:
                state = 1 // header
            case 1:
                seqLine = line
                state = 2
            case 2:
                state = 3 // '+'
            default:
                if let s = seqLine {
                    sum += s.count
                    seen += 1
                }
                seqLine = nil
                state = 0
            }
        }

        if startIndex > 0 {
            compactBuffer(&buffer, consumed: startIndex)
            startIndex = 0
        }
    }

    if seen == 0 { return 0 }
    return max(0, Int(Double(sum) / Double(seen)))
}

// =====================================================
// MARK: - Pass 2: Multicore processing
// =====================================================

func processFastqFileMulticore(
    file: URL,
    opt: UltraQ30Options,
    progress: ((Double) -> Void)? = nil
) throws -> FastqStats {
    let h = try FileHandle(forReadingFrom: file)
    defer { try? h.close() }

    let attrs = try FileManager.default.attributesOfItem(atPath: file.path)
    let fileSize = (attrs[.size] as? NSNumber)?.doubleValue ?? 0

    let workers = max(1, opt.threads)
    let workerQ = DispatchQueue(label: "ultraQ30.work", attributes: .concurrent)
    let group = DispatchGroup()

    // do we need seq bytes at all?
    let needSeq = opt.doGC || opt.doDup || opt.validateSeq

    var workerStats = (0..<workers).map { _ -> FastqStats in
        var s = FastqStats()
        if opt.doDup { s.enableDupSet() }
        return s
    }

    let mergeQueues = (0..<workers).map { DispatchQueue(label: "ultraQ30.merge.\($0)") }
    var nextWorker = 0

    func submit(_ batch: PackedBatch) {
        if batch.count == 0 { return }

        let idx = nextWorker
        nextWorker = (nextWorker + 1) % workers

        group.enter()
        workerQ.async {
            var local = FastqStats()
            if opt.doDup { local.enableDupSet() }

            batch.data.withUnsafeBytes { raw in
                let base = raw.bindMemory(to: UInt8.self).baseAddress!

                for i in 0..<batch.count {
                    let qo = batch.qualOff[i]
                    let ql = batch.qualLen[i]
                    if ql <= 0 { continue }

                    let qPtr = UnsafeBufferPointer(start: base.advanced(by: qo), count: ql)
                    let readLen = ql

                    var sPtr: UnsafeBufferPointer<UInt8>? = nil
                    if batch.hasSeq {
                        let so = batch.seqOff[i]
                        let sl = batch.seqLen[i]
                        if sl != ql || sl <= 0 { continue }
                        sPtr = UnsafeBufferPointer(start: base.advanced(by: so), count: sl)
                    }

                    if opt.validateSeq, let s = sPtr {
                        if !isValidSequence(s) { continue }
                    }

                    local.addRead(len: readLen)
                    local.q30Bases += countQ30UnalignedSIMD(qPtr)

                    if opt.doGC, let s = sPtr {
                        // GC full-length
                        local.gcBases += countGCUnalignedSIMD(s)
                    }

                    if opt.doDup, let s = sPtr {
                        // DUP full-length (hash de toute la séquence)
                        let h128 = hashFullRead128(s)
                        local.uniques?.insert(h128)
                    }
                }
            }

            mergeQueues[idx].sync {
                workerStats[idx].merge(local)
            }

            group.leave()
        }
    }

    // Producer parsing
    var buffer: [UInt8] = []
    buffer.reserveCapacity(opt.chunkSize)

    var startIndex = 0
    var state = 0
    var seqLine: ArraySlice<UInt8>? = nil
    var qualLine: ArraySlice<UInt8>? = nil

    var bytesRead: Double = 0

    var batch = PackedBatch()
    let approxBytesPerRead = needSeq ? 2 * 160 : 160
    batch.reserve(reads: opt.batchReads, bytes: opt.batchReads * approxBytesPerRead, hasSeq: needSeq)

    func flushBatch() {
        if batch.count > 0 {
            submit(batch)
            batch = PackedBatch()
            batch.reserve(reads: opt.batchReads, bytes: opt.batchReads * approxBytesPerRead, hasSeq: needSeq)
        }
    }

    while true {
        let chunk = try h.read(upToCount: opt.chunkSize) ?? Data()
        if chunk.isEmpty { break }

        bytesRead += Double(chunk.count)
        buffer.append(contentsOf: chunk)

        while true {
            guard let line = readLine(from: &buffer, startIndex: &startIndex) else { break }

            switch state {
            case 0:
                state = 1 // header
            case 1:
                seqLine = line
                state = 2
            case 2:
                state = 3 // '+'
            default:
                qualLine = line
                state = 0

                guard let q = qualLine else {
                    seqLine = nil
                    qualLine = nil
                    continue
                }

                if needSeq {
                    guard let s = seqLine else {
                        seqLine = nil
                        qualLine = nil
                        continue
                    }

                    // pack seq
                    let so = batch.data.count
                    batch.data.append(contentsOf: s)
                    let sl = s.count

                    // pack qual
                    let qo = batch.data.count
                    batch.data.append(contentsOf: q)
                    let ql = q.count

                    batch.seqOff.append(so)
                    batch.seqLen.append(sl)

                    batch.qualOff.append(qo)
                    batch.qualLen.append(ql)
                } else {
                    // pack qual only
                    let qo = batch.data.count
                    batch.data.append(contentsOf: q)
                    let ql = q.count

                    batch.qualOff.append(qo)
                    batch.qualLen.append(ql)
                }

                if batch.count >= opt.batchReads {
                    flushBatch()
                }

                seqLine = nil
                qualLine = nil
            }
        }

        if startIndex > 0 {
            compactBuffer(&buffer, consumed: startIndex)
            startIndex = 0
        }

        if fileSize > 0 {
            progress?(min(1.0, bytesRead / fileSize))
        }
    }

    flushBatch()
    group.wait()
    progress?(1.0)

    var total = FastqStats()
    if opt.doDup { total.enableDupSet() }
    for s in workerStats { total.merge(s) }

    if total.reads == 0 { total.minReadLen = 0 }
    return total
}

// =====================================================
// MARK: - IO helpers
// =====================================================

func listFastqFiles(in input: URL) throws -> [URL] {
    var isDir: ObjCBool = false
    guard FileManager.default.fileExists(atPath: input.path, isDirectory: &isDir) else { return [] }

    if isDir.boolValue {
        return try FileManager.default.contentsOfDirectory(
            at: input,
            includingPropertiesForKeys: nil,
            options: [.skipsHiddenFiles]
        )
        .filter { ["fastq", "fq"].contains($0.pathExtension.lowercased()) }
        .sorted { $0.lastPathComponent < $1.lastPathComponent }
    } else {
        return [input]
    }
}

func safeBase(_ url: URL) -> String {
    url.deletingPathExtension().lastPathComponent
}

func writeText(_ text: String, to url: URL) throws {
    try FileManager.default.createDirectory(
        at: url.deletingLastPathComponent(),
        withIntermediateDirectories: true
    )
    try text.write(to: url, atomically: true, encoding: .utf8)
}

func reportHeader(extra: Bool) -> String {
    if extra {
        return "file\treads\tbases\tq30_bases\tq30_percent\tgc_bases\tgc_percent\tmin_read_len\tmax_read_len\tmean_read_len\tdup_percent\tunique_exact\tmean_len_seed\tthreads\tno_gc\tno_dup\tno_validate"
    } else {
        return "file\treads\tbases\tq30_bases\tq30_percent\tgc_bases\tgc_percent\tmin_read_len\tmax_read_len\tmean_read_len\tdup_percent\tunique_exact"
    }
}

func fmtNA(_ x: Double?) -> String {
    guard let v = x else { return "NA" }
    return String(format: "%.2f", v)
}

func reportLine(fileName: String, s: FastqStats, meanSeed: Int, opt: UltraQ30Options) -> String {
    let gcPctStr: String = opt.doGC ? String(format: "%.2f", s.gcPercent) : "NA"
    let gcBasesStr: String = opt.doGC ? "\(s.gcBases)" : "NA"

    let dupStr: String = fmtNA(s.dupPercentExact)
    let uniqStr: String = s.uniqueExact.map(String.init) ?? "NA"

    if opt.extraColumns {
        return "\(fileName)\t\(s.reads)\t\(s.bases)\t\(s.q30Bases)\t\(String(format: "%.2f", s.q30Percent))\t\(gcBasesStr)\t\(gcPctStr)\t\(s.minReadLen)\t\(s.maxReadLen)\t\(String(format: "%.2f", s.meanReadLength))\t\(dupStr)\t\(uniqStr)\t\(meanSeed)\t\(opt.threads)\t\((opt.doGC ? 0 : 1))\t\((opt.doDup ? 0 : 1))\t\((opt.validateSeq ? 0 : 1))"
    } else {
        return "\(fileName)\t\(s.reads)\t\(s.bases)\t\(s.q30Bases)\t\(String(format: "%.2f", s.q30Percent))\t\(gcBasesStr)\t\(gcPctStr)\t\(s.minReadLen)\t\(s.maxReadLen)\t\(String(format: "%.2f", s.meanReadLength))\t\(dupStr)\t\(uniqStr)"
    }
}

func defaultOutDir() -> URL {
    FileManager.default.temporaryDirectory.appendingPathComponent("ultraQ30_results", isDirectory: true)
}

func printProgress(file: String, p: Double) {
    print(String(format: "\r[%@] %.2f%%", file, p * 100.0), terminator: "")
    fflush(stdout)
}

func stderr(_ s: String) {
    fputs(s + "\n", stderr)
}

func usage() {
    print("""
    Usage: ultraQ30 <folder|fastq> [--out-dir DIR] [--threads N] [--batch READS] [--chunk BYTES]
                  [--seed-reads N] [--no-gc] [--no-dup] [--no-validate] [--extra]

    Defaults (agressif):
      --chunk = 8388608 (8MB)
      --batch = 32768
    """)
}

// =====================================================
// MARK: - Main
// =====================================================
let argv = CommandLine.arguments

// ✅ help accessible sans input
if argv.contains("--help") || argv.contains("-h") {
    usage()
    exit(0)
}

guard argv.count >= 2 else {
    usage()
    exit(1)
}

let inputURL = URL(fileURLWithPath: argv[1])

var opt = UltraQ30Options()
var outDir: URL? = nil


var i = 2
while i < argv.count {
    switch argv[i] {
    case "--out-dir":
        i += 1
        if i < argv.count {
            outDir = URL(fileURLWithPath: argv[i], isDirectory: true)
        }

    case "--threads":
        i += 1
        if i < argv.count {
            opt.threads = max(1, Int(argv[i]) ?? opt.threads)
        }

    case "--batch":
        i += 1
        if i < argv.count {
            opt.batchReads = max(1024, Int(argv[i]) ?? opt.batchReads)
        }

    case "--chunk":
        i += 1
        if i < argv.count {
            opt.chunkSize = max(256 * 1024, Int(argv[i]) ?? opt.chunkSize)
        }

    case "--seed-reads":
        i += 1
        if i < argv.count {
            opt.seedReadsForMean = max(100, Int(argv[i]) ?? opt.seedReadsForMean)
        }

    case "--no-gc":
        opt.doGC = false

    case "--no-dup":
        opt.doDup = false

    case "--no-validate":
        opt.validateSeq = false

    case "--extra":
        opt.extraColumns = true

    case "--help", "-h":
        usage()
        exit(0)

    default:
        break
    }
    i += 1
}

let outputDir = outDir ?? defaultOutDir()

do {
    let files = try listFastqFiles(in: inputURL)
    if files.isEmpty {
        stderr("❌ No .fastq/.fq found at: \(inputURL.path)")
        exit(2)
    }

    try FileManager.default.createDirectory(at: outputDir, withIntermediateDirectories: true)

    var summary: [String] = [reportHeader(extra: opt.extraColumns)]
    summary.reserveCapacity(files.count + 1)

    for (idx, f) in files.enumerated() {
        print("\n(\(idx + 1)/\(files.count)) \(f.lastPathComponent)")

        let meanSeed = try estimateMeanReadLength(file: f, opt: opt)
        print(" meanLen(seed=\(opt.seedReadsForMean))=\(meanSeed), threads=\(opt.threads), noGC=\(!opt.doGC), noDup=\(!opt.doDup), noValidate=\(!opt.validateSeq)")

        let s = try processFastqFileMulticore(file: f, opt: opt) { p in
            printProgress(file: f.lastPathComponent, p: p)
        }
        print("")

        let outFile = outputDir.appendingPathComponent("\(safeBase(f)).ultraQ30.tsv")
        let content = reportHeader(extra: opt.extraColumns) + "\n"
            + reportLine(fileName: f.lastPathComponent, s: s, meanSeed: meanSeed, opt: opt)
            + "\n"
        try writeText(content, to: outFile)

        summary.append(reportLine(fileName: f.lastPathComponent, s: s, meanSeed: meanSeed, opt: opt))
    }

    let summaryURL = outputDir.appendingPathComponent("summary.tsv")
    try writeText(summary.joined(separator: "\n") + "\n", to: summaryURL)

    print("\nDone.")
    print("Results folder: \(outputDir.path)")
    print("Summary: \(summaryURL.path)")
} catch {
    stderr("❌ Error: \(error)")
    exit(3)
}

