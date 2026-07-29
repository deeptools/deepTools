use crate::filtering::Alignmentfilters;
use core::panic;
use ndarray::Array1;
use rust_htslib::bam::ext::BamRecordExtensions;
use rust_htslib::bam::{IndexedReader, Read};
use std::cmp::{max, min};
use std::collections::HashMap;
use std::fmt;
use std::fs::OpenOptions;
use std::io::{BufWriter, Write};
use tempfile::{Builder, TempPath};

pub fn parse_regions(region: &str, bam_ifile: Vec<&str>) -> (Vec<Region>, HashMap<String, u32>) {
    // Takes a vector of regions, and a bam reference
    // returns a vector of regions, with all chromosomes and full lengths if original regions was empty
    // Else it validates the regions against the information from the bam header
    // Finally, a Vec with chromsizes is returned as well.
    let mut found_chroms: HashMap<String, usize> = HashMap::new();
    for bam in bam_ifile.iter() {
        let bam = IndexedReader::from_path(bam).unwrap();
        let chroms: Vec<String> = bam
            .header()
            .target_names()
            .iter()
            .map(|x| String::from_utf8(x.to_vec()).unwrap())
            .collect();
        for chrom in chroms.iter() {
            // if it's not in the hashmap, add it, else increment count
            if !found_chroms.contains_key(chrom) {
                found_chroms.insert(chrom.clone(), 1);
            } else {
                let count = found_chroms.get_mut(chrom).unwrap();
                *count += 1;
            }
        }
    }
    let mut validchroms: Vec<String> = Vec::new();
    // loop over all chroms in the hashmap, if the count is expected, include them
    for (chrom, count) in found_chroms.iter() {
        if *count == bam_ifile.len() {
            validchroms.push(chrom.clone());
        } else {
            println!(
                "Chromosome {} is missing in at least one bam file, and thus ignored!",
                chrom
            );
        }
    }
    // Crash if validchroms is empty.
    assert!(
        !validchroms.is_empty(),
        "No chromosomes found that are present in all bam files. Did you mix references ?"
    );
    // Read header from first bam file
    let bam = IndexedReader::from_path(bam_ifile[0]).unwrap();
    let header = bam.header().clone();
    let mut chromregions: Vec<Region> = Vec::new();
    let mut chromsizes = HashMap::new();
    if region == "None" {
        // if regions is empty, we default to all chromosomes, full length
        for tid in 0..header.target_count() {
            let chromname = String::from_utf8(header.tid2name(tid).to_vec())
                .expect("Invalid UTF-8 in chromosome name");
            let chromlen = header
                .target_len(tid)
                .expect("Error retrieving length for chromosome");
            // If chromname is not in validchroms, skip it.
            if !validchroms.contains(&chromname) {
                continue;
            }
            let _reg = Region {
                chrom: chromname.clone(),
                start: Revalue::U(0),
                end: Revalue::U(chromlen as u32),
                score: String::from("."),
                strand: String::from("."),
                name: format!("{}:{}-{}", chromname, 0, chromlen),
                regionlength: chromlen as u32,
            };
            chromregions.push(_reg);
            chromsizes.insert(chromname.to_string(), chromlen as u32);
        }
    } else {
        // populate chromsizes
        for tid in 0..header.target_count() {
            let chromname = String::from_utf8(header.tid2name(tid).to_vec())
                .expect("Invalid UTF-8 in chromosome name");
            let chromlen = header
                .target_len(tid)
                .expect("Error retrieving length for chromosome");
            if validchroms.contains(&chromname) {
                chromsizes.insert(chromname, chromlen as u32);
            }
        }
        // Either the region is just a chromosome, or it is a 'true region' i.e. chr:start:end
        // first test if the region is a chromosome alone
        let parts = region.split(":").collect::<Vec<&str>>();
        if parts.len() == 1 {
            let chromname = parts[0].to_string();
            // Check if chromname is in validchroms
            assert!(
                validchroms.contains(&chromname),
                "Supplied chromosome {} is not found.",
                chromname
            );
            let chromlen = chromsizes.get(&chromname).unwrap();
            let _reg = Region {
                chrom: chromname.to_string(),
                start: Revalue::U(0),
                end: Revalue::U(*chromlen),
                score: String::from("."),
                strand: String::from("."),
                name: format!("{}:{}-{}", chromname, 0, chromlen),
                regionlength: *chromlen,
            };
            chromregions.push(_reg);
        } else {
            // We have a region, split it into chrom, start, end
            let chromname = parts[0].to_string();
            let start = parts[1]
                .parse::<u32>()
                .expect("Error reading supplied start position.");
            let end = parts[2]
                .parse::<u32>()
                .expect("Error reading supplied end position.");
            // Check if chromname is in validchroms
            assert!(
                validchroms.contains(&chromname),
                "Supplied chromosome {} is not found.",
                chromname
            );
            let chromlen = chromsizes.get(&chromname).unwrap();
            assert!(
                end <= *chromlen,
                "Suplied region end goes beyond chromosome boundary. {} > {}",
                end,
                chromlen
            );
            let _reg = Region {
                chrom: chromname.to_string(),
                start: Revalue::U(start),
                end: Revalue::U(end),
                score: String::from("."),
                strand: String::from("."),
                name: format!("{}:{}-{}", chromname, start, end),
                regionlength: end - start,
            };
            chromregions.push(_reg);
        }
    }
    // Sort regions to make our live easier down the line (and to have valid bigwigs written.)
    // Sort Vec of Regions per chromosome, and then by start.
    chromregions.sort_by(|a, b| {
        a.chrom
            .cmp(&b.chrom)
            .then(a.get_startu().cmp(&b.get_startu()))
    });
    return (chromregions, chromsizes);
}

/// Main workhorse for bamCoverage and bamCompare
/// Calculates coverage either per bp (bs = 1) or over bins (bs > 1)
#[allow(unused_assignments)]
#[allow(unused_variables)]
#[allow(unused_mut)]
pub fn bam_pileup<'a>(
    bam_ifile: &str,            // the bam file to consider.
    regionvec: &'a Vec<Region>, // Vector with regions to calculate coverage from.
    provbs: &u32, // Provisional bin size, in case we have gene mode later on this gets omitted.
    ispe: &bool,  // Wether or not the bam file is paired end or single end.
    ignorechr: &Vec<String>, // chromosomes to ignore for normalization calculation.
    filters: &Alignmentfilters, // a struct with filtering settings.
    collapse: bool, // collapse bins with same coverage, should be avoided in bamCompare scenario for example.
    gene_mode: bool, // Gene mode (BED / GTF) or not (bins)
    gather_lengths: bool, // Wether or not to collect frag/read lengths (blows up memory for big bam files in MBS, for example).
    smoothlength: u32,    // Distance in bp for center-weighted smoothing; 0 = no smoothing
) -> (
    Vec<TempPath>, // temp bedgraph file.
    u32,           // mapped reads
    u32,           // unmapped reads
    Vec<u32>,      // read lengths
    Vec<u32>,      // fragment lengths
) {
    let mut bs = *provbs;
    let mut binsize = &bs;
    // constant to check if read is first in pair (relevant later)
    const FREAD: u16 = 0x40;

    // init variables for mapping statistics and lengths
    let mut mapped_reads: u32 = 0;
    let mut unmapped_reads: u32 = 0;
    let mut readlens: Vec<u32> = Vec::new();
    let mut fraglens: Vec<u32> = Vec::new();

    // Create the output vector
    let bg = Builder::new()
        .prefix("deeptoolstmp_")
        .suffix(".bedgraph")
        .rand_bytes(12)
        .tempfile()
        .expect("Failed to create temporary file.");

    // Two cases: either the binsize is 1, or it is > 1.
    // Counting between the two modes is different. In binsize == 1 we compute pileups
    // for binsize > 1, we count the number of reads that overlap a bin.

    // Open BAM file once per thread, reuse across all regions
    let mut bam = IndexedReader::from_path(&bam_ifile).unwrap();

    for regstruct in regionvec.iter() {
        // There are two options here:
        // either we are supposed to calculate coverage over regions (variable binsize required) gene_mode = true
        // or we have a regular bin setting, gene_mode = false
        let mut region: (String, u32, u32);
        if gene_mode {
            region = (
                regstruct.chrom.clone(),
                regstruct.get_startu(),
                regstruct.get_endu(),
            );
            bs = region.2 - region.1;
            binsize = &bs;
        } else {
            region = (
                regstruct.chrom.clone(),
                regstruct.get_startu(),
                regstruct.get_endu(),
            );
        }
        bam.fetch((region.0.as_str(), region.1, region.2))
            .expect(&format!("Error fetching region: {:?}", region));
        let mut counts: Vec<u32>;
        let mut startstr: String = region.1.to_string();
        let mut endstr: String = region.2.to_string();
        if gene_mode {
            // It could be that we are in metagene mode, i.e. we only  want counts over exons
            // In this we need another iter - fetch per regstruct
            match (regstruct.start.clone(), regstruct.end.clone()) {
                (Revalue::U(start), Revalue::U(end)) => {
                    counts = vec![0u32; 1];
                    for record in bam.records() {
                        let mut record = record.expect("Error parsing record.");
                        if filters.filter {
                            if filters.filter_record(&record, &region.0.as_str()) {
                                continue;
                            }
                        }
                        if !ignorechr.contains(&region.0) {
                            if record.is_unmapped() {
                                unmapped_reads += 1;
                            } else {
                                mapped_reads += 1;
                                if *ispe {
                                    if record.is_paired()
                                        && record.is_proper_pair()
                                        && (record.flags() & FREAD != 0)
                                    {
                                        if gather_lengths {
                                            fraglens.push(record.insert_size().abs() as u32);
                                        }
                                    }
                                }
                                if gather_lengths {
                                    readlens.push(record.seq_len() as u32);
                                }
                            }
                        }
                        counts[0] += 1;
                    }
                }
                (Revalue::V(starts), Revalue::V(ends)) => {
                    // Make a string with the start values comma separated
                    startstr = starts
                        .iter()
                        .map(|x| x.to_string())
                        .collect::<Vec<String>>()
                        .join(",");
                    endstr = ends
                        .iter()
                        .map(|x| x.to_string())
                        .collect::<Vec<String>>()
                        .join(",");

                    counts = vec![0u32; 1];
                    let exons: Vec<(u32, u32)> = starts
                        .iter()
                        .zip(ends.iter())
                        .map(|(&s, &e)| (s, e))
                        .collect();
                    for exon in exons {
                        bam.fetch((regstruct.chrom.as_str(), exon.0, exon.1))
                            .expect(&format!(
                                "Error fetching region: {}:{},{}",
                                regstruct.chrom, exon.0, exon.1
                            ));
                        for record in bam.records() {
                            let mut record = record.expect("Error parsing record.");
                            if filters.filter {
                                if filters.filter_record(&record, region.0.as_str()) {
                                    continue;
                                }
                            }
                            if !ignorechr.contains(&region.0) {
                                if record.is_unmapped() {
                                    unmapped_reads += 1;
                                } else {
                                    mapped_reads += 1;
                                    if *ispe {
                                        if record.is_paired()
                                            && record.is_proper_pair()
                                            && (record.flags() & FREAD != 0)
                                        {
                                            if gather_lengths {
                                                fraglens.push(record.insert_size().abs() as u32);
                                            }
                                        }
                                    }
                                    if gather_lengths {
                                        readlens.push(record.seq_len() as u32);
                                    }
                                }
                            }
                            counts[0] += 1;
                        }
                    }
                }
                _ => panic!(
                    "Start and End are not either both u32, or Vecs. This means your regions file is ill-defined. Fix {}.",
                    regstruct.name
                ),
            }
        } else {
            // populate the bg vector with 0 counts over all bins
            counts = vec![0u32; (region.2 - region.1).div_ceil(*binsize) as usize];
            let mut bin_indices: Vec<usize> = Vec::with_capacity(512);

            for record in bam.records() {
                let mut record = record.expect("Error parsing record.");

                if filters.filter {
                    if filters.filter_record(&record, region.0.as_str()) {
                        continue;
                    }
                }
                bin_indices.clear();

                if filters.manipulate {
                    let manipulated_blockpos = filters.manipulate_record(&mut record);
                    if manipulated_blockpos.is_none() {
                        continue;
                    }
                    for &x in &manipulated_blockpos.unwrap() {
                        let ix = ((x - region.1) / binsize) as usize;
                        if ix < counts.len() {
                            bin_indices.push(ix);
                        }
                    }
                } else {
                    for x in record.aligned_blocks() {
                        if (x[1] as u32) >= region.1
                            && (x[1] as u32) <= region.2
                            && (x[0] as u32) >= region.1
                            && (x[0] as u32) <= region.2
                        {
                            for pos in x[0] as u32..x[1] as u32 {
                                let ix = ((pos - region.1) / binsize) as usize;
                                if ix < counts.len() {
                                    bin_indices.push(ix);
                                }
                            }
                        }
                    }
                }

                if *binsize == 1 {
                    // For binsize=1, each bp is its own bin — positions are unique, no dedup needed
                    for &ix in &bin_indices {
                        counts[ix] += 1;
                    }
                } else {
                    bin_indices.sort_unstable();
                    bin_indices.dedup();
                    for &ix in &bin_indices {
                        counts[ix] += 1;
                    }
                }

                if !ignorechr.contains(&region.0) && gather_lengths {
                    if record.is_unmapped() {
                        unmapped_reads += 1;
                    } else {
                        mapped_reads += 1;
                        if *ispe {
                            if record.is_paired()
                                && record.is_proper_pair()
                                && (record.flags() & FREAD != 0)
                            {
                                fraglens.push(record.insert_size().abs() as u32);
                            }
                        }
                        readlens.push(record.seq_len() as u32);
                    }
                }
            }
        }

        // Convert integer counts to f32, then apply smoothing if requested
        let fcounts: Vec<f32> = counts.into_iter().map(|c| c as f32).collect();
        let smoothed = if smoothlength > 0 && fcounts.len() > 1 {
            let n = fcounts.len();
            let smooth_tiles = (smoothlength / *binsize).max(1) as usize;
            let mut out = Vec::<f32>::with_capacity(n);
            let (tiles_left, tiles_right) = if smooth_tiles == 1 {
                (0, 0)
            } else {
                let side = (smooth_tiles - 1) as f64 / 2.0f64;
                let left = side.ceil() as usize;
                let right = side.floor() as usize + 1;
                (left, right)
            };
            for i in 0..n {
                let s = i.saturating_sub(tiles_left);
                let e = min(i + tiles_right, n);
                let sum: f32 = fcounts[s..e].iter().sum();
                out.push(sum / (e - s) as f32);
            }
            out
        } else {
            fcounts
        };

        // There are two scenarios:
        // bamCoverage mode -> we can collapse bins with same coverage (collapse = true)
        // bamCompare & others -> We cannot collapse the bins, yet. (collapse = false)
        // Note that collapse can also be passed as a CLI, for those that want that.
        let mut outbuf = Vec::with_capacity(smoothed.len().max(1) * 64);
        let mut push_line = |chrom: &str, s: u32, e: u32, v: f32| {
            use std::io::Write;
            writeln!(outbuf, "{}\t{}\t{}\t{}", chrom, s, e, v).unwrap();
        };
        if smoothed.len() == 1 {
            push_line(&region.0, 0, 0, smoothed[0]); // start/end overridden below for gene_mode
            outbuf.clear();
            writeln!(
                outbuf,
                "{}\t{}\t{}\t{}",
                region.0, startstr, endstr, smoothed[0]
            )
            .unwrap();
        } else {
            if collapse {
                let mut lcov = smoothed[0];
                let mut lstart = region.1;
                let mut lend = region.1 + binsize;
                let mut start = lstart;
                let mut end = lend;
                let mut bin: u32 = 0;

                for (ix, count) in smoothed.into_iter().skip(1).enumerate() {
                    bin = (ix + 1) as u32; // offset of 1 due to skip(1)
                    start = (bin * binsize) + region.1;
                    end = min(start + binsize, region.2);
                    if count != lcov {
                        push_line(&region.0, lstart, lend, lcov);
                        lstart = lend;
                        lcov = count;
                    }
                    lend = end;
                }
                push_line(&region.0, lstart, lend, lcov);
            } else {
                let mut start = region.1;
                let mut end = region.1 + binsize;
                push_line(&region.0, start, end, smoothed[0]);
                for (ix, count) in smoothed.into_iter().skip(1).enumerate() {
                    let bin = (ix + 1) as u32;
                    start = (bin * binsize) + region.1;
                    end = min(start + binsize, region.2);
                    push_line(&region.0, start, end, count);
                }
            }
        }
        if !outbuf.is_empty() {
            use std::io::Write;
            let file = OpenOptions::new()
                .append(true)
                .create(true)
                .open(&bg)
                .expect("Error opening tmp file.");
            let mut writer = BufWriter::new(file);
            writer.write_all(&outbuf).unwrap();
        }
    }
    let bgpath = bg.into_temp_path();
    let tmpvec = vec![bgpath];

    return (tmpvec, mapped_reads, unmapped_reads, readlens, fraglens);
}

#[derive(Clone, Debug)]
pub struct Region {
    pub chrom: String,
    pub start: Revalue,
    pub end: Revalue,
    pub score: String,
    pub strand: String,
    pub name: String,
    pub regionlength: u32,
}

impl Region {
    pub fn assert_end(&self, chromend: u32) {
        match &self.end {
            Revalue::U(end) => {
                assert!(
                    *end <= chromend,
                    "Region end goes beyond chromosome boundary. Fix {}. {} {} {} (chr end = {})",
                    self.name,
                    self.chrom,
                    self.start,
                    self.end,
                    chromend
                );
            }
            Revalue::V(ends) => {
                for end in ends.iter() {
                    assert!(
                        *end <= chromend,
                        "Region end goes beyond chromosome boundary. Fix {}. {} {} {} (chr end = {})",
                        self.name,
                        self.chrom,
                        self.start,
                        end,
                        chromend
                    );
                }
            }
        }
    }

    pub fn get_anchorpoint(&self, referencepoint: &str) -> u32 {
        // reference-point mode.
        // What is exactly returned depends on a couple of parameters
        // what is the referencepoint : TSS, center, TES
        // what is the strand: +, -, . (note . is assumed to be +)
        // depending on if we have exon blocks (start / end are Revalue V == Vectors) or not (start / end are Revalue U == u32's)
        match referencepoint {
            "TSS" => match self.strand.as_str() {
                "+" | "." => match &self.start {
                    Revalue::U(start) => *start,
                    Revalue::V(starts) => starts[0],
                },
                "-" => match &self.end {
                    Revalue::U(end) => *end,
                    Revalue::V(ends) => *ends.last().unwrap(),
                },
                _ => panic!(
                    "Strand should either be + or - or . {:?} is not supported.",
                    self.strand
                ),
            },
            "TES" => match self.strand.as_str() {
                "+" | "." => match &self.end {
                    Revalue::U(end) => *end,
                    Revalue::V(ends) => *ends.last().unwrap(),
                },
                "-" => match &self.start {
                    Revalue::U(start) => *start,
                    Revalue::V(starts) => starts[0],
                },
                _ => panic!(
                    "Strand should either be + or - or . {:?} is not supported.",
                    self.strand
                ),
            },
            "center" => {
                // Here + or - doesn't matter. It is important though if we have 'metagenes' or not.
                match (&self.start, &self.end) {
                    (Revalue::U(start), Revalue::U(end)) => (*start + *end) / 2,
                    (Revalue::V(starts), Revalue::V(ends)) => {
                        let exonlength: u32 =
                            starts.iter().zip(ends.iter()).map(|(s, e)| e - s).sum();
                        let middle = exonlength / 2;
                        let mut cumsum: u32 = 0;
                        for (s, e) in starts.iter().zip(ends.iter()) {
                            cumsum += e - s;
                            if cumsum >= middle {
                                return s + (middle - (cumsum - (e - s)));
                            }
                        }
                        panic!(
                            "Middle of region not found. Fix {}. {}:{}-{}",
                            self.name, self.chrom, self.start, self.end
                        )
                    }
                    _ => panic!(
                        "Start and End are not either both u32, or Vecs. This means your regions file is ill-defined. Fix {}. {}:{}-{}",
                        self.name, self.chrom, self.start, self.end
                    ),
                }
            }
            _ => panic!(
                "Reference should either be TSS, TES or center. {:?} is not supported.",
                referencepoint
            ),
        }
    }

    #[allow(unused_assignments)]
    #[allow(unused_mut)]
    pub fn get_anchor_bins(&self, scale_regions: &Scalingregions, chromend: u32) -> Vec<Bin> {
        // Given an anchorpoint, return a vector, start, end , middle
        // The order of the vector is always 5' -> 3', meaning 'increasing' for +/. regions, and 'decreasing' for - regions.
        // At this stage, two situations are possible:
        // - self.start / self.end are Revalue::U, meaning we are in 'non metagene' mode.
        // - self.start / self.end are Revalue::V, meaning we are in 'metagene' mode, and the bins returned are exon-aware.
        // We need a notion of bins that don't make sense (i.e. beyond chromosome boundaries). These are encoded as (0,0)

        let mut bins: Vec<Bin> = Vec::new();
        let mut bodybins: Vec<Bin> = Vec::new();

        // Define anchorpoints
        let anchorstart;
        let anchorstop;
        match scale_regions.mode.as_str() {
            "reference-point" => {
                anchorstart = self.get_anchorpoint(&scale_regions.referencepoint);
                anchorstop = anchorstart;
            }
            "scale-regions" => match (&self.start, &self.end) {
                (Revalue::U(start), Revalue::U(end)) => {
                    anchorstart = *start;
                    anchorstop = *end;
                }
                (Revalue::V(start), Revalue::V(end)) => {
                    anchorstart = *start.first().unwrap();
                    anchorstop = *end.last().unwrap();
                }
                _ => panic!(
                    "Start and End are not either both u32, or Vecs. This means your regions file is ill-defined. Fix {}.",
                    self.name
                ),
            },
            _ => panic!(
                "Mode should either be reference-point or scale-regions. {} is not supported.",
                scale_regions.mode
            ),
        }
        if scale_regions.mode != "reference-point" {
            // scale-regions mode. Assert
            assert!(
                scale_regions.regionbodylength != 0,
                "scale-regions mode, but regionbodylength is 0."
            );
            if self.regionlength - (scale_regions.unscaled5prime + scale_regions.unscaled3prime)
                < scale_regions.binsize
            {
                println!(
                    "Warning ! Region {} is shorter than the binsize (potentially after unscaled regions taken into account. Whole region encoded as 0 or NA",
                    self.name
                );
                let nbin = scale_regions.cols_expected / scale_regions.bwfiles;
                for _ in 0..nbin {
                    bins.push(Bin::Conbin(0, 0));
                }
                return bins;
            } else {
                bodybins.extend(self.scale_regionbody(scale_regions, chromend));
            }
        }

        // Get flanking regions.
        // Note that we still need to deal with exon - non-exon as reference-point mode could require metagene walks.
        match self.strand.as_str() {
            "+" | "." => {
                match (&self.start, &self.end) {
                    (Revalue::U(_start), Revalue::U(end)) => {
                        let mut leftbins: Vec<Bin> = Vec::new();
                        let mut rightbins: Vec<Bin> = Vec::new();

                        let mut absstart: i64 = anchorstart as i64 - scale_regions.upstream as i64;
                        let absstop: i64 = anchorstop as i64 + scale_regions.downstream as i64;
                        for binix in
                            (absstart..anchorstart as i64).step_by(scale_regions.binsize as usize)
                        {
                            if binix + scale_regions.binsize as i64 <= 0 {
                                // Entirely upstream of chromosome start -> invalid
                                leftbins.push(Bin::Conbin(0, 0));
                            } else if binix >= 0 && binix as u32 > chromend {
                                // Entirely downstream of chromosome end -> invalid
                                leftbins.push(Bin::Conbin(0, 0));
                            } else if (binix + scale_regions.binsize as i64) as u32 > chromend {
                                // bin is partially downstream of chromosome end -> keep valid but truncate
                                let end = min(binix as u32 + scale_regions.binsize, chromend);
                                let start = max(binix, 0) as u32;
                                leftbins.push(Bin::Conbin(start, end));
                            } else {
                                // bin is valid -> keep as is
                                let start = max(binix, 0) as u32;
                                leftbins.push(Bin::Conbin(
                                    start,
                                    (binix as u32) + scale_regions.binsize,
                                ));
                            }
                        }

                        for binix in
                            (anchorstop as i64..absstop).step_by(scale_regions.binsize as usize)
                        {
                            let upper_bound: u32 = if scale_regions.nan_after_end {
                                min(chromend, *end)
                            } else {
                                chromend
                            };

                            if binix < 0 || binix as u32 >= upper_bound {
                                // entire bin at or past the valid region -> fully invalid
                                rightbins.push(Bin::Conbin(0, 0));
                            } else if (binix as u32 + scale_regions.binsize) > upper_bound {
                                // bin straddles the valid boundary -> keep only the valid portion
                                rightbins.push(Bin::Conbin(binix as u32, upper_bound));
                            } else {
                                // fully within bounds
                                rightbins.push(Bin::Conbin(
                                    binix as u32,
                                    binix as u32 + scale_regions.binsize,
                                ));
                            }
                        }

                        for bin in leftbins.into_iter() {
                            bins.push(bin);
                        }
                        // If we have bodybins, they should be squeezed in here.
                        for bin in bodybins.into_iter() {
                            bins.push(bin);
                        }
                        // Reset bodybins, as they are consumed.
                        bodybins = Vec::new();
                        for bin in rightbins.into_iter() {
                            bins.push(bin);
                        }
                    }
                    (Revalue::V(start), Revalue::V(end)) => {
                        let exons: Vec<(u32, u32)> = start
                            .iter()
                            .zip(end.iter())
                            .map(|(&s, &e)| (s, e))
                            .collect();
                        // Right side.
                        let mut rightbins: Vec<Bin> = Vec::new();
                        let mut lastanchor: u32 = anchorstop;
                        let mut walked_bps: u32 = 0;
                        while walked_bps < scale_regions.downstream {
                            if lastanchor > chromend {
                                rightbins.push(Bin::Conbin(0, 0));
                                walked_bps += scale_regions.binsize;
                            } else {
                                let (bin, retanch) = refpoint_exonwalker(
                                    &exons,
                                    lastanchor,
                                    scale_regions.binsize,
                                    chromend,
                                    scale_regions.nan_after_end,
                                    true,
                                );
                                rightbins.push(bin);
                                walked_bps += scale_regions.binsize;
                                lastanchor = retanch;
                            }
                        }
                        // Left side.
                        let mut leftbins: Vec<Bin> = Vec::new();
                        let mut lastanchor: u32 = anchorstart;
                        let mut walked_bps: u32 = 0;

                        while walked_bps < scale_regions.upstream {
                            if lastanchor == 0 {
                                leftbins.push(Bin::Conbin(0, 0));
                                walked_bps += scale_regions.binsize;
                            } else {
                                let (bin, retanch) = refpoint_exonwalker(
                                    &exons,
                                    lastanchor,
                                    scale_regions.binsize,
                                    chromend,
                                    scale_regions.nan_after_end,
                                    false,
                                );
                                leftbins.push(bin);
                                walked_bps += scale_regions.binsize;
                                lastanchor = retanch;
                            }
                        }
                        // Now we need to reverse the leftbins, as we walked backwards.
                        leftbins.reverse();
                        for bin in leftbins.into_iter() {
                            bins.push(bin);
                        }
                        // If we have bodybins, they should be squeezed in here.
                        for bin in bodybins.into_iter() {
                            bins.push(bin);
                        }
                        // Reset bodybins, as they are consumed.
                        bodybins = Vec::new();
                        for bin in rightbins.into_iter() {
                            bins.push(bin);
                        }
                    }
                    _ => panic!(
                        "Start and End are not either both u32, or Vecs. This means your regions file is ill-defined. Fix {}. {}:{}-{}",
                        self.name, self.chrom, self.start, self.end
                    ),
                }
            }
            "-" => {
                match (&self.start, &self.end) {
                    (Revalue::U(start), Revalue::U(_end)) => {
                        let mut leftbins: Vec<Bin> = Vec::new();
                        let mut rightbins: Vec<Bin> = Vec::new();

                        let mut absstart: i64 = anchorstop as i64 + scale_regions.upstream as i64;
                        let absstop: i64 = anchorstart as i64 - scale_regions.downstream as i64;

                        let steps: Vec<_> = (anchorstop as i64..absstart)
                            .step_by(scale_regions.binsize as usize)
                            .collect();
                        for binix in steps.into_iter().rev() {
                            if binix < 0 || binix as u32 > chromend {
                                // entire bin off the chromosome -> invalid
                                rightbins.push(Bin::Conbin(0, 0));
                            } else if (binix + scale_regions.binsize as i64) as u32 > chromend {
                                // bin straddles chromend -> truncate
                                let end = min(binix as u32 + scale_regions.binsize, chromend);
                                rightbins.push(Bin::Conbin(binix as u32, end));
                            } else {
                                rightbins.push(Bin::Conbin(
                                    binix as u32,
                                    binix as u32 + scale_regions.binsize,
                                ));
                            }
                        }
                        let steps: Vec<_> = (absstop..anchorstart as i64)
                            .step_by(scale_regions.binsize as usize)
                            .collect();
                        for binix in steps.into_iter().rev() {
                            let lower_bound: i64 = if scale_regions.nan_after_end {
                                max(*start as i64, 0)
                            } else {
                                0
                            };

                            if binix + scale_regions.binsize as i64 <= lower_bound {
                                // entire bin below the valid range -> invalid
                                leftbins.push(Bin::Conbin(0, 0));
                            } else if binix < lower_bound {
                                // bin straddles the lower boundary -> keep only the valid portion
                                leftbins.push(Bin::Conbin(
                                    lower_bound as u32,
                                    (binix + scale_regions.binsize as i64) as u32,
                                ));
                            } else {
                                // fully within bounds
                                leftbins.push(Bin::Conbin(
                                    binix as u32,
                                    (binix + scale_regions.binsize as i64) as u32,
                                ));
                            }
                        }

                        for bin in rightbins.into_iter() {
                            bins.push(bin);
                        }
                        // If we have bodybins, they should be squeezed in here.
                        bodybins.reverse();
                        for bin in bodybins.into_iter() {
                            bins.push(bin);
                        }
                        // Reset bodybins, as they are consumed.
                        bodybins = Vec::new();
                        for bin in leftbins.into_iter() {
                            bins.push(bin);
                        }
                    }
                    (Revalue::V(start), Revalue::V(end)) => {
                        let exons: Vec<(u32, u32)> = start
                            .iter()
                            .zip(end.iter())
                            .map(|(&s, &e)| (s, e))
                            .collect();
                        // Right side.
                        let mut rightbins: Vec<Bin> = Vec::new();
                        let mut lastanchor: u32 = anchorstop;
                        let mut walked_bps: u32 = 0;
                        while walked_bps < scale_regions.upstream {
                            if lastanchor > chromend {
                                rightbins.push(Bin::Conbin(0, 0));
                                walked_bps += scale_regions.binsize;
                            } else {
                                let (bin, retanch) = refpoint_exonwalker(
                                    &exons,
                                    lastanchor,
                                    scale_regions.binsize,
                                    chromend,
                                    scale_regions.nan_after_end,
                                    true,
                                );
                                rightbins.push(bin);
                                walked_bps += scale_regions.binsize;
                                lastanchor = retanch;
                            }
                        }
                        // Left side.
                        let mut leftbins: Vec<Bin> = Vec::new();
                        let mut lastanchor: u32 = anchorstart;
                        let mut walked_bps: u32 = 0;
                        while walked_bps < scale_regions.downstream {
                            if lastanchor == 0 {
                                leftbins.push(Bin::Conbin(0, 0));
                                walked_bps += scale_regions.binsize;
                            } else {
                                let (bin, retanch) = refpoint_exonwalker(
                                    &exons,
                                    lastanchor,
                                    scale_regions.binsize,
                                    chromend,
                                    scale_regions.nan_after_end,
                                    false,
                                );
                                leftbins.push(bin);
                                walked_bps += scale_regions.binsize;
                                lastanchor = retanch;
                            }
                        }

                        // Note that now we need to go the exact opposite way as for the + strand as the 'highest position' is the 'starting point'.
                        rightbins.reverse();
                        for bin in rightbins.into_iter() {
                            bins.push(bin);
                        }
                        bodybins.reverse();
                        for bin in bodybins.into_iter() {
                            bins.push(bin);
                        }
                        bodybins = Vec::new();
                        for bin in leftbins.into_iter() {
                            bins.push(bin);
                        }
                    }
                    _ => panic!(
                        "Start and End are not either both u32, or Vecs. This means your regions file is ill-defined. Fix {}. {}:{}-{}",
                        self.name, self.chrom, self.start, self.end
                    ),
                }
            }
            _ => panic!(
                "Strand should either be + or - or . {:?} is not supported.",
                self.strand
            ),
        };
        assert_eq!(
            bins.len(),
            scale_regions.cols_expected / scale_regions.bwfiles,
            "Number of bins does not match expected number of columns: {} != {}",
            bins.len(),
            scale_regions.cols_expected / scale_regions.bwfiles
        );
        bins
    }

    pub fn get_endu(&self) -> u32 {
        match &self.end {
            Revalue::U(end) => *end,
            Revalue::V(ends) => *ends.last().unwrap(),
        }
    }

    pub fn get_startu(&self) -> u32 {
        match &self.start {
            Revalue::U(start) => *start,
            Revalue::V(starts) => *starts.first().unwrap(),
        }
    }

    pub fn scale_regionbody(&self, scale_regions: &Scalingregions, chromend: u32) -> Vec<Bin> {
        // A vector of bins needs to be constructed for regionbody.
        // Since this is scaling mode, 'linspace' functionality is reproduced.
        match self.strand.as_str() {
            "+" | "." => {
                match (&self.start, &self.end) {
                    (Revalue::U(start), Revalue::U(end)) => {
                        // No exons, forward strand. divide start  - end as such:
                        // |---un5prime---|---bodylength---|---un3prime---|
                        let mut un5bins: Vec<Bin> = Vec::new();
                        let mut un3bins: Vec<Bin> = Vec::new();
                        let mut innerbins: Vec<Bin> = Vec::new();
                        if scale_regions.unscaled5prime > 0 {
                            un5bins.extend(
                                (0..scale_regions.unscaled5prime)
                                    .step_by(scale_regions.binsize as usize)
                                    .map(|i| {
                                        Bin::Conbin(*start + i, *start + i + scale_regions.binsize)
                                    })
                                    .collect::<Vec<Bin>>(),
                            );
                        }

                        if scale_regions.unscaled3prime > 0 {
                            un3bins.extend(
                                (0..scale_regions.unscaled3prime)
                                    .step_by(scale_regions.binsize as usize)
                                    .rev()
                                    .map(|i| {
                                        Bin::Conbin(*end - i - scale_regions.binsize, *end - i)
                                    })
                                    .collect::<Vec<Bin>>(),
                            );
                        }
                        let bodystart = *start + scale_regions.unscaled5prime;
                        let bodyend = *end - scale_regions.unscaled3prime;

                        // Get the bins over the body length. These need to be scaled, so similar to deeptools < 4, linspace is used.
                        let neededbins =
                            (scale_regions.regionbodylength / scale_regions.binsize) as usize;
                        // There's multiple options here:
                        // transcriptlength >= regionbodylength -> linspace
                        // regionbodylength / binsize > transcriptlength <= regionbodylength -> 1 >= binsize > binsize.
                        // transcriptlength <= regionbodylength / binsize -> index repetitions with binsize of one.
                        let scaledbinsize = min(
                            max((bodyend - bodystart) / neededbins as u32, 1),
                            scale_regions.binsize,
                        );
                        innerbins.extend(
                            Array1::linspace(
                                bodystart as f32,
                                (bodyend - scaledbinsize) as f32,
                                neededbins,
                            )
                            .mapv(|x| x.floor() as u32)
                            .map(|x| Bin::Conbin(*x, *x + scaledbinsize))
                            .into_iter()
                            .collect::<Vec<_>>(),
                        );
                        // Combine the vectors and return
                        let mut combined_bins = Vec::new();
                        if scale_regions.unscaled5prime > 0 {
                            combined_bins.extend(un5bins.into_iter());
                        }
                        combined_bins.extend(innerbins.into_iter());
                        if scale_regions.unscaled3prime > 0 {
                            combined_bins.extend(un3bins.into_iter());
                        }
                        return combined_bins;
                    }
                    (Revalue::V(start), Revalue::V(end)) => {
                        let exons: Vec<(u32, u32)> = start
                            .iter()
                            .zip(end.iter())
                            .map(|(&s, &e)| (s, e))
                            .collect();
                        let mut un5bins: Vec<Bin> = Vec::new();
                        let mut un3bins: Vec<Bin> = Vec::new();
                        if scale_regions.unscaled5prime > 0 {
                            let mut walked_bps: u32 = 0;
                            let mut lastanchor: u32 = start[0];

                            while walked_bps < scale_regions.unscaled5prime {
                                let (bin, retanch) = refpoint_exonwalker(
                                    &exons,
                                    lastanchor,
                                    scale_regions.binsize,
                                    chromend,
                                    scale_regions.nan_after_end,
                                    true,
                                );
                                un5bins.push(bin);
                                walked_bps += scale_regions.binsize;
                                lastanchor = retanch;
                            }
                        }
                        if scale_regions.unscaled3prime > 0 {
                            let mut walked_bps: u32 = 0;
                            let mut lastanchor: u32 = *end.last().unwrap();

                            while walked_bps < scale_regions.unscaled3prime {
                                let (bin, retanch) = refpoint_exonwalker(
                                    &exons,
                                    lastanchor,
                                    scale_regions.binsize,
                                    chromend,
                                    scale_regions.nan_after_end,
                                    false,
                                );
                                un3bins.push(bin);
                                walked_bps += scale_regions.binsize;
                                lastanchor = retanch;
                            }
                        }
                        un3bins.reverse();

                        let bodystart: u32;
                        let bodyend: u32;
                        if scale_regions.unscaled5prime > 0 {
                            bodystart = un5bins.last().unwrap().get_end() - 1;
                        } else {
                            bodystart = *start.first().unwrap();
                        }
                        if scale_regions.unscaled3prime > 0 {
                            bodyend = un3bins.first().unwrap().get_start();
                        } else {
                            bodyend = *end.last().unwrap();
                        }
                        let truebodylength = self.regionlength
                            - scale_regions.unscaled5prime
                            - scale_regions.unscaled3prime;
                        let neededbins =
                            (scale_regions.regionbodylength / scale_regions.binsize) as usize;
                        let scaledbinsize = min(
                            max(truebodylength / neededbins as u32, 1),
                            scale_regions.binsize,
                        );
                        // Things are a bit tricky now, as we can do a linspace over the region, but we don't have a notion of the exons.
                        // I think easiest is to just pull a hashmap over the entire region, get linspace from hashmap to vec, and be done with it.
                        // technically we fetch a bunch of regions we don't need, but this operation is not too expensive.

                        let mut binmap: HashMap<u32, Bin> = HashMap::new();
                        let mut lastanchor: u32 = bodystart;

                        for ix in 0..((truebodylength / scaledbinsize) + 1) {
                            let (bin, anchor) = refpoint_exonwalker(
                                &exons,
                                lastanchor,
                                scaledbinsize,
                                chromend,
                                scale_regions.nan_after_end,
                                true,
                            );
                            lastanchor = anchor;
                            match bin {
                                Bin::Conbin(start, end) => {
                                    if end > bodyend {
                                        binmap.insert(ix, Bin::Conbin(start, bodyend));
                                    } else {
                                        binmap.insert(ix, Bin::Conbin(start, end));
                                    }
                                }
                                Bin::Catbin(bins) => {
                                    if bins.last().unwrap().1 > bodyend {
                                        let mut newbins: Vec<(u32, u32)> = Vec::new();
                                        for (s, e) in bins.iter() {
                                            if *e > bodyend {
                                                newbins.push((*s, bodyend));
                                            } else {
                                                newbins.push((*s, *e));
                                            }
                                        }
                                        binmap.insert(ix, Bin::Catbin(newbins));
                                    } else {
                                        binmap.insert(ix, Bin::Catbin(bins));
                                    }
                                }
                            }
                        }

                        let innerbins = Array1::linspace(
                            0 as f32,
                            ((truebodylength) / scaledbinsize) as f32,
                            neededbins,
                        )
                        .mapv(|x| x.floor() as u32)
                        .map(|x| binmap.get(&x).unwrap().clone())
                        .into_iter()
                        .collect::<Vec<Bin>>();

                        // Combine the vectors and return
                        let mut combined_bins = Vec::new();
                        if scale_regions.unscaled5prime > 0 {
                            combined_bins.extend(un5bins.into_iter());
                        }
                        combined_bins.extend(innerbins.into_iter());
                        if scale_regions.unscaled3prime > 0 {
                            combined_bins.extend(un3bins.into_iter());
                        }
                        return combined_bins;
                    }
                    _ => panic!(
                        "Start and End are not either both u32, or Vecs. This means your regions file is ill-defined."
                    ),
                }
            }
            "-" => {
                match (&self.start, &self.end) {
                    (Revalue::U(start), Revalue::U(end)) => {
                        // No exons, negative strand. divide start  - end as such:
                        // |---un3prime---|---bodylength---|---un5prime---|
                        let mut un5bins: Vec<Bin> = Vec::new();
                        let mut un3bins: Vec<Bin> = Vec::new();
                        let mut innerbins: Vec<Bin> = Vec::new();
                        if scale_regions.unscaled3prime > 0 {
                            un3bins.extend(
                                (0..scale_regions.unscaled3prime)
                                    .step_by(scale_regions.binsize as usize)
                                    .map(|i| {
                                        Bin::Conbin(*start + i, *start + i + scale_regions.binsize)
                                    })
                                    .collect::<Vec<Bin>>(),
                            );
                        }

                        if scale_regions.unscaled5prime > 0 {
                            un5bins.extend(
                                (0..scale_regions.unscaled5prime)
                                    .step_by(scale_regions.binsize as usize)
                                    .rev()
                                    .map(|i| {
                                        Bin::Conbin(*end - i - scale_regions.binsize, *end - i)
                                    })
                                    .collect::<Vec<Bin>>(),
                            );
                        }
                        let bodystart = *start + scale_regions.unscaled3prime;
                        let bodyend = *end - scale_regions.unscaled5prime;

                        // Get the bins over the body length. These need to be scaled, so similar to deeptools < 4, linspace is used.
                        let neededbins =
                            (scale_regions.regionbodylength / scale_regions.binsize) as usize;
                        // There's multiple options here:
                        // transcriptlength >= regionbodylength -> linspace
                        // regionbodylength / binsize > transcriptlength <= regionbodylength -> 1 >= binsize > binsize.
                        // transcriptlength <= regionbodylength / binsize -> index repetitions with binsize of one.
                        let scaledbinsize = min(
                            max((bodyend - bodystart) / neededbins as u32, 1),
                            scale_regions.binsize,
                        );
                        innerbins.extend(
                            Array1::linspace(
                                bodystart as f32,
                                (bodyend - scaledbinsize) as f32,
                                neededbins,
                            )
                            .mapv(|x| x.floor() as u32)
                            .map(|x| Bin::Conbin(*x, *x + scaledbinsize))
                            .into_iter()
                            .collect::<Vec<_>>(),
                        );
                        // Combine the vectors and return
                        let mut combined_bins = Vec::new();
                        if scale_regions.unscaled3prime > 0 {
                            combined_bins.extend(un3bins.into_iter());
                        }
                        combined_bins.extend(innerbins.into_iter());
                        if scale_regions.unscaled5prime > 0 {
                            combined_bins.extend(un5bins.into_iter());
                        }
                        return combined_bins;
                    }
                    (Revalue::V(start), Revalue::V(end)) => {
                        let exons: Vec<(u32, u32)> = start
                            .iter()
                            .zip(end.iter())
                            .map(|(&s, &e)| (s, e))
                            .collect();
                        let mut un5bins: Vec<Bin> = Vec::new();
                        let mut un3bins: Vec<Bin> = Vec::new();
                        if scale_regions.unscaled5prime > 0 {
                            let mut walked_bps: u32 = 0;
                            let mut lastanchor: u32 = *end.last().unwrap();

                            while walked_bps < scale_regions.unscaled5prime {
                                let (bin, retanch) = refpoint_exonwalker(
                                    &exons,
                                    lastanchor,
                                    scale_regions.binsize,
                                    chromend,
                                    scale_regions.nan_after_end,
                                    false,
                                );
                                un5bins.push(bin);
                                walked_bps += scale_regions.binsize;
                                lastanchor = retanch;
                            }
                        }
                        un5bins.reverse();
                        if scale_regions.unscaled3prime > 0 {
                            let mut walked_bps: u32 = 0;
                            let mut lastanchor: u32 = start[0];

                            while walked_bps < scale_regions.unscaled3prime {
                                let (bin, retanch) = refpoint_exonwalker(
                                    &exons,
                                    lastanchor,
                                    scale_regions.binsize,
                                    chromend,
                                    scale_regions.nan_after_end,
                                    true,
                                );
                                un3bins.push(bin);
                                walked_bps += scale_regions.binsize;
                                lastanchor = retanch;
                            }
                        }

                        let bodystart: u32;
                        let bodyend: u32;
                        if scale_regions.unscaled5prime > 0 {
                            bodyend = un5bins.first().unwrap().get_start();
                        } else {
                            bodyend = *start.first().unwrap();
                        }
                        if scale_regions.unscaled3prime > 0 {
                            bodystart = un3bins.last().unwrap().get_end() - 1;
                        } else {
                            bodystart = *end.last().unwrap();
                        }
                        let truebodylength = self.regionlength
                            - scale_regions.unscaled5prime
                            - scale_regions.unscaled3prime;
                        let neededbins =
                            (scale_regions.regionbodylength / scale_regions.binsize) as usize;
                        let scaledbinsize = min(
                            max(truebodylength / neededbins as u32, 1),
                            scale_regions.binsize,
                        );
                        // Things are a bit tricky now, as we can do a linspace over the region, but we don't have a notion of the exons.
                        // I think easiest is to just pull a hashmap over the entire region, get linspace from hashmap to vec, and be done with it.
                        // technically we fetch a bunch of regions we don't need, but this operation is not too expensive.

                        let mut binmap: HashMap<u32, Bin> = HashMap::new();
                        let mut lastanchor: u32 = bodystart;

                        for ix in 0..((truebodylength / scaledbinsize) + 1) {
                            let (bin, anchor) = refpoint_exonwalker(
                                &exons,
                                lastanchor,
                                scaledbinsize,
                                chromend,
                                scale_regions.nan_after_end,
                                true,
                            );
                            lastanchor = anchor;
                            match bin {
                                Bin::Conbin(start, end) => {
                                    if end > bodyend {
                                        binmap.insert(ix, Bin::Conbin(start, bodyend));
                                    } else {
                                        binmap.insert(ix, Bin::Conbin(start, end));
                                    }
                                }
                                Bin::Catbin(bins) => {
                                    if bins.last().unwrap().1 > bodyend {
                                        let mut newbins: Vec<(u32, u32)> = Vec::new();
                                        for (s, e) in bins.iter() {
                                            if *e > bodyend {
                                                newbins.push((*s, bodyend));
                                            } else {
                                                newbins.push((*s, *e));
                                            }
                                        }
                                        binmap.insert(ix, Bin::Catbin(newbins));
                                    } else {
                                        binmap.insert(ix, Bin::Catbin(bins));
                                    }
                                }
                            }
                        }

                        let innerbins = Array1::linspace(
                            0 as f32,
                            ((truebodylength) / scaledbinsize) as f32,
                            neededbins,
                        )
                        .mapv(|x| x.floor() as u32)
                        .map(|x| binmap.get(&x).unwrap().clone())
                        .into_iter()
                        .collect::<Vec<Bin>>();

                        // Combine the vectors and return
                        let mut combined_bins = Vec::new();
                        if scale_regions.unscaled3prime > 0 {
                            combined_bins.extend(un3bins.into_iter());
                        }
                        combined_bins.extend(innerbins.into_iter());
                        if scale_regions.unscaled5prime > 0 {
                            combined_bins.extend(un5bins.into_iter());
                        }
                        return combined_bins;
                    }
                    _ => panic!(
                        "Start and End are not either both u32, or Vecs. This means your regions file is ill-defined."
                    ),
                }
            }
            _ => panic!(
                "Strand should either be + or - or . {:?} is not supported.",
                self.strand
            ),
        }
    }
}

fn refpoint_exonwalker(
    exons: &Vec<(u32, u32)>,
    anchor: u32,
    binsize: u32,
    chromend: u32,
    nan_after_end: bool,
    forward: bool,
) -> (Bin, u32) {
    // Basic function that walks over exons, and returns a Bin (Either Conbin or Catbin) and the last anchorpoint.
    let mut anchorix: Option<usize> = None;

    for (ix, exon) in exons.iter().enumerate() {
        if anchor >= exon.0 && anchor <= exon.1 {
            anchorix = Some(ix);
        }
    }
    if forward {
        // Walk forward (downstream, towards chromosome end)
        match anchorix {
            Some(i) => {
                // anchor sits in an exon. Check if anchor + binsize is also in same exon.
                if anchor + binsize > chromend {
                    return (Bin::Conbin(0, 0), chromend);
                }
                if anchor + binsize <= exons[i].1 {
                    (Bin::Conbin(anchor, anchor + binsize), anchor + binsize)
                } else {
                    // anchor + binsize is not in same exon. We need a Catbin.
                    // Things are a bit more difficult here as well, as we need to walk exons.
                    let mut start_end_vec: Vec<(u32, u32)> = Vec::new();
                    start_end_vec.push((anchor, exons[i].1));

                    let mut remainingbin: u32 = binsize - (exons[i].1 - anchor);
                    let mut lastix: usize = i;
                    let mut lastanchor: u32 = exons[i].1;

                    while remainingbin != 0 {
                        if lastix + 1 < exons.len() {
                            // next exon is available.
                            // Two options here:
                            // the remainder fits in the lastix + 1 exon. We are done.
                            // the remainder doesn't fit in the lastix + 1 exon. We need to walk further.
                            if exons[lastix + 1].1 - exons[lastix + 1].0 >= remainingbin {
                                // remainder fits in next exon.
                                start_end_vec.push((
                                    exons[lastix + 1].0,
                                    exons[lastix + 1].0 + remainingbin,
                                ));
                                lastanchor = exons[lastix + 1].0 + remainingbin;
                                remainingbin = 0;
                            } else {
                                // Remainder is larger then our exon. We need another walk.
                                start_end_vec.push((exons[lastix + 1].0, exons[lastix + 1].1));
                                remainingbin -= exons[lastix + 1].1 - exons[lastix + 1].0;
                                lastix += 1;
                                lastanchor = exons[lastix].1;
                            }
                        } else {
                            // No next exon available. Remainder can just be genomic.
                            // The last entry here can be changed to include the last part.
                            if nan_after_end {
                                start_end_vec.push((0, 0));
                            } else {
                                let last = start_end_vec.last_mut().unwrap();
                                assert_eq!(
                                    last.1, lastanchor,
                                    "In the exon - genomic walk, our coordinates are not contiguous"
                                );
                                // Check we don't fall of the chromosome.
                                if lastanchor + remainingbin > chromend {
                                    last.1 = chromend;
                                    lastanchor = chromend;
                                } else {
                                    last.1 = lastanchor + remainingbin;
                                    lastanchor = lastanchor + remainingbin;
                                }
                            }
                            remainingbin = 0;
                        }
                    }
                    // We now have a Vec of start - end, we can construct a CatBin.
                    // Note that CatBins are (absstart, absstop, ((intstart1, intstart2), ...))
                    // This seems weird, but makes sure we need to slice the bigwig file only once per bin.
                    if start_end_vec.len() == 1 {
                        (Bin::Conbin(anchor, lastanchor), lastanchor)
                    } else {
                        (Bin::Catbin(start_end_vec), lastanchor)
                    }
                }
            }
            None => {
                // our anchor doesn't sit in exons. We just return the anchor + binsize as Bin
                if anchor + binsize > chromend {
                    (Bin::Conbin(0, 0), chromend)
                } else {
                    (Bin::Conbin(anchor, anchor + binsize), anchor + binsize)
                }
            }
        }
    } else {
        // Walk backwards (upstream, towards chromosome start)
        match anchorix {
            Some(i) => {
                // Run a check to see if binsize > anchor.
                if anchor < binsize {
                    // We are at the start of the chromosome. We need to return a Conbin.
                    return (Bin::Conbin(0, 0), 0);
                }
                // anchor sits in an exon. Check if anchor - binsize is also in same exon.
                if anchor - binsize >= exons[i].0 {
                    (Bin::Conbin(anchor - binsize, anchor), anchor - binsize)
                } else {
                    // anchor + binsize is not in same exon. We need a Catbin.
                    // Things are a bit more difficult here as well, as we need to walk exons.
                    let mut start_end_vec: Vec<(u32, u32)> = Vec::new();
                    start_end_vec.push((exons[i].0, anchor));

                    let mut remainingbin: u32 = binsize - (anchor - exons[i].0);
                    let mut lastix: usize = i;
                    let mut lastanchor: u32 = exons[i].0;

                    while remainingbin != 0 {
                        if lastix >= 1 {
                            // previous exon is available.
                            // Two options here:
                            // the remainder fits in the previous exon. We are done.
                            // the remainder doesn't fit in the previous exon. We need to walk further.
                            if exons[lastix - 1].1 - exons[lastix - 1].0 >= remainingbin {
                                // remainder fits in next exon.
                                start_end_vec.push((
                                    exons[lastix - 1].1 - remainingbin,
                                    exons[lastix - 1].1,
                                ));
                                lastanchor = exons[lastix - 1].1 - remainingbin;
                                remainingbin = 0;
                            } else {
                                // Remainder is larger then our exon. We need another walk.
                                start_end_vec.push((exons[lastix - 1].0, exons[lastix - 1].1));
                                remainingbin -= exons[lastix - 1].1 - exons[lastix - 1].0;
                                lastix -= 1;
                                lastanchor = exons[lastix].0;
                            }
                        } else {
                            // No previous exon available. Remainder can just be genomic.
                            // The last entry here can be changed to include the last part.
                            if nan_after_end {
                                start_end_vec.push((0, 0));
                            } else {
                                let last = start_end_vec.last_mut().unwrap();
                                assert_eq!(
                                    last.0, lastanchor,
                                    "In the exon - genomic walk (reverse), our coordinates are not contiguous"
                                );
                                // Check we don't go in the negative.
                                if lastanchor < remainingbin {
                                    last.0 = 0;
                                    lastanchor = 0;
                                } else {
                                    last.0 = lastanchor - remainingbin;
                                    lastanchor = lastanchor - remainingbin;
                                }
                            }
                            remainingbin = 0;
                        }
                    }
                    // We now have a Vec of start - end, we can construct a CatBin.
                    // Note that CatBins are (absstart, absstop, ((intstart1, intstart2), ...))
                    // This seems weird, but makes sure we need to slice the bigwig file only once per bin.
                    if start_end_vec.len() == 1 {
                        (Bin::Conbin(lastanchor, anchor), lastanchor)
                    } else {
                        start_end_vec.reverse();
                        (Bin::Catbin(start_end_vec), lastanchor)
                    }
                }
            }
            None => {
                // our anchor doesn't sit in exons. We just return the anchor - binsize as Bin
                if anchor < binsize {
                    (Bin::Conbin(0, 0), 0)
                } else {
                    (Bin::Conbin(anchor - binsize, anchor), anchor - binsize)
                }
            }
        }
    }
}

pub fn region_divider(regs: &Vec<Region>) -> Vec<Vec<Region>> {
    // This function decides on how regions are divided to process in parallel.
    // Since the Vec<Region> could be chromosomes on one hand, or genes on the other, we can decide on how to chop them up.
    // Note that this somehow relates to 'genomeChunkLength', but not entirely.
    let mut bplen: u32 = 0;
    let mut blocks: Vec<Vec<Region>> = Vec::new();
    let mut tempregionvec: Vec<Region> = Vec::new();

    for reg in regs.iter() {
        if reg.regionlength > 10000000 {
            if tempregionvec.len() > 0 {
                blocks.push(tempregionvec.clone());
                tempregionvec = Vec::new();
                bplen = 0
            }

            blocks.push(vec![reg.clone()]);
        } else {
            tempregionvec.push(reg.clone());
            bplen += reg.regionlength;
        }
        if bplen > 1000000 {
            blocks.push(tempregionvec.clone());
            tempregionvec = Vec::new();
            bplen = 0;
        }
    }
    if tempregionvec.len() > 0 {
        blocks.push(tempregionvec.clone());
    }
    blocks
}

#[derive(Debug)]
#[allow(dead_code)]
pub struct Scalingregions {
    pub upstream: u32,
    pub downstream: u32,
    pub unscaled5prime: u32,
    pub unscaled3prime: u32,
    pub regionbodylength: u32,
    pub binsize: u32,
    pub cols_expected: usize,
    pub bpsum: u32,
    pub missingdata_as_zero: bool,
    pub scale: f32,
    pub nan_after_end: bool,
    pub skipzero: bool,
    pub minthresh: f32,
    pub maxthresh: f32,
    pub referencepoint: String,
    pub mode: String,
    pub bwfiles: usize,
    pub avgtype: String,
    pub verbose: bool,
    pub proc_number: usize,
    pub regionlabels: Vec<String>,
    pub bwlabels: Vec<String>,
    pub startlabel: String,
    pub endlabel: String,
}

#[derive(Clone)]
pub struct Gtfparse {
    pub metagene: bool,
    pub txnid: String,
    pub exonid: String,
    pub txniddesignator: String,
}

#[derive(Clone)]
pub enum Revalue {
    U(u32),
    V(Vec<u32>),
}

impl Revalue {
    pub fn rewrites(&self) -> String {
        match self {
            Revalue::U(v) => format!("{}", v),
            Revalue::V(vs) => vs
                .iter()
                .map(|v| v.to_string())
                .collect::<Vec<_>>()
                .join(","),
        }
    }
    /// Sort key for a "start"-type value: first element if V, else the U value.
    pub fn start_key(&self) -> u32 {
        match self {
            Revalue::U(v) => *v,
            Revalue::V(vs) => *vs.first().expect("Revalue::V should not be empty"),
        }
    }

    /// Sort key for an "end"-type value: last element if V, else the U value.
    pub fn end_key(&self) -> u32 {
        match self {
            Revalue::U(v) => *v,
            Revalue::V(vs) => *vs.last().expect("Revalue::V should not be empty"),
        }
    }
}

impl fmt::Debug for Revalue {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Revalue::U(value) => write!(f, "U({})", value),
            Revalue::V(values) => write!(f, "V({:?})", values),
        }
    }
}

impl fmt::Display for Revalue {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Revalue::U(value) => write!(f, "U({})", value),
            Revalue::V(values) => write!(
                f,
                "V({})",
                values
                    .iter()
                    .map(|v| v.to_string())
                    .collect::<Vec<String>>()
                    .join(", ")
            ),
        }
    }
}

#[derive(Clone, Debug)]
pub enum Bin {
    Conbin(u32, u32),
    Catbin(Vec<(u32, u32)>),
}

impl Bin {
    pub fn get_start(&self) -> u32 {
        match self {
            Bin::Conbin(start, _) => *start,
            Bin::Catbin(starts) => starts.first().unwrap().0,
        }
    }
    pub fn get_end(&self) -> u32 {
        match self {
            Bin::Conbin(_, end) => *end,
            Bin::Catbin(ends) => ends.last().unwrap().1,
        }
    }
}

pub struct TempZip<I>
where
    I: Iterator,
{
    pub iterators: Vec<I>,
}

impl<I, T> Iterator for TempZip<I>
where
    I: Iterator<Item = T>,
{
    type Item = Vec<T>;
    fn next(&mut self) -> Option<Self::Item> {
        let o: Option<Vec<T>> = self.iterators.iter_mut().map(|x| x.next()).collect();
        o
    }
}
