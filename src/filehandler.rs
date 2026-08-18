use crate::calc::{max_float, mean_float, median_float, min_float, std_float, sum_float};
use crate::covcalc::{Bin, Gtfparse, Region, Revalue, Scalingregions};
use bigtools::beddata::BedParserStreamingIterator;
use bigtools::{BigWigRead, BigWigWrite, Value};
use flate2::Compression;
use flate2::read::MultiGzDecoder;
use flate2::write::GzEncoder;
use itertools::Itertools;
use rust_htslib::bam::{IndexedReader, Read, Reader};
use std::collections::HashMap;
use std::fs::File;
use std::io::prelude::*;
use std::io::{BufReader, BufWriter, Write};
use std::path::Path;

pub fn bam_ispaired(bam_ifile: &str) -> bool {
    let mut bam = Reader::from_path(bam_ifile)
        .unwrap_or_else(|e| panic!("Failed to open BAM file '{}': {}", bam_ifile, e));
    let mut count = 0;
    const MAX_READS: usize = 1000;
    for record in bam.records() {
        let record = record.expect("Error parsing record.");
        if record.is_paired() {
            return true;
        }
        count += 1;
        if count >= MAX_READS {
            break;
        }
    }
    return false;
}

pub fn write_covfile<LI>(
    lines: LI,
    ofile: &str,
    filetype: &str,
    chromsizes: HashMap<String, u32>,
    nproc: usize,
)
where
    LI: Iterator<Item = (String, Value)>,
{
    if filetype == "bedgraph" {
        // write output file, bedgraph
        let mut writer = BufWriter::new(
            File::create(ofile)
                .unwrap_or_else(|e| panic!("Failed to create output file '{}': {}", ofile, e)),
        );
        for (chrom, val) in lines {
            writeln!(
                writer,
                "{}\t{}\t{}\t{}",
                chrom, val.start, val.end, val.value
            )
            .unwrap_or_else(|e| panic!("Failed to write bedgraph line to '{}': {}", ofile, e));
        }
    } else {
        let vals = BedParserStreamingIterator::wrap_infallible_iter(lines, false);
        let mut writer = BigWigWrite::create_file(ofile, chromsizes)
            .unwrap_or_else(|e| panic!("Failed to create output bigwig file '{}': {}", ofile, e));
        // bigtools spawns one async task per chromosome (zoom aggregation + section
        // encoding) and funnels their output through a single sequential writer task
        // that owns the actual file handle; the crate's own channels serialize that
        // hand-off, so raising worker count only widens the per-chromosome encode
        // stage, it does not add any concurrent writers to `ofile`.
        // Mirrors bigtools' own bedgraphtobigwig CLI: skip the multi-thread runtime
        // entirely at nproc == 1, since it would only ever have one consumer.
        let runtime = if nproc <= 1 {
            writer.options.channel_size = 0;
            tokio::runtime::Builder::new_current_thread()
                .build()
                .expect("Unable to create tokio runtime for bw writing.")
        } else {
            tokio::runtime::Builder::new_multi_thread()
                .worker_threads(nproc)
                .build()
                .expect("Unable to create tokio runtime for bw writing.")
        };
        let _ = writer.write(vals, runtime);
    }
}

fn open_bed_or_gtf_reader(fp: &str) -> BufReader<Box<dyn std::io::Read>> {
    let mut file = File::open(fp).expect(format!("Failed to open file: {}", fp).as_str());
    let mut magic = [0u8; 2];
    let n = std::io::Read::read(&mut file, &mut magic).unwrap_or(0);
    std::io::Seek::seek(&mut file, std::io::SeekFrom::Start(0))
        .expect(format!("Failed to seek in file: {}", fp).as_str());
    if n == 2 && magic[0] == 0x1f && magic[1] == 0x8b {
        BufReader::new(Box::new(MultiGzDecoder::new(file)) as Box<dyn std::io::Read>)
    } else {
        BufReader::new(Box::new(file) as Box<dyn std::io::Read>)
    }
}

pub fn is_bed_or_gtf(fp: &str) -> String {
    // Check if file is a bed or gtf file.
    let reader = open_bed_or_gtf_reader(fp);
    // Get the first line that doesn't start with #
    for line in reader.lines() {
        let line = line.unwrap_or_else(|e| panic!("Failed to read a line from '{}': {}", fp, e));
        if !line.starts_with('#') {
            let fields: Vec<&str> = line.split('\t').collect();
            if fields.len() == 9 {
                return "gtf".to_string();
            } else {
                return "bed".to_string();
            }
        }
    }
    "Unknown".to_string()
}

pub fn read_gtffile(
    gtf_file: &String,
    gtfparse: &Gtfparse,
    chroms: Vec<&String>,
) -> (Vec<Region>, (String, u32)) {
    // At some point this zoo of String clones should be refactored. Not now though, We have a deadline.
    let mut regions: Vec<Region> = Vec::new();
    let mut names: HashMap<String, u32> = HashMap::new();
    let mut entries: u32 = 0;
    let mut txnids: Vec<String> = Vec::new();

    let gtffile = open_bed_or_gtf_reader(gtf_file);

    if gtfparse.metagene {
        // metagene implementation - more work here.
        let mut txn_hash: HashMap<String, Vec<(u32, u32)>> = HashMap::new();
        let mut txn_strand: HashMap<String, String> = HashMap::new();
        let mut txn_chrom: HashMap<String, String> = HashMap::new();
        // Store transcript-level coordinates as fallback when no exon entries exist.
        let mut txn_transcript: HashMap<String, (u32, u32)> = HashMap::new();

        for line in gtffile.lines() {
            let line =
                line.unwrap_or_else(|e| panic!("Failed to read a line from '{}': {}", gtf_file, e));
            // skip comments
            if line.starts_with('#') {
                continue;
            }
            let fields: Vec<&str> = line.split('\t').collect();
            let feature = fields[2].to_string();
            let mut start: u32 = fields[3].parse().unwrap_or_else(|e| {
                panic!(
                    "Failed to parse start position in GTF line '{}': {}",
                    line, e
                )
            });
            if start >= 1 {
                start -= 1;
            }
            let end: u32 = fields[4].parse().unwrap_or_else(|e| {
                panic!("Failed to parse end position in GTF line '{}': {}", line, e)
            });
            let txnid: Option<String> = fields[8]
                .split(';')
                .find(|x| x.trim().starts_with(gtfparse.txniddesignator.as_str()))
                .and_then(|x| x.split('"').nth(1))
                .map(|s| s.to_string());

            if feature == gtfparse.exonid {
                let tid = txnid.unwrap_or_else(|| {
                    panic!(
                        "GTF line has an exon feature but no transcript id (attribute '{}') was found: '{}'",
                        gtfparse.txniddesignator, line
                    )
                });
                if !txnids.contains(&tid) {
                    txnids.push(tid.clone());
                }

                let txnentry = txn_hash.entry(tid.clone()).or_insert(Vec::new());
                if txn_strand.contains_key(&tid) {
                    assert_eq!(
                        txn_strand.get(&tid).expect(
                            "Transcript id vanished from txn_strand map between check and get"
                        ),
                        fields[6]
                    );
                } else {
                    txn_strand.insert(tid.clone(), fields[6].to_string());
                }
                if txn_chrom.contains_key(&tid) {
                    assert_eq!(
                        txn_chrom.get(&tid).expect(
                            "Transcript id vanished from txn_chrom map between check and get"
                        ),
                        fields[0]
                    );
                } else {
                    txn_chrom.insert(tid.clone(), fields[0].to_string());
                }
                txnentry.push((start, end));
            } else if feature == gtfparse.txnid {
                // Store transcript-level coordinates as fallback.
                if let Some(tid) = txnid {
                    if !txn_strand.contains_key(&tid) {
                        txn_strand.insert(tid.clone(), fields[6].to_string());
                    }
                    if !txn_chrom.contains_key(&tid) {
                        txn_chrom.insert(tid.clone(), fields[0].to_string());
                    }
                    txn_transcript.insert(tid, (start, end));
                }
            }
        }

        // Also pull in transcript IDs that have no exon entries.
        for tid in txn_transcript.keys() {
            if !txnids.contains(tid) {
                txnids.push(tid.clone());
            }
        }

        for txnid in txnids.into_iter() {
            let has_exons = txn_hash.contains_key(&txnid);
            let has_transcript = txn_transcript.contains_key(&txnid);

            if has_exons && !has_transcript {
                continue;
            }

            let (starts, ends) = if has_exons {
                let txnentry = txn_hash
                    .get_mut(&txnid)
                    .expect("Transcript id vanished from txn_hash map between check and get");
                txnentry.sort_by(|a, b| a.0.cmp(&b.0));
                txnentry.iter().map(|(s, e)| (*s, *e)).unzip()
            } else {
                // No exons, fall back to transcript-level coordinates.
                let (s, e) = txn_transcript
                    .get(&txnid)
                    .expect("Transcript id vanished from txn_transcript map between check and get");
                (vec![*s], vec![*e])
            };
            let length: u32 = starts.iter().zip(ends.iter()).map(|(s, e)| e - s).sum();
            let chrom = txn_chrom
                .get(&txnid)
                .unwrap_or_else(|| panic!("Transcript '{}' has no chromosome recorded", txnid))
                .to_string();

            if !chroms.contains(&&chrom.to_string()) {
                println!(
                    "Warning, region {} not found in at least one of the bigwig/bam files. Skipping {}.",
                    chrom, txnid
                );
            } else {
                regions.push(Region {
                    chrom,
                    start: Revalue::V(starts),
                    end: Revalue::V(ends),
                    score: ".".to_string(),
                    strand: txn_strand
                        .get(&txnid)
                        .unwrap_or_else(|| panic!("Transcript '{}' has no strand recorded", txnid))
                        .to_string(),
                    name: txnid.to_string(),
                    regionlength: length,
                });
                entries += 1;
            }
        }
    } else {
        // Take fields with col 3 == gtfparse.txnid, start, end
        for line in gtffile.lines() {
            let line =
                line.unwrap_or_else(|e| panic!("Failed to read a line from '{}': {}", gtf_file, e));
            // skip comments
            if line.starts_with('#') {
                continue;
            }

            let fields: Vec<&str> = line.split('\t').collect();
            if fields[2].to_string() == gtfparse.txnid {
                let mut start: u32 = fields[3].parse().unwrap_or_else(|e| {
                    panic!(
                        "Failed to parse start position in GTF line '{}': {}",
                        line, e
                    )
                });
                if start >= 1 {
                    start -= 1;
                }
                let end: u32 = fields[4].parse().unwrap_or_else(|e| {
                    panic!("Failed to parse end position in GTF line '{}': {}", line, e)
                });
                let mut entryname = fields[8]
                    .split(';')
                    .find(|x| x.trim().starts_with(gtfparse.txniddesignator.as_str()))
                    .and_then(|x| x.split('"').nth(1))
                    .map(|s| s.to_string())
                    .unwrap_or_else(|| format!("{}:{}-{}", fields[0], fields[1], fields[2]));

                if names.contains_key(&entryname) {
                    let count = names
                        .get_mut(&entryname)
                        .expect("Entry name vanished from names map between check and get");
                    *count += 1;
                    entryname = format!("{}_r{}", entryname, count);
                } else {
                    names.insert(entryname.clone(), 0);
                }
                if !chroms.contains(&&fields[0].to_string()) {
                    println!(
                        "Warning, region {} not found in at least one of the bigwig/bam files. Skipping {}.",
                        fields[0], entryname
                    );
                } else {
                    regions.push(Region {
                        chrom: fields[0].to_string(),  //chrom
                        start: Revalue::U(start),      //start
                        end: Revalue::U(end),          //end
                        score: fields[5].to_string(),  //score
                        strand: fields[6].to_string(), //strand
                        name: entryname,               //region name
                        regionlength: end - start,     // regionlength
                    });
                    entries += 1;
                }
            }
        }
    }
    let filename = Path::new(gtf_file)
        .file_stem()
        .unwrap_or_else(|| {
            panic!(
                "Could not determine a file stem/label for GTF file '{}'",
                gtf_file
            )
        })
        .to_string_lossy()
        .into_owned();

    return (regions, (filename, entries));
}

pub fn read_bedfile(
    bed_file: &String,
    metagene: bool,
    chroms: &HashMap<String, u32>,
) -> (Vec<Region>, (String, u32)) {
    // read a provided bed_file into a vec of Region
    // Additional return is the filename and the number of entries (for sorting later on if needed).

    let mut regions: Vec<Region> = Vec::new();
    let mut names: HashMap<String, u32> = HashMap::new();
    let mut nonbed12: bool = false;
    let mut entries: u32 = 0;

    let bedfile = open_bed_or_gtf_reader(bed_file);

    for line in bedfile.lines() {
        let line =
            line.unwrap_or_else(|e| panic!("Failed to read a line from '{}': {}", bed_file, e));
        let fields: Vec<&str> = line.split('\t').collect();
        // Depending on bedfile, we have either BED3, BED6 or BED12
        // Note that this approach could allow somebody to have a 'mixed' bedfile, why not.
        match fields.len() {
            3 => {
                if !nonbed12 {
                    nonbed12 = true;
                }
                let chrom = fields[0];
                let mut entryname = format!("{}:{}-{}", fields[0], fields[1], fields[2]);
                // Only check for valid chroms, if chroms vec is not empty.

                let chromlen = match chroms.get(chrom) {
                    Some(len) => *len,
                    None => {
                        println!(
                            "Warning, region {} not found in at least one of the bigwig/bam files. Skipping {}.",
                            chrom, entryname
                        );
                        continue;
                    }
                };
                let start: u32 = fields[1].parse().unwrap_or_else(|e| {
                    panic!(
                        "Failed to parse start position in BED line '{}': {}",
                        line, e
                    )
                });
                let mut end: u32 = fields[2].parse().unwrap_or_else(|e| {
                    panic!("Failed to parse end position in BED line '{}': {}", line, e)
                });
                if start >= chromlen {
                    println!(
                        "Warning, region {} lies entirely beyond the end of {} ({}bp). Skipping.",
                        entryname, chrom, chromlen
                    );
                    continue;
                }
                end = end.min(chromlen);
                if names.contains_key(&entryname) {
                    let count = names
                        .get_mut(&entryname)
                        .expect("Entry name vanished from names map between check and get");
                    *count += 1;
                    entryname = format!("{}_r{}", entryname, count);
                } else {
                    names.insert(entryname.clone(), 0);
                }
                regions.push(Region {
                    chrom: fields[0].to_string(), //chrom
                    start: Revalue::U(start),     //start
                    end: Revalue::U(end),         //end
                    score: ".".to_string(),       //score
                    strand: ".".to_string(),      //score
                    name: entryname,              //region name
                    regionlength: end - start,    // regionlength
                });
                entries += 1;
            }
            6 => {
                if !nonbed12 {
                    nonbed12 = true;
                }
                let chrom = fields[0];
                let mut entryname = fields[3].to_string();
                let chromlen = match chroms.get(chrom) {
                    Some(len) => *len,
                    None => {
                        println!(
                            "Warning, region {} not found in at least one of the bigwig/bam files. Skipping {}.",
                            chrom, entryname
                        );
                        continue;
                    }
                };
                let start: u32 = fields[1].parse().unwrap_or_else(|e| {
                    panic!(
                        "Failed to parse start position in BED line '{}': {}",
                        line, e
                    )
                });
                let mut end: u32 = fields[2].parse().unwrap_or_else(|e| {
                    panic!("Failed to parse end position in BED line '{}': {}", line, e)
                });
                if start >= chromlen {
                    println!(
                        "Warning, region {} lies entirely beyond the end of {} ({}bp). Skipping.",
                        entryname, chrom, chromlen
                    );
                    continue;
                }
                end = end.min(chromlen);
                if names.contains_key(&entryname) {
                    let count = names
                        .get_mut(&entryname)
                        .expect("Entry name vanished from names map between check and get");
                    *count += 1;
                    entryname = format!("{}_r{}", entryname, count);
                } else {
                    names.insert(entryname.clone(), 0);
                }
                regions.push(Region {
                    chrom: fields[0].to_string(),  //chrom
                    start: Revalue::U(start),      //start
                    end: Revalue::U(end),          //end
                    score: fields[4].to_string(),  //score
                    strand: fields[5].to_string(), //strand
                    name: entryname,               //region name
                    regionlength: end - start,     // regionlength
                });
                entries += 1;
            }
            12 => {
                let chrom = fields[0];
                let mut entryname = fields[3].to_string();
                let chromlen = match chroms.get(chrom) {
                    Some(len) => *len,
                    None => {
                        println!(
                            "Warning, region {} not found in at least one of the bigwig/bam files. Skipping {}.",
                            chrom, entryname
                        );
                        continue;
                    }
                };
                let feat_start: u32 = fields[1].parse().unwrap_or_else(|e| {
                    panic!(
                        "Failed to parse start position in BED12 line '{}': {}",
                        line, e
                    )
                });
                if feat_start >= chromlen {
                    println!(
                        "Warning, region {} lies entirely beyond the end of {} ({}bp). Skipping.",
                        entryname, chrom, chromlen
                    );
                    continue;
                }
                if names.contains_key(&entryname) {
                    let count = names
                        .get_mut(&entryname)
                        .expect("Entry name vanished from names map between check and get");
                    *count += 1;
                    entryname = format!("{}_r{}", entryname, count);
                } else {
                    names.insert(entryname.clone(), 0);
                }
                if metagene {
                    let start: u32 = feat_start;
                    let blocksizes: Vec<u32> = fields[10]
                        .split(',')
                        .filter(|x| !x.is_empty())
                        .map(|x| {
                            x.parse().unwrap_or_else(|e| {
                                panic!("Failed to parse blockSizes in BED12 line '{}': {}", line, e)
                            })
                        })
                        .collect();
                    let blockstarts: Vec<u32> = fields[11]
                        .split(',')
                        .filter(|x| !x.is_empty())
                        .map(|x| {
                            x.parse::<u32>().unwrap_or_else(|e| {
                                panic!(
                                    "Failed to parse blockStarts in BED12 line '{}': {}",
                                    line, e
                                )
                            }) + start
                        })
                        .collect();

                    let (starts, ends): (Vec<u32>, Vec<u32>) = blocksizes
                        .into_iter()
                        .zip(blockstarts.into_iter())
                        .map(|(s, start)| (start, start + s))
                        .into_iter()
                        .unzip();
                    // Clip exons to the chromosome length: drop any exon that starts
                    // beyond it, and truncate one that only partially overflows.
                    let (starts, ends): (Vec<u32>, Vec<u32>) = starts
                        .into_iter()
                        .zip(ends.into_iter())
                        .filter(|&(s, _)| s < chromlen)
                        .map(|(s, e)| (s, e.min(chromlen)))
                        .unzip();
                    let length: u32 = starts.iter().zip(ends.iter()).map(|(s, e)| e - s).sum();
                    regions.push(Region {
                        chrom: fields[0].to_string(),  //chrom
                        start: Revalue::V(starts),     //start
                        end: Revalue::V(ends),         //end
                        score: fields[4].to_string(),  //score
                        strand: fields[5].to_string(), //strand
                        name: entryname,               //region name
                        regionlength: length,          // regionlength
                    });
                    entries += 1;
                } else {
                    let start = feat_start;
                    let end: u32 = fields[2]
                        .parse::<u32>()
                        .unwrap_or_else(|e| {
                            panic!(
                                "Failed to parse end position in BED12 line '{}': {}",
                                line, e
                            )
                        })
                        .min(chromlen);
                    regions.push(Region {
                        chrom: fields[0].to_string(),  //chrom
                        start: Revalue::U(start),      //start
                        end: Revalue::U(end),          //end
                        score: fields[4].to_string(),  //score
                        strand: fields[5].to_string(), //strand
                        name: entryname,               //region name
                        regionlength: end - start,     // regionlength
                    });
                    entries += 1;
                }
            }
            _ => panic!("Invalid BED format. BED file doesn't have 3, 6 or 12 fields."),
        }
    }

    let filename = Path::new(bed_file)
        .file_stem()
        .unwrap_or_else(|| {
            panic!(
                "Could not determine a file stem/label for BED file '{}'",
                bed_file
            )
        })
        .to_string_lossy()
        .into_owned();

    if metagene && nonbed12 {
        println!(
            "Warning: Metagene analysis is requested, but not all bedfiles and/or bedfile entries are in BED12 format. Proceed at your own risk."
        );
    }
    return (regions, (filename, entries));
}

pub fn chrombounds_from_bw(bwfile: &str) -> HashMap<String, u32> {
    // define chromsizes hashmap
    let mut chromsizes: HashMap<String, u32> = HashMap::new();
    let bwf = File::open(bwfile)
        .unwrap_or_else(|e| panic!("Failed to open bigwig file '{}': {}", bwfile, e));
    let reader = BigWigRead::open(bwf)
        .unwrap_or_else(|e| panic!("Failed to parse bigwig file '{}': {}", bwfile, e));
    for chrom in reader.chroms() {
        chromsizes.insert(chrom.name.clone(), chrom.length);
    }
    chromsizes
}

pub fn chrombounds_from_bam(bamfiles: Vec<&str>) -> HashMap<String, u32> {
    let mut found_chroms: HashMap<String, usize> = HashMap::new();
    for bam in bamfiles.iter() {
        let bamreader = IndexedReader::from_path(bam)
            .unwrap_or_else(|e| panic!("Failed to open indexed BAM file '{}': {}", bam, e));
        let chroms: Vec<String> = bamreader
            .header()
            .target_names()
            .iter()
            .map(|x| {
                String::from_utf8(x.to_vec()).unwrap_or_else(|e| {
                    panic!(
                        "BAM header for '{}' has a non-UTF-8 chromosome name: {}",
                        bam, e
                    )
                })
            })
            .collect();
        for chrom in chroms.iter() {
            // if it's not in the hashmap, add it, else increment count
            if !found_chroms.contains_key(chrom) {
                found_chroms.insert(chrom.clone(), 1);
            } else {
                let count = found_chroms.get_mut(chrom).expect(
                    "Chromosome key vanished from found_chroms map between check and update",
                );
                *count += 1;
            }
        }
    }
    let mut validchroms: Vec<String> = Vec::new();
    // loop over all chroms in the hashmap, if the count is expected, include them
    for (chrom, count) in found_chroms.iter() {
        if *count == bamfiles.len() {
            validchroms.push(chrom.clone());
        } else {
            println!(
                "Chromosome {} is missing in at least one bam file, and thus ignored!",
                chrom
            );
        }
    }

    let mut chromsizes: HashMap<String, u32> = HashMap::new();
    for bamfile in bamfiles.iter() {
        let bamreader = IndexedReader::from_path(bamfile)
            .unwrap_or_else(|e| panic!("Failed to open indexed BAM file '{}': {}", bamfile, e));
        let header = bamreader.header().clone();
        for tid in 0..header.target_count() {
            let chromname = String::from_utf8(header.tid2name(tid).to_vec())
                .expect("Invalid UTF-8 in chromosome name");
            if !validchroms.contains(&chromname) {
                continue;
            }
            let chromlen = header
                .target_len(tid)
                .expect("Error retrieving length for chromosome") as u32;
            chromsizes
                .entry(chromname)
                .and_modify(|len| *len = (*len).min(chromlen))
                .or_insert(chromlen);
        }
    }
    chromsizes
}

pub fn bwintervals(
    bwfile: &str,
    regions: &[Region],
    slopregions: &[Vec<Bin>],
    scale_regions: &Scalingregions,
    blacklist_index: Option<&crate::filtering::BlacklistIndex>,
) -> Vec<Vec<f32>> {
    // For a given bw file, a vector of slopregions (Bin enum))
    // return a vector with for every region a vector of f64.

    // Make sure regions and slopregions are of equal length
    assert_eq!(
        regions.len(),
        slopregions.len(),
        "Regions from bed file and parsed regions (slopped) do not have equal length. Something went wrong during computation."
    );

    // Define return vector, set up bw reader.
    let mut bwvals: Vec<Vec<f32>> = Vec::new();
    let bwf = File::open(bwfile)
        .unwrap_or_else(|e| panic!("Failed to open bigwig file '{}': {}", bwfile, e));
    let mut reader = BigWigRead::open(bwf)
        .unwrap_or_else(|e| panic!("Failed to parse bigwig file '{}': {}", bwfile, e));

    // Iterate over regions and slopregions synchronously.
    for (sls, region) in slopregions.iter().zip(regions.iter()) {
        let mut bwval: Vec<f32> = Vec::new();
        // get 'min and max' to query
        // at some point this should become an impl but man am I tired.
        let (min, max) = sls
            .iter()
            .flat_map(|bin| match bin {
                Bin::Conbin(a, b) => vec![*a, *b],
                Bin::PaddedConbin(a, b, _) => vec![*a, *b],
                Bin::Catbin(pairs) => pairs
                    .iter()
                    .flat_map(|(x, y)| vec![*x, *y])
                    .collect::<Vec<u32>>(),
                Bin::PaddedCatbin(pairs, _) => pairs
                    .iter()
                    .flat_map(|(x, y)| vec![*x, *y])
                    .collect::<Vec<u32>>(),
            })
            .fold((u32::MAX, u32::MIN), |(min, max), x| {
                (min.min(x), max.max(x))
            });
        let binvals = reader
            .get_interval(&region.chrom, min as u32, max as u32)
            .unwrap_or_else(|e| {
                panic!(
                    "Failed to query bigwig '{}' for interval {}:{}-{}: {}",
                    bwfile, region.chrom, min, max, e
                )
            });

        let span = max.saturating_sub(min) as usize;
        let mut bwvec: Vec<f32> = vec![f32::NAN; span];
        for interval in binvals {
            let interval = interval.unwrap_or_else(|e| {
                panic!("Failed to read an interval from bigwig '{}': {}", bwfile, e)
            });
            let start = (interval.start as u32).max(min);
            let end = (interval.end as u32).min(max);
            let val = interval.value as f32;
            if start < end {
                bwvec[(start - min) as usize..(end - min) as usize].fill(val);
            }
        }
        let gather_vals = |a: u32, b: u32, vals: &mut Vec<f32>| {
            for bp in a..b {
                if blacklist_index.is_some()
                    && blacklist_index
                        .expect("blacklist_index vanished between is_some() check and use")
                        .contains(&region.chrom, bp)
                {
                    continue;
                }
                let v = bwvec[(bp - min) as usize];
                if v.is_nan() {
                    if scale_regions.missingdata_as_zero {
                        vals.push(0.0);
                    }
                } else {
                    vals.push(v);
                }
            }
        };

        // Now we can iterate over the slopped regions, and get the values from the hashmap.
        for bin in sls {
            match bin {
                Bin::Conbin(a, b) => {
                    if a == b && *b == 0 {
                        if scale_regions.missingdata_as_zero {
                            bwval.push(0.0);
                        } else {
                            bwval.push(std::f32::NAN);
                        }
                    } else {
                        let mut vals: Vec<f32> = Vec::new();
                        gather_vals(*a, *b, &mut vals);
                        if vals.is_empty() {
                            if scale_regions.missingdata_as_zero {
                                bwval.push(0.0);
                            } else {
                                bwval.push(std::f32::NAN);
                            }
                        } else {
                            let valrefs: Vec<&f32> = vals.iter().collect();
                            let val = match scale_regions.avgtype.as_str() {
                                "mean" => mean_float(&valrefs),
                                "median" => median_float(&valrefs),
                                "min" => min_float(&valrefs),
                                "max" => max_float(&valrefs),
                                "std" => std_float(&valrefs),
                                "sum" => sum_float(&valrefs),
                                _ => panic!("Unknown avgtype."),
                            };
                            bwval.push(val);
                        }
                    }
                }
                Bin::PaddedConbin(a, b, pad) => {
                    let mut vals: Vec<f32> = Vec::new();
                    gather_vals(*a, *b, &mut vals);
                    if scale_regions.missingdata_as_zero && *pad > 0 {
                        vals.extend(std::iter::repeat(0.0f32).take(*pad as usize));
                    }
                    if vals.is_empty() {
                        if scale_regions.missingdata_as_zero {
                            bwval.push(0.0);
                        } else {
                            bwval.push(std::f32::NAN);
                        }
                    } else {
                        let valrefs: Vec<&f32> = vals.iter().collect();
                        let val = match scale_regions.avgtype.as_str() {
                            "mean" => mean_float(&valrefs),
                            "median" => median_float(&valrefs),
                            "min" => min_float(&valrefs),
                            "max" => max_float(&valrefs),
                            "std" => std_float(&valrefs),
                            "sum" => sum_float(&valrefs),
                            _ => panic!("Unknown avgtype."),
                        };
                        bwval.push(val);
                    }
                }
                Bin::Catbin(pairs) => {
                    let mut vals: Vec<f32> = Vec::new();

                    for (start, end) in pairs {
                        if start == end && *end == 0 {
                            continue;
                        }
                        gather_vals(*start, *end, &mut vals);
                    }

                    // Handle case where no valid values exist
                    if vals.is_empty() {
                        if scale_regions.missingdata_as_zero {
                            bwval.push(0.0);
                        } else {
                            bwval.push(std::f32::NAN);
                        }
                    } else {
                        let valrefs: Vec<&f32> = vals.iter().collect();
                        let val = match scale_regions.avgtype.as_str() {
                            "mean" => mean_float(&valrefs),
                            "median" => median_float(&valrefs),
                            "min" => min_float(&valrefs),
                            "max" => max_float(&valrefs),
                            "std" => std_float(&valrefs),
                            "sum" => sum_float(&valrefs),
                            _ => panic!("Unknown avgtype."),
                        };
                        bwval.push(val);
                    }
                }
                Bin::PaddedCatbin(pairs, pad) => {
                    let mut vals: Vec<f32> = Vec::new();

                    for (start, end) in pairs {
                        if start == end && *end == 0 {
                            continue;
                        }
                        gather_vals(*start, *end, &mut vals);
                    }
                    if scale_regions.missingdata_as_zero && *pad > 0 {
                        vals.extend(std::iter::repeat(0.0f32).take(*pad as usize));
                    }

                    if vals.is_empty() {
                        if scale_regions.missingdata_as_zero {
                            bwval.push(0.0);
                        } else {
                            bwval.push(std::f32::NAN);
                        }
                    } else {
                        let valrefs: Vec<&f32> = vals.iter().collect();
                        let val = match scale_regions.avgtype.as_str() {
                            "mean" => mean_float(&valrefs),
                            "median" => median_float(&valrefs),
                            "min" => min_float(&valrefs),
                            "max" => max_float(&valrefs),
                            "std" => std_float(&valrefs),
                            "sum" => sum_float(&valrefs),
                            _ => panic!("Unknown avgtype."),
                        };
                        bwval.push(val);
                    }
                }
            }
        }
        assert_eq!(
            bwval.len(),
            scale_regions.cols_expected / scale_regions.bwfiles
        );
        bwvals.push(bwval);
    }
    bwvals
}

pub fn header_matrix(
    scale_regions: &Scalingregions,
    regionsizes: &HashMap<String, u32>,
    sortregions: &str,
    sortusing: &str,
) -> String {
    // Create the header for the matrix.
    // This is quite ugly, but this is mainly because we need to accomodate delta bwfiles.
    let mut headstr = String::new();
    headstr.push_str("@{");
    headstr.push_str(&format!(
        "\"upstream\":[{}],",
        (0..scale_regions.bwfiles)
            .map(|_| scale_regions.upstream)
            .collect::<Vec<_>>()
            .into_iter()
            .join(",")
    ));
    headstr.push_str(&format!(
        "\"downstream\":[{}],",
        (0..scale_regions.bwfiles)
            .map(|_| scale_regions.downstream)
            .collect::<Vec<_>>()
            .into_iter()
            .join(",")
    ));
    headstr.push_str(&format!(
        "\"body\":[{}],",
        (0..scale_regions.bwfiles)
            .map(|_| scale_regions.regionbodylength)
            .collect::<Vec<_>>()
            .into_iter()
            .join(",")
    ));
    headstr.push_str(&format!(
        "\"bin size\":[{}],",
        (0..scale_regions.bwfiles)
            .map(|_| scale_regions.binsize)
            .collect::<Vec<_>>()
            .into_iter()
            .join(",")
    ));
    // ref point can be empty (for scale_regions, for example).
    // To keep compatibility with deepTools 3 it should be written as null
    let refpointstring = (0..scale_regions.bwfiles)
        .map(|_| {
            if scale_regions.referencepoint.is_empty() {
                "null".to_string()
            } else {
                format!("\"{}\"", scale_regions.referencepoint)
            }
        })
        .join(",");
    headstr.push_str(&format!("\"ref point\":[{}],", refpointstring));

    headstr.push_str(&format!("\"verbose\":{},", scale_regions.verbose));
    headstr.push_str(&format!("\"bin avg type\":\"{}\",", scale_regions.avgtype));
    headstr.push_str(&format!(
        "\"missing data as zero\":{},",
        scale_regions.missingdata_as_zero
    ));
    // Unimplemented arguments, but they need to be present in header for now anyway.
    let minthresh_str = if scale_regions.minthresh != 0.0 {
        format!("{}", scale_regions.minthresh)
    } else {
        "null".to_string()
    };
    let maxthresh_str = if scale_regions.maxthresh != 0.0 {
        format!("{}", scale_regions.maxthresh)
    } else {
        "null".to_string()
    };
    headstr.push_str(
        &format!("\"min threshold\":{},\"max threshold\":{},\"scale\":{},\"skip zeros\":{},\"nan after end\":{},", minthresh_str, maxthresh_str, scale_regions.scale, scale_regions.skipzero, scale_regions.nan_after_end)
    );
    headstr.push_str(&format!("\"proc number\":{},", scale_regions.proc_number));
    headstr.push_str(&format!(
        "\"sort regions\":\"{}\",\"sort using\":\"{}\",",
        sortregions, sortusing
    ));
    if scale_regions.referencepoint.is_empty() {
        headstr.push_str(&format!(
            "\"startLabel\":\"{}\",\"endLabel\":\"{}\",",
            scale_regions.startlabel, scale_regions.endlabel
        ));
    }

    headstr.push_str(&format!(
        "\"unscaled 5 prime\":[{}],",
        (0..scale_regions.bwfiles)
            .map(|_| scale_regions.unscaled5prime)
            .collect::<Vec<_>>()
            .into_iter()
            .join(",")
    ));
    headstr.push_str(&format!(
        "\"unscaled 3 prime\":[{}],",
        (0..scale_regions.bwfiles)
            .map(|_| scale_regions.unscaled3prime)
            .collect::<Vec<_>>()
            .into_iter()
            .join(",")
    ));
    headstr.push_str(&format!(
        "\"group_labels\":[\"{}\"],",
        scale_regions.regionlabels.join("\",\"")
    ));
    // Get cumulative sizes of regions
    let mut groupbounds: Vec<u32> = Vec::new();
    groupbounds.push(0);
    let mut cumsum: u32 = 0;
    for regionlabel in scale_regions.regionlabels.iter() {
        cumsum += regionsizes.get(regionlabel).unwrap_or_else(|| {
            panic!(
                "Region label '{}' not found in regionsizes map",
                regionlabel
            )
        });
        groupbounds.push(cumsum);
    }
    let groupbounds = format!(
        "{}",
        groupbounds
            .iter()
            .map(|&x| x.to_string())
            .collect::<Vec<String>>()
            .join(",")
    );

    headstr.push_str(&format!("\"group_boundaries\":[{}],", groupbounds));
    // Sample labels
    headstr.push_str(&format!(
        "\"sample_labels\":[\"{}\"],",
        scale_regions.bwlabels.join("\",\"")
    ));
    // Get cumulative sizes for sample boundaries
    let colsize_per_sample = scale_regions.cols_expected / scale_regions.bwfiles;
    let mut samplebounds: Vec<u64> = Vec::new();
    samplebounds.push(0);
    let mut cumsum: u64 = 0;
    for _ in 0..scale_regions.bwfiles {
        cumsum += colsize_per_sample as u64;
        samplebounds.push(cumsum);
    }
    let samplebounds = format!(
        "{}",
        samplebounds
            .iter()
            .map(|&x| x.to_string())
            .collect::<Vec<String>>()
            .join(",")
    );

    headstr.push_str(&format!("\"sample_boundaries\":[{}]", samplebounds));
    headstr.push_str("}\n");
    headstr
}

pub fn write_matrix(
    header: String,
    mat: &[Vec<f32>],
    ofile: &str,
    regions: &[Region],
    scale_regions: &Scalingregions,
) {
    // Write out the matrix to a compressed file.
    let omat = File::create(ofile)
        .unwrap_or_else(|e| panic!("Failed to create output matrix file '{}': {}", ofile, e));
    let mut encoder = GzEncoder::new(omat, Compression::default());
    encoder
        .write_all(header.as_bytes())
        .unwrap_or_else(|e| panic!("Failed to write header to matrix file '{}': {}", ofile, e));
    assert_eq!(regions.len(), mat.len());
    for (region, row) in regions.iter().zip(mat.iter()) {
        let mut writerow = format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t",
            region.chrom,            // Chromosome
            region.start.rewrites(), // Revalue for start (either u32, or Vec<u32>)
            region.end.rewrites(),   // Revalue for end (either u32, or Vec<u32>)
            region.name, // String field for name. Duplicates taken care of in read_bedfiles.
            region.score, // Score field persisted from bedfile
            region.strand, // Strand field persisted from bedfile
        );
        writerow.push_str(
            &row.iter()
                .map(|x| {
                    if x.is_nan() {
                        // as deepTools 3 encoded nan in matrix.
                        "nan".to_string()
                    } else {
                        format!("{:.6}", scale_regions.scale * x)
                    }
                })
                .collect::<Vec<String>>()
                .join("\t"),
        );
        writerow.push_str("\n");
        encoder
            .write_all(writerow.as_bytes())
            .unwrap_or_else(|e| panic!("Failed to write row to matrix file '{}': {}", ofile, e));
    }
}

pub fn write_matrix_values(
    file_name: &str,
    mat: &[Vec<f32>],
    scale_regions: &Scalingregions,
    regionsizes: &HashMap<String, u32>,
) {
    use std::io::Write;
    let mut fh = File::create(file_name)
        .unwrap_or_else(|e| panic!("Failed to create matrix values file '{}': {}", file_name, e));

    // Header line 1: group labels with region counts
    let info: Vec<String> = scale_regions
        .regionlabels
        .iter()
        .map(|label| format!("{}:{}", label, regionsizes.get(label).unwrap_or(&0)))
        .collect();
    fh.write_all(format!("#{}\n", info.join("\t")).as_bytes())
        .unwrap_or_else(|e| {
            panic!(
                "Failed to write header to matrix values file '{}': {}",
                file_name, e
            )
        });

    // Header line 2: region dimension parameters
    let header2 = format!(
        "#downstream:{}\tupstream:{}\tbody:{}\tbin size:{}\tunscaled 5 prime:{}\tunscaled 3 prime:{}\n",
        scale_regions.downstream,
        scale_regions.upstream,
        scale_regions.regionbodylength,
        scale_regions.binsize,
        scale_regions.unscaled5prime,
        scale_regions.unscaled3prime,
    );
    fh.write_all(header2.as_bytes()).unwrap_or_else(|e| {
        panic!(
            "Failed to write header line 2 to matrix values file '{}': {}",
            file_name, e
        )
    });

    // Header line 3: sample labels repeated per column
    let cols_per_sample = scale_regions.cols_expected / scale_regions.bwfiles;
    let sample_info: Vec<String> = scale_regions
        .bwlabels
        .iter()
        .flat_map(|label| (0..cols_per_sample).map(move |_| label.clone()))
        .collect();
    fh.write_all(format!("{}\t{}\n", info.join("\t"), sample_info.join("\t")).as_bytes())
        .unwrap_or_else(|e| {
            panic!(
                "Failed to write header line 3 to matrix values file '{}': {}",
                file_name, e
            )
        });
    fh.flush()
        .unwrap_or_else(|e| panic!("Failed to flush matrix values file '{}': {}", file_name, e));

    // Reopen in append mode and write matrix data
    let mut fh = std::fs::OpenOptions::new()
        .create(false)
        .append(true)
        .open(file_name)
        .unwrap_or_else(|e| {
            panic!(
                "Failed to reopen matrix values file '{}' for appending: {}",
                file_name, e
            )
        });
    for row in mat.iter() {
        let line = row
            .iter()
            .map(|x| {
                if x.is_nan() {
                    "nan".to_string()
                } else {
                    format!("{:.6}", scale_regions.scale * x)
                }
            })
            .collect::<Vec<_>>()
            .join("\t");
        writeln!(fh, "{}", line)
            .unwrap_or_else(|e| panic!("Failed to write matrix row to '{}': {}", file_name, e));
    }
}

pub fn write_sorted_regions_bed(
    file_name: &str,
    regions: &[Region],
    scale_regions: &Scalingregions,
    regionsizes: &HashMap<String, u32>,
) {
    use std::io::Write;
    let mut fh = File::create(file_name).unwrap_or_else(|e| {
        panic!(
            "Failed to create sorted regions BED file '{}': {}",
            file_name, e
        )
    });
    let header = "#chrom\tstart\tend\tname\tscore\tstrand\tthickStart\tthickEnd\titemRGB\tblockCount\tblockSizes\tblockStarts\tdeepTools_group\n";
    fh.write_all(header.as_bytes()).unwrap_or_else(|e| {
        panic!(
            "Failed to write header to sorted regions BED file '{}': {}",
            file_name, e
        )
    });
    // Build group boundaries from regionsizes (cumulative region counts per group label)
    let mut group_boundaries: Vec<u32> = vec![0];
    let mut cumsum: u32 = 0;
    for label in &scale_regions.regionlabels {
        cumsum += *regionsizes.get(label).unwrap_or(&0);
        group_boundaries.push(cumsum);
    }
    for (idx, region) in regions.iter().enumerate() {
        let label_idx = group_boundaries
            .iter()
            .take_while(|&&b| b <= idx as u32)
            .count()
            .saturating_sub(1)
            .min(scale_regions.regionlabels.len() - 1);
        let start_first = match &region.start {
            Revalue::U(v) => *v,
            Revalue::V(vs) => *vs.first().unwrap_or_else(|| {
                panic!("Region '{}' has an empty exon-start vector", region.name)
            }),
        };
        let end_last = match &region.end {
            Revalue::U(v) => *v,
            Revalue::V(vs) => *vs
                .last()
                .unwrap_or_else(|| panic!("Region '{}' has an empty exon-end vector", region.name)),
        };
        let (block_count, block_sizes, block_starts) = match (&region.start, &region.end) {
            (Revalue::U(_), _) => (
                "1".to_string(),
                (end_last - start_first).to_string(),
                "0".to_string(),
            ),
            (Revalue::V(starts), Revalue::V(ends)) => {
                let bc = starts.len().to_string();
                let sz: Vec<String> = starts
                    .iter()
                    .zip(ends.iter())
                    .map(|(s, e)| (e - s).to_string())
                    .collect();
                let st: Vec<String> = starts
                    .iter()
                    .map(|s| (s - start_first).to_string())
                    .collect();
                (bc, sz.join(","), st.join(","))
            }
            _ => (
                "1".to_string(),
                (end_last - start_first).to_string(),
                "0".to_string(),
            ),
        };
        let group_label = &scale_regions.regionlabels[label_idx];
        let score_str = region
            .score
            .parse::<f64>()
            .map(|v| format!("{:.1}", v))
            .unwrap_or_else(|_| region.score.clone());

        let line = format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t0\t{}\t{}\t{}\t{}\n",
            region.chrom,
            start_first,
            end_last,
            region.name,
            score_str,
            region.strand,
            start_first,
            end_last,
            block_count,
            block_sizes,
            block_starts,
            group_label,
        );
        fh.write_all(line.as_bytes()).unwrap_or_else(|e| {
            panic!(
                "Failed to write line to sorted regions BED file '{}': {}",
                file_name, e
            )
        });
    }
}
