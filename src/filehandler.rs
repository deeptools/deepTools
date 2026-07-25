use rust_htslib::bam::{Read, Reader, IndexedReader};
use itertools::Itertools;
use std::io::{BufReader, BufWriter, Write};
use std::io::prelude::*;
use std::fs::File;
use std::path::Path;
use bigtools::{BigWigRead, BigWigWrite, Value};
use bigtools::beddata::BedParserStreamingIterator;
use std::collections::HashMap;
use flate2::write::GzEncoder;
use flate2::Compression;
use crate::covcalc::{Scalingregions, Gtfparse, Revalue, Region, Bin};
use crate::calc::{mean_float, median_float, min_float, max_float, sum_float, std_float};

/// Format f32 with 4 significant digits, matching Python's `%.4g`.
fn format_g4(x: f32) -> String {
    if x == 0.0 {
        return "0".to_string();
    }
    let abs = x.abs();
    let exp = abs.log10().floor() as i32;
    // Decide which format to use
    if exp >= -4 && exp < 4 {
        let prec = 4 - exp as usize - 1; // significant digits after decimal = 4 - (digits before decimal)
        if prec < 0 { return format!("{}", x); }
        let s = format!("{:.*}", prec, x);
        // Strip trailing zeros after decimal point
        if s.contains('.') {
            s.trim_end_matches('0').trim_end_matches('.').to_string()
        } else {
            s
        }
    } else {
        // Use scientific notation
        let prec = 3; // 4 sig figs = 1 before dot + 3 after
        let s = format!("{:.*e}", prec, x);
        // Rust formats as 1.234e+01, normalize to 1.234e+1
        if let Some(pos) = s.find('e') {
            let mantissa = &s[..pos];
            let exp_str = &s[pos + 1..];
            let exp_val: i32 = exp_str.parse().unwrap_or(0);
            let sign = if exp_val >= 0 { '+' } else { '-' };
            format!("{}{}e{}{}", mantissa.trim_end_matches('0').trim_end_matches('.'), sign, exp_val.abs(), "")
        } else {
            s
        }
    }
}

pub fn bam_ispaired(bam_ifile: &str) -> bool {
    let mut bam = Reader::from_path(bam_ifile).unwrap();
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

pub fn write_covfile<LI>(lines: LI, ofile: &str, filetype: &str, chromsizes: HashMap<String, u32>)
where
    LI: Iterator<Item = (String, Value)>,
 {
    if filetype == "bedgraph" {
        // write output file, bedgraph
        let mut writer = BufWriter::new(File::create(ofile).unwrap());
        for (chrom, val) in lines {
            writeln!(writer, "{}\t{}\t{}\t{}", chrom, val.start, val.end, val.value).unwrap();
        }
    } else {
        let vals = BedParserStreamingIterator::wrap_infallible_iter(
            lines,
            false
        );
        let runtime = tokio::runtime::Builder::new_multi_thread()
            .worker_threads(1)
            .build()
            .expect("Unable to create tokio runtime for bw writing.");
        let writer = BigWigWrite::create_file(ofile, chromsizes).unwrap();
        let _ = writer.write(vals, runtime);
    }
}

pub fn is_bed_or_gtf(fp: &str) -> String {
    // Check if file is a bed or gtf file.
    let file = File::open(fp).expect(format!("Failed to open file: {}", fp).as_str());
    let reader = BufReader::new(file);
    // Get the first line that doesn't start with #
    for line in reader.lines() {
        let line = line.unwrap();
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

pub fn read_gtffile(gtf_file: &String, gtfparse: &Gtfparse, chroms: Vec<&String>) -> (Vec<Region>, (String, u32)) {
    // At some point this zoo of String clones should be refactored. Not now though, We have a deadline.
    let mut regions: Vec<Region> = Vec::new();
    let mut names: HashMap<String, u32> = HashMap::new();
    let mut entries: u32 = 0;
    let mut txnids: Vec<String> = Vec::new();

    let gtffile = BufReader::new(File::open(gtf_file).unwrap());

    if gtfparse.metagene {
        // metagene implementation - more work here.
        let mut txn_hash: HashMap<String, Vec<(u32, u32)>> = HashMap::new();
        let mut txn_strand: HashMap<String, String> = HashMap::new();
        let mut txn_chrom: HashMap<String, String> = HashMap::new();

        for line in gtffile.lines() {
            let line = line.unwrap();
            // skip comments
            if line.starts_with('#') {
                continue;
            }
            let fields: Vec<&str> = line.split('\t').collect();
            if fields[2].to_string() == gtfparse.exonid {
                let start: u32 = fields[3].parse().unwrap();
                let end: u32 = fields[4].parse().unwrap();
                let txnid = fields[8]
                    .split(';')
                    .find(|x| x.trim().starts_with(gtfparse.txniddesignator.as_str()))
                    .and_then(|x| x.split('"').nth(1))
                    .map(|s| s.to_string())
                    .unwrap();
                if !txnids.contains(&txnid) {
                    txnids.push(txnid.clone());
                }

                let txnentry = txn_hash.entry(txnid.clone()).or_insert(Vec::new());
                // Just to verify all exons are on the same strand.
                if txn_strand.contains_key(&txnid) {
                    assert_eq!(txn_strand.get(&txnid).unwrap(), fields[6]);
                } else {
                    txn_strand.insert(txnid.clone(), fields[6].to_string());
                }
                // Same for chromosome
                if txn_chrom.contains_key(&txnid) {
                    assert_eq!(txn_chrom.get(&txnid).unwrap(), fields[0]);
                } else {
                    txn_chrom.insert(txnid, fields[0].to_string());
                }
                txnentry.push((start, end));
            }
        }

        for txnid in txnids.into_iter() {
            
            let txnentry = txn_hash.get_mut(&txnid).unwrap();
            let length: u32 = txnentry.iter().map(|(s, e)| e - s).sum();
            txnentry.sort_by(|a, b| a.0.cmp(&b.0));
            let (starts, ends): (Vec<u32>, Vec<u32>) = txnentry.iter().map(|(s, e)| (*s, *e)).unzip();
            let chrom = txn_chrom.get(&txnid).unwrap().to_string();

            if !chroms.contains(&&chrom.to_string()) {
                println!("Warning, region {} not found in at least one of the bigwig/bam files. Skipping {}.", chrom, txnid);
            } else {
                regions.push(
                    Region {
                        chrom: txn_chrom.get(&txnid).unwrap().to_string(), //chrom
                        start: Revalue::V(starts), //start
                        end: Revalue::V(ends), //end
                        score: ".".to_string(), //score
                        strand: txn_strand.get(&txnid).unwrap().to_string(), // strand
                        name: txnid.to_string(), 
                        regionlength: length // regionlength
                    }
                );
                entries += 1;
            }
        }
    } else {
        // Take fields with col 3 == gtfparse.txnid, start, end
        for line in gtffile.lines() {
            let line = line.unwrap();
            // skip comments
            if line.starts_with('#') {
                continue;
            }

            let fields: Vec<&str> = line.split('\t').collect();
            if fields[2].to_string() == gtfparse.txnid {
                let start: u32 = fields[3].parse().unwrap();
                let end: u32 = fields[4].parse().unwrap();
                let mut entryname = fields[8]
                    .split(';')
                    .find(|x| x.trim().starts_with(gtfparse.txniddesignator.as_str()))
                    .and_then(|x| x.split('"').nth(1))
                    .map(|s| s.to_string())
                    .unwrap_or_else(|| format!("{}:{}-{}", fields[0], fields[1], fields[2]));
                
                if names.contains_key(&entryname) {
                    let count = names.get_mut(&entryname).unwrap();
                    *count += 1;
                    entryname = format!("{}_r{}", entryname, count);
                } else {
                    names.insert(entryname.clone(), 0);
                }
                if !chroms.contains(&&fields[0].to_string()) {
                    println!("Warning, region {} not found in at least one of the bigwig/bam files. Skipping {}.", fields[0], entryname);
                } else {
                    regions.push(
                        Region {
                            chrom: fields[0].to_string(), //chrom
                            start: Revalue::U(start), //start
                            end: Revalue::U(end), //end
                            score: fields[5].to_string(), //score
                            strand: fields[6].to_string(), //strand
                            name: entryname, //region name
                            regionlength: end - start // regionlength
                        }
                    );
                    entries += 1;
                }
            }
        }
    }
    let filename = Path::new(gtf_file)
        .file_stem()
        .unwrap()
        .to_string_lossy()
        .into_owned();

    return (regions, (filename, entries));
}

pub fn read_bedfile(bed_file: &String, metagene: bool, chroms: Vec<&String>) -> (Vec<Region>, (String, u32)) {
    // read a provided bed_file into a vec of Region
    // Additional return is the filename and the number of entries (for sorting later on if needed).

    let mut regions: Vec<Region> = Vec::new();
    let mut names: HashMap<String, u32> = HashMap::new();
    let mut nonbed12: bool = false;
    let mut entries: u32 = 0;

    let bedfile = BufReader::new(File::open(bed_file).unwrap());
    
    for line in bedfile.lines() {
        let line = line.unwrap();
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

                if !chroms.contains(&&chrom.to_string()) {
                    println!("Warning, region {} not found in at least one of the bigwig/bam files. Skipping {}.", chrom, entryname);
                    continue;
                }
                if names.contains_key(&entryname) {
                    let count = names.get_mut(&entryname).unwrap();
                    *count += 1;
                    entryname = format!("{}_r{}", entryname, count);
                } else {
                    names.insert(entryname.clone(), 0);
                }
                let start = fields[1].parse().unwrap();
                let end = fields[2].parse().unwrap();
                regions.push(
                    Region {
                        chrom: fields[0].to_string(), //chrom
                        start: Revalue::U(start), //start
                        end: Revalue::U(end), //end
                        score: ".".to_string(), //score
                        strand: ".".to_string(), //score
                        name: entryname, //region name
                        regionlength: end - start // regionlength
                    }
                );
                entries += 1;
            },
            6 => {
                if !nonbed12 {
                    nonbed12 = true;
                }
                let chrom = fields[0];
                let mut entryname = fields[3].to_string();
                if !chroms.contains(&&chrom.to_string()) {
                    println!("Warning, region {} not found in at least one of the bigwig/bam files. Skipping {}.", chrom, entryname);
                    continue;
                }
                let start = fields[1].parse().unwrap();
                let end = fields[2].parse().unwrap();
                if names.contains_key(&entryname) {
                    let count = names.get_mut(&entryname).unwrap();
                    *count += 1;
                    entryname = format!("{}_r{}", entryname, count);
                } else {
                    names.insert(entryname.clone(), 0);
                }
                regions.push(
                    Region {
                        chrom: fields[0].to_string(), //chrom
                        start: Revalue::U(start), //start
                        end: Revalue::U(end), //end
                        score: fields[4].to_string(), //score
                        strand: fields[5].to_string(), //strand
                        name: entryname, //region name
                        regionlength: end - start // regionlength
                    }
                );
                entries += 1;
            },
            12 => {
                let chrom = fields[0];
                let mut entryname = fields[3].to_string();
                if !chroms.contains(&&chrom.to_string()) {
                    println!("Warning, region {} not found in at least one of the bigwig/bam files. Skipping {}.", chrom, entryname);
                    continue;
                }
                if names.contains_key(&entryname) {
                    let count = names.get_mut(&entryname).unwrap();
                    *count += 1;
                    entryname = format!("{}_r{}", entryname, count);
                } else {
                    names.insert(entryname.clone(), 0);
                }
                if metagene {                        
                    let start: u32 = fields[1].parse().unwrap();
                    let blocksizes: Vec<u32> = fields[10]
                        .split(',')
                        .filter(|x| !x.is_empty())
                        .map(|x| x.parse().unwrap())
                        .collect();
                    let length: u32 = blocksizes.iter().sum();
                    let blockstarts: Vec<u32> = fields[11]
                        .split(',')
                        .filter(|x| !x.is_empty())
                        .map(|x| x.parse::<u32>().unwrap() + start)
                        .collect();

                    let (starts, ends) = blocksizes
                        .into_iter()
                        .zip(blockstarts.into_iter())
                        .map(|(s, start)| (start, start + s))
                        .into_iter()
                        .unzip();
                    regions.push(
                        Region {
                            chrom: fields[0].to_string(), //chrom
                            start: Revalue::V(starts), //start
                            end: Revalue::V(ends), //end
                            score: fields[4].to_string(), //score
                            strand: fields[5].to_string(), //strand
                            name: entryname, //region name
                            regionlength: length // regionlength
                        }
                    );
                    entries += 1;
                } else {
                    let start = fields[1].parse().unwrap();
                    let end = fields[2].parse().unwrap();
                    regions.push(
                        Region {
                            chrom: fields[0].to_string(), //chrom
                            start: Revalue::U(start), //start
                            end: Revalue::U(end), //end
                            score: fields[4].to_string(), //score
                            strand: fields[5].to_string(), //strand
                            name: entryname, //region name
                            regionlength: end - start // regionlength
                        }
                    );
                    entries += 1;
                }
            },
            _ => panic!("Invalid BED format. BED file doesn't have 3, 6 or 12 fields."),
        }
    }

    let filename = Path::new(bed_file)
        .file_stem()
        .unwrap()
        .to_string_lossy()
        .into_owned();

    if metagene && nonbed12 {
        println!("Warning: Metagene analysis is requested, but not all bedfiles and/or bedfile entries are in BED12 format. Proceed at your own risk.");
    }
    return (regions, (filename, entries));
}

pub fn chrombounds_from_bw(bwfile: &str) -> HashMap<String, u32> {
    // define chromsizes hashmap
    let mut chromsizes: HashMap<String, u32> = HashMap::new();
    let bwf = File::open(bwfile).expect("Failed to open bw file.");
    let reader = BigWigRead::open(bwf).unwrap();
    for chrom in reader.chroms() {
        chromsizes.insert(chrom.name.clone(), chrom.length);
    }
    chromsizes
}

pub fn chrombounds_from_bam(bamfiles: Vec<&str>) -> HashMap<String, u32> {
    let mut found_chroms: HashMap<String, usize> = HashMap::new();
    for bam in bamfiles.iter() {
        let bam = IndexedReader::from_path(bam).unwrap();
        let chroms: Vec<String> = bam.header().target_names().iter().map(|x| String::from_utf8(x.to_vec()).unwrap()).collect();
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
        if *count == bamfiles.len() {
            validchroms.push(chrom.clone());
        } else {
            println!("Chromosome {} is missing in at least one bam file, and thus ignored!", chrom);
        }
    }
    let bam = IndexedReader::from_path(bamfiles[0]).unwrap();
    let header = bam.header().clone();
    let mut chromsizes = HashMap::new();
    for tid in 0..header.target_count() {
        let chromname = String::from_utf8(header.tid2name(tid).to_vec())
            .expect("Invalid UTF-8 in chromosome name");
        let chromlen = header.target_len(tid)
            .expect("Error retrieving length for chromosome");
        if validchroms.contains(&chromname) {
            chromsizes.insert(chromname, chromlen as u32);
        }
    }
    chromsizes
}

pub fn bwintervals(
    bwfile: &str,
    regions: &Vec<Region>,
    slopregions: &Vec<Vec<Bin>>,
    scale_regions: &Scalingregions
) -> Vec<Vec<f32>> {
    // For a given bw file, a vector of slopregions (Bin enum))
    // return a vector with for every region a vector of f64.

    // Make sure regions and slopregions are of equal length
    assert_eq!(
        regions.len(), slopregions.len(),
        "Regions from bed file and parsed regions (slopped) do not have equal length. Something went wrong during computation."
    );
    
    // Define return vector, set up bw reader.
    let mut bwvals: Vec<Vec<f32>> = Vec::new();
    let bwf = File::open(bwfile).expect("Failed to open bw file.");
    let mut reader = BigWigRead::open(bwf).unwrap();

    // Iterate over regions and slopregions synchronously.
    for (sls, region) in slopregions.iter().zip(regions.iter()) {
        let mut bwval: Vec<f32> = Vec::new();
        // get 'min and max' to query
        // at some point this should become an impl but man am I tired.
        let (min, max) = sls
            .iter()
            .flat_map(|bin| match bin {
                Bin::Conbin(a, b) => vec![*a, *b],
                Bin::Catbin(pairs) => pairs.iter().flat_map(|(x, y)| vec![*x, *y]).collect::<Vec<u32>>(),
            })
            .fold((u32::MAX, u32::MIN), |(min, max), x| (min.min(x), max.max(x)));
        let binvals = reader.get_interval(&region.chrom, min as u32, max as u32).unwrap();
        // since binvals (can) be over binsizes, we expand them to bp and push them to a hashmap
        let mut bwhash: HashMap<u32, f32> = HashMap::new();
        for interval in binvals {
            let interval = interval.unwrap();
            let start = interval.start as u32;
            let end = interval.end as u32;
            let val = interval.value as f32;
            bwhash.extend((start..end).map(|bp| (bp, val)));
        }
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
                        // Get values from the hashmap
                        let vals: Vec<&f32> = (*a..*b)
                            .filter_map(|bp| bwhash.get(&bp))
                            .collect();
                        let val = match scale_regions.avgtype.as_str() {
                            "mean" => mean_float(vals),
                            "median" => median_float(vals),
                            "min" => min_float(vals),
                            "max" => max_float(vals),
                            "std" => std_float(vals),
                            "sum" => sum_float(vals),
                            _ => panic!("Unknown avgtype."),
                        };
                        bwval.push(val);
                    }
                },
                Bin::Catbin(pairs) => {
                    let mut vals: Vec<&f32> = Vec::new();

                    for (start, end) in pairs {
                        if start == end && *end == 0 {
                            if scale_regions.missingdata_as_zero {
                                vals.push(&0.0);
                            } else {
                                vals.push(&std::f32::NAN);
                            }
                        } else {
                            // Get values from the hashmap
                            (*start..*end)
                                .filter_map(|bp| bwhash.get(&bp))
                                .for_each(|v| vals.push(v));
                        }
                    }

                    let val = match scale_regions.avgtype.as_str() {
                        "mean" => mean_float(vals),
                        "median" => median_float(vals),
                        "min" => min_float(vals),
                        "max" => max_float(vals),
                        "std" => std_float(vals),
                        "sum" => sum_float(vals),
                        _ => panic!("Unknown avgtype."),
                    };
                    bwval.push(val);
                }
            }
        }
        assert_eq!(bwval.len(), scale_regions.cols_expected / scale_regions.bwfiles);
        bwvals.push(bwval);
    }
    bwvals
}

pub fn header_matrix(scale_regions: &Scalingregions, regionsizes: &HashMap<String, u32>, sortregions: &str, sortusing: &str) -> String {
    // Create the header for the matrix.
    // This is quite ugly, but this is mainly because we need to accomodate delta bwfiles.
    let mut headstr = String::new();
    headstr.push_str("@{");
    headstr.push_str(
        &format!("\"upstream\":[{}],", (0..scale_regions.bwfiles).map(|_| scale_regions.upstream).collect::<Vec<_>>().into_iter().join(","))
    );
    headstr.push_str(
        &format!("\"downstream\":[{}],", (0..scale_regions.bwfiles).map(|_| scale_regions.downstream).collect::<Vec<_>>().into_iter().join(","))
    );
    headstr.push_str(
        &format!("\"body\":[{}],", (0..scale_regions.bwfiles).map(|_| scale_regions.regionbodylength).collect::<Vec<_>>().into_iter().join(","))
    );
    headstr.push_str(
        &format!("\"bin size\":[{}],", (0..scale_regions.bwfiles).map(|_| scale_regions.binsize).collect::<Vec<_>>().into_iter().join(","))
    );
    // ref point can be empty (for scale_regions, for example).
    // To keep compatibility with deepTools 3 it should be written as null
    let refpointstring = (0..scale_regions.bwfiles)
        .map(|_| scale_regions.referencepoint.clone())
        .collect::<Vec<_>>()
        .into_iter()
        .join("\",\"");

    if refpointstring.is_empty() {
        headstr.push_str(
            &format!("\"ref point\":[null],")
        );
    } else {
        headstr.push_str(
            &format!("\"ref point\":[\"{}\"],", refpointstring)
        );
    }
    
    headstr.push_str(
        &format!("\"verbose\":{},", scale_regions.verbose)
    );
    headstr.push_str(
        &format!("\"bin avg type\":\"{}\",", scale_regions.avgtype)
    );
    headstr.push_str(
        &format!("\"missing data as zero\":{},", scale_regions.missingdata_as_zero)
    );
    // Unimplemented arguments, but they need to be present in header for now anyway.
    let minthresh_str = if scale_regions.minthresh != 0.0 { format!("{}", scale_regions.minthresh) } else { "null".to_string() };
    let maxthresh_str = if scale_regions.maxthresh != 0.0 { format!("{}", scale_regions.maxthresh) } else { "null".to_string() };
    headstr.push_str(
        &format!("\"min threshold\":{},\"max threshold\":{},\"scale\":{},\"skip zeros\":{},\"nan after end\":{},", minthresh_str, maxthresh_str, scale_regions.scale, scale_regions.skipzero, scale_regions.nan_after_end)
    );
    headstr.push_str(
        &format!("\"proc number\":{},", scale_regions.proc_number)
    );
    headstr.push_str(
        &format!("\"sort regions\":\"{}\",\"sort using\":\"{}\",", sortregions, sortusing)
    );
    headstr.push_str(
        &format!("\"startLabel\":\"{}\",\"endLabel\":\"{}\",", scale_regions.startlabel, scale_regions.endlabel)
    );
    headstr.push_str(
        &format!("\"unscaled 5 prime\":[{}],", (0..scale_regions.bwfiles).map(|_| scale_regions.unscaled5prime).collect::<Vec<_>>().into_iter().join(","))
    );
    headstr.push_str(
        &format!("\"unscaled 3 prime\":[{}],", (0..scale_regions.bwfiles).map(|_| scale_regions.unscaled3prime).collect::<Vec<_>>().into_iter().join(","))
    );
    headstr.push_str(
        &format!("\"group_labels\":[\"{}\"],", scale_regions.regionlabels.join("\",\""))
    );
    // Get cumulative sizes of regions
    let mut groupbounds: Vec<u32> = Vec::new();
    groupbounds.push(0);
    let mut cumsum: u32 = 0;
    for regionlabel in scale_regions.regionlabels.iter() {
        cumsum += regionsizes.get(regionlabel).unwrap();
        groupbounds.push(cumsum);
    }
    let groupbounds = format!("{}", groupbounds.iter()
        .map(|&x| x.to_string())
        .collect::<Vec<String>>()
        .join(","));

    headstr.push_str(
        &format!("\"group_boundaries\":[{}],", groupbounds)
    );
    // Sample labels
    headstr.push_str(
        &format!("\"sample_labels\":[\"{}\"],", scale_regions.bwlabels.join("\",\""))
    );
    // Get cumulative sizes for sample boundaries
    let colsize_per_sample = scale_regions.cols_expected / scale_regions.bwfiles;
    let mut samplebounds: Vec<u64> = Vec::new();
    samplebounds.push(0);
    let mut cumsum: u64 = 0;
    for _ in 0..scale_regions.bwfiles {
        cumsum += colsize_per_sample as u64;
        samplebounds.push(cumsum);
    }
    let samplebounds = format!("{}", samplebounds.iter()
        .map(|&x| x.to_string())
        .collect::<Vec<String>>()
        .join(","));

    headstr.push_str(
        &format!("\"sample_boundaries\":[{}]", samplebounds)
    );
    headstr.push_str(
        "}\n"
    );
    headstr
}

pub fn write_matrix(
    header: String,
    mat: Vec<Vec<f32>>,
    ofile: &str,
    regions: Vec<Region>,
    scale_regions: &Scalingregions
) {
    // Write out the matrix to a compressed file.
    let omat = File::create(ofile).unwrap();
    let mut encoder = GzEncoder::new(omat, Compression::default());
    encoder.write_all(header.as_bytes()).unwrap();
    // Final check to make sure our regions and mat iter are of equal length.
    assert_eq!(regions.len(), mat.len());
    for (region, row) in regions.into_iter().zip(mat.into_iter()) {
        // Skipping rules.
        // skip_zeros
        if scale_regions.skipzero && row.iter().all(|&x| x == 0.0) {
            continue;
        }
        // min threshold (Python: coverage.min() <= minThreshold → skip)
        // Python: if any element is NaN, coverage.min() returns NaN, NaN <= threshold is False → row passes
        // Rust: NaN <= threshold is False, so NaN is ignored in .any(). To match Python, skip threshold
        // check entirely if the row contains any NaN values.
        if scale_regions.minthresh != 0.0 && !row.iter().any(|&x| x.is_nan()) && row.iter().any(|&x| x <= scale_regions.minthresh) {
            continue;
        }
        // max threshold (Python: if max(coverage) >= maxThreshold → skip)
        // NaN matching: Python's np.min()/np.max() propagate NaN, so NaN comparisons are False. If row contains NaN, skip filtering.
        if scale_regions.maxthresh != 0.0 && !row.iter().any(|&x| x.is_nan()) && row.iter().any(|&x| x >= scale_regions.maxthresh) {
            continue;
        }
        let mut writerow = format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t",
            region.chrom,                // Chromosome
            region.start.rewrites(),     // Revalue for start (either u32, or Vec<u32>)
            region.end.rewrites(),       // Revalue for end (either u32, or Vec<u32>)
            region.name,                 // String field for name. Duplicates taken care of in read_bedfiles.
            region.score,                // Score field persisted from bedfile
            region.strand,               // Strand field persisted from bedfile
        );
        writerow.push_str(
            &row
                .iter()
                .map(
                    |x| {
                        if x.is_nan() {
                            // as deepTools 3 encoded nan in matrix.
                            "nan".to_string()
                        } else {
                            ((scale_regions.scale * x * 100.0).round() / 100.0).to_string()
                        }
                    }
                )
                .collect::<Vec<String>>().join("\t")
        );
        writerow.push_str("\n");
        encoder.write_all(writerow.as_bytes()).unwrap();
    }
}

pub fn write_matrix_values(file_name: &str, mat: &[Vec<f32>], scale_regions: &Scalingregions, regionsizes: &HashMap<String, u32>) {
    use std::io::Write;
    let mut fh = File::create(file_name).unwrap();
    // Header line 1: group labels with region counts
    let info: Vec<String> = scale_regions.regionlabels.iter()
        .map(|label| format!("{}:{}", label, regionsizes.get(label).unwrap_or(&0)))
        .collect();
    fh.write_all(format!("#{}\n", info.join("\t")).as_bytes()).unwrap();
    // Header line 2: region dimension parameters
    let header2 = format!(
        "#downstream:{}\tupstream:{}\tbody:{}\tbin size:{}\tunscaled 5 prime:{}\tunscaled 3 prime:{}\n",
        scale_regions.downstream, scale_regions.upstream, scale_regions.regionbodylength,
        scale_regions.binsize, scale_regions.unscaled5prime, scale_regions.unscaled3prime,
    );
    fh.write_all(header2.as_bytes()).unwrap();
    // Header line 3: sample labels repeated per column
    let cols_per_sample = scale_regions.cols_expected / scale_regions.bwfiles;
    let sample_info: Vec<String> = scale_regions.bwlabels.iter()
        .flat_map(|label| (0..cols_per_sample).map(move |_| label.clone()))
        .collect();
    fh.write_all(format!("{}\n", sample_info.join("\t")).as_bytes()).unwrap();
    fh.flush().unwrap();

    // Reopen in append mode and write matrix data (matches Python np.savetxt fmt="%.4g")
    let mut fh = std::fs::OpenOptions::new().create(false).append(true).open(file_name).unwrap();
    for row in mat.iter() {
        let line = row.iter()
            .map(|x| if x.is_nan() { "nan".to_string() } else { format_g4(*x) })
            .collect::<Vec<_>>()
            .join("\t");
        writeln!(fh, "{}", line).unwrap();
    }
}

pub fn write_sorted_regions_bed(file_name: &str, regions: &[Region], scale_regions: &Scalingregions, regionsizes: &HashMap<String, u32>) {
    use std::io::Write;
    let mut fh = File::create(file_name).unwrap();
    let header = "#chrom\tstart\tend\tname\tscore\tstrand\tthickStart\tthickEnd\titemRGB\tblockCount\tblockSizes\tblockStarts\tdeepTools_group\n";
    fh.write_all(header.as_bytes()).unwrap();
    // Build group boundaries from regionsizes (cumulative region counts per group label)
    let mut group_boundaries: Vec<u32> = vec![0];
    let mut cumsum: u32 = 0;
    for label in &scale_regions.regionlabels {
        cumsum += *regionsizes.get(label).unwrap_or(&0);
        group_boundaries.push(cumsum);
    }
    for (idx, region) in regions.iter().enumerate() {
        // Find label_idx: last boundary <= idx, matching Python: np.flatnonzero(boundaries <= idx)[-1]
        let label_idx = group_boundaries.iter()
            .take_while(|&&b| b <= idx as u32)
            .count().saturating_sub(1)
            .min(scale_regions.regionlabels.len() - 1);
        let start_first = match &region.start {
            Revalue::U(v) => *v,
            Revalue::V(vs) => *vs.first().unwrap(),
        };
        let end_last = match &region.end {
            Revalue::U(v) => *v,
            Revalue::V(vs) => *vs.last().unwrap(),
        };
        let (block_count, block_sizes, block_starts) = match (&region.start, &region.end) {
            (Revalue::U(_), _) => ("1".to_string(), (end_last - start_first).to_string(), "0".to_string()),
            (Revalue::V(starts), Revalue::V(ends)) => {
                let bc = starts.len().to_string();
                let sz: Vec<String> = starts.iter().zip(ends.iter()).map(|(s, e)| (e - s).to_string()).collect();
                let st: Vec<String> = starts.iter().map(|s| (s - start_first).to_string()).collect();
                (bc, sz.join(","), st.join(","))
            }
            _ => ("1".to_string(), (end_last - start_first).to_string(), "0".to_string()),
        };
        let group_label = &scale_regions.regionlabels[label_idx];
        let line = format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t0\t{}\t{}\t{}\t{}\n",
            region.chrom, start_first, end_last, region.name, region.score,
            region.strand, start_first, end_last, block_count, block_sizes,
            block_starts, group_label,
        );
        fh.write_all(line.as_bytes()).unwrap();
    }
}