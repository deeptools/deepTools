//! Scale factor calculation for coverage normalization.
//!
//! All functions are pure: no I/O, no printing, no panics on bad input. Invalid
//! input is reported through [`NormalizationError`].

use std::cmp;
use std::str::FromStr;

/// Errors returned by the scale factor functions.
#[derive(Debug, Clone, PartialEq, thiserror::Error)]
pub enum NormalizationError {
    /// A method name that is not one of the accepted values.
    #[error("unknown {kind} '{value}', expected one of: {expected}")]
    UnknownMethod {
        /// What was being parsed (e.g. `normalization method`).
        kind: &'static str,
        /// The rejected input.
        value: String,
        /// Comma separated list of accepted values.
        expected: &'static str,
    },
    /// RPGC divides by the effective genome size, which must be non-zero.
    #[error("RPGC normalization requires a non-zero effective genome size")]
    MissingEffectiveGenomeSize,
}

/// Coverage normalization method (`--normalizeUsing`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Normalization {
    /// Reads per kilobase per million mapped reads.
    Rpkm,
    /// Counts per million mapped reads.
    Cpm,
    /// Bins per million mapped reads.
    Bpm,
    /// Reads per genomic content (1x coverage).
    Rpgc,
    /// No normalization; scale factor is 1.
    None,
}

impl FromStr for Normalization {
    type Err = NormalizationError;

    /// Accepts `RPKM`, `CPM`, `BPM`, `RPGC` and `None` (case-sensitive).
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "RPKM" => Ok(Self::Rpkm),
            "CPM" => Ok(Self::Cpm),
            "BPM" => Ok(Self::Bpm),
            "RPGC" => Ok(Self::Rpgc),
            "None" => Ok(Self::None),
            _ => Err(NormalizationError::UnknownMethod {
                kind: "normalization method",
                value: s.to_string(),
                expected: "RPKM, CPM, BPM, RPGC, None",
            }),
        }
    }
}

/// How two samples are scaled against each other (`--scaleFactorsMethod`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ScaleFactorsMethod {
    /// Scale the larger sample down to the read count of the smaller one.
    ReadCount,
    /// Scale each sample independently with a [`Normalization`] method.
    None,
}

impl FromStr for ScaleFactorsMethod {
    type Err = NormalizationError;

    /// Accepts `readCount` and `None` (case-sensitive).
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "readCount" => Ok(Self::ReadCount),
            "None" => Ok(Self::None),
            _ => Err(NormalizationError::UnknownMethod {
                kind: "scale factors method",
                value: s.to_string(),
                expected: "readCount, None",
            }),
        }
    }
}

/// Scale factor for one sample.
///
/// * `mapped`: number of mapped reads.
/// * `binsize`: bin size in bp.
/// * `effective_genome_size`: only used (and required non-zero) for RPGC.
/// * `readlen`: median read length in bp, only used for RPGC.
/// * `preset`: user supplied scale factor. When `Some`, it is returned as is
///   and `method` is ignored.
pub fn scale_factor(
    method: Normalization,
    mapped: u32,
    binsize: u32,
    effective_genome_size: u64,
    readlen: f32,
    preset: Option<f32>,
) -> Result<f32, NormalizationError> {
    if let Some(sf) = preset {
        return Ok(sf);
    }
    let mmr = mapped as f32 / 1e6;
    let bs_kb = binsize as f32 / 1000.0;
    Ok(match method {
        // reads per tile / total reads (millions) / tile length (kb)
        Normalization::Rpkm => 1.0 / (mmr * bs_kb),
        // reads per tile / total reads (millions)
        Normalization::Cpm => 1.0 / mmr,
        // bins per million mapped reads
        Normalization::Bpm => {
            let tmp = (mapped as f32 / bs_kb) / 1e6;
            1.0 / (tmp * bs_kb)
        }
        // 1x coverage: mapped reads * read length / effective genome size
        Normalization::Rpgc => {
            if effective_genome_size == 0 {
                return Err(NormalizationError::MissingEffectiveGenomeSize);
            }
            1.0 / ((mapped as f32 * readlen) / effective_genome_size as f32)
        }
        Normalization::None => 1.0,
    })
}

/// Per-sample read statistics used by [`scale_factor_bamcompare`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SampleStats {
    /// Number of mapped reads.
    pub mapped: u32,
    /// Median read length in bp.
    pub readlen: f32,
}

/// Scale factors `(first, second)` for a pair of samples.
///
/// With [`ScaleFactorsMethod::ReadCount`] the sample with more reads is scaled
/// down to the read count of the other. With [`ScaleFactorsMethod::None`] each
/// sample is scaled independently using `norm`.
pub fn scale_factor_bamcompare(
    method: ScaleFactorsMethod,
    samples: [SampleStats; 2],
    binsize: u32,
    effective_genome_size: u64,
    norm: Normalization,
) -> Result<(f32, f32), NormalizationError> {
    match method {
        ScaleFactorsMethod::ReadCount => {
            let min = cmp::min(samples[0].mapped, samples[1].mapped);
            Ok((
                min as f32 / samples[0].mapped as f32,
                min as f32 / samples[1].mapped as f32,
            ))
        }
        ScaleFactorsMethod::None => {
            let sf = |s: SampleStats| {
                scale_factor(norm, s.mapped, binsize, effective_genome_size, s.readlen, None)
            };
            Ok((sf(samples[0])?, sf(samples[1])?))
        }
    }
}
