//! Scale factors for coverage normalization.

use std::cmp;
use std::str::FromStr;

/// Errors from the scale factor functions.
#[derive(Debug, Clone, PartialEq, thiserror::Error)]
pub enum NormalizationError {
    /// Unknown method name.
    #[error("unknown {kind} '{value}', expected one of: {expected}")]
    UnknownMethod {
        /// What was parsed.
        kind: &'static str,
        /// Rejected input.
        value: String,
        /// Accepted values.
        expected: &'static str,
    },
    /// RPGC without an effective genome size.
    #[error("RPGC normalization requires a non-zero effective genome size")]
    MissingEffectiveGenomeSize,
}

/// Normalization method (`--normalizeUsing`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Normalization {
    /// Reads per kilobase per million.
    Rpkm,
    /// Counts per million.
    Cpm,
    /// Bins per million.
    Bpm,
    /// Reads per genomic content (1x).
    Rpgc,
    /// No normalization.
    None,
}

impl FromStr for Normalization {
    type Err = NormalizationError;

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

/// Scaling between two samples (`--scaleFactorsMethod`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ScaleFactorsMethod {
    /// Scale the larger sample down to the smaller.
    ReadCount,
    /// Normalize each sample on its own.
    None,
}

impl FromStr for ScaleFactorsMethod {
    type Err = NormalizationError;

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

/// Scale factor for one sample; `preset` (`--scaleFactor`) overrides `method`. Zero mapped reads gives `inf`.
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

/// Read statistics of one sample.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SampleStats {
    /// Mapped reads.
    pub mapped: u32,
    /// Median read length.
    pub readlen: f32,
}

/// Scale factors for a pair of samples. An empty sample gives `NaN` with `ReadCount`.
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
