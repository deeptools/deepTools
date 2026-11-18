//! The `@{...}` JSON header line of the computeMatrix matrix file.

use crate::covcalc::Scalingregions;
use itertools::Itertools;
use std::collections::HashMap;

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
    let threshold_json = |threshold: Option<f32>| {
        threshold.map_or_else(|| "null".to_string(), |value| value.to_string())
    };
    let minthresh_str = threshold_json(scale_regions.minthresh);
    let maxthresh_str = threshold_json(scale_regions.maxthresh);
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
