use crate::fixedfmt::push_f32_fixed6;

fn reference(v: f32) -> String {
    format!("{:.6}", v)
}

fn mine(v: f32) -> String {
    let mut o = Vec::new();
    push_f32_fixed6(&mut o, v);
    String::from_utf8(o).unwrap()
}

/// Small xorshift PRNG so the tests need no extra dependency.
pub struct Rng(pub u64);
impl Rng {
    pub fn next(&mut self) -> u64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        self.0
    }
}

#[test]
fn test_fixedfmt_special_values() {
    let vals = [
        0.0f32,
        -0.0,
        1.0,
        -1.0,
        0.5,
        1e-7,
        -1e-7,
        4.9999997e-7,
        5e-7,
        5.0000006e-7,
        0.0078125, // exact tie at the 7th decimal: round half to even
        0.0234375,
        -0.0078125,
        1.0000005,
        123456.79,
        16777216.0,
        3.4028235e38,
        -3.4028235e38,
        f32::MIN_POSITIVE,
        1e-45,
        f32::INFINITY,
        f32::NEG_INFINITY,
        f32::NAN,
        -f32::NAN,
    ];
    for v in vals {
        assert_eq!(mine(v), reference(v), "value bits {:#010x}", v.to_bits());
    }
}

#[test]
fn test_fixedfmt_random_bits() {
    let mut rng = Rng(0x9E3779B97F4A7C15);
    for _ in 0..2_000_000 {
        let v = f32::from_bits(rng.next() as u32);
        assert_eq!(mine(v), reference(v), "value bits {:#010x}", v.to_bits());
    }
    // Typical coverage values: small magnitudes near rounding boundaries.
    for i in 0..2_000_000u32 {
        let v = (i as f32) * 1e-7 - 0.1;
        assert_eq!(mine(v), reference(v), "value bits {:#010x}", v.to_bits());
    }
}

/// Exhaustive check over all 2^32 bit patterns (~25 CPU-minutes).
/// Run with `cargo test --release -- --ignored test_fixedfmt_exhaustive`.
#[test]
#[ignore]
fn test_fixedfmt_exhaustive() {
    use rayon::prelude::*;
    let bad: u64 = (0u64..(1 << 16))
        .into_par_iter()
        .map(|hi| {
            let mut bad = 0;
            for lo in 0u64..(1 << 16) {
                let v = f32::from_bits(((hi << 16) | lo) as u32);
                if mine(v) != reference(v) {
                    bad += 1;
                }
            }
            bad
        })
        .sum();
    assert_eq!(bad, 0);
}
