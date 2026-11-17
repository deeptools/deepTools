//! Allocation-free fixed-point formatting of `f32` values.
//!
//! `push_f32_fixed6(out, v)` appends exactly the bytes that `format!("{:.6}", v)`
//! produces, for every `f32` bit pattern (including `NaN`, `inf`, `-inf`, `-0.0`).
//! Equivalence was checked exhaustively over all 2^32 bit patterns with
//! rustc 1.98 (see `tests::test_fixedfmt`; the exhaustive test is `#[ignore]`d).
//!
//! An `f32` is `m * 2^e` with an integer mantissa `m < 2^24`, so its exact decimal
//! value scaled by 10^6 is `m * 10^6 / 2^-e`; this fits in 64 bits whenever
//! `e < 0`. The quotient is rounded half-to-even on the exact value, which is what
//! core's `{:.N}` formatting does (it formats the exact binary value, not a
//! shortest representation such as `ryu`'s). Non-negative exponents are integers
//! and are printed via `u128`.

const SCALE: u64 = 1_000_000;

#[inline]
fn push_u64(out: &mut Vec<u8>, mut n: u64) {
    let mut buf = [0u8; 20];
    let mut i = buf.len();
    loop {
        i -= 1;
        buf[i] = b'0' + (n % 10) as u8;
        n /= 10;
        if n == 0 {
            break;
        }
    }
    out.extend_from_slice(&buf[i..]);
}

#[inline]
fn push_u128(out: &mut Vec<u8>, mut n: u128) {
    if n <= u64::MAX as u128 {
        return push_u64(out, n as u64);
    }
    let mut buf = [0u8; 40];
    let mut i = buf.len();
    while n > 0 {
        i -= 1;
        buf[i] = b'0' + (n % 10) as u8;
        n /= 10;
    }
    out.extend_from_slice(&buf[i..]);
}

/// Append the 6 fractional digits of `frac` (< 10^6), zero padded.
#[inline]
fn push_frac6(out: &mut Vec<u8>, frac: u64) {
    let f = frac as u32;
    out.extend_from_slice(&[
        b'.',
        b'0' + (f / 100_000) as u8,
        b'0' + (f / 10_000 % 10) as u8,
        b'0' + (f / 1_000 % 10) as u8,
        b'0' + (f / 100 % 10) as u8,
        b'0' + (f / 10 % 10) as u8,
        b'0' + (f % 10) as u8,
    ]);
}

/// Append `format!("{:.6}", v)` to `out`.
#[inline]
pub fn push_f32_fixed6(out: &mut Vec<u8>, v: f32) {
    let bits = v.to_bits();
    let neg = (bits >> 31) != 0;
    let exp_bits = ((bits >> 23) & 0xff) as i32;
    let frac_bits = (bits & 0x7f_ffff) as u64;
    if exp_bits == 0xff {
        if frac_bits != 0 {
            out.extend_from_slice(b"NaN");
        } else if neg {
            out.extend_from_slice(b"-inf");
        } else {
            out.extend_from_slice(b"inf");
        }
        return;
    }
    if neg {
        // core prints the sign of negative values, including -0.0 and values
        // that round to zero (e.g. "-0.000000").
        out.push(b'-');
    }
    let (m, e) = if exp_bits == 0 {
        (frac_bits, -149i32)
    } else {
        (frac_bits | 0x80_0000, exp_bits - 150)
    };
    if e >= 0 {
        // Integer value m * 2^e < 2^128.
        push_u128(out, (m as u128) << e);
        out.extend_from_slice(b".000000");
        return;
    }
    let shift = (-e) as u32; // 1..=149
    let num = m * SCALE; // < 2^44
    let q = if shift >= 63 {
        // num < 2^44 < 2^(shift-1): always rounds down to 0.
        0
    } else {
        let mut q = num >> shift;
        let rem = num & ((1u64 << shift) - 1);
        let half = 1u64 << (shift - 1);
        if rem > half || (rem == half && (q & 1) == 1) {
            q += 1;
        }
        q
    };
    push_u64(out, q / SCALE);
    push_frac6(out, q % SCALE);
}
