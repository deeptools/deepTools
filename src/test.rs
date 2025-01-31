use crate::calc::median;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_median() {
        let v: Vec<u32> = vec![1,2,3,4,5];
        assert_eq!(median(v), 3.0);
    }
}