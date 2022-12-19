use ark_ff::FftField;
use ark_poly::{univariate::{DensePolynomial, DenseOrSparsePolynomial}, UVPolynomial, domain, Polynomial};

/*
    Fast evaluation algorithm based on reminders tree
 */
pub struct FastEval<F: FftField> {
    pub(crate) domain_name: String,
    pub(crate) subproduct_tree: Vec<Vec<DensePolynomial<F>>>
}

impl<F: FftField> FastEval<F> {
    pub fn build_tree(domain_name: &str, roots: &[F]) -> Self {
        let n = roots.len(); 

        if n == 0 {
            panic!("FastEval: Empty roots");
        }

        if n & n - 1 != 0 {
            panic!("FastEval: Roots are not power of 2");
        }

        let k: usize = n.trailing_zeros().try_into().unwrap();
        let mut subproduct_tree = vec![vec![]; k + 1]; 

        for &root in roots {
            subproduct_tree[0].push(DensePolynomial::from_coefficients_slice(&[-root, F::one()]));
        }

        let mut nodes_on_layer = n;
        for i in 1..=k {
            nodes_on_layer /= 2;
            for j in 0..nodes_on_layer {
                let lhs_node = subproduct_tree[i-1][2*j].clone();
                let rhs_node = subproduct_tree[i-1][2*j+1].clone();

                subproduct_tree[i].push(&lhs_node * &rhs_node);
            }
        }

        Self { 
            domain_name: domain_name.into(),
            subproduct_tree 
        }

    }

    pub fn multipoint_eval(&self, n: usize, root: (usize, usize), f: &DensePolynomial<F>) -> Vec<F> {
        assert!(f.degree() < n);

        if n == 1 {
            // it's just constant since f.degree() < 1
            return vec![f.coeffs[0]] 
        }

        let f_ds = DenseOrSparsePolynomial::from(f);
        let lhs_divisor = DenseOrSparsePolynomial::from(&self.subproduct_tree[root.0 - 1][2*root.1]);
        let rhs_divisor = DenseOrSparsePolynomial::from(&self.subproduct_tree[root.0 - 1][2*root.1 + 1]);

        let (_, r0) = f_ds.divide_with_q_and_r(&lhs_divisor).unwrap();
        let (_, r1) = f_ds.divide_with_q_and_r(&rhs_divisor).unwrap();

        let mut lhs_evals = self.multipoint_eval(n/2, (root.0 - 1, 2*root.1), &r0);
        let rhs_evals = self.multipoint_eval(n/2, (root.0 - 1, 2*root.1 + 1), &r1);

        lhs_evals.extend_from_slice(&rhs_evals);
        lhs_evals

    }
}

#[cfg(test)]
mod fast_eval_tests {
    use ark_bn254::Fr;
    use ark_poly::{GeneralEvaluationDomain, EvaluationDomain, univariate::DensePolynomial, UVPolynomial};
    use ark_std::test_rng;

    use super::FastEval; 
    #[test]
    fn test_fft() {
        let mut rng = test_rng();
        let n: usize = 1024; 
        let k: usize = n.trailing_zeros().try_into().unwrap();
        let domain = GeneralEvaluationDomain::<Fr>::new(n).unwrap();

        let roots: Vec<Fr> = domain.elements().collect();
        let fast_eval = FastEval::build_tree("omegas", &roots);

        let f = DensePolynomial::<Fr>::rand(n - 1, &mut rng);
        let f_fft_evals = domain.fft(&f); 

        let f_evals = fast_eval.multipoint_eval(n, (k, 0), &f);
        assert_eq!(f_evals, f_fft_evals);
    }
}