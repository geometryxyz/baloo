use ark_ff::{FftField, BigInteger, Zero, batch_inversion};
use ark_poly::{univariate::{DensePolynomial, DenseOrSparsePolynomial}, UVPolynomial, Polynomial};

/*
    Fast evaluation algorithm based on reminders tree
    We require FftField to have multiplication in O(nlog(n))
 */
pub struct FastEval<F: FftField> {
    pub(crate) subproduct_tree: Vec<Vec<DensePolynomial<F>>>
}

impl<F: FftField> FastEval<F> {
    pub fn build_tree(roots: &[F]) -> Self {
        let n = roots.len(); 

        if n == 0 {
            panic!("FastEval: Empty roots");
        }

        if n & n - 1 != 0 {
            panic!("FastEval: Roots are not power of 2");
        }

        let k: usize = n.trailing_zeros().try_into().unwrap();
        let mut subproduct_tree = vec![vec![]; k + 1]; 

        subproduct_tree[0] = Vec::with_capacity(n);
        for &root in roots {
            subproduct_tree[0].push(DensePolynomial::from_coefficients_slice(&[-root, F::one()]));
        }

        let mut nodes_on_layer = n;
        for i in 1..=k {
            nodes_on_layer /= 2;
            subproduct_tree[i] = Vec::with_capacity(nodes_on_layer);
            for j in 0..nodes_on_layer {
                let lhs_node = subproduct_tree[i-1][2*j].clone();
                let rhs_node = subproduct_tree[i-1][2*j+1].clone();

                subproduct_tree[i].push(&lhs_node * &rhs_node);
            }
        }

        Self { 
            subproduct_tree 
        }

    }

    pub fn eval_over_domain(&self, f: &DensePolynomial<F>) -> Vec<F> {
        let n = self.subproduct_tree[0].len(); 
        let k = self.subproduct_tree.len() - 1; 
        self.multipoint_eval(n, (k, 0), f)
    }

    pub fn evaluate_lagrange_polys(&self, vanishing: &DensePolynomial<F>, point: &F) -> Vec<F> {
        // TODO: store vanishing and vanishing_derivative in fast_eval
        let mut vanishing_derivative = DensePolynomial::default();
        let mut monomials_evals = Vec::with_capacity(self.subproduct_tree[0].len());
        let zi_d_or_s = DenseOrSparsePolynomial::from(vanishing);
        for root_monomial in &self.subproduct_tree[0] {
            monomials_evals.push(root_monomial.evaluate(point));
            let (q, r) = zi_d_or_s.divide_with_q_and_r(&(root_monomial.into())).unwrap(); 
            assert!(r.is_zero());
            vanishing_derivative += &q;
        }

        let zi_eval = vanishing.evaluate(&point); 

        let ri_inverses = self.eval_over_domain(&vanishing_derivative);
        let mut ri_monomial_evals = ri_inverses.iter().zip(monomials_evals.iter()).map(|(&ri, &beta_minus_omega_i)| ri * beta_minus_omega_i).collect::<Vec<_>>();
        batch_inversion(&mut ri_monomial_evals);

        ri_monomial_evals.iter().map(|&l_part_i| l_part_i * zi_eval).collect()
    }

    fn multipoint_eval(&self, n: usize, root: (usize, usize), f: &DensePolynomial<F>) -> Vec<F> {
        assert!(f.degree() < n);

        if n == 1 {
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
    use ark_ff::{UniformRand, One};
    use ark_poly::{GeneralEvaluationDomain, EvaluationDomain, univariate::DensePolynomial, UVPolynomial, Polynomial};
    use ark_std::test_rng;

    use crate::utils::construct_lagrange_basis;

    use super::FastEval; 
    #[test]
    fn test_fft() {
        let mut rng = test_rng();
        let n: usize = 1024; 
        let domain = GeneralEvaluationDomain::<Fr>::new(n).unwrap();

        let roots: Vec<Fr> = domain.elements().collect();
        let fast_eval = FastEval::build_tree(&roots);

        let f = DensePolynomial::<Fr>::rand(n - 1, &mut rng);
        let f_fft_evals = domain.fft(&f); 

        let f_evals = fast_eval.eval_over_domain(&f);
        assert_eq!(f_evals, f_fft_evals);
    }

    #[test]
    fn test_non_multiplicative_subgroup() {
        let mut rng = test_rng();
        let n: usize = 8; 

        let roots: Vec<_> = (0..n).map(|_| Fr::rand(&mut rng)).collect();
        let fast_eval = FastEval::build_tree(&roots);
        let basis = construct_lagrange_basis(&roots);


        let mut vanishing = DensePolynomial::from_coefficients_slice(&[Fr::one()]);
        for ui in roots {
            vanishing = &vanishing * &DensePolynomial::from_coefficients_slice(&[-ui, Fr::one()]);
        }

        let beta = Fr::rand(&mut rng); 
        let lagrange_evals = fast_eval.evaluate_lagrange_polys(&vanishing, &beta);

        let lagrange_evals_slow = basis.iter().map(|li| li.evaluate(&beta)).collect::<Vec<_>>();
        assert_eq!(lagrange_evals, lagrange_evals_slow);
    }
}