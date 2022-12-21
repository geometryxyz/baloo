use ark_ff::{FftField, Zero, batch_inversion};
use ark_poly::{univariate::{DensePolynomial, DenseOrSparsePolynomial}, UVPolynomial, Polynomial};

/*
    - Fast evaluation and interpolation algorithm based on reminders tree
    - FftField is required in order to have multiplication in O(nlog(n))
    - Additional elements are computed for maximal baloo optimization
 */
pub struct FastEval<F: FftField> {
    pub(crate) vanishing: DensePolynomial<F>, 
    pub(crate) vanishing_derivative: DensePolynomial<F>,
    pub(crate) ri_evals: Option<Vec<F>>,
    pub(crate) subproduct_tree: Vec<Vec<DensePolynomial<F>>>
}

impl<F: FftField> FastEval<F> {
    pub fn prepare(roots: &[F]) -> Self {
        let n = roots.len(); 

        if n == 0 {
            panic!("FastEval: Empty roots");
        }

        if n & n - 1 != 0 {
            panic!("FastEval: Roots are not power of 2");
        }

        let k: usize = n.trailing_zeros().try_into().unwrap();
        let mut subproduct_tree = vec![vec![]; k + 1]; 
        let mut vanishing = DensePolynomial::from_coefficients_slice(&[F::one()]);

        subproduct_tree[0] = Vec::with_capacity(n);
        for &root in roots {
            let root_monomial = DensePolynomial::from_coefficients_slice(&[-root, F::one()]);
            vanishing = &vanishing * &root_monomial;
            subproduct_tree[0].push(root_monomial);
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

        let mut vanishing_derivative = DensePolynomial::default();
        let vanishing_into = DenseOrSparsePolynomial::from(&vanishing);
        for root_monomial in &subproduct_tree[0] {
            let (q, r) = vanishing_into.divide_with_q_and_r(&(root_monomial.into())).unwrap(); 
            assert!(r.is_zero());
            vanishing_derivative += &q;
        }

        Self { 
            vanishing, 
            vanishing_derivative,
            ri_evals: None,
            subproduct_tree 
        }

    }

    pub fn compute_ri_evals(&mut self) {
        if self.ri_evals.is_some() {}

        let mut ri_evals = self.eval_over_domain(&self.vanishing_derivative);
        batch_inversion(&mut ri_evals);

        self.ri_evals = Some(ri_evals);
    }

    pub fn eval_over_domain(&self, f: &DensePolynomial<F>) -> Vec<F> {
        let n = self.subproduct_tree[0].len(); 
        let k = self.subproduct_tree.len() - 1; 
        self.multipoint_eval(n, (k, 0), f)
    }

    pub fn interpolate(&self, evals: &Vec<F>) -> DensePolynomial<F> {
        if self.ri_evals.is_none() {
            panic!("Compute ri evals")
        }

        let ri_evals = self.ri_evals.as_ref().unwrap();
        let k = self.subproduct_tree.len() - 1; 
        let evals = evals.iter().zip(ri_evals.iter()).map(|(&vi, &ri)| vi * ri).collect::<Vec<_>>();
        self.interpolate_from_subtree((0, evals.len() - 1), (k, 0), &evals)
    }

    pub fn reconstruct_vanishing_derivative(&self) -> DensePolynomial<F> {
        let n = self.subproduct_tree[0].len(); 
        let k = self.subproduct_tree.len() - 1; 
        let evals = vec![F::one(); n];
        self.interpolate_from_subtree((0, evals.len() - 1), (k, 0), &evals)
    }

    pub fn interpolate_from_subtree(&self, index_bounds: (usize, usize), root: (usize, usize), evals: &Vec<F>) -> DensePolynomial<F> {
        if index_bounds.1 - index_bounds.0 == 0 {
            return DensePolynomial::from_coefficients_slice(&[evals[index_bounds.0]])
        }

        let len = (index_bounds.1 - index_bounds.0) / 2;
        let lhs_bounds = (index_bounds.0, index_bounds.0 + len);
        let rhs_bounds = (lhs_bounds.1 + 1, index_bounds.1);    

        let r0 = self.interpolate_from_subtree(lhs_bounds, (root.0 - 1, 2*root.1), evals);
        let r1 = self.interpolate_from_subtree(rhs_bounds, (root.0 - 1, 2*root.1 + 1), evals);

        let lhs = &self.subproduct_tree[root.0 - 1][2*root.1];
        let rhs = &self.subproduct_tree[root.0 - 1][2*root.1 + 1];
        &r0 * rhs  + &r1 * lhs
    }

    pub fn evaluate_lagrange_polys(&self, point: &F) -> Vec<F> {
        let mut monomials_evals = Vec::with_capacity(self.subproduct_tree[0].len());
        for root_monomial in &self.subproduct_tree[0] {
            monomials_evals.push(root_monomial.evaluate(point));
        }

        let zi_eval = self.vanishing.evaluate(&point); 


        // If ri_evals are already computed, we use them, if not, we combine it's inverses with monomial evals and then batch invert it
        let li_part_evals = if let Some(ri_evals) = self.ri_evals.as_ref() {
            batch_inversion(&mut monomials_evals);
            ri_evals.iter().zip(monomials_evals.iter()).map(|(&ri, &beta_minus_omega_i)| ri * beta_minus_omega_i).collect::<Vec<_>>()
        } else {
            let ri_inverses = self.eval_over_domain(&self.vanishing_derivative);
            let mut ri_monomial_evals = ri_inverses.iter().zip(monomials_evals.iter()).map(|(&ri, &beta_minus_omega_i)| ri * beta_minus_omega_i).collect::<Vec<_>>();
            batch_inversion(&mut ri_monomial_evals);
            ri_monomial_evals
        };

        li_part_evals.iter().map(|&l_part_i| l_part_i * zi_eval).collect()
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
        let fast_eval = FastEval::prepare(&roots);

        let f = DensePolynomial::<Fr>::rand(n - 1, &mut rng);
        let f_fft_evals = domain.fft(&f); 

        let f_evals = fast_eval.eval_over_domain(&f);
        assert_eq!(f_evals, f_fft_evals);
    }

    #[test]
    fn test_non_multiplicative_subgroup() {
        let mut rng = test_rng();
        let n: usize = 128; 

        let roots: Vec<_> = (0..n).map(|_| Fr::rand(&mut rng)).collect();
        let fast_eval = FastEval::prepare(&roots);
        let basis = construct_lagrange_basis(&roots);

        let beta = Fr::rand(&mut rng); 
        let lagrange_evals = fast_eval.evaluate_lagrange_polys(&beta);

        let lagrange_evals_slow = basis.iter().map(|li| li.evaluate(&beta)).collect::<Vec<_>>();
        assert_eq!(lagrange_evals, lagrange_evals_slow);
    }

    #[test]
    fn test_interpolation() {
        let mut rng = test_rng();
        let n: usize = 128; 

        let roots: Vec<_> = (0..n).map(|_| Fr::rand(&mut rng)).collect();
        let mut fast_eval = FastEval::prepare(&roots);
        fast_eval.compute_ri_evals();

        let basis = construct_lagrange_basis(&roots);
        let evals: Vec<_> = (0..n).map(|_| Fr::rand(&mut rng)).collect();

        let poly_fast = fast_eval.interpolate(&evals);

        let mut poly_slow = DensePolynomial::default(); 
        for (li, &vi) in basis.iter().zip(evals.iter()) {
            poly_slow += (vi, li);
        }

        assert_eq!(poly_fast, poly_slow);
    }

    #[test]
    fn test_derivative() {
        let mut rng = test_rng();
        let n: usize = 128; 

        let roots: Vec<_> = (0..n).map(|_| Fr::rand(&mut rng)).collect();
        let fast_eval = FastEval::prepare(&roots);

        let vh_derivative = fast_eval.reconstruct_vanishing_derivative();
        assert_eq!(vh_derivative, fast_eval.vanishing_derivative);
    }

    #[test]
    fn test_monomial_multiplication_in_smaller_fft() {
        let mut rng = test_rng(); 
        
        let u0 = Fr::rand(&mut rng); 
        let u1 = Fr::rand(&mut rng);

        let m0 = DensePolynomial::from_coefficients_slice(&[-u0, Fr::one()]);
        let m1 = DensePolynomial::from_coefficients_slice(&[-u1, Fr::one()]);
        let prod_poly = &m0 * &m1;
        // println!("prod poly: {:?}", prod_poly);
        for coeff in prod_poly.coeffs() {
            println!("ci: {}", coeff);
        }

        println!("===========================");

        let n = 2; 
        let domain_2 = GeneralEvaluationDomain::<Fr>::new(n).unwrap(); 
        let domain_2n = GeneralEvaluationDomain::<Fr>::new(2*n).unwrap(); 

        let m0_evals = domain_2n.fft(&m0); 
        let m1_evals = domain_2n.fft(&m1);

        let prod_evals: Vec<Fr> = m0_evals.iter().zip(m1_evals.iter()).map(|(&m0, &m1)| m0 * m1).collect();
        let prod = DensePolynomial::from_coefficients_slice(&domain_2n.ifft(&prod_evals));   
        assert_eq!(prod, prod_poly);
        
        let m0_evals = domain_2.fft(&m0); 
        let m1_evals = domain_2.fft(&m1);
        let prod_evals: Vec<Fr> = m0_evals.iter().zip(m1_evals.iter()).map(|(&m0, &m1)| m0 * m1).collect();
        let mut prod = DensePolynomial::from_coefficients_slice(&domain_2.ifft(&prod_evals));   
        for coeff in prod.coeffs() {
            println!("ci: {}", coeff);
        }

        prod[0] -= Fr::one(); 
        prod.coeffs.push(Fr::one());
        assert_eq!(prod, prod_poly);

    }
}