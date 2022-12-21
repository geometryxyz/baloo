use std::marker::PhantomData;

use ark_ff::PrimeField;
use ark_poly::{
    univariate::{DenseOrSparsePolynomial, DensePolynomial},
    UVPolynomial,
};

pub struct GeneralizedInnerProduct<F: PrimeField> {
    _f: PhantomData<F>,
}

/*
   - zI doesn't have to be vanishing polynomial of multiplicative subgroup.
   - Caller of this submodule has to make sure that one of {a, b} is in normalized encoding, namely for lagrange basis [Li], normalized form of p(X):
            ----
   p(X) =   \      pi * Li(X)
            /           -----
            ----        Li(0)
*/
impl<F: PrimeField> GeneralizedInnerProduct<F> {
    /*
       a(X)b(X) - sigma = XR(X) + zI(X)Q(X)
    */
    pub fn prove_with_xrx(
        a: &DensePolynomial<F>,
        b: &DensePolynomial<F>,
        sigma: F,
        zi: &DensePolynomial<F>,
    ) -> (DensePolynomial<F>, DensePolynomial<F>) {
        let mut lhs = a * b;
        lhs[0] -= sigma;
        let (q, r) = DenseOrSparsePolynomial::from(lhs)
            .divide_with_q_and_r(&DenseOrSparsePolynomial::from(zi))
            .unwrap();
        (q, DensePolynomial::from_coefficients_slice(&r.coeffs[1..]))
    }

    /*
       a(X)b(X) - sigma = R(X) + zI(X)Q(X)
       Where additionally it has to be proven that R(0) = 0. Look at https://eprint.iacr.org/2022/1565.pdf, page 17, part: (iii)
    */
    pub fn prove_with_rx(
        a: &DensePolynomial<F>,
        b: &DensePolynomial<F>,
        sigma: F,
        zi: &DensePolynomial<F>,
    ) -> (DensePolynomial<F>, DensePolynomial<F>) {
        let mut lhs = a * b;
        lhs[0] -= sigma;
        let (q, r) = DenseOrSparsePolynomial::from(lhs)
            .divide_with_q_and_r(&DenseOrSparsePolynomial::from(zi))
            .unwrap();
        assert_eq!(r[0], F::zero());
        (q, r)
    }
}

#[cfg(test)]
mod generalized_inner_product_tests {
    use ark_bn254::Fr;
    use ark_ff::{batch_inversion, FftField, UniformRand};
    use ark_poly::{univariate::DensePolynomial, Polynomial, UVPolynomial};
    use ark_std::{
        rand::{rngs::StdRng, RngCore},
        test_rng,
    };

    use super::GeneralizedInnerProduct;


    /// given x coords construct Li polynomials
    pub fn construct_lagrange_basis<F: FftField>(evaluation_domain: &[F]) -> Vec<DensePolynomial<F>> {
        let mut bases = Vec::with_capacity(evaluation_domain.len());
        for i in 0..evaluation_domain.len() {
            let mut l_i = DensePolynomial::from_coefficients_slice(&[F::one()]);
            let x_i = evaluation_domain[i];
            for j in 0..evaluation_domain.len() {
                if j != i {
                    let xi_minus_xj_inv = (x_i - evaluation_domain[j]).inverse().unwrap();
                    l_i = &l_i
                        * &DensePolynomial::from_coefficients_slice(&[
                            -evaluation_domain[j] * xi_minus_xj_inv,
                            xi_minus_xj_inv,
                        ]);
                }
            }

            bases.push(l_i);
        }

        bases
    }

    pub fn setup<F: FftField, R: RngCore>(
        n: usize,
        rng: &mut R,
    ) -> (
        DensePolynomial<F>,
        DensePolynomial<F>,
        DensePolynomial<F>,
        F,
    ) {
        let zero = F::zero();

        let roots = (0..n).map(|_| F::rand(rng)).collect::<Vec<F>>();
        let ls = construct_lagrange_basis(&roots);
        let mut ls_zero_inverses = ls.iter().map(|li| li.evaluate(&zero)).collect::<Vec<F>>();
        batch_inversion(&mut ls_zero_inverses);

        let mut zI = DensePolynomial::from_coefficients_slice(&[F::one()]);
        for &omega_i in roots.iter() {
            let current_root = DensePolynomial::from_coefficients_slice(&[omega_i, -F::one()]);
            zI = &zI * &current_root;
        }

        let a_evals = (0..n).map(|_| F::rand(rng)).collect::<Vec<F>>();
        let b_evals = (0..n).map(|_| F::rand(rng)).collect::<Vec<F>>();

        let inner_product: F = a_evals
            .iter()
            .zip(b_evals.iter())
            .map(|(&a, &b)| a * b)
            .sum();

        let mut a = DensePolynomial::default();
        for ((&a_eval, &li_0_inv), li) in a_evals.iter().zip(ls_zero_inverses.iter()).zip(ls.iter())
        {
            a += (a_eval * li_0_inv, li);
        }

        let mut b = DensePolynomial::default();
        for (&b_eval, li) in b_evals.iter().zip(ls.iter()) {
            b += (b_eval, li);
        }

        (zI, a, b, inner_product)
    }

    #[test]
    pub fn test_xrx() {
        let mut rng = test_rng();
        let (zi, a, b, sigma) = setup::<Fr, StdRng>(20, &mut rng);

        let (q, r) = GeneralizedInnerProduct::prove_with_xrx(&a, &b, sigma, &zi);

        let alpha = Fr::rand(&mut rng);

        let lhs = a.evaluate(&alpha) * b.evaluate(&alpha) - sigma;
        let rhs = alpha * r.evaluate(&alpha) + q.evaluate(&alpha) * zi.evaluate(&alpha);
        assert_eq!(lhs, rhs);
    }

    #[test]
    pub fn test_rx() {
        let mut rng = test_rng();
        let (zi, a, b, sigma) = setup::<Fr, StdRng>(20, &mut rng);

        let (q, r) = GeneralizedInnerProduct::prove_with_rx(&a, &b, sigma, &zi);

        let alpha = Fr::rand(&mut rng);

        let lhs = a.evaluate(&alpha) * b.evaluate(&alpha) - sigma;
        let rhs = r.evaluate(&alpha) + q.evaluate(&alpha) * zi.evaluate(&alpha);
        assert_eq!(lhs, rhs);
    }
}
