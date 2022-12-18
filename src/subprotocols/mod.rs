pub mod generalized_inner_product;
pub mod prepare_subvector;
pub mod well_formation;

#[cfg(test)]
mod subprotocols_tests {
    use std::collections::BTreeMap;

    use ark_bn254::Fr;
    use ark_ff::{FftField, Field, One, UniformRand, Zero};
    use ark_poly::{
        domain::general::GeneralElements, univariate::DensePolynomial, EvaluationDomain,
        GeneralEvaluationDomain, Polynomial, UVPolynomial,
    };
    use ark_std::{
        rand::{rngs::StdRng, RngCore},
        test_rng,
    };

    use super::prepare_subvector::SubvectorPreprocessor;
    use crate::subprotocols::{
        generalized_inner_product::GeneralizedInnerProduct, well_formation::WellFormation,
    };

    fn prepare<F: FftField, R: RngCore>(
        h: usize,
        m: usize,
        subvector_positions: &[usize],
        rng: &mut R,
    ) -> (Vec<F>, Vec<F>) {
        assert_eq!(subvector_positions.len(), m);

        let c_evals: Vec<_> = (0..h).map(|_| F::rand(rng)).collect();
        let a_evals: Vec<_> = subvector_positions.iter().map(|&i| c_evals[i]).collect();

        (c_evals, a_evals)
    }

    #[test]
    fn test_subprotocols() {
        let mut rng = test_rng();
        let h = 128;
        let m = 8;

        let domain_v = GeneralEvaluationDomain::<Fr>::new(m).unwrap();

        let subvector_positions = [95usize, 43, 16, 100, 4, 4, 12, 43];
        let (c_evals, a_evals) = prepare::<Fr, StdRng>(h, m, &subvector_positions, &mut rng);

        let c_mapping: BTreeMap<Fr, usize> = c_evals
            .iter()
            .enumerate()
            .map(|(i, ci)| (ci.clone(), i))
            .collect();

        let (zi, v, t, tau_col_j_hat) =
            SubvectorPreprocessor::compute_subvector_related_oracles(&a_evals, &c_mapping).unwrap();

        let alpha = Fr::rand(&mut rng);

        let mu_alphas = domain_v.evaluate_all_lagrange_coefficients(alpha);

        let mut d = DensePolynomial::default();
        for (&mu_alpha, tau_hat_col_j) in mu_alphas.iter().zip(tau_col_j_hat.iter()) {
            d += (mu_alpha, tau_hat_col_j)
        }

        let phi = DensePolynomial::from_coefficients_slice(&domain_v.ifft(&a_evals));
        let phi_at_alpha = phi.evaluate(&alpha);
        let (q_2, r) = GeneralizedInnerProduct::prove_with_rx(&d, &t, phi_at_alpha, &zi);

        let beta = Fr::rand(&mut rng);

        let lhs = d.evaluate(&beta) * t.evaluate(&beta) - phi_at_alpha;
        let rhs = r.evaluate(&beta) + zi.evaluate(&beta) * q_2.evaluate(&beta);

        // check inner product
        assert_eq!(lhs, rhs);

        let e_evals: Vec<_> = tau_col_j_hat
            .iter()
            .map(|tau| tau.evaluate(&beta))
            .collect();
        let e = DensePolynomial::from_coefficients_slice(&domain_v.ifft(&e_evals));

        // check that E and D are colspace and rowspace encodings of same matrix
        assert_eq!(e.evaluate(&alpha), d.evaluate(&beta));

        let q_1 = WellFormation::prove(&e, &v, &zi, domain_v, beta);

        // check well formation of lookup matrix
        let lhs = e.evaluate(&alpha) * (beta * v.evaluate(&alpha) - Fr::one())
            + zi.evaluate(&beta) * zi.evaluate(&Fr::zero()).inverse().unwrap();
        let rhs = domain_v.evaluate_vanishing_polynomial(alpha) * q_1.evaluate(&alpha);
        assert_eq!(lhs, rhs);
    }
}
