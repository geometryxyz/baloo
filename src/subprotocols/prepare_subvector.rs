use std::{collections::BTreeMap, marker::PhantomData};

use crate::{error::Error, utils::construct_lagrange_basis};
use ark_ff::{batch_inversion, FftField};
use ark_poly::{
    domain, univariate::DensePolynomial, EvaluationDomain, GeneralEvaluationDomain, Polynomial,
    UVPolynomial,
};

/*
   Given public low degree extension of vector c and witness vector a, find subvector t of c such that for each a_j there exists some root w_i and value t_i, s.t. a_j = t_i, and c(w_i) = t_i
   Then define:
       - zI(X): vanishing polynomial over roots w_i
       - v(X): low degree extension of w_i^-1
       - t(X): low degree extension of values t_i

   Suppose values in c are distinct and that len is power of 2
*/

pub struct SubvectorPreprocessor<F: FftField> {
    _f: PhantomData<F>,
}

impl<F: FftField> SubvectorPreprocessor<F> {
    pub fn compute_subvector_related_oracles(
        a: &[F],
        c: &BTreeMap<F, usize>,
    ) -> Result<
        (
            DensePolynomial<F>,
            DensePolynomial<F>,
            DensePolynomial<F>,
            Vec<DensePolynomial<F>>,
        ),
        Error,
    > {
        let domain_m = GeneralEvaluationDomain::<F>::new(a.len()).unwrap();
        let domain_h = GeneralEvaluationDomain::<F>::new(c.len()).unwrap();

        let mut roots_mapping = BTreeMap::<F, (F, Vec<usize>)>::default();
        let mut v_evals = Vec::with_capacity(domain_m.size());

        /*
           In order to optimize subvector search we preserve quite complex structure: mapping: root_of_unity => (value, indices of a)
           - BTreeMap will provide us sorted roots and corresponding values, from which we construct t(X)
           - col: [m] -> [k] is given as inverse of mapping: index of t => index of a that we also keep in map
           - To optimize further we keep array of subvector roots as evals of v(X)
        */
        for (j, aj) in a.iter().enumerate() {
            match c.get(aj) {
                Some(index) => {
                    let eta_i = domain_h.element(*index);
                    v_evals.push(eta_i);
                    let (_, a_indices) = roots_mapping.entry(eta_i).or_insert((*aj, vec![]));
                    a_indices.push(j);
                }
                None => return Err(Error::ValueNotInTable),
            }
        }

        batch_inversion(&mut v_evals);
        let v = DensePolynomial::from_coefficients_slice(&domain_m.ifft(&v_evals));

        let roots: Vec<F> = roots_mapping.keys().map(|eta_i| eta_i.clone()).collect();
        let t_evals: Vec<F> = roots_mapping.values().map(|(ti, _)| ti.clone()).collect();

        let mut zi = DensePolynomial::from_coefficients_slice(&[F::one()]);
        for &eta_i in roots.iter() {
            let current_root = DensePolynomial::from_coefficients_slice(&[eta_i, -F::one()]);
            zi = &zi * &current_root;
        }

        let tau_basis = construct_lagrange_basis(&roots);
        let mut tau_normalizers: Vec<_> = tau_basis
            .iter()
            .map(|tau_i| tau_i.evaluate(&F::zero()))
            .collect();
        batch_inversion(&mut tau_normalizers);

        let tau_normalized_basis: Vec<_> = tau_basis
            .iter()
            .zip(tau_normalizers.iter())
            .map(|(tau_i, &tau_i_zero_inverse)| tau_i * tau_i_zero_inverse)
            .collect();

        let mut t = DensePolynomial::default();
        for (&t_eval, tau_i) in t_evals.iter().zip(tau_basis.iter()) {
            t += (t_eval, tau_i);
        }

        let mut tau_col_j_hat = vec![DensePolynomial::default(); domain_m.size()];
        for (j, (_, indices)) in roots_mapping.values().enumerate() {
            for &i in indices {
                tau_col_j_hat[i] += &tau_normalized_basis[j];
            }
        }

        Ok((zi, v, t, tau_col_j_hat))
    }
}
