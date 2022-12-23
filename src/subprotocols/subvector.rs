use std::{collections::BTreeMap, marker::PhantomData};

use crate::error::Error;
use ark_ff::{batch_inversion, FftField};
use ark_poly::{
    univariate::DensePolynomial, EvaluationDomain, GeneralEvaluationDomain, UVPolynomial,
};
use fast_eval::{PolyProcessor, PolyProcessorStrategy};
// use fast_eval::*;

/*
   Given public low degree extension of vector c and witness vector a, find subvector t of c such that for each a_j there exists some root w_i and value t_i, s.t. a_j = t_i, and c(w_i) = t_i
   Then define:
       - zI(X): vanishing polynomial over roots w_i
       - v(X): low degree extension of w_i^-1
       - t(X): low degree extension of values t_i

   Suppose values in c are distinct and that len is power of 2
*/

pub struct SubvectorExtractor<F: FftField> {
    _f: PhantomData<F>,
}

impl<F: FftField> SubvectorExtractor<F> {
    pub fn compute_subvector_related_oracles(
        a: &[F],
        c: &BTreeMap<F, usize>,
    ) -> Result<
        (
            DensePolynomial<F>,
            DensePolynomial<F>,
            Vec<usize>,
            Vec<usize>,
            Box<dyn PolyProcessor<F>>,
        ),
        Error,
    > {
        let domain_h = GeneralEvaluationDomain::<F>::new(c.len()).unwrap();
        let domain_m = GeneralEvaluationDomain::<F>::new(a.len()).unwrap();

        let mut roots_mapping = BTreeMap::<usize, (F, F, Vec<usize>)>::default();
        let mut col = vec![0usize; domain_m.size()];
        let mut v_evals = Vec::with_capacity(domain_m.size());

        /*
           In order to optimize subvector search we preserve quite complex structure: mapping: root_index => (root, eval, indices of a)
           - BTreeMap will provide us with roots and evals from which we construct t(X) and btree for fast polynomial processing
           - col: [m] -> [k] is given as inverse of mapping: root_i => a_indices, which we also keep in btree map
           - To optimize further we keep array of subvector roots as evals of v(X)
        */
        for (j, aj) in a.iter().enumerate() {
            match c.get(aj) {
                Some(index) => {
                    let eta_i = domain_h.element(*index);
                    v_evals.push(eta_i);
                    let (_, _, a_indices) = roots_mapping.entry(*index).or_insert((eta_i, *aj, vec![]));
                    a_indices.push(j);
                }
                None => return Err(Error::ValueNotInTable),
            }
        }

        batch_inversion(&mut v_evals);
        let v = DensePolynomial::from_coefficients_slice(&domain_m.ifft(&v_evals));

        let mut t_evals = Vec::with_capacity(roots_mapping.keys().len());
        let mut roots = Vec::with_capacity(roots_mapping.keys().len());

        for (j, (eta_i, t_eval, indices)) in roots_mapping.values().enumerate() {
            for &i in indices {
                col[i] = j;
            }

            t_evals.push(*t_eval);
            roots.push(*eta_i);
        }

        let subvector_indices: Vec<usize> = roots_mapping.keys().map(|i| i.clone()).collect();

        let poly_processor = PolyProcessorStrategy::resolve(&roots).unwrap();
        let t = poly_processor.interpolate(&t_evals);

        Ok((v, t, col, subvector_indices, poly_processor))
    }
}
