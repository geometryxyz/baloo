use std::collections::BTreeMap;

use ark_ec::PairingEngine;
use ark_ff::One;
use ark_poly::{
    univariate::DensePolynomial, EvaluationDomain, GeneralEvaluationDomain, Polynomial,
    UVPolynomial,
};
// use ark_serialize::{CanonicalDeserialize, CanonicalSerialize, Read, SerializationError, Write};

use crate::kzg::Kzg;

/*
   Precomputed data will be stored in <key, value> map for key = index, value = [w_{1,2}^i]_2

   We can precompute all data, but it's very possible that just some indices will be needed,
   so we optimize precomputed data needed to store
*/
// #[derive(CanonicalSerialize, CanonicalDeserialize)]
pub struct Precomputed<E: PairingEngine> {
    w1_mapping: BTreeMap<usize, E::G1Affine>,
    w2_mapping: BTreeMap<usize, E::G1Affine>,
}

impl<E: PairingEngine> Precomputed<E> {
    pub fn empty() -> Self {
        Self {
            w1_mapping: BTreeMap::default(),
            w2_mapping: BTreeMap::default(),
        }
    }
    pub fn get_w1_i(&self, index: &usize) -> E::G1Affine {
        match self.w1_mapping.get(index) {
            Some(element) => *element,
            None => panic!("Element on index: {} is not precomputed", index),
        }
    }

    pub fn get_w2_i(&self, index: &usize) -> E::G1Affine {
        match self.w2_mapping.get(index) {
            Some(element) => *element,
            None => panic!("Element on index: {} is not precomputed", index),
        }
    }

    pub fn precompute_w1(
        &mut self,
        srs: &[E::G1Affine],
        indices: &[usize],
        c: &DensePolynomial<E::Fr>,
        domain: &GeneralEvaluationDomain<E::Fr>,
    ) {
        for index in indices {
            let w_i = domain.element(*index);
            let mut num = c.clone();
            num[0] -= c.evaluate(&w_i);

            let denom = DensePolynomial::from_coefficients_slice(&[-w_i, E::Fr::one()]);
            let w1_i = &num / &denom;
            let w1_i = Kzg::<E>::commit_g1(srs, &w1_i);
            self.w1_mapping.insert(*index, w1_i.into());
        }
    }

    pub fn precompute_w2(
        &mut self,
        srs: &[E::G1Affine],
        indices: &[usize],
        domain: &GeneralEvaluationDomain<E::Fr>,
    ) {
        let zh: DensePolynomial<_> = domain.vanishing_polynomial().into();
        for index in indices {
            let w2_i = &zh
                / &DensePolynomial::from_coefficients_slice(&[
                    -domain.element(*index),
                    E::Fr::one(),
                ]);
            let w2_i = Kzg::<E>::commit_g1(srs, &w2_i);
            self.w2_mapping.insert(*index, w2_i.into());
        }
    }
}

#[cfg(test)]
pub mod precomputed_tests {
    use ark_bn254::{Bn254, Fr, G1Affine, Fq12};
    use ark_ec::{AffineCurve, PairingEngine, ProjectiveCurve};
    use ark_ff::{One, Zero};
    use ark_poly::{EvaluationDomain, GeneralEvaluationDomain, univariate::DensePolynomial};
    use ark_std::{rand::rngs::StdRng, test_rng};
    use fast_eval::PolyProcessorStrategy;

    use crate::{kzg::Kzg, utils::unsafe_setup_from_rng};

    use super::Precomputed;

    #[test]
    fn test_zh() {
        let mut rng = test_rng();
        let n = 64;
        let domain = GeneralEvaluationDomain::<Fr>::new(n).unwrap();
        let zh: DensePolynomial<_> = domain.vanishing_polynomial().into();

        let (srs_g1, srs_g2) = unsafe_setup_from_rng::<Bn254, StdRng>(n + 1, n + 1, &mut rng);

        let mut precomputed = Precomputed::<Bn254>::empty();

        // precompute all
        let indices: Vec<usize> = (0..n).collect();
        precomputed.precompute_w2(&srs_g1, &indices, &domain);

        let subvector_indices = [4usize, 12, 19, 25, 31, 45, 61, 60];
        let roots: Vec<Fr> = subvector_indices
            .iter()
            .map(|index| domain.element(*index))
            .collect();

        let poly_processor = PolyProcessorStrategy::resolve(&roots).unwrap();
        let zi = poly_processor.get_vanishing();

        let ri = poly_processor.get_ri();
        assert_ne!(ri, vec![Fr::one()]); // tmp check while it's not implemented for FftBased Processor

        // commitment phase
        let mut w2 = G1Affine::zero().into_projective();
        for (i, index) in subvector_indices.iter().enumerate() {
            w2 += &(precomputed.get_w2_i(index).mul(ri[i]));
        }

        let zi = Kzg::<Bn254>::commit_g2(&srs_g2, &zi);
        let zh = Kzg::<Bn254>::commit_g1(&srs_g1, &zh);

        let res = Bn254::product_of_pairings(&[
            ((-zh).into_affine().into(), srs_g2[0].into()), 
            (w2.into_affine().into(), zi.into_affine().into())
        ]);
        assert_eq!(res, Fq12::one());
    }

    #[test]
    fn test_c() {}
}
