use ark_ec::PairingEngine;
use ark_ff::FftField;
use ark_poly::{univariate::DensePolynomial, EvaluationDomain, GeneralEvaluationDomain};
use std::collections::BTreeMap;

/*
    Note: We store CommitterKey independently from TableProvingKey since in general there will be multiple TableProvingKeys and CommitterKey will be shared across them
*/

pub struct CommitterKey<E: PairingEngine> {
    pub(crate) srs_g1: Vec<E::G1Affine>,
    pub(crate) srs_g2: Vec<E::G2Affine>,
}

pub struct CommonInput<E: PairingEngine> {
    pub(crate) c_commit: E::G1Affine,
    pub(crate) zh_commit: E::G1Affine,
}

pub struct TableProvingKey<F: FftField> {
    pub(crate) domain_size: usize,
    pub(crate) table_values: Vec<F>,
    pub(crate) table_index_mapping: BTreeMap<F, usize>,
}

impl<F: FftField> TableProvingKey<F> {
    pub fn from_table_evals(c_evals: &[F]) -> Self {
        Self {
            domain_size: c_evals.len(),
            table_values: c_evals.to_vec(),
            table_index_mapping: c_evals
                .iter()
                .enumerate()
                .map(|(i, ci)| (ci.clone(), i))
                .collect(),
        }
    }

    pub fn from_table_poly_and_domain(
        domain_h: GeneralEvaluationDomain<F>,
        c: &DensePolynomial<F>,
    ) -> Self {
        let table_evals = domain_h.fft(&c);
        Self::from_table_evals(&table_evals)
    }
}

pub struct Proof<E: PairingEngine> {
    // Commitments

    // pi_1
    pub(crate) v: E::G1Affine,
    pub(crate) zi: E::G2Affine,
    pub(crate) t: E::G1Affine,

    // pi_2
    pub(crate) d: E::G1Affine,
    pub(crate) r: E::G1Affine,
    pub(crate) q_2: E::G1Affine,

    // pi_3
    pub(crate) e: E::G1Affine,
    pub(crate) q_1: E::G1Affine,

    // Openings
    pub(crate) u1: E::Fr, // e at alpha
    pub(crate) u2: E::Fr, // phi at alpha
    pub(crate) u3: E::Fr, // zi at zero
    pub(crate) u4: E::Fr, // zi at beta
    pub(crate) u5: E::Fr, // e at rho

    // Opening proofs
    pub(crate) a: E::G1Affine,
    pub(crate) w1: E::G1Affine,
    pub(crate) w2: E::G1Affine,
    pub(crate) w3: E::G1Affine,
    pub(crate) w4: E::G1Affine,
}

pub struct TableWitness<F: FftField> {
    pub(crate) domain_size: usize,
    pub(crate) phi: DensePolynomial<F>,
    pub(crate) phi_evals: Option<Vec<F>>,
}
