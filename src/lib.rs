pub mod error;
pub mod kzg;
pub mod subprotocols;
pub mod utils;

pub mod data_structures;
pub mod pairings;
pub mod prover;
pub mod verifier;

pub mod fast_eval;

#[cfg(test)]
mod lib_tests {
    use ark_bn254::{Bn254, Fr};
    use ark_ec::PairingEngine;
    use ark_poly::univariate::DensePolynomial;
    use ark_poly::{EvaluationDomain, GeneralEvaluationDomain, UVPolynomial};
    use ark_std::rand::rngs::StdRng;
    use ark_std::rand::RngCore;
    use ark_std::{test_rng, UniformRand};

    use crate::data_structures::{CommitterKey, TableProvingKey, TableWitness};
    use crate::kzg::Kzg;
    use crate::prover::Prover;
    use crate::utils::unsafe_setup_from_rng;
    use crate::verifier::Verifier;

    fn prepare<E: PairingEngine, R: RngCore>(
        h: usize,
        m: usize,
        subvector_positions: &[usize],
        rng: &mut R,
    ) -> (CommitterKey<E>, TableProvingKey<E::Fr>, TableWitness<E::Fr>) {
        assert_eq!(subvector_positions.len(), m);
        let domain_v = GeneralEvaluationDomain::<E::Fr>::new(m).unwrap();

        let (srs_g1, srs_g2) = unsafe_setup_from_rng::<E, R>(h - 1, h - 1, rng);
        let ck = CommitterKey::<E> { srs_g1, srs_g2 };

        let c_evals: Vec<_> = (0..h).map(|_| E::Fr::rand(rng)).collect();
        let phi_evals: Vec<_> = subvector_positions.iter().map(|&i| c_evals[i]).collect();
        let phi = DensePolynomial::from_coefficients_slice(&domain_v.ifft(&phi_evals));

        let table_key = TableProvingKey::from_table_evals(&c_evals);
        let table_witness = TableWitness {
            domain_size: m,
            phi,
            phi_evals: Some(phi_evals),
        };

        (ck, table_key, table_witness)
    }

    #[test]
    fn full_test() {
        let mut rng = test_rng();
        let h = 128;
        let m = 8;
        let domain_v = GeneralEvaluationDomain::<Fr>::new(m).unwrap();

        let subvector_positions = [95usize, 43, 16, 100, 4, 26, 12, 43];

        let (ck, table_key, table_witness) =
            prepare::<Bn254, StdRng>(h, m, &subvector_positions, &mut rng);

        let proof = Prover::prove(&ck, &table_key, &table_witness).unwrap();

        let cm = Kzg::<Bn254>::commit_g1(&ck.srs_g1, &table_witness.phi).into();
        Verifier::verify(&ck.srs_g1, &ck.srs_g2, &proof, &cm, domain_v);
    }
}
