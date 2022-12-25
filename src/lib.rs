pub mod error;
pub mod kzg;
pub mod subprotocols;
pub mod utils;

pub mod data_structures;
pub mod precomputed;
pub mod prover;
pub mod verifier;

// pub mod fast_eval;

#[cfg(test)]
mod lib_tests {
    use ark_bn254::{Bn254, Fr};
    use ark_ec::PairingEngine;
    use ark_poly::univariate::DensePolynomial;
    use ark_poly::{EvaluationDomain, GeneralEvaluationDomain, UVPolynomial};
    use ark_std::rand::rngs::StdRng;
    use ark_std::rand::RngCore;
    use ark_std::{test_rng, UniformRand};

    use crate::data_structures::{CommitterKey, CommonInput, TableProvingKey, TableWitness};
    use crate::kzg::Kzg;
    use crate::precomputed::Precomputed;
    use crate::prover::Prover;
    use crate::utils::unsafe_setup_from_rng;
    use crate::verifier::Verifier;

    fn prepare<E: PairingEngine, R: RngCore>(
        h: usize,
        m: usize,
        subvector_positions: &[usize],
        rng: &mut R,
    ) -> (
        CommitterKey<E>,
        TableProvingKey<E::Fr>,
        TableWitness<E::Fr>,
        Precomputed<E>,
        CommonInput<E>,
        E::G1Affine,
    ) {
        assert_eq!(subvector_positions.len(), m);
        let domain_h = GeneralEvaluationDomain::<E::Fr>::new(h).unwrap();
        let domain_v = GeneralEvaluationDomain::<E::Fr>::new(m).unwrap();

        let (srs_g1, srs_g2) = unsafe_setup_from_rng::<E, R>(h, h, rng);

        let zh: DensePolynomial<E::Fr> = domain_h.vanishing_polynomial().into();
        let zh_commit = Kzg::<E>::commit_g1(&srs_g1, &zh);

        let c_evals: Vec<_> = (0..h).map(|_| E::Fr::rand(rng)).collect();
        let c = DensePolynomial::from_coefficients_slice(&domain_h.ifft(&c_evals));
        let c_cm = Kzg::<E>::commit_g1(&srs_g1, &c);

        let phi_evals: Vec<_> = subvector_positions.iter().map(|&i| c_evals[i]).collect();
        let phi = DensePolynomial::from_coefficients_slice(&domain_v.ifft(&phi_evals));
        let cm = Kzg::<E>::commit_g1(&srs_g1, &phi);

        let table_key = TableProvingKey::from_table_evals(&c_evals);
        let table_witness = TableWitness {
            domain_size: m,
            phi,
            phi_evals: Some(phi_evals),
        };

        let mut precomputed = Precomputed::<E>::empty();
        let indices: Vec<usize> = (0..h).collect();

        // precompute all
        precomputed.precompute_w2(&srs_g1, &indices, &domain_h);
        precomputed.precompute_w1(&srs_g1, &indices, &c, &domain_h);

        let ck = CommitterKey::<E> { srs_g1, srs_g2 };

        let common_input = CommonInput::<E> {
            c_commit: c_cm.into(),
            zh_commit: zh_commit.into(),
        };

        (
            ck,
            table_key,
            table_witness,
            precomputed,
            common_input,
            cm.into(),
        )
    }

    #[test]
    fn full_test() {
        let mut rng = test_rng();
        let h = 128;
        let m = 8;
        let domain_v = GeneralEvaluationDomain::<Fr>::new(m).unwrap();

        let subvector_positions = [95usize, 43, 16, 100, 4, 26, 12, 43];

        let (ck, table_key, table_witness, precomputed, common_input, cm) =
            prepare::<Bn254, StdRng>(h, m, &subvector_positions, &mut rng);

        let proof = Prover::prove(&ck, &precomputed, &table_key, &table_witness).unwrap();

        Verifier::verify(&ck.srs_g1, &ck.srs_g2, &common_input, &proof, &cm, domain_v);
    }
}
