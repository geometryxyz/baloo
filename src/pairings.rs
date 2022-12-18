use std::iter;

use ark_ec::{AffineCurve, PairingEngine, ProjectiveCurve};
use ark_ff::{One, PrimeField, Zero};
use ark_std::{rand::RngCore, UniformRand};

pub struct EvaluationProof<E: PairingEngine> {
    p: E::G1Affine,
    q: E::G1Affine,
    opening_challenge: E::Fr,
    opening: E::Fr,
}

pub fn batch_kzg_opening_proofs<E: PairingEngine, R: RngCore>(
    evaluation_proofs: &[EvaluationProof<E>],
    rng: &mut R,
) -> (E::G1Affine, E::G1Affine) {
    let u = E::Fr::rand(rng);
    let powers_of_u = iter::successors(Some(E::Fr::one()), |u_pow| Some(u_pow.clone() * u));

    let mut lhs = E::G1Projective::zero();
    let mut rhs = E::G1Projective::zero();

    for (eval_proof, u_pow) in evaluation_proofs.iter().zip(powers_of_u) {
        let u_pow_rep = u_pow.into_repr();
        lhs = lhs + eval_proof.q.mul(u_pow_rep.clone());

        let rhs_i = {
            let q_part = eval_proof
                .q
                .mul((eval_proof.opening_challenge * u_pow).into_repr());
            let p_part = eval_proof.p.mul(u_pow_rep.clone());
            let y_part = E::G1Affine::prime_subgroup_generator()
                .mul((eval_proof.opening * u_pow).into_repr());

            q_part + p_part - y_part
        };

        rhs = rhs + rhs_i
    }

    (lhs.into_affine(), rhs.into_affine())
}
