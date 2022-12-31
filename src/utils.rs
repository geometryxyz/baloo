use ark_ec::{AffineCurve, PairingEngine};
use ark_ff::{Field, One, PrimeField};
use ark_poly::{univariate::DensePolynomial, UVPolynomial};
use ark_std::rand::RngCore;
use ark_std::UniformRand;
use std::{cmp::max, iter};

/// Create srs from rng
pub fn unsafe_setup_from_rng<E: PairingEngine, R: RngCore>(
    max_power_g1: usize,
    max_power_g2: usize,
    rng: &mut R,
) -> (Vec<E::G1Affine>, Vec<E::G2Affine>) {
    let tau = E::Fr::rand(rng);
    let size = max(max_power_g1 + 1, max_power_g2 + 1);
    let powers_of_tau: Vec<E::Fr> = iter::successors(Some(E::Fr::one()), |p| Some(*p * tau))
        .take(size)
        .collect();

    let g1_gen = E::G1Affine::prime_subgroup_generator();
    let g2_gen = E::G2Affine::prime_subgroup_generator();

    let srs_g1: Vec<E::G1Affine> = powers_of_tau
        .iter()
        .take(max_power_g1 + 1)
        .map(|tp| g1_gen.mul(tp.into_repr()).into())
        .collect();

    let srs_g2: Vec<E::G2Affine> = powers_of_tau
        .iter()
        .take(max_power_g2 + 1)
        .map(|tp| g2_gen.mul(tp.into_repr()).into())
        .collect();
    (srs_g1, srs_g2)
}

/// Create srs from specific tau
pub fn unsafe_setup_from_tau<E: PairingEngine, R: RngCore>(
    max_power_g1: usize,
    max_power_g2: usize,
    tau: E::Fr,
) -> (Vec<E::G1Affine>, Vec<E::G2Affine>) {
    let size = max(max_power_g1 + 1, max_power_g2 + 1);
    let powers_of_tau: Vec<E::Fr> = iter::successors(Some(E::Fr::one()), |p| Some(*p * tau))
        .take(size)
        .collect();

    let g1_gen = E::G1Affine::prime_subgroup_generator();
    let g2_gen = E::G2Affine::prime_subgroup_generator();

    let srs_g1: Vec<E::G1Affine> = powers_of_tau
        .iter()
        .take(max_power_g1 + 1)
        .map(|tp| g1_gen.mul(tp.into_repr()).into())
        .collect();

    let srs_g2: Vec<E::G2Affine> = powers_of_tau
        .iter()
        .take(max_power_g2 + 1)
        .map(|tp| g2_gen.mul(tp.into_repr()).into())
        .collect();
    (srs_g1, srs_g2)
}

pub fn x_pow_d<F: Field>(d: usize) -> DensePolynomial<F> {
    let mut coeffs = vec![F::zero(); d];
    coeffs.push(F::one());
    DensePolynomial::from_coefficients_slice(&coeffs)
}
