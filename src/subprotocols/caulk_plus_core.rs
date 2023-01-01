use std::marker::PhantomData;

use ark_ec::{AffineCurve, PairingEngine};
use ark_ff::Zero;

use crate::precomputed::Precomputed;

pub struct CaulkPlusCore<E: PairingEngine> {
    _e: PhantomData<E>,
}

impl<E: PairingEngine> CaulkPlusCore<E> {
    pub fn compute_quotients(
        ri: &[E::Fr],
        subvector_indices: &[usize],
        precomputed: &Precomputed<E>,
    ) -> (E::G1Affine, E::G1Affine) {
        let mut w1 = E::G1Projective::zero();
        for (i, index) in subvector_indices.iter().enumerate() {
            w1 += &(precomputed.get_w1_i(index).mul(ri[i]));
        }

        let mut w2 = E::G1Projective::zero();
        for (i, index) in subvector_indices.iter().enumerate() {
            w2 += &(precomputed.get_w2_i(index).mul(ri[i]));
        }

        (w1.into(), w2.into())
    }
}
