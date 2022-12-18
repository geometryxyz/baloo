use std::marker::PhantomData;

use ark_ec::{AffineCurve, PairingEngine, ProjectiveCurve};
use ark_ff::{Field, One, PrimeField};
use ark_poly::{EvaluationDomain, GeneralEvaluationDomain};

use crate::data_structures::Proof;

pub struct Verifier<E: PairingEngine> {
    _e: PhantomData<E>,
}

impl<E: PairingEngine> Verifier<E> {
    pub fn verify(
        srs_g1: &[E::G1Affine],
        srs_g2: &[E::G2Affine],
        proof: &Proof<E>,
        cm: &E::G1Affine,
        domain_v: GeneralEvaluationDomain<E::Fr>,
    ) {
        assert_eq!(srs_g1.len(), srs_g2.len());
        let alpha = E::Fr::from(9321381u64);
        let beta = E::Fr::from(231984213152u64);

        let rho = E::Fr::from(1293139123131u64);
        let gamma = E::Fr::from(858374289743829u64);

        // TODO: squeeze pairing separation challenge for batching further

        let gamma_squared = gamma * gamma;
        let gamma_cube = gamma_squared * gamma;

        let g1_gen = E::G1Affine::prime_subgroup_generator();
        let g2_gen = E::G2Affine::prime_subgroup_generator();

        let bound = domain_v.size();
        let d = srs_g1.len();
        let s = d - bound + 1;
        let (x_pow_s_g1, x_pow_s_plus_1_g1, x_pow_m_g1) =
            { (srs_g1[s], srs_g1[s + 1], srs_g1[bound]) };
        let (x_pow_s_g2, x_pow_s_plus_1_g2) = { (srs_g2[s], srs_g2[s + 1]) };

        let p1 = proof.t.mul(proof.u1)
        + g1_gen.mul(-proof.u2)
        + proof.q_2.mul(-proof.u4)
        + proof.r.mul(-E::Fr::one());

        let p2 = proof.v.mul(proof.u5 * beta)
            + g1_gen.mul(-proof.u5 + proof.u3.inverse().unwrap() * proof.u4)
            + proof.q_1.mul(-domain_v.evaluate_vanishing_polynomial(rho));
        

        // TODO - add caulk+ check

        let lhs_x_pow_s: E::G1Affine = {
            (cm.mul(gamma.into_repr())
                + g1_gen
                    .mul(-(proof.u1 + gamma * proof.u2))
                    .add_mixed(&proof.e))
            .into()
        };
        
        // pairing 2
        {
            let lhs_1 = proof.w1.mul(alpha).into_affine();
            let lhs_x = -proof.w1;

            let res = E::product_of_pairings(&[
                (lhs_1.into(), g2_gen.into()),
                (lhs_x.into(), srs_g2[1].into()),
                (lhs_x_pow_s.into(), x_pow_s_g2.into())
            ]);

            assert_eq!(res, E::Fqk::one());
        }

        // pairing 3 
        {
            let lhs_1 = g1_gen.mul(-proof.u3) + proof.r.mul(gamma); 
            let lhs_x = -proof.w2; 

            let mut lhs_zi = x_pow_s_plus_1_g1.mul(gamma_squared);
            lhs_zi.add_assign_mixed(&g1_gen);

            let lhs_x_pow_s_plus_1 = proof.r.mul(gamma_cube);

            let res = E::product_of_pairings(&[
                (lhs_1.into_affine().into(), g2_gen.into()),
                (lhs_x.into(), srs_g2[1].into()),
                (lhs_zi.into_affine().into(), proof.zi.into()),
                (lhs_x_pow_s_plus_1.into_affine().into(), x_pow_s_plus_1_g2.into())
            ]);

            assert_eq!(res, E::Fqk::one());
        }

        // pairing 4
        {
            let mut lhs_1 = proof.w3.mul(beta) + p1.mul(gamma_squared.into_repr()) + g1_gen.mul(-(proof.u1 + gamma * proof.u4));
            lhs_1.add_assign_mixed(&proof.d);

            let lhs_x = -proof.w3;

            let lhs_zi = g1_gen.mul(gamma);

            let res = E::product_of_pairings(&[
                (lhs_1.into_affine().into(), g2_gen.into()),
                (lhs_x.into(), srs_g2[1].into()),
                (lhs_zi.into_affine().into(), proof.zi.into()),
            ]);
    
            assert_eq!(res, E::Fqk::one());
        }

        // pairing 5
        {
            let lhs_1 = proof.w4.mul(rho) + p2.mul(gamma.into_repr()) + g1_gen.mul(-proof.u5);
            let lhs_1 = lhs_1.add_mixed(&proof.e).into_affine();

            let lhs_x = -proof.w4; 
            let res = E::product_of_pairings(&[
                (lhs_1.into(), g2_gen.into()),
                (lhs_x.into(), srs_g2[1].into()),
            ]);
    
            assert_eq!(res, E::Fqk::one());
        }
    }
}
