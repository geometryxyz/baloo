use std::marker::PhantomData;

use ark_ec::{AffineCurve, PairingEngine, ProjectiveCurve};
use ark_ff::{batch_inversion, Field, One, Zero};
use ark_poly::{
    univariate::DensePolynomial, EvaluationDomain, GeneralEvaluationDomain, Polynomial,
    UVPolynomial,
};

use crate::{
    data_structures::{CommitterKey, Proof, TableProvingKey, TableWitness},
    error::Error,
    kzg::{DegreeBound, Kzg},
    precomputed::{self, Precomputed},
    subprotocols::{
        caulk_plus_core::CaulkPlusCore, generalized_inner_product::GeneralizedInnerProduct,
        subvector::SubvectorExtractor, well_formation::WellFormation,
    },
};

pub struct Prover<E: PairingEngine> {
    _e: PhantomData<E>,
}

impl<E: PairingEngine> Prover<E> {
    pub fn prove(
        ck: &CommitterKey<E>,
        precomputed: &Precomputed<E>,
        // c_cm: &E::G1Affine,
        // cm: &E::G1Affine,
        table_key: &TableProvingKey<E::Fr>,
        table_witness: &TableWitness<E::Fr>,
    ) -> Result<Proof<E>, Error> {
        //TODO: commit to public data

        let domain_v = GeneralEvaluationDomain::<E::Fr>::new(table_witness.domain_size).unwrap();

        let phi_evals = match &table_witness.phi_evals {
            Some(evals) => evals.clone(),
            None => domain_v.fft(&table_witness.phi),
        };

        let (v, t, col, subvector_indices, poly_processor) =
            SubvectorExtractor::compute_subvector_related_oracles(
                &phi_evals,
                &table_key,
            )
            .unwrap();

        let zi = poly_processor.get_vanishing();
        let mut tau_normalizers = poly_processor.batch_evaluate_lagrange_basis(&E::Fr::zero());
        batch_inversion(&mut tau_normalizers);

        let zi_commit = Kzg::<E>::commit_g2(&ck.srs_g2, &zi).into();
        let v_commit = Kzg::<E>::commit_g1(&ck.srs_g1, &v).into();
        let t_commit = Kzg::<E>::commit_g1(&ck.srs_g1, &t).into();

        let alpha = E::Fr::from(9321381u64);

        let mu_alphas = domain_v.evaluate_all_lagrange_coefficients(alpha);

        let mut d_evals = vec![E::Fr::zero(); domain_v.size()];
        for i in 0..domain_v.size() {
            d_evals[col[i]] += tau_normalizers[col[i]] * mu_alphas[i];
        }

        let d = poly_processor.interpolate(&d_evals);

        // let phi = DensePolynomial::from_coefficients_slice(&domain_v.ifft(&a_evals));
        let phi_at_alpha = table_witness.phi.evaluate(&alpha);
        let (q_2, r) = GeneralizedInnerProduct::prove_with_rx(&d, &t, phi_at_alpha, &zi);

        let d_commit = Kzg::<E>::commit_g1(&ck.srs_g1, &d).into();
        let r_commit = Kzg::<E>::commit_g1(&ck.srs_g1, &r).into();
        let q_2_commit = Kzg::<E>::commit_g1(&ck.srs_g1, &q_2).into();

        let beta = E::Fr::from(231984213152u64);

        let lhs = d.evaluate(&beta) * t.evaluate(&beta) - phi_at_alpha;
        let rhs = r.evaluate(&beta) + zi.evaluate(&beta) * q_2.evaluate(&beta);

        // check inner product
        assert_eq!(lhs, rhs);

        let tau_beta = poly_processor.batch_evaluate_lagrange_basis(&beta);

        let e_evals: Vec<_> = (0..domain_v.size())
            .map(|i| tau_normalizers[col[i]] * tau_beta[col[i]])
            .collect();
        let e = DensePolynomial::from_coefficients_slice(&domain_v.ifft(&e_evals));

        // check that E and D are colspace and rowspace encodings of same matrix
        assert_eq!(e.evaluate(&alpha), d.evaluate(&beta));

        let zi_at_zero = zi.evaluate(&E::Fr::zero());
        let zi_at_beta = zi.evaluate(&beta);
        let q_1 = WellFormation::prove(&e, &v, zi_at_zero, zi_at_beta, domain_v, beta);

        let e_commit = Kzg::<E>::commit_g1(&ck.srs_g1, &e).into();
        let q_1_commit = Kzg::<E>::commit_g1(&ck.srs_g1, &q_1).into();

        let rho = E::Fr::from(1293139123131u64);
        let gamma = E::Fr::from(858374289743829u64);

        // Begin with KZG proofs and openings
        let (a1, a2) = CaulkPlusCore::compute_quotients(
            &poly_processor.get_ri(),
            &subvector_indices,
            precomputed,
        );
        let a = a2.mul(gamma).add_mixed(&a1);

        let e_at_alpha = e.evaluate(&alpha);
        let e_at_rho = e.evaluate(&rho);

        let mut p1 = &(&t * e_at_alpha) - &r;
        p1 += (-zi_at_beta, &q_2);
        p1[0] -= phi_at_alpha;
        // sanity
        assert_eq!(E::Fr::zero(), p1.evaluate(&beta));

        let mut p2 = &v * beta;
        p2[0] -= E::Fr::one();
        p2 = &p2 * e_at_rho;
        p2[0] += zi_at_beta * zi_at_zero.inverse().unwrap();
        p2 += (-domain_v.evaluate_vanishing_polynomial(rho), &q_1);
        // sanity
        assert_eq!(E::Fr::zero(), p2.evaluate(&rho));

        let w1 = Kzg::<E>::batch_open_g1(
            &ck.srs_g1,
            &[e.clone(), table_witness.phi.clone()],
            alpha,
            gamma,
            Some(domain_v.size()),
        );
        // TODO: Check why deg(zI) = exactly X^m. Imo it should be < m or to check that it's exactly X^k
        let w2 = Kzg::<E>::batch_open_g1_at_zero_with_bounds(
            &ck.srs_g1,
            &[zi.clone(), r.clone()],
            &[DegreeBound::Exact, DegreeBound::LessThan],
            gamma,
            domain_v.size(),
        );

        let w3 = Kzg::<E>::batch_open_g1(
            &ck.srs_g1,
            &[d.clone(), zi.clone(), p1.clone()],
            beta,
            gamma,
            None,
        );
        let w4 = Kzg::<E>::batch_open_g1(&ck.srs_g1, &[e, p2], rho, gamma, None);

        let proof = Proof {
            zi: zi_commit,
            v: v_commit,
            t: t_commit,
            d: d_commit,
            r: r_commit,
            q_2: q_2_commit,
            e: e_commit,
            q_1: q_1_commit,
            u1: e_at_alpha,
            u2: phi_at_alpha,
            u3: zi_at_zero,
            u4: zi_at_beta,
            u5: e_at_rho,
            a: a.into(),
            w1,
            w2,
            w3,
            w4,
        };

        Ok(proof)
    }
}
