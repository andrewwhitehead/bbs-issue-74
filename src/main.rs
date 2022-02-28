use bls12_381::*;
use ff::Field;
use group::{Curve, Group};
use rand::thread_rng;

fn main() {
    // check forged proof issue 74
    let mut rng = thread_rng();
    let k = Scalar::random(&mut rng);
    let w = (G2Projective::generator() * k).to_affine();

    let h0 = G1Projective::random(&mut rng).to_affine();
    let h1 = G1Projective::random(&mut rng).to_affine();
    let h2 = G1Projective::random(&mut rng).to_affine();
    let m = Scalar::random(&mut rng);
    let mh = Scalar::random(&mut rng);
    let s = Scalar::random(&mut rng);
    let e = Scalar::random(&mut rng);
    let b = (G1Projective::generator() + h0 * s + h1 * m + h2 * mh).to_affine();
    let a = (b * (k + e).invert().unwrap()).to_affine();

    assert_eq!(
        pairing(&b, &G2Affine::generator()),
        pairing(&a, &(w + G2Projective::generator() * e).to_affine())
    );

    let r1 = Scalar::random(&mut rng);
    let r2 = Scalar::random(&mut rng);
    let r3 = r1.invert().unwrap();
    let e_rand = Scalar::random(&mut rng);
    let r2_rand = Scalar::random(&mut rng);
    let r3_rand = Scalar::random(&mut rng);
    let s_rand = Scalar::random(&mut rng);
    let mh_rand = Scalar::random(&mut rng);

    let a_prime = (a * r1).to_affine();
    let a_bar = (a_prime * (-e) + b * r1).to_affine();
    let d = (b * r1 + h0 * r2).to_affine();
    let s_prime = s + r2 * r3;

    let c1 = (a_prime * e_rand + h0 * r2_rand).to_affine();
    let c2 = (d * (-r3_rand) + h0 * s_rand + h2 * mh_rand).to_affine();

    let c = Scalar::random(&mut rng);

    let e_resp = e_rand + c * e;
    let r2_resp = r2_rand + c * r2;
    let r3_resp = r3_rand + c * r3;
    let s_resp = s_rand + c * s_prime;
    let mh_resp = mh_rand + c * mh;
    let m_adv = mh_resp * c.invert().unwrap();

    let c1_v = ((G1Projective::from(a_bar) - d) * c + a_prime * e_resp + h0 * r2_resp).to_affine();
    let c2_v =
        ((G1Projective::generator() + h1 * m + h2 * m_adv) * c + d * (-r3_resp) + h0 * s_resp)
            .to_affine();

    assert_eq!(c1, c1_v);
    assert_eq!(c2, c2_v);
    assert_eq!(
        pairing(&a_bar, &G2Affine::generator()),
        pairing(&a_prime, &w)
    );
    println!("done");
}
