use ark_bn254::Fr;
use ark_std::{Zero,One,UniformRand};
use ark_poly::univariate::DensePolynomial;
use ark_poly::DenseUVPolynomial;
use ark_poly::Polynomial;

fn precompute(r:&Vec<Fr>) -> Vec<Fr>{
  let mut G = vec![Fr::one();1<<r.len()];
  for i in 0..r.len(){
    for b in (0..1<<i).rev(){
      G[(b<<1)|1] = G[b] * r[i];
      G[b<<1] = G[b] * (Fr::one() - r[i]);
    }
  }
  G
}

fn prove_ABC(res_com:Fr, A:&[Fr], B:&[Fr], C:&[Fr]) -> (Fr, Fr, Fr, Vec<Fr>){
  let mut rng = rand::thread_rng();
  let mut A = A.to_vec();
  let mut B = B.to_vec();
  let mut C = C.to_vec();
  let mut last_msg = DensePolynomial::zero();
  let mut A_com=Fr::zero();
  let mut B_com=Fr::zero();
  let mut C_com=Fr::zero();
  let mut rs=Vec::new();
  for i in (0..6+6+6).rev(){
    println!("{i}");
    let r=Fr::rand(&mut rng);
    let mut msg=DensePolynomial::zero();
    for b in 0..1<<i{
      let A_poly = DensePolynomial::from_coefficients_vec(vec![A[2 * b],A[2 * b + 1]-A[2 * b]]);
      let B_poly = DensePolynomial::from_coefficients_vec(vec![B[2 * b],B[2 * b + 1]-B[2 * b]]);
      let C_poly = if i<6{
        DensePolynomial::from_coefficients_vec(vec![C[2 * b],C[2 * b + 1]-C[2 * b]])
      }else{
        DensePolynomial::from_coefficients_vec(vec![C[2 * b / (1<<(i-5))],C[(2 * b + 1)/(1<<(i-5))]-C[2 * b/(1<<(i-5))]])
      };
      if i==0{
        A_com=A_poly.evaluate(&r);
        B_com=B_poly.evaluate(&r);
        C_com=C_poly.evaluate(&r);
      }
      msg = msg + A_poly.naive_mul(&B_poly).naive_mul(&C_poly);
      A[b] = A[2 * b]*(Fr::one()-r) + A[2 * b + 1]*r;
      B[b] = B[2 * b]*(Fr::one()-r) + B[2 * b + 1]*r;
    }
    if i<6{
      for b in 0..1<<i{
        C[b] = C[2 * b]*(Fr::one()-r) + C[2 * b + 1]*r;
      }
    }
    if i==6+6+6-1{
      assert!(res_com == msg.evaluate(&Fr::zero()) + msg.evaluate(&Fr::one()));
    }else{
      assert!(last_msg.evaluate(rs.last().unwrap()) == msg.evaluate(&Fr::zero()) + msg.evaluate(&Fr::one()));
    }
    last_msg = msg;
    rs.push(r)
  }
  assert!(last_msg.evaluate(rs.last().unwrap()) == A_com * B_com * C_com);
  (A_com, B_com, C_com, rs)
}

fn commit(A:&[Fr], r:&[Fr]) -> Fr{
  let mut A = A.to_vec();
  let n = r.len();
  for i in (1..n).rev(){
    for b in 0..1<<i{
      A[b] = A[b]*(Fr::one()-r[n-1-i]) + A[b+(1<<i)]*r[n-1-i];
    }
  }
  let A_poly = DensePolynomial::from_coefficients_vec(vec![A[0],A[1]-A[0]]);
  A_poly.evaluate(&r[n-1])
}

fn main() {
  let mut rng = rand::thread_rng();
  let A:Vec<_> = (0..64*64*64).map(|_|Fr::rand(&mut rng)).collect();
  let B:Vec<_> = (0..64).map(|_|Fr::rand(&mut rng)).collect();
  let C:Vec<_> = (0..64*64*64).map(|i|A[i] * B[i/64/64]).collect();
  let u:Vec<_> = (0..6+6+6).map(|_|Fr::rand(&mut rng)).collect();
  let circuit_eval=precompute(&u);
  let com_C=commit(&C,&u);
  let start=std::time::Instant::now();
  let (com_circuit, com_A, com_B, v) = prove_ABC(com_C,&circuit_eval,&A,&B);
  println!("{:?}",start.elapsed());
}
