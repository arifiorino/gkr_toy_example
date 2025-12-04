use ark_bn254::Fr;
use ark_std::{Zero,One,UniformRand};
use ark_poly::univariate::DensePolynomial;
use ark_poly::DenseUVPolynomial;
use ark_poly::Polynomial;

const logN:usize = 10;
const N:usize = 1<<logN;

fn precompute(r:&Vec<Fr>) -> Vec<Fr>{
  let mut G = vec![Fr::one()];
  for i in 0..logN{
    let mut G2 = vec![Fr::zero();1<<(i+1)];
    for b in 0..1<<i{
      G2[b<<1]=G[b] * (Fr::one() - r[i]);
      G2[(b<<1)|1]=G[b] * r[i];
    }
    G = G2;
  }
  G
}

fn bookeeping(A:&Vec<Fr>, r:&Vec<Fr>) -> (Vec<Vec<DensePolynomial<Fr>>>,Fr){
  let mut A = A.clone();
  let mut res:Vec<Vec<_>> = (0..logN).map(|i|vec![DensePolynomial::zero();1<<i]).collect();
  for i in 0..logN{
    for b in 0..1<<(logN-1-i){
      res[logN-1-i][b]=DensePolynomial::from_coefficients_vec(vec![A[b],A[b+(1<<(logN-1-i))]-A[b]]);
      A[b] = A[b]*(Fr::one()-r[i]) + A[b+(1<<(logN-1-i))]*r[i];
    }
  }
  let sum = res[0][0].evaluate(&r[logN-1]);
  (res,sum)
}

fn sumcheck_1(book_1:Vec<Vec<DensePolynomial<Fr>>>,
              book_2:Vec<Vec<DensePolynomial<Fr>>>,
              book_3:Vec<Vec<DensePolynomial<Fr>>>)->Vec<DensePolynomial<Fr>>{
  book_1.iter().zip(book_2.iter()).zip(book_3.iter()).map(|((a,b),c)|{
    a.iter().zip(b.iter()).zip(c.iter()).map(|((a,b),c)|{
      a * b * c
    }).reduce(|a,b|a+b).unwrap()
  }).collect()
}

fn verify_sumcheck(sum:Fr, messages:&Vec<DensePolynomial<Fr>>, r:&Vec<Fr>, test:Fr){
  assert!(messages[0].evaluate(&r[logN-1]) == test);
  for i in 0..logN-1{
    assert!(messages[i+1].evaluate(&r[logN-i-2]) ==
            messages[i].evaluate(&Fr::zero()) + messages[i].evaluate(&Fr::one()));
  }
  assert!(messages[logN-1].evaluate(&Fr::zero()) + messages[logN-1].evaluate(&Fr::one()) == sum);
}

fn main() {

  let mut proof_time = std::time::Duration::ZERO;

  let mut rng = rand::thread_rng();
  let a:Vec<_> = (0..N).map(|_|Fr::rand(&mut rng)).collect();
  let b:Vec<_> = (0..N).map(|_|Fr::rand(&mut rng)).collect();
  let d:Vec<_> = (0..N).map(|_|Fr::rand(&mut rng)).collect();
  let c:Vec<_> = a.iter().zip(b.iter()).map(|(x,y)|x * y).collect();
  let e:Vec<_> = c.iter().zip(d.iter()).map(|(x,y)|x * y).collect();

  let g_0:Vec<_> = (0..logN).map(|_|Fr::rand(&mut rng)).collect();
  let (_,e_g_0) = bookeeping(&e, &g_0);

  let g_1:Vec<_> = (0..logN).map(|_|Fr::rand(&mut rng)).collect();

  let start = std::time::Instant::now();
  let pre_g_0 = precompute(&g_0);
  let (book_pre_g_0,pre_g_0_g_1) = bookeeping(&pre_g_0, &g_1);
  let (book_c,c_g_1) = bookeeping(&c, &g_1);
  let (book_d,d_g_1) = bookeeping(&d, &g_1);
  let messages = sumcheck_1(book_pre_g_0,book_c,book_d);
  proof_time += start.elapsed();

  verify_sumcheck(e_g_0, &messages, &g_1, pre_g_0_g_1 * c_g_1 * d_g_1);

  let g_2:Vec<_> = (0..logN).map(|_|Fr::rand(&mut rng)).collect();

  let start = std::time::Instant::now();
  let pre_g_1 = precompute(&g_1);
  let (book_pre_g_1,pre_g_1_g_2) = bookeeping(&pre_g_1, &g_2);
  let (book_a,a_g_2) = bookeeping(&a, &g_2);
  let (book_b,b_g_2) = bookeeping(&b, &g_2);
  let messages = sumcheck_1(book_pre_g_1,book_a,book_b);
  proof_time += start.elapsed();

  verify_sumcheck(c_g_1, &messages, &g_2, pre_g_1_g_2 * a_g_2 * b_g_2);

  println!("{proof_time:?}");
}
