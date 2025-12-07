use ark_bn254::Fr;
use ark_std::{Zero,One,UniformRand};
use ark_poly::univariate::DensePolynomial;
use ark_poly::DenseUVPolynomial;
use ark_poly::Polynomial;

const logN:usize = 12;
const N:usize = 1<<logN;
const M:usize = 100;

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
      (a.naive_mul(b)).naive_mul(c)
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
  let inputs:Vec<Vec<_>> = (0..M).map(|_|(0..N).map(|_|Fr::rand(&mut rng)).collect()).collect();
  let mut intermediates:Vec<Vec<_>> = vec![(0..N).map(|j|inputs[0][j] * inputs[1][j]).collect()];
  for i in 2..M{
    intermediates.push((0..N).map(|j|intermediates[i-2][j] * inputs[i][j]).collect());
  }
  let mut queries = vec![(0..logN).map(|_|Fr::rand(&mut rng)).collect()];

  let (_,mut com_last) = bookeeping(&intermediates[M-2], &queries[0]);

  for i in 0..M-2{
    queries.push((0..logN).map(|_|Fr::rand(&mut rng)).collect());
    let start = std::time::Instant::now();
    let (book_circuit,com_circuit) = bookeeping(&precompute(&queries[i]), &queries[i+1]);
    let (book_first,com_first) = bookeeping(&intermediates[M-3-i], &queries[i+1]);
    let (book_second,com_second) = bookeeping(&inputs[M-1-i], &queries[i+1]);
    let messages = sumcheck_1(book_circuit,book_first,book_second);
    proof_time += start.elapsed();

    verify_sumcheck(com_last, &messages, &queries[i+1], com_circuit * com_first * com_second);
    com_last = com_first
  }

  queries.push((0..logN).map(|_|Fr::rand(&mut rng)).collect());
  let start = std::time::Instant::now();
  let (book_circuit,com_circuit) = bookeeping(&precompute(&queries[M-2]), &queries[M-1]);
  let (book_first,com_first) = bookeeping(&inputs[0], &queries[M-1]);
  let (book_second,com_second) = bookeeping(&inputs[1], &queries[M-1]);
  let messages = sumcheck_1(book_circuit,book_first,book_second);
  proof_time += start.elapsed();

  verify_sumcheck(com_last, &messages, &queries[M-1], com_circuit * com_first * com_second);

  println!("{proof_time:?}");
}
