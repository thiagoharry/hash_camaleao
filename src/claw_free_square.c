#include <stdio.h>
#include <gmp.h>

 #define mpz_rshift(A,B,l) mpz_tdiv_q_2exp(A, B, l)

typedef mzp_t PT;
typedef mzp_t PPK;
typedef mzp_t PSK[2];

int _root(mpz_t result, const mpz_t arg, const mpz_t prime){
  mpz_t y, b, t;
  unsigned int r, m;
  if (mpz_divisible_p(arg, prime)) {
    mpz_set_ui(result, 0);
    return 1;
  }
  if (mpz_legendre(arg, prime) == -1)
    return -1;
  mpz_init(b);
  mpz_init(t);     
  mpz_init_set_ui(y, 2);
  while(mpz_legendre(y, prime) != -1)
    mpz_add_ui(y, y, 1);
  mpz_sub_ui(result, prime, 1);
  r = mpz_scan1(result, 0);
  mpz_rshift(result, result, r); 
  mpz_powm(y, y, result, prime);   
  mpz_rshift(result, result, 1);
  mpz_powm(b, arg, result, prime); 
  mpz_mul(result, arg, b);
  mpz_mod(result, result, prime);  
  mpz_mul(b, result, b);
  mpz_mod(b, b, prime);  
  while(mpz_cmp_ui(b, 1)){   
    mpz_mul(t, b, b);
    mpz_mod(t, t, prime);
    for(m = 1; mpz_cmp_ui(t, 1); m++){
      mpz_mul(t, t, t);
      mpz_mod(t, t, prime);
    }
    mpz_set_ui(t, 0);
    mpz_setbit(t, r - m - 1);
    mpz_powm(t, y, t, prime); 
    mpz_mul(y, t, t);
    r = m;
    mpz_mul(result, result, t);
    mpz_mod(result, result, prime);
    mpz_mul(b, b, y);
    mpz_mod(b, b, prime);
  }
  mpz_clear(y);
  mpz_clear(b);
  mpz_clear(t);
  return 1;
}


// Permutation 0: x^2 mod q
PT P0(PPK mod, PT x){
  mpz_t result;
  mpz_init(result);
  mpz_mul(result, x, x);
  mpz_tdiv_q(result, result, mod);
  return result;
}

// Permutation 1: 4x^2 mod q
PT P1(PPK mod, PT x){
  mpz_t result, four;
  mpz_init(result);
  mpz_init(four);
  mpz_mul(result, x, x);
  mpz_mul(result, result, four);
  mpz_tdiv_q(result, result, mod);
  mpz_clear(four);
  return result;
}

// Inverse of P0: sqrt(x)
PT iP0(PSK psk, PT x){
  mpz_t f0 = psk[0]; // Factor 0
  mpz_t f1 = psk[1]; // Factor 1
  mpz_t root0, root1;
  // Obtaining two roots module p and q:
  _root(root0, x, f0);
  _root(root1, x, f1);
  // Combining the results with the chinese remainder theorem
  return result;
}


#include "claw_free_permutations_chameleon_hash.h"

