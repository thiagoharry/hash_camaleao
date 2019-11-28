#include <gmp.h>
#include "mod_math.h"

 #define mpz_rshift(A,B,l) mpz_tdiv_q_2exp(A, B, l)

int root_mod(mpz_t result, const mpz_t arg, const mpz_t prime){
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
