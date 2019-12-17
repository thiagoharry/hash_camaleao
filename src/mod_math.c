#include <gmp.h>
#include <stdlib.h>
#include <bsd/stdlib.h>
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

void root_mod_pq(mpz_t result, const mpz_t arg, const mpz_t n, const mpz_t p,
		const mpz_t q){
  mpz_t exp, tmp;
  mpz_t root0, root1;
  mpz_init(exp);
  mpz_init(tmp);
  mpz_init(root0);
  mpz_init(root1);
  root_mod(root0, arg, p);
  // If root0 is a quadratic residue, ok. Otherwise, choose -root0:
  {
    mpz_set(exp, p);
    mpz_sub_ui(exp, exp, 1);
    mpz_divexact_ui(exp, exp, 2);
    mpz_powm(tmp, root0, exp, p);
    if(mpz_cmp_ui(tmp, 1) != 0){
      mpz_sub(root0, p, root0);
    }
  }
  root_mod(root1, arg, q);
  // If root1 is a quadratic residue, ok. Otherwise, choose -root1:
  {
    mpz_set(exp, q);
    mpz_sub_ui(exp, exp, 1);
    mpz_divexact_ui(exp, exp, 2);
    mpz_powm(tmp, root1, exp, q);
    if(mpz_cmp_ui(tmp, 1) != 0){
      mpz_sub(root1, q, root1);
      // set root1 = -root1
    }
    mpz_clear(exp);
    mpz_clear(tmp);
  }
  // Combining the results with the chinese remainder theorem
  {
    mpz_t prod, sum, r, inv;
    mpz_init(sum);
    mpz_init(r);
    mpz_init(inv);
    mpz_init_set(prod, n);
    
    mpz_divexact(r, prod, p);
    mpz_invert(inv, r, p);
    mpz_mul(r, r, inv);
    mpz_mul(r, r, root0);
    mpz_add(sum, sum, r);

    mpz_divexact(r, prod, q);
    mpz_invert(inv, r, q);
    mpz_mul(r, r, inv);
    mpz_mul(r, r, root1);
    mpz_add(sum, sum, r);

    mpz_mod(result, sum, prod);
    mpz_clear(r);
    mpz_clear(sum);
    mpz_clear(prod);
    mpz_clear(inv);
    mpz_clear(root0);
    mpz_clear(root1);
  }
}

void mod_random_prime(mpz_t p, unsigned n){
  int i;
  char *string = (char *) malloc(n + 1);
  do{ // Generating x positive lesser than q
    string[0] = '1';
    for(i = 1; i < n - 1; i ++)
      if(arc4random_uniform(2))
	string[i] = '1';
      else
	string[i] = '0';
    string[n - 1] = '1';
    string[n] = '\0';
    mpz_set_str(p, string, 2);
  }while(mpz_probab_prime_p(p, 50) <= 0);
  free(string);
}

bool is_quadratic_residue(mpz_t a, mpz_t n){
  // Is 'a' a quadratic residue molulo n?
  mpz_t exp, tmp;
  mpz_init(exp);
  mpz_init(tmp);
  mpz_sub_ui(exp, n, 1);
  mpz_divexact_ui(exp, exp, 2);
  mpz_powm(tmp, a, exp, n);
  if(mpz_cmp_ui(tmp, 1) != 0){
    mpz_clear(exp);
    mpz_clear(tmp);
    return false;
  }
  else{
    mpz_clear(exp);
    mpz_clear(tmp);
    return true;
  }
}
