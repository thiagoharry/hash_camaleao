#include <stdio.h>
#include <gmp.h>
#include <bsd/stdlib.h>
#include <string.h>
#include "mod_math.h"

typedef mpz_t PT;
typedef mpz_t PPK;
typedef mpz_t PSK[2];

#define copyPT(dst, src) mpz_init_set(*dst, src)

void PermutationKeyGen(unsigned n, PPK ppk, PSK psk){
  // Private key: p, q primes such as p = 3 (mod 8) and q = 7 (mod 8).
  // Public key: pq
  // n is the number of bits of pq.
  unsigned long size_p, size_q;
  char *string;
  int i;
  mpz_t mod;
  mpz_init(psk[0]);
  mpz_init(psk[1]);
  mpz_init(ppk);
  mpz_init(mod);
  size_p = n / 2;
  size_q = (n+1) / 2;
  string = (char *) malloc(size_p + 1);
  do{ // Generating p
    string[0] = '1';
    for(i = 1; i < size_p; i ++)
      if(arc4random_uniform(2))
	string[i] = '1';
      else
	string[i] = '0';
    string[size_p] = '\0';
    mpz_set_str(psk[0], string, 2);
    mpz_mod_ui(mod, psk[0], 8);
  }while(mpz_probab_prime_p(psk[0], 50) <= 0 && mpz_get_ui(mod) != 3);
  free(string);
  string = (char *) malloc(size_q + 1);
  do{ // Generating q
    string[0] = '1';
    for(i = 1; i < size_q; i ++)
      if(arc4random_uniform(2))
	string[i] = '1';
      else
	string[i] = '0';
    string[size_q] = '\0';
    mpz_set_str(psk[0], string, 2);
    mpz_mod_ui(mod, psk[0], 8);
  }while(mpz_probab_prime_p(psk[0], 50) <= 0 && mpz_get_ui(mod) != 7);
  free(string);
  // Getting public key:
  mpz_mul(ppk, psk[0], psk[1]);
  mpz_set(psk[2], ppk);
}

// Permutation 0: x^2 mod q
void P0(PPK mod, PT x, PT *result){
  mpz_mul(*result, x, x);
  mpz_tdiv_q(*result, *result, mod);
}

// Permutation 1: 4x^2 mod q
void P1(PPK mod, PT x, PT *result){
  mpz_t four;
  mpz_init(four);
  mpz_mul(*result, x, x);
  mpz_mul(*result, *result, four);
  mpz_tdiv_q(*result, *result, mod);
  mpz_clear(four);
}

// Inverse of P0: sqrt(x)
void iP0(PSK psk, PT x, PT *result){
  mpz_t f0, f1, n;
  copyPT(&f0, psk[0]);
  copyPT(&f1, psk[1]);
  mpz_t root0, root1;
  // Obtaining two roots module p and q:
  root_mod(root0, x, f0);
  root_mod(root1, x, f1);
  // Combining the results with the chinese remainder theorem
  {
    mpz_t prod, sum, p, inv;
    mpz_init(prod);
    mpz_init(sum);
    mpz_init(p);
    mpz_init(inv);
    mpz_mul(prod, f0, f1);
    
    mpz_divexact(p, prod, f0);
    mpz_invert(inv, p, f0);
    mpz_mul(p, p, inv);
    mpz_mul(p, p, root0);
    mpz_add(sum, sum, p);

    mpz_divexact(p, prod, f1);
    mpz_invert(inv, p, f1);
    mpz_mul(p, p, inv);
    mpz_mul(p, p, root1);
    mpz_add(sum, sum, p);

    //mpz_divexact(p, prod, n);
    //mpz_invert(inv, p, n);
    //mpz_mul(p, p, inv);
    //mpz_mul(p, p, x);
    //mpz_sum(sum, sum, p);
    mpz_tdiv_q(*result, sum, prod);
    mpz_clear(p);
    mpz_clear(sum);
    mpz_clear(prod);
    mpz_clear(inv);
  }
}

// Inverse of P1: sqrt(x)/2
void iP1(PSK psk, PT x, PT *result){
  mpz_t f0, f1;
  iP0(psk, x, result);
  mpz_divexact_ui(*result, *result, 2);
}


#include "claw_free_permutations_chameleon_hash.h"

int main(int argc, char **argv){
  PAIR_OF_KEYS pksk;
  PK pk;
  SK sk;
  char *string1, *string2;
  struct chameleon_hash_scheme *CH = new_chameleon_hash_scheme();
  pksk = CH -> KeyGen(10);
  pk = pksk -> pk;
  sk = pksk -> sk;
  free(pksk);
  string1 = mpz_get_str(NULL, 10, sk -> psk[0]);
  string2 = mpz_get_str(NULL, 10, sk -> psk[1]);
  printf("SK = (%s, %s)\n", string1, string2);
  free(string1);
  free(string2);
  string1 = mpz_get_str(NULL, 10, pk -> ppk);
  printf("PK = %s\n", string1);
  free(string1);
  // ...
  free(pk);
  free(sk);
  free(CH);
  return 0;
}
