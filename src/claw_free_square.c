#include <stdio.h>
#include <gmp.h>
#include <bsd/stdlib.h>
#include <string.h>
#include <stdint.h>
#include "mod_math.h"

typedef mpz_t PT;
typedef mpz_t PPK;
typedef mpz_t PSK[4];

#define copyPT(dst, src) mpz_init_set(*dst, src)

void FreePT(PT *pt){
  mpz_clear(*pt);
  free(pt);
}

// Permutation 0: x^2 mod q
void P0(PPK mod, PT x, PT *result){
  //printf("P0: %lu -> ", mpz_get_ui(x));
  mpz_mul(*result, x, x);
  mpz_mod(*result, *result, mod);
  //printf("%lu\n", mpz_get_ui(*result));
}

// Permutation 1: 4x^2 mod q
void P1(PPK mod, PT x, PT *result){
  //printf("P1: %lu -> ", mpz_get_ui(x));
  mpz_mul(*result, x, x);
  mpz_mul_ui(*result, *result, 4);
  mpz_mod(*result, *result, mod);
  //printf("%lu\n", mpz_get_ui(*result));
}

// Inverse of P0: sqrt(x)
void iP0(PSK psk, PT x, PT *result){
  mpz_t n;
  mpz_t root0, root1;
  mpz_t exp, tmp;
  //printf("iP0: %lu -> ", mpz_get_ui(x));
  mpz_init(n);
  mpz_init(root0);
  mpz_init(root1);
  mpz_init(tmp);
  mpz_init(exp);
  // Obtaining two roots module p and q:
  root_mod(root0, x, psk[0]);
  // If root0 is a quadratic residue, ok. Otherwise, choose -root0:
  {
    mpz_set(exp, psk[0]);
    mpz_sub_ui(exp, exp, 1);
    mpz_divexact_ui(exp, exp, 2);
    mpz_powm(tmp, root0, exp, psk[0]);
    if(mpz_cmp_ui(tmp, 1) != 0){
      mpz_sub(root0, psk[0], root0);
    }
  }
  root_mod(root1, x, psk[1]);
  // If root1 is a quadratic residue, ok. Otherwise, choose -root1:
  {
    mpz_set(exp, psk[1]);
    mpz_sub_ui(exp, exp, 1);
    mpz_divexact_ui(exp, exp, 2);
    mpz_powm(tmp, root1, exp, psk[1]);
    if(mpz_cmp_ui(tmp, 1) != 0){
      mpz_sub(root1, psk[1], root1);
      // set root1 = -root1
    }
    mpz_clear(exp);
    mpz_clear(tmp);
  }

  // Combining the results with the chinese remainder theorem
  {
    mpz_t prod, sum, p, inv;
    mpz_init(sum);
    mpz_init(p);
    mpz_init(inv);
    mpz_init_set(prod, psk[2]);
    
    mpz_divexact(p, prod, psk[0]);
    mpz_invert(inv, p, psk[0]);
    mpz_mul(p, p, inv);
    mpz_mul(p, p, root0);
    mpz_add(sum, sum, p);

    mpz_divexact(p, prod, psk[1]);
    mpz_invert(inv, p, psk[1]);
    mpz_mul(p, p, inv);
    mpz_mul(p, p, root1);
    mpz_add(sum, sum, p);

    mpz_mod(*result, sum, prod);
    mpz_clear(p);
    mpz_clear(sum);
    mpz_clear(prod);
    mpz_clear(inv);
    mpz_clear(n);
    mpz_clear(root0);
    mpz_clear(root1);
  }
  //printf("%lu\n", mpz_get_ui(*result));
}

// 4x^2 = y       2x = sqrt(y)      x = sqrt(y)/2

// Inverse of P1: sqrt(x)/2
void iP1(PSK psk, PT x, PT *result){
  //printf("iP1: %lu -> ", mpz_get_ui(x));
  mpz_t f0, f1;
  iP0(psk, x, result);
  mpz_mul(*result, *result, psk[3]);
  mpz_mod(*result, *result, psk[2]);
  //printf("%lu\n", mpz_get_ui(*result));
}

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
  // 3 e 7: 2 bits e 3 bits
  do{ // Generating p
    string[0] = '1';
    for(i = 1; i < size_p - 1; i ++)
      if(arc4random_uniform(2))
	string[i] = '1';
      else
	string[i] = '0';
    string[size_p - 1] = '1';
    string[size_p] = '\0';
    mpz_set_str(psk[0], string, 2);
    mpz_mod_ui(mod, psk[0], 8);
  }while((mpz_probab_prime_p(psk[0], 50) <= 0) || (mpz_get_ui(mod) != 3));
  free(string);
  string = (char *) malloc(size_q + 1);
  do{ // Generating q
    string[0] = '1';
    for(i = 1; i < size_q - 1; i ++)
      if(arc4random_uniform(2))
	string[i] = '1';
      else
	string[i] = '0';
    string[size_q - 1] = '1';
    string[size_q] = '\0';
    mpz_set_str(psk[1], string, 2);
    mpz_mod_ui(mod, psk[1], 8);
  }while((mpz_probab_prime_p(psk[1], 50) <= 0) || (mpz_get_ui(mod) != 7));
  free(string);
  mpz_clear(mod);
  // Getting public key:
  mpz_mul(ppk, psk[0], psk[1]);
  mpz_init_set(psk[2], ppk);
  // Storing in psk[3] the inverse of square root of 4 module ppk:
  {
    mpz_t four;
    mpz_init_set_ui(four, 4);
    mpz_init(psk[3]);
    iP0(psk, four, &psk[3]);
    mpz_invert(psk[3], psk[3], ppk);
    mpz_clear(four);
  }
}


void FreePermutationKeys(PSK psk, PPK ppk){
  mpz_clear(ppk);
  mpz_clear(psk[0]);
  mpz_clear(psk[1]);
  mpz_clear(psk[2]);
  mpz_clear(psk[3]);
}

void _RandomPT(PPK ppk, PT *pt){
  char *num;
  size_t size;
  int i;
  mpz_init(*pt);
  size = mpz_sizeinbase(ppk, 2);
  num = (char *) malloc(size + 1);
  // Generate random value uniformly distributed between 0 and ppk
  // Whats the size of ppk in base 32767? (int minimal size)
  do{
    for(i = 0; i < size; i ++){
      if(arc4random_uniform(2))
	num[i] = '1';
      else
	num[i] = '0';
    }
    num[size] = '\0';
    mpz_set_str(*pt, num, 2);
  }while(mpz_cmp(*pt, ppk) >= 0);
  free(num);
  // Now return a quadratic residue:
  mpz_mul(*pt, *pt, *pt);
  mpz_mod(*pt, *pt, ppk);
}

#include "claw_free_permutations_chameleon_hash.h"

int main(int argc, char **argv){
  int i;
  PAIR_OF_KEYS pksk;
  PK pk;
  SK sk;
  char *string1, *string2;
  DIGEST *digest, *digest2;
  RND r, r2;
  char *m =  "0";
  char *m2 = "2";
  /*{// iP1(4) -> 1
    PT result;
    PSK psk;
    PT x;
    mpz_init(result);
    mpz_init_set_ui(psk[0], 3);
    mpz_init_set_ui(psk[1], 7);
    mpz_init_set_ui(psk[2], 21);
    mpz_init_set_ui(psk[3], 4);
    mpz_init_set_ui(x, 4);
    iP1(psk, x, &result);
    printf("sqrt(%lu)/2 mod 21 = %lu \n", mpz_get_ui(x),
	   mpz_get_ui(result));
    mpz_clear(x);
    mpz_clear(result);
    mpz_clear(psk[0]);
    mpz_clear(psk[1]);
    mpz_clear(psk[2]);
    mpz_clear(psk[3]);
    exit(0);
    }*/
  struct chameleon_hash_scheme *CH = new_chameleon_hash_scheme();  
  pksk = CH -> KeyGen(5);
  pk = pksk -> pk;
  sk = pksk -> sk;
  string1 = mpz_get_str(NULL, 10, sk -> psk[0]);
  string2 = mpz_get_str(NULL, 10, sk -> psk[1]);
  printf("SK = (%s, %s)\n", string1, string2);
  free(string1);
  free(string2);
  string1 = mpz_get_str(NULL, 10, pk -> ppk);
  printf("PK = %s\n", string1);
  free(string1);
  CH -> RandomR(pk, &r);
  string1 = mpz_get_str(NULL, 10, r);
  printf("r = %s\n", string1);
  digest = CH -> Hash(pk, m, strlen(m), r);
  string2 = mpz_get_str(NULL, 10, *digest);
  printf("Hash(\"%s\", %s) = %s\n", m, string1, string2);
  free(string1);
  free(string2);
  CH -> Collision(sk, m, strlen(m), r, m2, &r2);
  string1 = mpz_get_str(NULL, 10, r2);
  digest2 = CH -> Hash(pk, m2, strlen(m2), r2);
  string2 = mpz_get_str(NULL, 10, *digest2);
  printf("Hash(\"%s\", %s) = %s\n", m2, string1, string2);
  free(string1);
  free(string2);
  FreePT(digest);
  FreePT(digest2);
  mpz_clear(r);
  mpz_clear(r2);
  CH -> FreePairOfKeys(pksk);
  free(CH);
  return 0;
}
