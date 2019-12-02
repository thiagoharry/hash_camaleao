#include <stdio.h>
#include <gmp.h>
#include <bsd/stdlib.h>
#include <string.h>
#include <stdint.h>
#include <unistd.h>
#include "mod_math.h"

#define MSG_SIZE 54

typedef mpz_t PT;
typedef mpz_t PPK[3]; // e0, e1, n
typedef mpz_t PSK[3]; // d0, d1, n

/********************* TIMER **************************************/
#include <sys/time.h>
#include <math.h>
#define N 10
unsigned long t_sum = 0;
unsigned long measures[N];
int _i = 0;
#define TIMER_BEGIN() { struct timeval _begin, _end;	\
  gettimeofday(&_begin, NULL);
#define TIMER_END() gettimeofday(&_end, NULL);		  \
  measures[_i] = 1000000 * (_end.tv_sec - _begin.tv_sec) +	\
    _end.tv_usec - _begin.tv_usec;				\
  t_sum += measures[_i];					\
  _i ++;}
#define TIMER_RESULT() {					\
    double mean = ((double) t_sum) / ((double) N);			\
    unsigned long _dif_squared = 0;					\
    for(_i = 0; _i < N; _i ++)						\
      _dif_squared += (measures[_i] - mean) * (measures[i] - mean);	\
    printf("Mean: %.6fs ± %.6fs\n", 0.000001 * mean,			\
	   0.000001 * (sqrt(((double) _dif_squared) / (double) (N-1)))); \
    _i = t_sum = 0;							\
  }
/********************* TIMER **************************************/

#define copyPT(dst, src) mpz_init_set(*dst, src)

void FreePT(PT *pt){
  mpz_clear(*pt);
  free(pt);
}

// Permutation 0: x^2 mod q
void P0(PPK ppk, PT x, PT *result){
  mpz_powm(*result, x, ppk[0], ppk[2]);
}

// Permutation 1: 4x^2 mod q
void P1(PPK ppk, PT x, PT *result){
  mpz_powm(*result, x, ppk[1], ppk[2]);
}

// Inverse of P0: sqrt(x)
void iP0(PSK psk, PT x, PT *result){
  mpz_powm(*result, x, psk[0], psk[2]);
}

// Inverse of P1: sqrt(x)/2
void iP1(PSK psk, PT x, PT *result){
  mpz_powm(*result, x, psk[1], psk[2]);
}

void PermutationKeyGen(unsigned n, PPK ppk, PSK psk){
  // Private key: p, q primes such as p = 3 (mod 8) and q = 7 (mod 8).
  // Public key: pq
  // n is the number of bits of pq.
  unsigned long size_p, size_q;
  char *string;
  mpz_t p, q, carm;
  int i;
  mpz_init(p);
  mpz_init(q);
  mpz_init(psk[0]);
  mpz_init(psk[1]);
  mpz_init(psk[2]);
  mpz_init(ppk[0]);
  mpz_init(ppk[1]);
  mpz_init(ppk[2]);
  size_p = n / 2;
  size_q = (n+1) / 2;
  string = (char *) malloc(size_p + 1);
  do{ // Generating p
    string[0] = '1';
    for(i = 1; i < size_p - 1; i ++)
      if(arc4random_uniform(2))
	string[i] = '1';
      else
	string[i] = '0';
    string[size_p - 1] = '1';
    string[size_p] = '\0';
    mpz_set_str(p, string, 2);
  }while(mpz_probab_prime_p(p, 50) <= 0);
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
    mpz_set_str(q, string, 2);
  }while(mpz_probab_prime_p(q, 50) <= 0);
  free(string);
  // Getting n:
  mpz_mul(ppk[2], p, q);
  mpz_set(psk[2], ppk[2]);
  // Computing λ(n):
  mpz_init(carm);
  mpz_sub_ui(p, p, 1);
  mpz_sub_ui(q, q, 1);
  mpz_lcm(carm, p, q);
  // Setting e0 and e1:
  mpz_set_ui(ppk[0], 65537);
  mpz_set_ui(ppk[1], 257);
  // Setting d0 and d1:
  mpz_invert(psk[0], ppk[0], carm);
  mpz_invert(psk[1], ppk[1], carm);
  // Cleaning
  mpz_clear(carm);
  mpz_clear(p);
  mpz_clear(q);
}


void FreePermutationKeys(PSK psk, PPK ppk){
  mpz_clear(ppk[0]);
  mpz_clear(ppk[1]);
  mpz_clear(ppk[2]);
  mpz_clear(psk[0]);
  mpz_clear(psk[1]);
  mpz_clear(psk[2]);
}

void _RandomPT(PPK ppk, PT *pt){
  char *num;
  size_t size;
  int i;
  mpz_init(*pt);
  size = mpz_sizeinbase(ppk[2], 2);
  num = (char *) malloc(size + 1);
  // Generate random value uniformly distributed between 0 and n
  do{
    for(i = 0; i < size; i ++){
      if(arc4random_uniform(2))
	num[i] = '1';
      else
	num[i] = '0';
    }
    num[size] = '\0';
    mpz_set_str(*pt, num, 2);
  }while(mpz_cmp(*pt, ppk[2]) >= 0);
  free(num);
}

#include "claw_free_permutations_chameleon_hash.h"

void benchmark(int security_parameter){
  int i;
  PAIR_OF_KEYS pksk;
  PK pk;
  SK sk;
  DIGEST *digest, *digest2;
  RND r, r2;
  uint8_t m[MSG_SIZE], m2[MSG_SIZE];
  struct chameleon_hash_scheme *CH = new_chameleon_hash_scheme();
  // Keygen
  for(i = 0; i < N; i ++){
    TIMER_BEGIN();
    pksk = CH -> KeyGen(security_parameter);
    TIMER_END();
    CH -> FreePairOfKeys(pksk);
  }
  printf("KeyGen: ");
  TIMER_RESULT();
  pksk = CH -> KeyGen(security_parameter);
  pk = pksk -> pk;
  sk = pksk -> sk;
  // Hash
  for(i = 0; i < N; i ++){
    int j;
    for(j = 0; j < MSG_SIZE; j ++)
      m[j] = arc4random_uniform(256);
    TIMER_BEGIN();
    digest = CH -> Hash(pk, (char *) m, MSG_SIZE, r);
    TIMER_END();
    FreePT(digest);
  }
  printf("Hash: ");
  TIMER_RESULT();
  digest = CH -> Hash(pk, (char *) m, MSG_SIZE, r);
  // Collision
  for(i = 0; i < N; i ++){
    int j;
    for(j = 0; j < MSG_SIZE; j ++)
      m2[j] = arc4random_uniform(256);
    TIMER_BEGIN();
    CH -> FirstPreImage(sk, (char *) m, MSG_SIZE, *digest, &r2);
    TIMER_END();
  }
  printf("Collision: ");
  TIMER_RESULT();
}

int main(int argc, char **argv){
  int i;
  PAIR_OF_KEYS pksk;
  PK pk;
  SK sk;
  char *string1, *string2, *string3;
  DIGEST *digest, *digest2;
  RND r, r2;
  char *m =  "ChameleonHash";
  char *m2 = "ApplePieApple";
  int security_parameter;
  if(argc >= 2)
    security_parameter = atoi(argv[1]);
  if(argc < 2 || (security_parameter < 8 && security_parameter != 5)){
    fprintf(stderr, "Usage: chamhash SECURITY_PARAMETER [--benchmark]\n");
    fprintf(stderr, "Where SECURITY_PARAMETER is 5 or bigger than 8.\n");
    exit(1);
  }
  if(argc >= 3 && !strcmp("--benchmark", argv[2])){
    benchmark(security_parameter);
    exit(0);
  }
  // 65537 e 17 mod 361
  // iP0: 103 -> 278
  // P0:  278 -> 46
  /*{
    PPK ppk;
    PSK psk;
    mpz_t x, result;
    mpz_init(x);
    mpz_init(result);
    mpz_init_set_ui(ppk[0], 7);
    mpz_init_set_ui(ppk[1], 7);
    mpz_init_set_ui(ppk[2], 55);
    mpz_init_set_ui(psk[0], 23);
    mpz_init_set_ui(psk[1], 23);
    mpz_init_set_ui(psk[2], 55);
    mpz_set_ui(x, 8);
    P0(ppk, x, &result);
    printf("8^7 mod 55 = %lu\n", mpz_get_ui(result));
    iP0(psk, result, &result);
    printf("(8^7)^23 mod 55 = %lu\n", mpz_get_ui(result));
    mpz_clear(x);
    mpz_clear(result);
    exit(0);
    }*/
  struct chameleon_hash_scheme *CH = new_chameleon_hash_scheme();
  pksk = CH -> KeyGen(security_parameter);
  pk = pksk -> pk;
  sk = pksk -> sk;
  string1 = mpz_get_str(NULL, 10, sk -> psk[0]);
  string2 = mpz_get_str(NULL, 10, sk -> psk[1]);
  string3 = mpz_get_str(NULL, 10, sk -> psk[2]);
  printf("SK = (%s, %s, %s)\n", string1, string2, string3);
  free(string1);
  free(string2);
  free(string3);
  string1 = mpz_get_str(NULL, 10, pk -> ppk[0]);
  string2 = mpz_get_str(NULL, 10, pk -> ppk[1]);
  string3 = mpz_get_str(NULL, 10, pk -> ppk[2]);
  printf("PK = (%s, %s, %s)\n", string1, string2, string3);
  free(string1);
  free(string2);
  free(string3);
  CH -> RandomR(pk, &r);
  string1 = mpz_get_str(NULL, 10, r);
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
