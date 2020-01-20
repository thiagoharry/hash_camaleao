#include <stdio.h>
#include <gmp.h>
#include <bsd/stdlib.h>
#include <string.h>
#include <stdint.h>
#include <unistd.h>
#include "mod_math.h"

#include "timer.c"

#define MSG_SIZE 64

// KeyGen:    0.43400s
// Hash:      0.00104s
// Collision: 1.79000s

typedef struct{
  // p e q são primos tais que p mod 8 = 3 e q mod 8 = 7
  // n = pq
  // sqrt4_1 é o inverso multiplicativo da raíz de 4 módulo n
  mpz_t p, q, n, sqrt4_1;
} SK;

typedef struct{
  // n = pq
  mpz_t n;
} PK;

typedef struct{
  char *msg;
  int size;
} MSG;

typedef struct{
  mpz_t rnd;
} RND;

typedef RND DIGEST;

void init_digest(DIGEST *digest){
  mpz_init(digest -> rnd);
}

void free_digest(DIGEST *digest){
  mpz_clear(digest -> rnd);
}

void print_keys(PK *pk, SK *sk){
  char *string1, *string2, *string3;
  string1 = mpz_get_str(NULL, 10, sk -> p);
  string2 = mpz_get_str(NULL, 10, sk -> q);
  string3 = mpz_get_str(NULL, 10, pk -> n);
  printf("SK=(p=%s, \n    q=%s)\n", string1, string2);
  printf("PK=(n=%s)\n", string3);
  free(string1);
  free(string2);
  free(string3);
}

void print_rnd(RND *r){
  char *string1;
  string1 = mpz_get_str(NULL, 10, r -> rnd);
  printf("[%s]\n", string1);
  free(string1);
}

void print_hash(RND *r, DIGEST *digest){
  char *string1, *string2;
  string1 = mpz_get_str(NULL, 10, r -> rnd);
  string2 = mpz_get_str(NULL, 10, digest -> rnd);
  printf("Hash(..., %s) = \n %s\n", string1, string2);
  free(string1);
  free(string2);
}

void random_rnd(PK *pk, RND *rnd){
  mod_random_number(&(rnd -> rnd), pk -> n);
  // Garantindo um resíduo quadrático:
  mpz_mul(rnd -> rnd, rnd -> rnd, rnd -> rnd);
  mpz_mod(rnd -> rnd, rnd -> rnd, pk -> n);
}

void random_msg(MSG *msg, int size){
  int i;
  msg -> msg = (char *) malloc(size);
  for(i = 0; i < size; i ++)
    arc4random_uniform(256);
  msg -> size = size;
}

void free_msg(MSG *msg){
  free(msg -> msg);
}

// Permutation 0: x^2 mod q
void P0(PK *pk, RND *x, RND *result){
  mpz_mul(result -> rnd, x -> rnd, x -> rnd);
  mpz_mod(result -> rnd, result -> rnd, pk -> n);
}

// Permutation 1: 4x^2 mod q
void P1(PK *pk, RND *x, RND *result){
  mpz_mul(result -> rnd, x -> rnd, x -> rnd);
  mpz_mul_ui(result -> rnd, result -> rnd, 4);
  mpz_mod(result -> rnd, result -> rnd, pk -> n);
}

// Inverse of P0: sqrt(x)
void iP0(SK *sk, RND *x, RND *result){
  root_mod_pq(result -> rnd, x -> rnd, sk -> n, sk -> p, sk -> q);
}

// Inverse of P1: sqrt(x)/2
void iP1(SK *sk, RND *x, RND *result){
  iP0(sk, x, result);
  mpz_mul(result -> rnd, result -> rnd, sk -> sqrt4_1);
  mpz_mod(result -> rnd, result -> rnd, sk -> n);
}


void keygen(unsigned n, PK *pk, SK *sk){
  unsigned long size_p, size_q;
  char *string;
  int i;
  mpz_t mod;
  mpz_init(sk -> p);
  mpz_init(sk -> q);
  mpz_init(sk -> n);
  mpz_init(sk -> sqrt4_1);
  mpz_init(pk -> n);
  mpz_init(mod);
  size_p = n / 2;
  size_q = (n+1) / 2;
  string = (char *) malloc(size_p + 1);
  do{ // Generating p
    string[0] = '1';
    for(i = 1; i < size_p - 3; i ++)
      if(arc4random_uniform(2))
	string[i] = '1';
      else
	string[i] = '0';
    // Setando os últimos bits desta forma, garantimos p \equiv 3 (mod 8)
    if(size_p - 3 >= 0)
      string[size_p - 3] = '0';
    string[size_p - 2] = '1';
    string[size_p - 1] = '1';
    string[size_p] = '\0';
    mpz_set_str(sk -> p, string, 2);
    mpz_mod_ui(mod, sk -> p, 8);
  }while(mpz_probab_prime_p(sk -> p, 50) <= 0);
  free(string);
  string = (char *) malloc(size_q + 1);
  do{ // Generating q
    string[0] = '1';
    for(i = 1; i < size_q - 3; i ++)
      if(arc4random_uniform(2))
	string[i] = '1';
      else
	string[i] = '0';
    // Setando os últimos bits desta forma, garantimos q \equiv 7 (mod 8)
    string[size_q - 3] = '1';
    string[size_q - 2] = '1';
    string[size_q - 1] = '1';
    string[size_q] = '\0';
    mpz_set_str(sk -> q, string, 2);
    mpz_mod_ui(mod, sk -> q, 8);
  }while(mpz_probab_prime_p(sk -> q, 50) <= 0);
  free(string);
  mpz_clear(mod);
  // Getting public key:
  mpz_mul(pk -> n, sk -> p, sk -> q);
  mpz_init_set(sk -> n, pk -> n);
  // Storing the inverse of square root of 4 module pk:
  {
    mpz_t four;
    mpz_init_set_ui(four, 4);
    root_mod_pq(sk -> sqrt4_1, four, sk -> n, sk -> p, sk -> q);
    mpz_invert(sk -> sqrt4_1, sk -> sqrt4_1, pk -> n);
    mpz_clear(four);
  }
}

void hash(PK *pk, MSG *msg, RND *rnd, DIGEST *digest){
  RND r;
  int i, j;
  mpz_init_set(r.rnd, rnd -> rnd);
  for(i = 0; i < msg -> size; i ++){
    uint8_t c = msg -> msg[i];
    for(j = 0; j < 8; j ++){
	(c / 128)?(P1(pk, &r, &r)):(P0(pk, &r, &r));
	c = c << 1;
    }
  }
  mpz_init_set(digest -> rnd, r.rnd);
  mpz_clear(r.rnd);
}

void firstpreimage(SK *sk, MSG *msg, DIGEST *digest, RND *result){
  int i, j;
  RND r;
  mpz_init_set(r.rnd, digest -> rnd);
  for(i = msg -> size - 1; i >= 0; i --){
    uint8_t c = msg -> msg[i];
    for(j = 0; j < 8; j ++){
      (c % 2)?(iP1(sk, &r, &r)):(iP0(sk, &r, &r));
      c = c >> 1;
    }
  }
  mpz_set(result -> rnd, r.rnd);
  mpz_clear(r.rnd);
}

int main(int argc, char **argv){
  MSG msg1, msg2;
  PK pk;
  SK sk;
  DIGEST digest;
  RND r;
  int security_parameter;
  if(argc >= 2)
    security_parameter = atoi(argv[1]);
  if(argc < 2 || security_parameter == 0){
    fprintf(stderr, "Usage: chamhash SECURITY_PARAMETER [--benchmark]\n");
    exit(1);
  }
  init_digest(&digest);
  keygen(security_parameter, &pk, &sk);
  random_msg(&msg1, MSG_SIZE);
  random_msg(&msg2, MSG_SIZE);
  print_keys(&pk, &sk);
  random_rnd(&pk, &r);
  
  hash(&pk, &msg1, &r, &digest);
  print_hash(&r, &digest);
  firstpreimage(&sk, &msg2, &digest, &r);
  hash(&pk, &msg2, &r, &digest);
  print_hash(&r, &digest);
  free_msg(&msg1);
  free_msg(&msg2);
  free_digest(&digest);
  return 0;
}

/*

#define copyPT(dst, src) mpz_init_set(*dst, src)

void FreePT(PT *pt){
  mpz_clear(*pt);
  free(pt);
}



// Inverse of P0: sqrt(x)
void iP0(PSK psk, PT x, PT *result){
  mpz_t n;
  mpz_t root0, root1;
  mpz_t exp, tmp;
  mpz_init(n);
  mpz_init(root0);
  mpz_init(root1);
  mpz_init(tmp);
  mpz_init(exp);
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
  }while(mpz_cmp(*pt, ppk) >= 0 || mpz_cmp_ui(*pt, 0) == 0);
  free(num);
  // Now return a quadratic residue:
  mpz_mul(*pt, *pt, *pt);
  mpz_mod(*pt, *pt, ppk);
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
    CH -> RandomR(pk, &r);
    TIMER_BEGIN();
    digest = CH -> Hash(pk, (char *) m, MSG_SIZE, r);
    TIMER_END();
    mpz_clear(r);
    FreePT(digest);
  }
  printf("Hash: ");
  TIMER_RESULT();
  // Collision
  for(i = 0; i < N; i ++){
    int j;
    for(j = 0; j < MSG_SIZE; j ++)
      m2[j] = arc4random_uniform(256);
    CH -> RandomR(pk, &r);
    TIMER_BEGIN();
    CH -> FirstPreImage(sk, (char *) m2, MSG_SIZE, r, &r2);
    TIMER_END();
    mpz_clear(r2);
    mpz_clear(r);
  }
  printf("Collision: ");
  TIMER_RESULT();
  CH -> FreePairOfKeys(pksk);
  free(CH);
}

int main(int argc, char **argv){
  int i;
  PAIR_OF_KEYS pksk;
  PK pk;
  SK sk;
  char *string1, *string2;
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
  struct chameleon_hash_scheme *CH = new_chameleon_hash_scheme();
  pksk = CH -> KeyGen(security_parameter);
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
  digest = CH -> Hash(pk, m, strlen(m), r);
  string2 = mpz_get_str(NULL, 10, *digest);
  printf("Hash(\"%s\", %s) = %s\n", m, string1, string2);
  free(string1);
  free(string2);  
  CH -> Collision(sk, (char *) m, strlen(m), r, m2, &r2);
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
*/
