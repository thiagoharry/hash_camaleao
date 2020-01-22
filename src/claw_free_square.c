#include <stdio.h>
#include <gmp.h>
#include <bsd/stdlib.h>
#include <string.h>
#include <stdint.h>
#include <unistd.h>
#include "mod_math.h"

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

bool has_first_pre_image = true;

void init_rnd(RND *r){
  mpz_init(r -> rnd);
}

#define init_digest(a) init_rnd(a)

void free_rnd(RND *r){
  mpz_clear(r -> rnd);
}

#define free_digest(a) free_rnd(a)

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

void print_msg(MSG *m){
  int i;
  for(i = 0; i < m -> size; i ++){
    int a;
    a = ((m -> msg[i] / 128)?(8):(0));
    a += (((m -> msg[i] / 64) % 2)?(4):(0));
    a += (((m -> msg[i] / 32) % 2)?(2):(0));
    a += (((m -> msg[i] / 16) % 2)?(1):(0));
    if(a < 10)
      printf("%d", a);
    else
      printf("%c", 'a' + (a - 10));
    a = (((m -> msg[i] / 8) % 2)?(8):(0));
    a += (((m -> msg[i] / 4) % 2)?(4):(0));
    a += (((m -> msg[i] / 2) % 2)?(2):(0));
    a += ((m -> msg[i] % 2)?(1):(0));
    if(a < 10)
      printf("%d", a);
    else
      printf("%c", 'a' + (a - 10));
  }
}

void print_hash(MSG *msg, RND *r, DIGEST *digest){
  char *string1, *string2;
  string1 = mpz_get_str(NULL, 10, r -> rnd);
  string2 = mpz_get_str(NULL, 10, digest -> rnd);
  printf("Hash(");
  print_msg(msg);
  printf(",\n     %s) = \n %s\n", string1, string2);
  free(string1);
  free(string2);
}

void random_rnd(PK *pk, RND *rnd){
  mod_random_number(&(rnd -> rnd), pk -> n);
  // Garantindo um resíduo quadrático:
  mpz_mul(rnd -> rnd, rnd -> rnd, rnd -> rnd);
  mpz_mod(rnd -> rnd, rnd -> rnd, pk -> n);
}

#define random_digest(a, b) random_rnd(a, b)

void random_msg(PK *pk, MSG *msg, int size){
  int i;
  msg -> msg = (char *) malloc(size);
  for(i = 0; i < size; i ++)
    msg -> msg[i] = arc4random_uniform(256);
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
  mpz_set(sk -> n, pk -> n);
  // Storing the inverse of square root of 4 module pk:
  {
    mpz_t four;
    mpz_init_set_ui(four, 4);
    root_mod_pq(sk -> sqrt4_1, four, sk -> n, sk -> p, sk -> q);
    mpz_invert(sk -> sqrt4_1, sk -> sqrt4_1, pk -> n);
    mpz_clear(four);
  }
}

void free_keys(PK *pk, SK *sk){
  mpz_clear(pk -> n);
  mpz_clear(sk -> p);
  mpz_clear(sk -> q);
  mpz_clear(sk -> n);
  mpz_clear(sk -> sqrt4_1);
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
  mpz_set(digest -> rnd, r.rnd);
  mpz_clear(r.rnd);
}

void collision(SK *sk, MSG *msg, RND *r, MSG *msg2, RND *r2){
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

#include "main.c"
