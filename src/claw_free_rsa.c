#include <stdio.h>
#include <gmp.h>
#include <bsd/stdlib.h>
#include <string.h>
#include <stdint.h>
#include <unistd.h>
#include "mod_math.h"

// KeyGen:    0.4040s
// Hash:      0.0146s
// Collision: 2.1200s

// x^e0^d0 = x^d0^e0 = x (mod n)
// x^e1^d1 = x^d1^e1 = x (mod n)
typedef struct{
  mpz_t d0, d1, n;
} SK;

typedef struct{
  mpz_t e0, e1, n;
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

void free_keys(PK *pk, SK *sk){
  mpz_clear(sk -> d0);
  mpz_clear(sk -> d1);
  mpz_clear(sk -> n);
  mpz_clear(pk -> e0);
  mpz_clear(pk -> e1);
  mpz_clear(pk -> n);
}

void init_digest(DIGEST *digest){
  mpz_init(digest -> rnd);
}

void init_rnd(RND *r){
  mpz_init(r -> rnd);
}

void free_digest(DIGEST *digest){
  mpz_clear(digest -> rnd);
}

void free_rnd(RND *r){
  mpz_clear(r -> rnd);
}

void print_keys(PK *pk, SK *sk){
  char *string1, *string2, *string3;
  string1 = mpz_get_str(NULL, 10, sk -> d0);
  string2 = mpz_get_str(NULL, 10, sk -> d1);
  string3 = mpz_get_str(NULL, 10, sk -> n);
  printf("SK=(%s,\n    %s,\n    %s)\n", string1, string2, string3);
  free(string1);
  free(string2);
  string1 = mpz_get_str(NULL, 10, pk -> e0);
  string2 = mpz_get_str(NULL, 10, pk -> e1);  
  printf("PK=(%s,\n    %s\n    %s)\n", string1, string2, string3);
  free(string1);
  free(string2);
  free(string3);
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

// Permutations
void P0(PK *pk, RND *x, RND *result){
  mpz_powm(result -> rnd, x -> rnd, pk -> e0, pk -> n);
}

void P1(PK *pk, RND *x, RND *result){
  mpz_powm(result -> rnd, x -> rnd, pk -> e1, pk -> n);
}

void iP0(SK *sk, RND *x, RND *result){
  mpz_powm(result -> rnd, x -> rnd, sk -> d0, sk -> n);
}

void iP1(SK *sk, RND *x, RND *result){
  mpz_powm(result -> rnd, x -> rnd, sk -> d1, sk -> n);
}

void keygen(unsigned n, PK *pk, SK *sk){
  // Gerando p, q e n=pq
  unsigned long size_p, size_q;
  mpz_t p, q, carm;
  mpz_init(p);
  mpz_init(q);
  mpz_init(carm);
  mpz_init(sk -> d0);
  mpz_init(sk -> d1);
  mpz_init(sk -> n);
  mpz_init(pk -> e0);
  mpz_init(pk -> e1);
  mpz_init(pk -> n);
  size_p = n / 2;
  size_q = (n+1) / 2;
  do{
    mod_random_prime(p, size_p);
    mod_random_prime(q, size_q);
    // Getting n:
    mpz_mul(pk -> n, p, q);
    mpz_set(sk -> n, pk -> n);
    // Computing λ(n):
    mpz_sub_ui(p, p, 1);
    mpz_sub_ui(q, q, 1);
    mpz_lcm(carm, p, q);
    // Setting e0 and e1:
    mpz_set_ui(pk -> e0, 65537);
    mpz_set_ui(pk -> e1, 257);
    // Checking if the numbers are valid
    mpz_gcd(sk -> d0, carm, pk -> e0);
    mpz_gcd(sk -> d1, carm, pk -> e1);    
  } while(mpz_cmp_ui(sk -> d0, 1) != 0 || mpz_cmp_ui(sk -> d1, 1) != 0);
  // Setting d0 and d1:
  mpz_invert(sk -> d0, pk -> e0, carm);
  mpz_invert(sk -> d1, pk -> e1, carm);
  // Cleaning
  mpz_clear(carm);
  mpz_clear(p);
  mpz_clear(q);
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
