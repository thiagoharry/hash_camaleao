#include <string.h>
#include <gmp.h>
#include <bsd/stdlib.h>
#include <stdbool.h>
#include "mod_math.h"

// KeyGen:    0.
// Hash:      0.
// Collision: 0.


typedef struct{
  mpz_t y, g, p, q;
} PK;

typedef struct{
  mpz_t x, p, q, g;
} SK;

typedef struct{
  mpz_t msg;
} MSG;

typedef struct{
  mpz_t r, s;
} RND;

typedef struct{
  mpz_t rnd;
} DIGEST;

bool has_first_pre_image = true;

void free_keys(PK *pk, SK *sk){
  mpz_clear(sk -> x);
  mpz_clear(sk -> p);
  mpz_clear(sk -> q);
  mpz_clear(sk -> g);
  mpz_clear(pk -> y);
  mpz_clear(pk -> g);
  mpz_clear(pk -> p);
  mpz_clear(pk -> q);
}

void init_digest(DIGEST *digest){
  mpz_init(digest -> rnd);
}

void init_rnd(RND *r){
  mpz_init(r -> r);
  mpz_init(r -> s);
}


void free_rnd(RND *r){
  mpz_clear(r -> r);
  mpz_clear(r -> s);
}

void free_digest(DIGEST *digest){
  mpz_clear(digest -> rnd);
}

void print_keys(PK *pk, SK *sk){
  char *string1, *string2;
  string1 = mpz_get_str(NULL, 10, sk -> x);
  printf("SK=(x = %s)\n", string1);
  free(string1);
  string1 = mpz_get_str(NULL, 10, pk -> y);
  string2 = mpz_get_str(NULL, 10, pk -> g);
  printf("PK=(y = %s, g= %s)\n", string1, string2);
  free(string1);
  free(string2);
}

void print_msg(MSG *m){
  char *string1;
  string1 = mpz_get_str(NULL, 10, m -> msg);
  printf("%s", string1);
  free(string1);
}

void print_hash(MSG *msg, RND *r, DIGEST *digest){
  char *string1, *string2, *string3;
  string1 = mpz_get_str(NULL, 10, r -> r);
  string3 = mpz_get_str(NULL, 10, r -> s);
  string2 = mpz_get_str(NULL, 10, digest -> rnd);
  printf("Hash(");
  print_msg(msg);
  printf(",\n     %s,\n     %s) = \n %s\n", string1, string3, string2);
  free(string1);
  free(string2);
  free(string3);
}

void random_rnd(PK *pk, RND *rnd){
  mod_random_number(&(rnd -> r), pk -> p);
  mod_random_number(&(rnd -> s), pk -> q);
}

void random_digest(PK *pk, DIGEST *digest){
  mod_random_number(&(digest -> rnd), pk -> p);
}

void random_msg(PK *pk, MSG *msg, int size){
  mpz_init(msg -> msg);
  mod_random_number(&(msg -> msg), pk -> q);
}

void free_msg(MSG *msg){
  mpz_clear(msg -> msg);
}

void keygen(unsigned n, PK *pk, SK *sk){
  mpz_init(sk -> x);
  mpz_init(sk -> p);
  mpz_init(sk -> q);
  mpz_init(sk -> g);
  mpz_init(pk -> y);
  mpz_init(pk -> p);
  mpz_init(pk -> q);
  mpz_init(pk -> g);
  // DH Group 14
  // Generating p:
  mpz_set_str(pk -> p, "FFFFFFFFFFFFFFFFC90FDAA22168C234C4C6628B80DC1CD129024E088A67CC74020BBEA63B139B22514A08798E3404DDEF9519B3CD3A431B302B0A6DF25F14374FE1356D6D51C245E485B576625E7EC6F44C42E9A637ED6B0BFF5CB6F406B7EDEE386BFB5A899FA5AE9F24117C4B1FE649286651ECE45B3DC2007CB8A163BF0598DA48361C55D39A69163FA8FD24CF5F83655D23DCA3AD961C62F356208552BB9ED529077096966D670C354E4ABC9804F1746C08CA18217C32905E462E36CE3BE39E772C180E86039B2783A2EC07A28FB5C55DF06F4C52C9DE2BCBF6955817183995497CEA956AE515D2261898FA051015728E5A8AACAA68FFFFFFFFFFFFFFFF", 16);
  mpz_set(sk -> p, pk -> p);
  // Generating q:
  mpz_sub_ui(sk -> q, pk -> p, 1);
  mpz_divexact_ui(sk -> q, sk -> p, 2);
  mpz_set(pk -> q, sk -> q);
  // Choosing h random between 2 and p-2
  {
    mpz_t h, _p, aux;
    mpz_init(h);
    mpz_init(_p);
    mpz_init(aux);
    mpz_sub_ui(_p, sk -> p, 2);
    do{
      mod_random_number(&h, _p);
      // Setting g:
      mpz_sub_ui(_p, sk -> p, 1);
      mpz_invert(aux, sk -> q, sk -> p);
      mpz_mul(aux, aux, _p);
      mpz_mod(aux, aux, sk -> p);
      mpz_powm(pk -> g, h, aux, sk -> p);
    } while(mpz_cmp_ui(pk -> g, 1) == 0);
    mpz_set(sk -> g, pk -> g);
    mpz_clear(h);
    mpz_clear(_p);
    mpz_clear(aux);
  }
  // Random x mod q
  mod_random_number(&(sk -> x), sk -> q);
  // Setting y = g^x mod p:
  mpz_powm(pk -> y, pk -> g, sk -> x, pk -> p);
}

void hash(PK *pk, MSG *msg, RND *rnd, DIGEST *digest){
  mpz_t aux, s_inv;
  mpz_init(aux);
  mpz_init(s_inv);
  // Computing aux = r/s mod q:
  mpz_invert(s_inv, rnd -> s, pk -> q);
  mpz_mul(aux, s_inv, rnd -> r);
  mpz_mod(aux, aux, pk -> q);
  // Computing digest = y^(r/s) mod p
  mpz_powm(digest -> rnd, pk -> y, aux, pk -> p);
  // Computing s_inv = m/s mod q
  mpz_mul(s_inv, s_inv, msg -> msg);
  mpz_mod(s_inv, s_inv, pk -> q);
  // Computing s_inv = g^(m/s) mod p
  mpz_powm(s_inv, pk -> g, s_inv, pk -> p);
  // Computing digest = r - ((digest)(s_inv) mod p) mod q
  mpz_mul(digest -> rnd, digest -> rnd, s_inv);
  mpz_mod(digest -> rnd, digest -> rnd, pk -> p);
  mpz_sub(digest -> rnd, rnd -> r, digest -> rnd);
  mpz_mod(digest -> rnd, digest -> rnd, pk -> q);
  mpz_clear(s_inv);
  mpz_clear(aux);
}

void firstpreimage(SK *sk, MSG *msg, DIGEST *digest, RND *result){
  mpz_t k;
  mpz_init(k);
  // Getting random k
  mod_random_number(&k, sk -> q);
  // r = c + (g^k mod p) mod q
  mpz_powm(result -> r, sk -> g, k, sk -> p);
  mpz_add(result -> r, result -> r, digest -> rnd);
  mpz_mod(result -> r, result -> r, sk -> q);
  // s = (k^-1 (m+xr)) mod q
  mpz_invert(k, k, sk -> q);
  mpz_mul(result -> s, sk -> x, result -> r);
  mpz_add(result -> s, result -> s, msg -> msg);
  //mpz_mod(result -> s, result -> s, sk -> q);
  mpz_mul(result -> s, result -> s, k);
  mpz_mod(result -> s, result -> s, sk -> q);
  mpz_clear(k);
}

void collision(SK *sk, MSG *msg, RND *r, MSG *msg2, RND *r2){
}

#include "main.c"
