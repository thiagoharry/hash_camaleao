#include <string.h>
#include <gmp.h>
#include <bsd/stdlib.h>
#include <stdbool.h>
#include "mod_math.h"

// KeyGen:    0.005490s
// Hash:      0.009530s
// Collision: 0.000032s

typedef struct{
  // x e q tal que q é primo
  mpz_t x, q;
} SK;

typedef struct{
  // y = g^x (mod p), p é primo e p=kq+1
  mpz_t y, g, p, q;
} PK;

typedef struct{
  mpz_t msg;
} MSG;

typedef struct{
  mpz_t rnd;
} RND;

typedef RND DIGEST;

bool has_first_pre_image = false;

void free_keys(PK *pk, SK *sk){
  mpz_clear(sk -> x);
  mpz_clear(sk -> q);
  mpz_clear(pk -> y);
  mpz_clear(pk -> g);
  mpz_clear(pk -> p);
  mpz_clear(pk -> q);
}

void init_digest(DIGEST *digest){
  mpz_init(digest -> rnd);
}

void init_rnd(RND *r){
  mpz_init(r -> rnd);
}

void free_rnd(RND *r){
  mpz_clear(r -> rnd);
}

void free_digest(DIGEST *digest){
  mpz_clear(digest -> rnd);
}

void print_keys(PK *pk, SK *sk){
  char *string1;
  string1 = mpz_get_str(NULL, 10, sk -> x);
  printf("SK=(x = %s)\n", string1);
  free(string1);
  string1 = mpz_get_str(NULL, 10, pk -> y);
  printf("PK=(y = %s)\n", string1);
  free(string1);
}

void print_msg(MSG *m){
  char *string1;
  string1 = mpz_get_str(NULL, 10, m -> msg);
  printf("%s", string1);
  free(string1);
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
  mod_random_number(&(rnd -> rnd), pk -> q);
}

void random_digest(PK *pk, RND *rnd){
  mod_random_number(&(rnd -> rnd), pk -> p);
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
  mpz_init(sk -> q);
  mpz_init(pk -> y);
  mpz_init(pk -> g);
  mpz_init(pk -> p);
  mpz_init(pk -> q);
  // DH Group 14
  // Generating p:
  mpz_set_str(pk -> p, "FFFFFFFFFFFFFFFFC90FDAA22168C234C4C6628B80DC1CD129024E088A67CC74020BBEA63B139B22514A08798E3404DDEF9519B3CD3A431B302B0A6DF25F14374FE1356D6D51C245E485B576625E7EC6F44C42E9A637ED6B0BFF5CB6F406B7EDEE386BFB5A899FA5AE9F24117C4B1FE649286651ECE45B3DC2007CB8A163BF0598DA48361C55D39A69163FA8FD24CF5F83655D23DCA3AD961C62F356208552BB9ED529077096966D670C354E4ABC9804F1746C08CA18217C32905E462E36CE3BE39E772C180E86039B2783A2EC07A28FB5C55DF06F4C52C9DE2BCBF6955817183995497CEA956AE515D2261898FA051015728E5A8AACAA68FFFFFFFFFFFFFFFF", 16);
  // Generating q (q=(p-1)/2):
  mpz_sub_ui(sk -> q, pk -> p, 1);
  mpz_divexact_ui(sk -> q, sk -> q, 2);
  mpz_set(pk -> q, sk -> q);
  // Setting generator g:
  mpz_set_ui(pk -> g, 2);
  // Generating x positive lesser than q
  mod_random_number(&(sk -> x), sk -> q);
  // Setting y = g^x mod p:
  mpz_powm(pk -> y, pk -> g, sk -> x, pk -> p);  
}

void chamhash(PK *pk, MSG *msg, RND *rnd, DIGEST *digest){
  mpz_t aux;
  mpz_init(aux);
  mpz_powm(digest -> rnd, pk -> g, msg -> msg, pk -> p);
  mpz_powm(aux, pk -> y, rnd -> rnd, pk -> p);
  mpz_mul(digest -> rnd, digest -> rnd, aux);
  mpz_mod(digest -> rnd, digest -> rnd, pk -> p);
  mpz_clear(aux);
}

void firstpreimage(SK *sk, MSG *msg, DIGEST *digest, RND *result){

}

void collision(SK *sk, MSG *msg, RND *r, MSG *msg2, RND *r2){
  // r2 = (m +xr - m2)x^(-1) mod q
  mpz_t inv;
  mpz_init(inv);
  // Getting x^(-1)
  mpz_invert(inv, sk -> x, sk -> q);
  // Calculating:
  mpz_mul(r2 -> rnd, sk -> x, r -> rnd); // xr
  mpz_add(r2 -> rnd, r2 -> rnd, msg -> msg); // xr+m
  mpz_sub(r2 -> rnd, r2 -> rnd, msg2 -> msg); // xr+m-m2
  mpz_mul(r2 -> rnd, r2 -> rnd, inv);
  mpz_mod(r2 -> rnd, r2 -> rnd, sk -> q);
  // Ending
  mpz_clear(inv);
}

#include "main.c"
