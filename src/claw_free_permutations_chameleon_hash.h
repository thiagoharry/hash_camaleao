#ifndef _CLAW_FREE_CHAMELEON_HASH_
#define _CLAW_FREE_CHAMELEON_HASH_

#include <stdlib.h>
#include <string.h>

// Please, define the following elements before including this file:
// PT:  The type of the permutations domain
// PPK: Permutation public key
// PSK: Permutation secret key
// P0:  (PPK, PT, PT *) a first permutation given a public key
// P1:  (PPK, PT, PT *) a second permutation claw-free with the first
// iP0: (PSK, PT, PT *) The inverse of P1
// iP1: (PSK, PT, PT *) The inverse of P2
// void copyPT(PT *, PT) Copy a PT value
// void PermutationKeyGen(unsigned, &PPK, &PSK): Creates permutation keys
// void FreePermutationKeys(PSK, PPK): frees the permutation keys
// PT *_RandomPT(PPK) A random PT given a PPK
// void FreePT(PT *) Frees a PT memory


typedef char *MSG;
typedef PT RND;
typedef PT DIGEST;

typedef struct _PK{
  unsigned n; // Security parameter
  PPK ppk; // Permutation public key
  void (*P0)(PPK, PT, PT *);
  void (*P1)(PPK, PT, PT *);
} *PK;

typedef struct _SK{
  unsigned n; // Security parameter
  PK pk; // Associated public key
  PSK psk; // Permutation secret key
  void (*iP0)(PSK, PT, PT *);
  void (*iP1)(PSK, PT, PT *);
} *SK;

typedef struct _PAIR_OF_KEYS{
  PK pk;
  SK sk;
} *PAIR_OF_KEYS;

#include "chameleon_hash.h"

PAIR_OF_KEYS _KeyGen(unsigned n){
  SK sk = NULL;
  PK pk = NULL;
  PAIR_OF_KEYS pair = (PAIR_OF_KEYS) malloc(sizeof(struct _PAIR_OF_KEYS));
  if(pair == NULL)
    return NULL;
  sk = (SK) malloc(sizeof(struct _SK));
  if(sk == NULL){
    free(pair);
    return NULL;
  }
  pk = (PK) malloc(sizeof(struct _PK));
  if(pk == NULL){
    free(pair);
    free(sk);
    return NULL;
  }  
  sk -> n = pk -> n = n;
  pk -> P0 = P0;
  pk -> P1 = P1;
  sk -> iP0 = iP0;
  sk -> iP1 = iP1;
  sk -> pk = pk;
  PermutationKeyGen(n, pk -> ppk, sk -> psk);
  pair -> pk = pk;
  pair -> sk = sk;
  return pair;
}

DIGEST *_Hash(PK pk, MSG msg, unsigned msg_size, RND rnd){
  int i, j;
  DIGEST *digest = (DIGEST *) malloc(sizeof(DIGEST));
  copyPT(digest, rnd);
  for(i = 0; i < msg_size; i ++){
    uint8_t c = msg[i];
    for(j = 0; j < 8; j ++){
      //printf("P%d: %lu -> ", (c/128), mpz_get_ui(*digest));
      (c / 128)?(pk -> P1(pk -> ppk, *digest, digest)):
	(pk -> P0(pk -> ppk, *digest, digest));
      //printf("%lu\n", mpz_get_ui(*digest));
      c = c << 1;
    }
  }
  return digest;
}

void _FirstPreImage(SK sk, MSG msg, unsigned msg_size, DIGEST digest,
		    RND *result){
  int i, j;
  copyPT(result, digest);
  for(i = msg_size - 1; i >= 0; i --){ 
    char c = msg[i];
    for(j = 0; j < 8; j ++){
      //printf("iP%d: %lu -> ", (c%2), mpz_get_ui(*result));
      (c % 2)?(sk -> iP1(sk -> psk, *result, result)):
	(sk -> iP0(sk -> psk, *result, result));
      //printf("%lu\n", mpz_get_ui(*result));
      c = c >> 1;
    }
  }
}

void _Collision(SK sk, MSG msg, unsigned msg_size, RND rnd, MSG msg2, RND *r){
  DIGEST *digest = _Hash(sk -> pk, msg, msg_size, rnd);
  _FirstPreImage(sk, msg2, msg_size, *digest, r);
  FreePT(digest);
}

void _FreePairOfKeys(PAIR_OF_KEYS pksk){
  FreePermutationKeys(pksk -> sk -> psk, pksk -> pk -> ppk);
  free(pksk -> sk);
  free(pksk -> pk);
  free(pksk);
}

void _RandomR(PK pk, RND *rnd){
  _RandomPT(pk -> ppk, rnd);
}

struct chameleon_hash_scheme *new_chameleon_hash_scheme(void){
  struct chameleon_hash_scheme *new;
  new = (struct chameleon_hash_scheme *)
    malloc(sizeof(struct chameleon_hash_scheme));
  if(new == NULL)
    return NULL;
  new -> KeyGen = _KeyGen;
  new -> Hash = _Hash;
  new -> FirstPreImage = _FirstPreImage;
  new -> FreePairOfKeys = _FreePairOfKeys;
  new -> RandomR = _RandomR;
  new -> Collision = _Collision;
  new -> IForge = NULL;
  return new;
}



#endif
