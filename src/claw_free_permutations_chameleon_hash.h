#ifndef _CLAW_FREE_CHAMELEON_HASH_
#define _CLAW_FREE_CHAMELEON_HASH_

#include <stdlib.h>

// Please, define the following elements before including this file:
// PT:  The type of the permutations domain
// PPK: Permutation public key
// PSK: Permutation secret key
// P0:  (PPK, PT) -> PT a first permutation given a public key
// P1:  (PPK, PT) -> PT a second permutation claw-free with the first
// iP0: (PSK, PT) -> PT The inverse of P1
// iP1: (PSK, PT) -> PT The inverse of P2

typedef char * MSG;
typedef PT RND;
typedef PT DIGEST;

typedef struct{
  unsigned n; // Security parameter
  PPK -> ppk; // Permutation public key
  RND (*P0)(PPK, RND);
  RND (*P1)(PPK, RND);
} _PK *PK;

typedef struct{
  unsigned n; // Security parameter
  PK pk; // Associated public key
  PSK psk; // Permutation secret key
  RND (*iP0)(SSK, RND);
  RND (*iP1)(SSK, RND);
} _SK *SK;

typedef struct{
  PK pk;
  SK sk;
} _PAIR_OF_KEYS *PAIR_OF_KEYS;

#include "chameleon_hash.h"

PAIR_OF_KEYS *_KeyGen(unsigned n){
  SK sk = NULL;
  PK pk = NULL;
  PAIR_OF_KEYS pair = (PAIR_OF_KEYS) malloc(sizeof(_PAIR_POF_KEYS));
  if(pair == NULL)
    return NULL;
  sk = (SK) malloc(sizeof(_SK));
  if(sk == NULL){
    free(pair);
    return NULL;
  }
  pk = (PK) malloc(sizeof(_PK));
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
  pair -> pk = pk;
  pair -> sk = sk;
  return pair;
}

DIGEST _Hash(PK pk, MSG msg, RND rnd){
  char *p = MSG;
  DIGEST digest = rnd;
  while(*p != '\0'){
    char c = *p;
    digest = ((c / 128)?(pk -> P1(digest)):(pk -> P0(pk -> ppk, digest)));
    c = c << 1;
    digest = ((c / 128)?(pk -> P1(digest)):(pk -> P0(pk -> ppk, digest)));
    c = c << 1;
    digest = ((c / 128)?(pk -> P1(digest)):(pk -> P0(pk -> ppk, digest)));
    c = c << 1;
    digest = ((c / 128)?(pk -> P1(digest)):(pk -> P0(pk -> ppk, digest)));
    c = c << 1;
    digest = ((c / 128)?(pk -> P1(digest)):(pk -> P0(pk -> ppk, digest)));
    c = c << 1;
    digest = ((c / 128)?(pk -> P1(digest)):(pk -> P0(pk -> ppk, digest)));
    c = c << 1;
    digest = ((c / 128)?(pk -> P1(digest)):(pk -> P0(pk -> ppk, digest)));
    c = c << 1;
    digest = ((c / 128)?(pk -> P1(digest)):(pk -> P0(pk -> ppk, digest)));
    p ++;
  }
  return digest;
}

RND _FirstPreImage(SK sk, MSG msg, DIGEST digest){
  char *p = MSG;
  while(*p != '\0') p ++;
  while(p != msg){
    char c = *p;
    digest = ((c % 2)?(sk -> iP1(digest)):(sk -> iP0(digest)));
    c = c >> 1;
    digest = ((c % 2)?(sk -> iP1(digest)):(sk -> iP0(digest)));
    c = c >> 1;
    digest = ((c % 2)?(sk -> iP1(digest)):(sk -> iP0(digest)));
    c = c >> 1;
    digest = ((c % 2)?(sk -> iP1(digest)):(sk -> iP0(digest)));
    c = c >> 1;
    digest = ((c % 2)?(sk -> iP1(digest)):(sk -> iP0(digest)));
    c = c >> 1;
    digest = ((c % 2)?(sk -> iP1(digest)):(sk -> iP0(digest)));
    c = c >> 1;
    digest = ((c % 2)?(sk -> iP1(digest)):(sk -> iP0(digest)));
    c = c >> 1;
    digest = ((c % 2)?(sk -> iP1(digest)):(sk -> iP0(digest)));
    p --;
  }
  return digest;
}

RND _Collision(SK sk, MSG msg, RND rnd, MSG msg2){
  DIGEST digest = _Hash(sk -> pk, msg, rnd);
  RND result = _FirstPreImage(sk, msg2, digest);
  return result;
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
  new -> Collision = _Collision;
  new -> IForge = NULL;
  return new;
}



#endif
