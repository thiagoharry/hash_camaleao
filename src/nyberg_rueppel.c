#include <string.h>
#include <gmp.h>
#include <bsd/stdlib.h>

#include "timer.c"

// KeyGen:    0.
// Hash:      0.
// Collision: 0.

#define X 0
#define Y 0
#define P 1
#define Q 2
#define G 3

#define MSG_SIZE 54

typedef struct _PK{
  mpz_t y, g, p, q;
} PK;

typedef struct _SK{
  mpz_t x, p, q, g;
  PK *pk;
} SK;

//typedef mpz_t SK[3]; // x, p e q tal que q e p são primos e p=2q+1
//typedef mpz_t PK[4]; // y, g, p e q tal que y=g^x mod p, p é primo e p=2q+1
typedef mpz_t MSG;
typedef mpz_t RND[2];
typedef mpz_t DIGEST;

typedef struct _PAIR_OF_KEYS{
  PK pk;
  SK sk;
} *PAIR_OF_KEYS;

#include "chameleon_hash.h"


void FreeDigest(DIGEST *digest){
  mpz_clear(*digest);
  free(digest);
}

PAIR_OF_KEYS _KeyGen(unsigned n){
  unsigned long size_p, size_q;
  char *string;
  int i;
  PAIR_OF_KEYS pair = (PAIR_OF_KEYS) malloc(sizeof(struct _PAIR_OF_KEYS));
  if(pair == NULL)
    return NULL;
  mpz_init(pair -> sk.x);
  mpz_init(pair -> sk.p);
  mpz_init(pair -> sk.q);
  mpz_init(pair -> sk.g);
  mpz_init(pair -> pk.y);
  mpz_init(pair -> pk.p);
  mpz_init(pair -> pk.q);
  mpz_init(pair -> pk.g);
  
  size_p = n;
  size_q = n - 1;
  // DH Group 14
  // Generating p:
  mpz_set_str(pair -> pk.p, "FFFFFFFFFFFFFFFFC90FDAA22168C234C4C6628B80DC1CD129024E088A67CC74020BBEA63B139B22514A08798E3404DDEF9519B3CD3A431B302B0A6DF25F14374FE1356D6D51C245E485B576625E7EC6F44C42E9A637ED6B0BFF5CB6F406B7EDEE386BFB5A899FA5AE9F24117C4B1FE649286651ECE45B3DC2007CB8A163BF0598DA48361C55D39A69163FA8FD24CF5F83655D23DCA3AD961C62F356208552BB9ED529077096966D670C354E4ABC9804F1746C08CA18217C32905E462E36CE3BE39E772C180E86039B2783A2EC07A28FB5C55DF06F4C52C9DE2BCBF6955817183995497CEA956AE515D2261898FA051015728E5A8AACAA68FFFFFFFFFFFFFFFF", 16);
  mpz_set(pair -> sk.p, pair -> pk.p);
  // Generating q:
  mpz_sub_ui(pair -> sk.q, pair -> pk.p, 1);
  mpz_divexact_ui(pair -> sk.q, pair -> sk.p, 2);
  mpz_set(pair -> pk.q, pair -> sk.q);
  // Setting generator g:
  mpz_set_ui(pair -> pk.g, 2);
  mpz_set_ui(pair -> sk.g, 2);
  string = (char *) malloc(size_q + 1);
  do{ // Generating x positive lesser than q
    string[0] = '1';
    for(i = 1; i < size_q - 1; i ++)
      if(arc4random_uniform(2))
	string[i] = '1';
      else
	string[i] = '0';
    string[size_q - 1] = '1';
    string[size_q] = '\0';
    mpz_set_str(pair -> sk.x, string, 2);
  }while(mpz_cmp(pair -> sk.x, pair -> sk.q) >= 0 ||
    mpz_cmp_ui(pair -> sk.x, 0) == 0);
  free(string);
  // Setting y = g^x mod p:
  mpz_powm(pair -> pk.y, pair -> pk.g, pair -> sk.x, pair -> pk.p);
  pair -> sk.pk = &(pair -> pk);
  return pair;
}

void _FreePairOfKeys(PAIR_OF_KEYS pksk){
  mpz_clear(pksk -> pk.y);
  mpz_clear(pksk -> pk.g);
  mpz_clear(pksk -> pk.p);
  mpz_clear(pksk -> pk.q);
  mpz_clear(pksk -> sk.x);
  mpz_clear(pksk -> sk.p);
  mpz_clear(pksk -> sk.q);
  mpz_clear(pksk -> sk.g);
  free(pksk);
}

void _RandomR(PK pk, RND *rnd){
  int i;
  mpz_init((*rnd)[0]);
  mpz_init((*rnd)[1]);
  char *string = (char*) malloc(2048);
  // Get a random number smaller than q (r1)
  do{
    for(i = 0; i < 2047; i ++)
      if(arc4random_uniform(2))
	string[i] = '1';
      else
	string[i] = '0';
    string[2047] = '\0';
    mpz_set_str((*rnd)[0], string, 2);
  }while(mpz_cmp((*rnd)[0], pk.q) >= 0 || mpz_cmp_ui((*rnd)[0], 0) == 0);
  // Get a random number smaller than q (r2)
  do{
    for(i = 0; i < 2047; i ++)
      if(arc4random_uniform(2))
	string[i] = '1';
      else
	string[i] = '0';
    string[2047] = '\0';
    mpz_set_str((*rnd)[1], string, 2);
  }while(mpz_cmp((*rnd)[1], pk.q) >= 0 || mpz_cmp_ui((*rnd)[1], 0) == 0);
  free(string);
}

DIGEST *_Hash(PK pk, MSG msg, unsigned msg_size, RND rnd){
  char hash[64];
  mpz_t aux;
  DIGEST *digest = (DIGEST *) malloc(sizeof(DIGEST));
  mpz_init(*digest);
  mpz_init(aux);
  // Computing y^e:
  mpz_powm(*digest, pk.y, msg, pk.p);
  // Computing g^r[1]
  mpz_powm(aux, pk.g, rnd[1], pk.p);
  // Multiplying the results, mod p
  mpz_mul(*digest, *digest, aux);
  mpz_mod(*digest, *digest, pk.p);
  // Making r[0] - result:
  mpz_sub(*digest, rnd[0], *digest);
  mpz_mod(*digest, *digest, pk.q);
  mpz_clear(aux);
  return digest;
}

void _FirstPreImage(SK sk, MSG msg, unsigned msg_size, DIGEST C, RND *r2){
  //r2[0] = C + g^k mod p mod q
  unsigned int K = 65537;
  mpz_powm_ui((*r2)[0], sk.g, K, sk.p);
  mpz_add((*r2)[0], (*r2)[0], C);
  mpz_mod((*r2)[0], (*r2)[0], sk.p);
  mpz_mod((*r2)[0], (*r2)[0], sk.q);
  // r2[1] = k-(msg)x mod q
  mpz_mul((*r2)[1], msg, sk.x);
  mpz_ui_sub((*r2)[1], K, (*r2)[1]);
  mpz_mod((*r2)[1], (*r2)[1], sk.q);
}

void _Collision(SK sk, MSG msg, unsigned msg_size, RND r, MSG m2, RND *r2){
  DIGEST *digest = _Hash(*(sk.pk), msg, msg_size, r);
  _FirstPreImage(sk, m2, msg_size, *digest, r2);
  mpz_clear(*digest);
  free(digest);
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

void benchmark(int security_parameter){
  int i;
  PAIR_OF_KEYS pksk;
  DIGEST *digest, *digest2;
  RND r, r2;
  RND m;
  char *string;
  mpz_init(r2[0]);
  mpz_init(r2[1]);
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
  // Hash
  string = (char *) malloc(security_parameter + 1);
  for(i = 0; i < N; i ++){
    int j;
    CH -> RandomR(pksk -> pk, &r);
    CH -> RandomR(pksk -> pk, &m);
    TIMER_BEGIN();
    digest = CH -> Hash(pksk -> pk, m[0], 0, r);
    TIMER_END();
    FreeDigest(digest);
    mpz_clear(r[0]);
    mpz_clear(r[1]);
    mpz_clear(m[0]);
    mpz_clear(m[1]);
  }
  printf("Hash: ");
  TIMER_RESULT();
  CH -> RandomR(pksk -> pk, &r);
  CH -> RandomR(pksk -> pk, &m);
  digest = CH -> Hash(pksk -> pk, m[0], 0, r);
  mpz_clear(m[0]);
  mpz_clear(m[1]);
  // Collision
  for(i = 0; i < N; i ++){
    int j;
    CH -> RandomR(pksk -> pk, &m);
    TIMER_BEGIN();
    CH -> FirstPreImage(pksk -> sk, m[0], 0, *digest, &r2);
    TIMER_END();
    mpz_clear(m[0]);
    mpz_clear(m[1]);
  }
  printf("Collision: ");
  TIMER_RESULT();
  mpz_clear(r[0]);
  mpz_clear(r[1]);
  free(string);
  FreeDigest(digest);
  mpz_clear(r2[0]);
  mpz_clear(r2[1]);
  CH -> FreePairOfKeys(pksk);
  free(CH);
}


int main(int argc, char **argv){
  int i;
  PAIR_OF_KEYS pair;
  char *string1 = NULL, *string2 = NULL, *string3 = NULL;
  DIGEST *digest, *digest2;
  RND r, r2;
  RND m;
  int security_parameter;
  if(argc >= 2)
    security_parameter = atoi(argv[1]);
  if(argc < 2 || (security_parameter != 2048)){
    fprintf(stderr, "Usage: chamhash SECURITY_PARAMETER [--benchmark]\n");
    fprintf(stderr, "Where SECURITY_PARAMETER is 2048.\n");
    exit(1);
  }
  if(argc >= 3 && !strcmp("--benchmark", argv[2])){
    benchmark(security_parameter);
    exit(0);
  }
  mpz_init(r2[0]);
  mpz_init(r2[1]);
  struct chameleon_hash_scheme *CH = new_chameleon_hash_scheme();
  pair = CH -> KeyGen(security_parameter);
  string1 = mpz_get_str(NULL, 10, pair -> sk.x);
  printf("SK = %s\n", string1);
  free(string1);
  string1 = mpz_get_str(NULL, 10, pair -> pk.y);
  printf("PK = %s\n", string1);
  free(string1);
  CH -> RandomR(pair -> pk, &r);
  CH -> RandomR(pair -> pk, &m);
  string1 = mpz_get_str(NULL, 10, r[0]);
  string3 = mpz_get_str(NULL, 10, r[1]);
  digest = CH -> Hash(pair -> pk, m[0], 0, r);
  string2 = mpz_get_str(NULL, 10, *digest);
  printf("Hash(M1, (%s, %s)) = %s\n", string1, string3, string2);
  free(string1);
  free(string2);
  free(string3);
  CH -> FirstPreImage(pair -> sk, m[1], 0, *digest, &r2);
  string1 = mpz_get_str(NULL, 10, r2[0]);
  string3 = mpz_get_str(NULL, 10, r2[1]);
  digest2 = CH -> Hash(pair -> pk, m[1], 0, r2);
  string2 = mpz_get_str(NULL, 10, *digest2);
  printf("Hash(M2, (%s, %s)) = %s\n", string1, string3, string2);
  free(string1);
  free(string2);
  free(string3);
  FreeDigest(digest);
  FreeDigest(digest2);
  mpz_clear(r[0]);
  mpz_clear(r[1]);
  mpz_clear(m[0]);
  mpz_clear(m[1]);
  mpz_clear(r2[0]);
  mpz_clear(r2[1]);
  CH -> FreePairOfKeys(pair);
  free(CH);
  return 0;
}
