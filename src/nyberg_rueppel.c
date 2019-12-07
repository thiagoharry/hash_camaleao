#include <string.h>
#include <gmp.h>
#include <bsd/stdlib.h>

// KeyGen:    0.
// Hash:      0.
// Collision: 0.

#define X 0
#define Y 0
#define P 1
#define Q 2
#define G 3

#define MSG_SIZE 54

typedef mpz_t SK[3]; // x, p e q tal que q e p são primos e p=2q+1
typedef mpz_t PK[4]; // y, g, p e q tal que y=g^x mod p, p é primo e p=2q+1
typedef char MSG[MSG_SIZE];
typedef mpz_t RND[2];
typedef mpz_t DIGEST;

typedef struct _PAIR_OF_KEYS{
  PK pk;
  SK sk;
} *PAIR_OF_KEYS;

#include "chameleon_hash.h"

/********************* TIMER **************************************/
#include <sys/time.h>
#include <math.h>
#define N 10000
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
  mpz_init(pair -> sk[X]);
  mpz_init(pair -> sk[P]);
  mpz_init(pair -> sk[Q]);
  mpz_init(pair -> pk[Y]);
  mpz_init(pair -> pk[P]);
  mpz_init(pair -> pk[Q]);
  mpz_init(pair -> pk[G]);
  
  size_p = n;
  size_q = n - 1;
  // DH Group 14
  // Generating p:
  mpz_set_str(pair -> pk[P], "FFFFFFFFFFFFFFFFC90FDAA22168C234C4C6628B80DC1CD129024E088A67CC74020BBEA63B139B22514A08798E3404DDEF9519B3CD3A431B302B0A6DF25F14374FE1356D6D51C245E485B576625E7EC6F44C42E9A637ED6B0BFF5CB6F406B7EDEE386BFB5A899FA5AE9F24117C4B1FE649286651ECE45B3DC2007CB8A163BF0598DA48361C55D39A69163FA8FD24CF5F83655D23DCA3AD961C62F356208552BB9ED529077096966D670C354E4ABC9804F1746C08CA18217C32905E462E36CE3BE39E772C180E86039B2783A2EC07A28FB5C55DF06F4C52C9DE2BCBF6955817183995497CEA956AE515D2261898FA051015728E5A8AACAA68FFFFFFFFFFFFFFFF", 16);
  mpz_set(pair -> sk[P], pair -> pk[P]);
  // Generating q:
  mpz_sub_ui(pair -> sk[Q], pair -> pk[P], 1);
  mpz_divexact_ui(pair -> sk[Q], pair -> sk[P], 2);
  mpz_set(pair -> pk[Q], pair -> sk[Q]);
  // Setting generator g:
  mpz_set_ui(pair -> pk[G], 2);
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
    mpz_set_str(pair -> sk[X], string, 2);
  }while(mpz_cmp(pair -> sk[X], pair -> sk[Q]) >= 0 ||
    mpz_cmp_ui(pair -> sk[X], 0) == 0);
  free(string);
  // Setting y = g^x mod p:
  mpz_powm(pair -> pk[Y], pair -> pk[G], pair -> sk[X], pair -> pk[P]);  
  return pair;
}

void _FreePairOfKeys(PAIR_OF_KEYS pksk){
  mpz_clear(pksk -> pk[Y]);
  mpz_clear(pksk -> pk[G]);
  mpz_clear(pksk -> pk[P]);
  mpz_clear(pksk -> pk[Q]);
  mpz_clear(pksk -> sk[X]);
  mpz_clear(pksk -> sk[P]);
  mpz_clear(pksk -> sk[Q]);
  free(pksk);
}

void _RandomR(PK pk, RND *rnd){
  int i;
  mpz_init(*rnd[0]);
  mpz_init(*rnd[1]);
  char *string = (char*) malloc(2048);
  // Get a random number smaller than q (r1)
  do{
    for(i = 0; i < 2047; i ++)
      if(arc4random_uniform(2))
	string[i] = '1';
      else
	string[i] = '0';
    string[2047] = '\0';
    mpz_set_str(*rnd[0], string, 2);
  }while(mpz_cmp(*rnd, pk[Q]) >= 0 || mpz_cmp_ui(*rnd, 0) == 0);
  // Get a random number smaller than q (r2)
  do{
    for(i = 0; i < 2047; i ++)
      if(arc4random_uniform(2))
	string[i] = '1';
      else
	string[i] = '0';
    string[2047] = '\0';
    mpz_set_str(*rnd[1], string, 2);
  }while(mpz_cmp(*rnd[1], pk[Q]) >= 0 || mpz_cmp_ui(*rnd, 0) == 0);
  free(string);
}

DIGEST *_Hash(PK pk, MSG msg, unsigned msg_size, RND rnd){
  char hash[64];
  mpz_t aux;
  DIGEST *digest = (DIGEST *) malloc(sizeof(DIGEST));
  mpz_init(*digest);
  mpz_init(aux);
  // Computing y^e:
  mpz_powm(*digest, pk[Y], msg, pk[P]);
  // Computing g^r[1]
  mpz_powm(aux, pk[G], rnd[1], pk[P]);
  // Multiplying the results, mod p
  mpz_mul(*digest, *digest, aux);
  mpz_mod(*digest, *digest, pk[P]);
  // Making r[0] - result:
  mpz_sub(*digest, rnd[0], *digest);
  mpz_mod(*digest, *digest, pk[Q]);
  mpz_clear(aux);
  return digest;
}

void _FirstPreImage(SK sk, MSG msg, unsigned msg_size, DIGEST C, RND *r2){
  //r2[0] = C + g^k mod p mod q
  unsigned int K = 65537;
  mpz_powm_ui(*r2[0], sk[G], K, sk[P]);
  mpz_add(*r2[0], *r2[0], C);
  mpz_mod(*r2[0], *r2[0], sk[P]);
  mpz_mod(*r2[0], *r2[0], sk[Q]);
  // r2[1] = k-(msg)x mod q
  mpz_mul(*r2[1], msg, sk[X]);
  mpz_ui_sub(*r2[1], K, *r2[1]);
  mpz_mod(*r2[1], *r2[1], sk[Q]);
}

void _Collision(SK sk, MSG msg, unsigned msg_size, RND r, MSG m2, RND *r2){
  DIGEST *digest = _Hash(sk -> pk, msg, msg_size, rnd);
  _FirstPreImage(sk, msg2, msg_size, *digest, r);
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
  mpz_init(r2);
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
    mpz_clear(r);
    mpz_clear(m);
  }
  printf("Hash: ");
  TIMER_RESULT();
  CH -> RandomR(pksk -> pk, &r);
  CH -> RandomR(pksk -> pk, &m);
  digest = CH -> Hash(pksk -> pk, m[0], 0, r);
  // Collision
  for(i = 0; i < N; i ++){
    int j;
    CH -> RandomR(pksk -> pk, &m);
    TIMER_BEGIN();
    CH -> FirstPreImage(pair -> sk, m[0], 0, r, m[1], &r2);
    TIMER_END();
    mpz_clear(m2);
  }
  printf("Collision: ");
  TIMER_RESULT();
  mpz_clear(r);
  mpz_clear(m);
  free(string);
  FreeDigest(digest);
  mpz_clear(r2);
  CH -> FreePairOfKeys(pksk);
  free(CH);
}


int main(int argc, char **argv){
  int i;
  PAIR_OF_KEYS pair;
  char *string1, *string2, *string3;
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
  mpz_init(r2);
  struct chameleon_hash_scheme *CH = new_chameleon_hash_scheme();
  pair = CH -> KeyGen(security_parameter);
  string1 = mpz_get_str(NULL, 10, pair -> sk[X]);
  printf("SK = %s\n", string1);
  free(string1);
  string1 = mpz_get_str(NULL, 10, pair -> pk[Y]);
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
  CH -> FirstPreImage(pair -> sk, m[0], 0, r, m[1], &r2);
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
  mpz_clear(r);
  mpz_clear(m);
  mpz_clear(r2);
  mpz_clear(m2);
  CH -> FreePairOfKeys(pair);
  free(CH);
  return 0;
}
