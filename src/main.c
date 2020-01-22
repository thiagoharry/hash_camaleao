#define MSG_SIZE 32

/********************* TIMER **************************************/
#ifndef __timer_h_
#define __timer_h_
#include <sys/time.h>
#include <math.h>
#define N 2
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
      _dif_squared += (measures[_i] - mean) * (measures[_i] - mean);	\
    printf("Mean: %.6fs ± %.6fs\n", 0.000001 * mean,			\
	   0.000001 * (sqrt(((double) _dif_squared) / (double) (N-1)))); \
    _i = t_sum = 0;							\
  }
#endif
/********************* TIMER **************************************/


void benchmark(int security_parameter){
  int i;
  PK pk;
  SK sk;
  DIGEST digest;
  RND r1, r2;
  MSG msg1, msg2;
  init_digest(&digest);
  init_rnd(&r1);
  init_rnd(&r2);
  // Keygen
  for(i = 0; i < N; i ++){
    TIMER_BEGIN();
    keygen(security_parameter, &pk, &sk);
    TIMER_END();
    free_keys(&pk, &sk);
  }
  printf("KeyGen: ");
  TIMER_RESULT();
  keygen(security_parameter, &pk, &sk);
  // Hash
  for(i = 0; i < N; i ++){
    random_msg(&pk, &msg1, MSG_SIZE);
    random_rnd(&pk, &r1);
    TIMER_BEGIN();
    hash(&pk, &msg1, &r1, &digest);
    TIMER_END();
    free_msg(&msg1);
  }
  printf("Hash: ");
  TIMER_RESULT();
  // Collision
  for(i = 0; i < N; i ++){
    if(has_first_pre_image){
      random_msg(&pk, &msg1, MSG_SIZE);
      random_digest(&pk, &digest);
      TIMER_BEGIN();
      firstpreimage(&sk, &msg1, &digest, &r1);
      TIMER_END();
      free_msg(&msg1);
    }
    else{
      random_msg(&pk, &msg1, MSG_SIZE);
      random_msg(&pk, &msg2, MSG_SIZE);
      random_rnd(&pk, &r1);
      TIMER_BEGIN();
      collision(&sk, &msg1, &r1, &msg2, &r2);
      TIMER_END();
      free_msg(&msg1);
      free_msg(&msg2);
    }
  }
  printf("Collision: ");
  TIMER_RESULT();
  free_keys(&pk, &sk);
  free_digest(&digest);
  free_rnd(&r1);
  free_rnd(&r2);
}

int main(int argc, char **argv){
  MSG msg1, msg2;
  PK pk;
  SK sk;
  DIGEST digest;
  RND r, r2;
  int security_parameter;
  if(argc >= 2)
    security_parameter = atoi(argv[1]);
  if(argc < 2 || security_parameter == 0){
    fprintf(stderr, "Usage: chamhash SECURITY_PARAMETER [--benchmark]\n");
    exit(1);
  }
  if(argc >= 3 && !strcmp("--benchmark", argv[2])){
    benchmark(security_parameter);
    exit(0);
  }
  init_digest(&digest);
  init_rnd(&r);
  init_rnd(&r2);
  keygen(security_parameter, &pk, &sk);
  random_msg(&pk, &msg1, MSG_SIZE);
  random_msg(&pk, &msg2, MSG_SIZE);
  print_keys(&pk, &sk);
  random_rnd(&pk, &r);
  hash(&pk, &msg1, &r, &digest);
  print_hash(&msg1, &r, &digest);
  if(has_first_pre_image)
    firstpreimage(&sk, &msg2, &digest, &r2);
  else
    collision(&sk, &msg1, &r, &msg2, &r2);
  hash(&pk, &msg2, &r2, &digest);
  print_hash(&msg2, &r2, &digest);
  free_msg(&msg1);
  free_msg(&msg2);
  free_digest(&digest);
  free_rnd(&r);
  free_rnd(&r2);
  free_keys(&pk, &sk);
  return 0;
}

