#define MSG_SIZE 64

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
  RND r;
  MSG msg;
  init_digest(&digest);
  // Keygen
  for(i = 0; i < N; i ++){
    TIMER_BEGIN();
    keygen(security_parameter, &pk, &sk);
    TIMER_END();
  }
  printf("KeyGen: ");
  TIMER_RESULT();
  keygen(security_parameter, &pk, &sk);
  // Hash
  for(i = 0; i < N; i ++){
    random_msg(&msg, MSG_SIZE);
    random_rnd(&pk, &r);
    TIMER_BEGIN();
    hash(&pk, &msg, &r, &digest);
    TIMER_END();
    free_msg(&msg);
  }
  printf("Hash: ");
  TIMER_RESULT();
  // Collision
  for(i = 0; i < N; i ++){
    random_msg(&msg, MSG_SIZE);
    random_digest(&pk, &digest);
    TIMER_BEGIN();
    firstpreimage(&sk, &msg, &digest, &r);
    TIMER_END();
    free_msg(&msg);
  }
  printf("Collision: ");
  TIMER_RESULT();
  free_digest(&digest);
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
  if(argc >= 3 && !strcmp("--benchmark", argv[2])){
    benchmark(security_parameter);
    exit(0);
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

