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
