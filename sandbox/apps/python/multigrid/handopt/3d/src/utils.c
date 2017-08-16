#include "utils.h"

#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>

void print_msg(char* c) {
  printf( "\n %s\n", c );
  exit(1);
}

double get_time() {
  struct timeval t;
  gettimeofday(&t, NULL);
  return (double)t.tv_sec + (double)t.tv_usec * 1.e-6;
}
