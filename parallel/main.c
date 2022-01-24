#include <stdio.h>
#include <math.h>
#define _GNU_SOURCE
#include <getopt.h>
#include <stdlib.h>
#include <string.h>

#include "./src/icolor.h"
#include "./src/icolorant.h"
#include "./src/tabucol.h"
#include "./src/util.h"

int main(int argc, char *argv[]) { /*{{{*/

  gcp_solution_t *results;

#if defined NRAND
  unsigned long int seed;
  int x;
#endif

  problem = malloc_(sizeof(gcp_t));
  init_flag(problem->flags);
  problem->nof_vertices = 0;
  problem->nof_edges = 0;
  problem->cycles = 0;
  problem->convergence_cycles = 0;
  problem->time = 0;
  problem->colors = 0;
  problem->flags = 0;
  problem->degree = 0;
  problem->adj_matrix = 0;
  problem->adj_list = 0;
  problem->fileout = stdout;
#if defined LRAND
  problem->seed = 0;
#endif
#if defined NRAND
  for (x = 0; x < 3; x++)
    problem->seed[x] = 0;
#endif

  colorant_malloc();
  tabucol_malloc();

  parseargs(argc, argv);

  if (!(get_flag(problem->flags, FLAG_CYCLE)) &&
      !(get_flag(problem->flags, FLAG_TIME))) {
    printf("You need to set the stop criterion.\n");
    exit(0);
  }

  set_flag(problem->flags, FLAG_S_ASSIGN);
  set_flag(problem->flags, FLAG_ADJ_MATRIX);

  initialization();

  colorant_initialization();
  tabucol_initialization();

  if (!(get_flag(problem->flags, FLAG_COLOR))) {
    problem->colors = problem->nof_vertices;
  }

  if (!(get_flag(problem->flags, FLAG_SEED))) {
#if defined LRAND
    problem->seed = create_seed();
#endif
#if defined NRAND
    seed = create_seed();
    memcpy(problem->seed, &seed, sizeof(unsigned short) * 3);
#endif
    set_flag(problem->flags, FLAG_SEED);
  }

#if defined LRAND
  srand48_r(problem->seed, &problem->buffer);
#endif
#if defined NRAND
  seed48_r(problem->seed, &problem->buffer);
#endif

  printbanner();

  time_initial = current_time_secs(TIME_INITIAL, 0);

  results = find_solution();

  problem->real_time = current_time_secs(TIME_FINAL, time_initial);

  show_solution(results);

  fclose(problem->fileout);

  return 0;

} /*}}}*/