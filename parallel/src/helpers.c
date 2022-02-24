#include <stdio.h>
#include "./icolorant/icolorant.h"
#include "./icolorant/icolorant.utils.h"
#include "tabucol.h"

void show_help(char *nameprog) {
  printf("Usage: %s [options] <file>\n\n", nameprog);
  printf("ColorAnt options\n");
  printf("\t[ -a, --alpha                            ] <value>\tDefine alpha "
         "parameter as <value>. Default: %d.\n",
         COLORANT_ALPHA);
  printf("\t[ -b, --beta                             ] <value>\tDefine beta "
         "parameter as <value>. Default: %d.\n",
         COLORANT_BETA);
  printf("\t[ -r, --rho                              ] <value>\tDefine rho "
         "parameter as <value>. Default: %f.\n",
         COLORANT_RHO);
  printf("\t[ -A, --ants                             ] <value>\tDefine number "
         "of ants as <value>. Default: %d.\n",
         COLORANT_ANTS);
  printf("\t[ -R, --use_ants_ratio                   ] \tDefine number of ants "
         "as ratio of vertices. Default: FALSE.\n");
  printf("\t[ -p, --pheromone-scheme                 ] <value>\tDefine the "
         "pheromone scheme. Default: %i.\n\t\t\t\t\t\t\t1: All ants + Best ant "
         "+ Best colony.\n\t\t\t\t\t\t\t2: Best ant + Best "
         "colony\n\t\t\t\t\t\t\t3: Best ant + Best colony (gap).\n",
         COLORANT_PHEROMONE_SCHEME);
  printf("\t[ -n, --change-phero-scheme-iterations   ] <value>\tDefine "
         "iterations as <value>. Default: %d without improvement.\n",
         COLORANT_CHANGE_PHERO_SCHEME_ITERATIONS);

  printf("   Memory Usage:\n\t[ -m, --memory-size                      ] "
         "<value>\tDefine memory size as <value>. Default: %d.\n",
         COLORANT_MEMORY_SIZE);
  printf("\t[ -M, --use-memory-ratio                 ] \tDefine memory size as "
         "ratio of ants. Default: FALSE.\n");
  printf("\t[ -d, --delta                            ] <value>\tDefine delta "
         "parameter as <value>. Default: %f.\n",
         COLORANT_DELTA);
  printf("   Pheromone Scheme 3:\n\t[ -g, --gap                              ] "
         "<value>\tDefine gap parameter as <value>. Default: %d.\n",
         COLORANT_GAP);
  printf("   Try to reuse color:\n\t[ -x, --x                                ] "
         "<value>\tDefine x parameter as <value>. Default: %f.\n",
         COLORANT_X);
  printf("\t[ -y, --y                                ] <value>\tDefine y "
         "parameter as <value>. Default: %f.\n",
         COLORANT_Y);

  printf("   Change alpha and beta:\n\t[ -G, --gamma                           "
         " ] <value>\tDefine gamma parameter as <value>. Default: %f.\n",
         COLORANT_GAMMA);
  printf("\t[ -o, --omega                            ] <value>\tDefine omega "
         "parameter as <value>. Default: %f.\n",
         COLORANT_OMEGA);
  printf("\t[ -i, --iterations-alpha-beta            ] <value>\tDefine "
         "iterations parameter as <value>. Default: %i.\n",
         COLORANT_ITERATIONS_ALPHA_BETA);

  printf("\nTabucol options\n");
  printf("\t[ -t, --tabucol-cycles                   ] <value>\tDefine maximum "
         "number of cycles in local search as <value>. Default: %d.\n",
         TABUCOL_CYCLES);
  printf("\t[ -T, --tabucol-convergence-cycles       ] <value>\tDefine maximum "
         "number of local search cycles without improvement.Default: %d.\n",
         TABUCOL_CONVERGENCE_CYCLES);
  printf("\t[ -s, --tabucol-scheme                   ] \t\tDefine the tabucol "
         "scheme for tabu tenure (0 - reactive, 1 - dynamic). Default: "
         "reactive scheme.\n");
  printf("\t[ -N, --change-tabucol-scheme-iterations ] <value>\tDefine "
         "iterations as <value>. Default: %d.\n",
         TABUCOL_CHANGE_SCHEME_ITERATIONS);
  printf("\t[ -F, --diff-tabucol-scheme-iterations   ] <value>\tDefine "
         "iterations as <value>. Default: %d.\n",
         TABUCOL_DIFF_SCHEME_ITERATIONS);

  printf("\t[ -u, --apply-tabucol-all-ants                 ] \tApply tabucol "
         "on all ants. Default: only on the best ant.\n");

  printf("\nCriterion Stopping\n");
  printf("\t[ -c, --cycles                           ] <value>\tDefine number "
         "of cycles.\n");
  printf("\t[ -E, --time                             ] <value> \tDefine time "
         "in seconds.\n");
  printf("\t[ -Y, --convergence-cycles               ] <value> \tDefine "
         "maximum number of local search cycles without improvement. \n");

  printf("\nGeneral options\n");
  printf("\t[ -k, --colors                           ] <value>\tDefine promote "
         "iterations as <value>. Default: vertices.\n");
  printf("\t[ -v, --verbose                          ] \t\tDisplay "
         "informations during execution.\n");
  printf("\t[ -V, --tabucol-verbose                  ] \t\tDisplay "
         "informations during execution about local search.\n");
  printf("\t[ -S, --seed                             ] <value>\tDefine <value> "
         "as the seed of rand function. Default: time\n");
  printf("\t[ -f, --output-filename                  ] <value>\tDefine the "
         "output filename. Default: stdout.\n");
  printf("\t[ -w, --threads                          ] \t\t Define number "
         "of threads.\n\n");
  printf("\t[ -h, --help                             ] \t\tDisplay this "
         "information.\n\n");
}

void test_map(gcp_solution_t *solution) {
  int i, j, n;
  int confs = 0;
  for (i = 0; i < problem->nof_vertices; i++) {
    // printf("color of %d: %d\n", i+1, solution->color_of[i]);
    if (get_flag(problem->flags, FLAG_ADJ_MATRIX)) {
      for (j = i; j < problem->nof_vertices; j++) {
        if (problem->adj_matrix[i][j] &&
            solution->color_of[i] == solution->color_of[j]) {
          //	printf("ERROR!! Conflicting edge %d--%d \n", i+1, j+1);
          confs++;
        }
      }
    } else {
      for (j = 1; j <= problem->adj_list[i][0]; j++) {
        n = problem->adj_list[i][j];
        if (solution->color_of[i] == solution->color_of[n]) {
          //	printf("ERROR!! Conflicting edge %d--%d \n", i+1, n+1);
          confs++;
        }
      }
    }
  }
  if (confs != solution->nof_confl_edges) {
    fprintf(problem->fileout, "ERROR!! Confl edges = %d; Calculated = %d\n",
            confs, solution->nof_confl_edges);
  }
}

void show_solution(gcp_solution_t *solution) {
  fprintf(problem->fileout,
          "\n-------------------------------------------------\n");
  fprintf(problem->fileout, "SOLUTION:\n");
  fprintf(problem->fileout,
          "-------------------------------------------------\n");
  fprintf(problem->fileout, "No. of colors utilized.............: %d\n",
          solution->nof_colors);
  fprintf(problem->fileout, "No. of conflicting edges...........: %d\n",
          solution->nof_confl_edges);
  fprintf(problem->fileout, "No. of conflicting vertices........: %d\n",
          solution->nof_confl_vertices);
  fprintf(problem->fileout, "Real time..........................: %lf\n",
          problem->real_time);
  fprintf(problem->fileout, "Spent Time.........................: %lf\n",
          solution->spent_time);

  fprintf(problem->fileout, "Spent time (ACO)...................: %lf\n",
          solution->spent_time - tabucol_info->spent_time);
  fprintf(problem->fileout, "Spent time (LS)....................: %lf\n",
          tabucol_info->spent_time);

  fprintf(problem->fileout, "Time to the best...................: %lf\n",
          solution->time_to_best);
  fprintf(problem->fileout, "Total of cycles....................: %d\n",
          solution->total_cycles);
  fprintf(problem->fileout, "Cycles to the best.................: %d\n",
          solution->cycles_to_best);
  fprintf(problem->fileout, "Stop criterion.....................: %d\n",
          solution->stop_criterion);

  colorant_show_solution();

  fprintf(problem->fileout,
          "-------------------------------------------------\n");
  test_map(solution);
}

void printbanner(void) {

  fprintf(problem->fileout,
          "-------------------------------------------------\n");
  colorant_printbanner();
  fprintf(problem->fileout,
          "-------------------------------------------------\n");
  tabucol_printbanner();

  fprintf(problem->fileout,
          "-------------------------------------------------\n");
  fprintf(problem->fileout, "GENERAL Options\n");
  fprintf(problem->fileout,
          "-------------------------------------------------\n");
  fprintf(problem->fileout, "  K......................................: %i\n",
          problem->colors);
  fprintf(problem->fileout, "  No. of threads ........................: %d\n",
          aco_info->threads);
  if (get_flag(problem->flags, FLAG_CYCLE))
    fprintf(problem->fileout, "  Cycles.................................: %d\n",
            problem->cycles);
  if (get_flag(problem->flags, FLAG_TIME))
    fprintf(problem->fileout,
            "  Time...................................: %.2lf\n",
            problem->time);

#if defined LRAND
  fprintf(problem->fileout,
          "  Seed...................................: %lu (lrand)\n",
          problem->seed);
#elif defined NRAND
  fprintf(problem->fileout,
          "  Seed...................................: %lu (nrand)\n",
          print_seed(problem->seed));
#endif
  if (get_flag(problem->flags, FLAG_CONV))
    fprintf(problem->fileout, "  Stop after %d cycles without improvement.\n",
            problem->convergence_cycles);

  if (problem->flags & FLAG_VERBOSE)
    fprintf(problem->fileout, "  Running on verbose mode.\n");
  if (problem->flags & FLAG_TABUCOL_VERBOSE)
    fprintf(problem->fileout, "  Running Tabu search on verbose mode.\n");
  fprintf(problem->fileout,
          "-------------------------------------------------\n");
}