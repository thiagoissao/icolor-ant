#include <stdio.h>
#include "./icolorant/icolorant.h"
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
