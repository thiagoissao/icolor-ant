
/*************************************************************************
 * Created: Seg 7 Fev 2011 16:36:24 BRST
 *
 * Author: Carla N. Lintzmayer, carla0negri@gmail.com
 *
 *************************************************************************
 *
 *
 *************************************************************************/
#include <stdio.h>
#include <math.h>
#define _GNU_SOURCE
#include <getopt.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>

#include "icolor.h"
#include "./icolorant/icolorant.h"
#include "./icolorant/icolorant.utils.h"
#include "tabucol.h"
#include "util.h"
#include "helpers.h"
#include "../main.h"

static char *namefilein;

void parseargs(int argc, char *argv[]) {

  extern char *optarg;
  char op;

#if defined NRAND
  unsigned long seed;
#endif

  /* Usando getopt para tratamento dos argumentos */
  struct option longopts[] = {
      {"alpha", 1, NULL, 'a'},
      {"beta", 1, NULL, 'b'},
      {"rho", 1, NULL, 'r'},
      {"ants", 1, NULL, 'A'},
      {"use-ants-ratio", 0, NULL, 'R'},
      {"pheromone-scheme", 1, NULL, 'p'},
      {"change-phero-scheme-iterations", 1, NULL, 'n'},
      {"memory-size", 1, NULL, 'm'},
      {"use-memory-mratio", 0, NULL, 'M'},
      {"delta", 1, NULL, 'd'},
      {"gap", 1, NULL, 'g'},
      {"x", 1, NULL, 'x'},
      {"y", 1, NULL, 'y'},
      {"gama", 1, NULL, 'G'},
      {"omega", 1, NULL, 'o'},
      {"iterations-alpha-beta", 1, NULL, 'i'},
      {"tabucol-cycles", 1, NULL, 't'},
      {"tabucol-convergence-cycles", 1, NULL, 'T'},
      {"tabucol-scheme", 1, NULL, 's'},
      {"change-tabucol-scheme-iterations", 1, NULL, 'N'},
      {"diff-tabucol-scheme-iterations", 1, NULL, 'F'},
      {"apply-tabucol-all-ants", 0, NULL, 'u'},
      {"cycles", 1, NULL, 'c'},
      {"time", 1, NULL, 'E'},
      {"convergence-cycles", 1, NULL, 'Y'},
      {"colors", 1, NULL, 'k'},
      {"threads", 1, NULL, 'w'},
      {"verbose", 0, NULL, 'v'},
      {"tabucol-verbose", 0, NULL, 'V'},
      {"seed", 1, NULL, 'S'},
      {"output-filename", 1, NULL, 'f'},
      {"help", 0, NULL, 'h'}};

  while ((op = getopt_long(
              argc, argv,
              "w:a:b:r:A:Rp:n:m:Md:g:x:y:G:o:i:t:T:s:N:F:uc:E:Y:k:vVS:f:h",
              longopts, NULL)) != -1) {

    switch (op) {
    case 'w':
      aco_info->threads = THREADS;
      int arg = atoi(optarg);
      if (arg > 0) {
        aco_info->threads = arg;
      }
      break;
    case 'a':
      aco_info->alpha = atof(optarg);
      if (aco_info->alpha <= 0.0)
        aco_info->alpha = COLORANT_ALPHA;
      break;
    case 'b':
      aco_info->beta = atof(optarg);
      if (aco_info->beta <= 0.0)
        aco_info->beta = COLORANT_BETA;
      break;
    case 'r':
      aco_info->rho = atof(optarg);
      if (aco_info->rho <= 0.0)
        aco_info->rho = COLORANT_RHO;
      break;
    case 'A':
      aco_info->ants = atoi(optarg);
      if (aco_info->ants < 1)
        aco_info->ants = COLORANT_ANTS;
      break;
    case 'R':
      set_flag(problem->flags, FLAG_ANTS_RATIO);
      break;
    case 'p':
      aco_info->pheromone_scheme = atoi(optarg);
      if ((aco_info->pheromone_scheme < PHEROMONE_SCHEME_1) ||
          (aco_info->pheromone_scheme > PHEROMONE_SCHEME_3))
        aco_info->pheromone_scheme = PHEROMONE_SCHEME_1;
      break;
    case 'n':
      aco_info->change_phero_scheme_iterations = atoi(optarg);
      if (aco_info->change_phero_scheme_iterations < 1)
        aco_info->change_phero_scheme_iterations =
            COLORANT_CHANGE_PHERO_SCHEME_ITERATIONS;
      set_flag(problem->flags, FLAG_CHANGE_PHEROMONE_SCHEME);
      break;
    case 'm':
      aco_info->memory_size = atoi(optarg);
      if (aco_info->memory_size < 1)
        aco_info->memory_size = COLORANT_MEMORY_SIZE;
      set_flag(problem->flags, FLAG_MEMORY);
      break;
    case 'M':
      set_flag(problem->flags, FLAG_MEMORY);
      set_flag(problem->flags, FLAG_MEMORY_RATIO);
      break;
    case 'd':
      aco_info->delta = atof(optarg);
      if (aco_info->delta < 0.0)
        aco_info->delta = COLORANT_DELTA;
      set_flag(problem->flags, FLAG_MEMORY);
      break;
    case 'g':
      aco_info->gap = atoi(optarg);
      if (aco_info->gap < 1)
        aco_info->gap = COLORANT_GAP;
      break;

    case 'x':
      aco_info->x = atof(optarg);
      if (aco_info->x < 1.0)
        aco_info->x = COLORANT_X;
      set_flag(problem->flags, FLAG_REUSE_COLOR);
      break;
    case 'y':
      aco_info->y = atof(optarg);
      if (aco_info->y < 1.0)
        aco_info->y = COLORANT_Y;
      set_flag(problem->flags, FLAG_REUSE_COLOR);
      break;
    case 'G':
      aco_info->gamma = atof(optarg);
      if (aco_info->gamma < 1.0)
        aco_info->gamma = COLORANT_GAMMA;
      set_flag(problem->flags, FLAG_CHANGE_ALPHA_BETA);
      break;
    case 'o':
      aco_info->omega = atof(optarg);
      if (aco_info->omega < 1.0)
        aco_info->omega = COLORANT_OMEGA;
      set_flag(problem->flags, FLAG_CHANGE_ALPHA_BETA);
      break;
    case 'i':
      aco_info->iterations_alpha_beta = atoi(optarg);
      if (aco_info->iterations_alpha_beta < 1)
        aco_info->iterations_alpha_beta = COLORANT_ITERATIONS_ALPHA_BETA;
      set_flag(problem->flags, FLAG_CHANGE_ALPHA_BETA);
      break;
    case 't':
      tabucol_info->cycles = atoi(optarg);
      // if (tabucol_info->cycles < 1)
      // tabucol_info->cycles = TABUCOL_CYCLES;
      break;
    case 'T':
      tabucol_info->convergence_cycles = atoi(optarg);
      if (tabucol_info->convergence_cycles < 1)
        tabucol_info->convergence_cycles = TABUCOL_CONVERGENCE_CYCLES;
      set_flag(problem->flags, FLAG_TABUCOL_CONV);
      break;
    case 's':
      tabucol_info->tl_style = atoi(optarg);
      if ((tabucol_info->tl_style != TABUCOL_REACTIVE) &&
          (tabucol_info->tl_style != TABUCOL_DYNAMIC))
        tabucol_info->tl_style = TABUCOL_REACTIVE;
      break;
    case 'N':
      tabucol_info->change_scheme_iterations = atoi(optarg);
      if (tabucol_info->change_scheme_iterations < 1)
        tabucol_info->change_scheme_iterations =
            TABUCOL_CHANGE_SCHEME_ITERATIONS;
      set_flag(problem->flags, FLAG_CHANGE_TABUCOL_SCHEME);
      break;
    case 'F':
      tabucol_info->diff_scheme_iterations = atoi(optarg);
      if (tabucol_info->diff_scheme_iterations < 1)
        tabucol_info->diff_scheme_iterations = TABUCOL_DIFF_SCHEME_ITERATIONS;
      set_flag(problem->flags, FLAG_DIFF_TABUCOL_SCHEME);
      break;
    case 'u':
      set_flag(problem->flags, FLAG_TABUCOL_ALL_ANTS);
      break;
    case 'c':
      problem->cycles = atoi(optarg);
      if (problem->cycles < 1)
        problem->cycles = DEFAULT_CYCLES;
      set_flag(problem->flags, FLAG_CYCLE);
      break;
    case 'E':
      problem->time = atof(optarg);
      if (problem->time < 1.0)
        problem->time = DEFAULT_TIME;
      set_flag(problem->flags, FLAG_TIME);
      break;
    case 'Y':
      problem->convergence_cycles = atoi(optarg);
      if (problem->convergence_cycles < 1)
        problem->convergence_cycles = DEFAULT_CONVERGENCE_CYCLES;
      set_flag(problem->flags, FLAG_CONV);
      break;
    case 'k':
      problem->colors = atoi(optarg);
      set_flag(problem->flags, FLAG_COLOR);
      break;

    case 'v':
      set_flag(problem->flags, FLAG_VERBOSE);
      break;
    case 'V':
      set_flag(problem->flags, FLAG_TABUCOL_VERBOSE);
      break;
    case 'S':
#if defined LRAND
      problem->seed = atol(optarg);
#elif defined NRAND
      seed = atol(optarg);
      memcpy(problem->seed, &seed, sizeof(unsigned short) * 3);
#endif
      set_flag(problem->flags, FLAG_SEED);
      fprintf(stdout,
              "Sem imprimir flags (%i), está gerando semente ao invés de pegar "
              "a passada por parâmetro!!\nSei lá o que está acontecendo!!\nCom "
              "este print funciona, então vai assim!!!!\n",
              problem->flags);
      break;
    case 'f':
      problem->fileout = fopen(optarg, "w");
      break;
    case 'h':
      show_help(argv[0]);
      exit(0);
    }
  }

  /* O único argumento não capturado acima é o nome do arquivo de entrada,
   * se existir */
  if (optind < argc) {
    namefilein = malloc_(sizeof(char) * strlen(argv[optind]) + 1);
    strcpy(namefilein, argv[optind++]);

    /* verificar se foi passado algum argumento a mais */
    if (optind < argc) {
      printf("error: invalid argument. Use '-h'\n");
      exit(0);
    }
  } else {
    printf("error: no input files\n");
    exit(0);
  }
}

void initialization(void) {

  FILE *in;

  int i, j, vi, vj;
  char f, t[50];

  in = fopen(namefilein, "r");
  if (!in) {
    printf("error: no input files\n");
    exit(0);
  }

  /* Ignoring initial informations */
  while ((j = fscanf(in, "%c", &f)) && f != 'p') {
    while (f != '\n') {
      j = fscanf(in, "%c", &f);
    }
  }

  j = fscanf(in, "%s %d %d\n", t, &problem->nof_vertices, &problem->nof_edges);
  problem->degree = malloc_(sizeof(int) * problem->nof_vertices);
  if (get_flag(problem->flags, FLAG_ADJ_MATRIX)) {
    problem->adj_matrix = malloc_(sizeof(int *) * problem->nof_vertices);
  }
  if (get_flag(problem->flags, FLAG_ADJ_LIST)) {
    problem->adj_list = malloc_(sizeof(int *) * problem->nof_vertices);
  }

  for (i = 0; i < problem->nof_vertices; i++) {

    if (get_flag(problem->flags, FLAG_ADJ_MATRIX)) {
      problem->adj_matrix[i] = malloc_(sizeof(int) * problem->nof_vertices);
    }
    if (get_flag(problem->flags, FLAG_ADJ_LIST)) {
      problem->adj_list[i] = malloc_(sizeof(int) * (problem->nof_edges + 1));
    }

    for (j = 0; j < problem->nof_vertices; j++) {

      if (get_flag(problem->flags, FLAG_ADJ_MATRIX)) {
        problem->adj_matrix[i][j] = 0;
      }
      if (get_flag(problem->flags, FLAG_ADJ_LIST)) {
        problem->adj_list[i][j] = 0;
      }
    }

    if (get_flag(problem->flags, FLAG_ADJ_LIST)) {
      problem->adj_list[i][problem->nof_vertices] = 0;
    }

    problem->degree[i] = 0;
  }

  for (i = 0; i < problem->nof_edges; i++) {
    j = fscanf(in, "%c %d %d\n", &f, &vi, &vj);

    if (get_flag(problem->flags, FLAG_ADJ_MATRIX)) {
      problem->adj_matrix[vi - 1][vj - 1] = 1;
      problem->adj_matrix[vj - 1][vi - 1] = 1;
    }

    if (get_flag(problem->flags, FLAG_ADJ_LIST)) {
      problem->adj_list[vi - 1][0]++;
      problem->adj_list[vi - 1][problem->adj_list[vi - 1][0]] = vj - 1;
      problem->adj_list[vj - 1][0]++;
      problem->adj_list[vj - 1][problem->adj_list[vj - 1][0]] = vi - 1;
    }

    problem->degree[vi - 1]++;
    problem->degree[vj - 1]++;
  }

  fclose(in);
}

void cpy_solution(gcp_solution_t *src, gcp_solution_t *dst) {

  int i, j;
  if (get_flag(problem->flags, FLAG_S_ASSIGN)) {
    for (i = 0; i < problem->nof_vertices; i++) {
      dst->color_of[i] = src->color_of[i];
    }
  }
  if (get_flag(problem->flags, FLAG_S_PARTITION)) {
    for (i = 0; i <= problem->nof_vertices; i++) {
      for (j = 0; j < problem->colors; j++)
        dst->class_color[j][i] = src->class_color[j][i];
    }
  }

  dst->spent_time = src->spent_time;
  dst->spent_time_ls = src->spent_time_ls;
  dst->time_to_best = src->time_to_best;
  dst->total_cycles = src->total_cycles;
  dst->cycles_to_best = src->cycles_to_best;
  dst->nof_colors = src->nof_colors;
  dst->nof_confl_edges = src->nof_confl_edges;
  dst->nof_confl_vertices = src->nof_confl_vertices;
  dst->stop_criterion = src->stop_criterion;
  dst->h1 = src->h1;
  dst->f1 = src->f1;
}

gcp_solution_t *init_solution(void) {
  int i;
  gcp_solution_t *solution;

  solution = malloc_(sizeof(gcp_solution_t));

  if (get_flag(problem->flags, FLAG_S_ASSIGN)) {
    solution->color_of = malloc_(sizeof(int) * problem->nof_vertices);
  }
  if (get_flag(problem->flags, FLAG_S_PARTITION)) {
    solution->class_color = malloc_(sizeof(int *) * problem->colors);
    for (i = 0; i < problem->colors; i++) {
      solution->class_color[i] =
          malloc_(sizeof(int) * (problem->nof_vertices + 1));
    }
  }

  solution->nof_colors = 0;
  solution->total_cycles = 0;
  solution->cycles_to_best = 0;
  solution->nof_confl_edges = 0;
  solution->nof_confl_vertices = 0;
  solution->stop_criterion = -1;

  return solution;
}

void find_global_best_ant() {

  ant_t *local_ant = malloc_(sizeof(ant_t *));
  pthread_mutex_lock(&global_best_ant_mutex);

  global_best_ant = execute_colorant(&local_ant);
  pthread_mutex_unlock(&global_best_ant_mutex);
}

int terminate_conditions(gcp_solution_t *solution, int cycle, int converg) {

  if (get_flag(problem->flags, FLAG_CONV) &&
      converg >= problem->convergence_cycles) {
    solution->stop_criterion = STOP_CONV;
    return TRUE;
  } else if (get_flag(problem->flags, FLAG_CYCLE) && cycle >= problem->cycles) {
    solution->stop_criterion = STOP_CYCLES;
    return TRUE;
  } else if (get_flag(problem->flags, FLAG_TIME) &&
             current_time_secs(TIME_FINAL, time_initial) >= problem->time) {
    solution->stop_criterion = STOP_TIME;
    return TRUE;
  }
  return FALSE;
}
