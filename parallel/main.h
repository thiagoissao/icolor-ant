/**********************************************************
 * Created: Wed 23 Fev 2022 06:45:20 BRT
 *
 * Author: Thiago I. Yasunaka, thiagoyasunaka@hotmail.com
 *
 *
 **********************************************************/
#include <pthread.h>
#include "./src/icolorant/icolorant.h"

#ifndef __MAIN_H
#define __MAIN_H

gcp_solution_t *global_best_ant;
pthread_mutex_t global_best_ant_mutex;
pthread_t *workers;

#endif /* __MAIN_H */
