/*
  Monte Carlo Simulation of a 2D Hard Sphere Gas
  ----------------------------------------------
  This program simulates the dynamics of hard colloids and depletants on a 2D grid
  using a basic Monte Carlo algorithm with periodic boundary conditions.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

// Grid dimensions and number of Monte Carlo steps
#define L 3.0
#define H 2.0
#define epochs 50000000

#ifndef step_multiplier
#define step_multiplier 1
#endif

// Utility functions
double absolute(double a) {
  return (a < 0) ? -a : a;
}

void warning(void* a) {
  if (a == NULL) {
    printf("Allocation failed. Exiting...\n");
    exit(3);
  }
}

double random_uniform() {
  return rand() / (RAND_MAX + 1.0);
}

double modulo_float(double a, double limit) {
  if (a > limit) return a - limit;
  else if (a < 0) return a + limit;
  else return a;
}

struct coord {
  double x;
  double y;
};

double compute_distance(struct coord a, struct coord b) {
  double dx = a.x - b.x;
  double dy = a.y - b.y;

  if (dx > L / 2) dx -= L;
  else if (dx < -L / 2) dx += L;

  if (dy > H / 2) dy -= H;
  else if (dy < -H / 2) dy += H;

  return sqrt(dx * dx + dy * dy);
}

void print_counters(int *counter) {
  printf("Depletant steps rejected due to colloid overlap: %d\n", counter[0]);
  printf("Depletant steps accepted: %d\n", counter[1]);
  printf("Colloid steps rejected due to exceeding max distance: %d\n", counter[2]);
  printf("Colloid steps rejected due to colloid overlap: %d\n", counter[3]);
  printf("Colloid steps rejected due to depletant overlap: %d\n", counter[4]);
  printf("Colloid steps accepted: %d\n\n", counter[5]);
}

void initialize_particles(struct coord *colloids, struct coord *depletants, double sigmac, double sigmap, double max_dist, int N) {
  colloids[0].x = sigmac / 2;
  colloids[0].y = H / 2;

  colloids[1].x = random_uniform() * L;
  colloids[1].y = H / 2;

  while (fabs(colloids[0].x - colloids[1].x) < sigmac || fabs(colloids[0].x - colloids[1].x) > max_dist) {
    colloids[1].x = random_uniform() * L;
  }

  for (int i = 0; i < N; i++) {
    int overlap = 1;

    while (overlap) {
      depletants[i].x = random_uniform() * L;
      depletants[i].y = random_uniform() * H;

      overlap = 0;
      for (int j = 0; j < 2; j++) {
        if (compute_distance(depletants[i], colloids[j]) < (sigmac + sigmap) / 2) {
          overlap = 1;
          break;
        }
      }
    }
  }
}

double* move_particles(struct coord *colloids, struct coord *depletants, double sigmac, double sigmap, double delta, double max_dist, int transient, int N) {
  int i, j, k, selector2, overlap, *counter;
  double selector1, reldist[2];
  struct coord temp_depletant, temp_colloid;

  double *distances = (double *)calloc(epochs - transient, sizeof(double));
  warning(distances);

  counter = (int *)calloc(6, sizeof(int));
  warning(counter);

  double prob = (double)N / (N + 1.0);

  for (i = 0; i < epochs; i++) {
    if (i >= transient) {
      distances[i - transient] = fabs(colloids[1].x - colloids[0].x);
    }

    for (j = 0; j < (N + 1) * step_multiplier; j++) {
      selector1 = random_uniform();

      if (selector1 < prob) {
        selector2 = (int)(random_uniform() * N);

        temp_depletant.x = depletants[selector2].x + (random_uniform() - 0.5) * delta;
        temp_depletant.y = depletants[selector2].y + (random_uniform() - 0.5) * delta;

        temp_depletant.x = modulo_float(temp_depletant.x, L);
        temp_depletant.y = modulo_float(temp_depletant.y, H);

        overlap = 0;
        for (k = 0; k < 2; k++) {
          reldist[k] = compute_distance(temp_depletant, colloids[k]);
          if (reldist[k] < (sigmac + sigmap) / 2) {
            overlap = 1;
            break;
          }
        }

        if (overlap)
          counter[0]++;
        else {
          counter[1]++;
          depletants[selector2] = temp_depletant;
        }
      } else {
        temp_colloid.x = colloids[1].x + (random_uniform() - 0.5) * delta;
        temp_colloid.y = colloids[1].y;
        temp_colloid.x = modulo_float(temp_colloid.x, L);

        overlap = 0;
        reldist[0] = fabs(temp_colloid.x - colloids[0].x);

        if (reldist[0] > max_dist) {
          counter[2]++;
          overlap = 1;
        }

        if (reldist[0] < sigmac) {
          counter[3]++;
          overlap = 1;
        }

        for (k = 0; k < N; k++) {
          reldist[1] = compute_distance(temp_colloid, depletants[k]);
          if (reldist[1] < (sigmac + sigmap) / 2) {
            counter[4]++;
            overlap = 1;
            break;
          }
        }

        if (!overlap) {
          counter[5]++;
          colloids[1] = temp_colloid;
        }
      }
    }

    for (j = 0; j < N; j++) {
      if (depletants[j].x < 0 || depletants[j].y < 0 || depletants[j].x > L || depletants[j].y > H) {
        printf("Error: a depletant escaped the lattice.\nPosition: %.4f %.4f\n", depletants[j].x, depletants[j].y);
        exit(1);
      }
    }

    if (fabs(colloids[1].x - colloids[0].x) > max_dist || fabs(colloids[1].x - colloids[0].x) < sigmac) {
      printf("Error: invalid colloid displacement.\n");
      exit(2);
    }
  }

  print_counters(counter);
  return distances;
}

void print_results(double *r, double sigmac, double sigmap, double max_dist, int transient, int N, FILE *fp1, FILE *fp2) {
  fprintf(fp2, "%.10f %.10f %d %.10lf %.10lf %.10lf\n", L, H, N, sigmac, sigmap, max_dist);
  for (int i = 0; i < epochs - transient; i++) {
    fprintf(fp1, "%d %.10lf\n", i, r[i]);
  }
}

int main() {
  clock_t begin = clock();
  double sigmap = 0.1, sigmac = 1.0, max_dist, delta;
  int transient = (int)(0.2 * epochs);
  int N[10];
  double *r;
  struct coord *colloids, *depletants;
  FILE *fp1, *fp2;
  char filename[50];

  for (int i = 0; i < 10; i++) N[i] = 50;

  delta = sigmac * 0.1;
  max_dist = sigmac + 2 * sigmap;

  fp2 = fopen("SimulationParameters_v3.dat", "w");

  for (int i = 0; i < 10; i++) {
    srand(time(NULL));

    colloids = (struct coord *)calloc(2, sizeof(struct coord));
    warning(colloids);
    depletants = (struct coord *)calloc(N[i], sizeof(struct coord));
    warning(depletants);
    r = (double *)calloc(epochs, sizeof(double));

    initialize_particles(colloids, depletants, sigmac, sigmap, max_dist, N[i]);
    r = move_particles(colloids, depletants, sigmac, sigmap, delta, max_dist, transient, N[i]);

    sprintf(filename, "Simulation%02d_v3.dat", i);
    fp1 = fopen(filename, "w");
    print_results(r, sigmac, sigmap, max_dist, transient, N[i], fp1, fp2);
    fclose(fp1);

    free(colloids);
    free(depletants);
    free(r);
  }

  fclose(fp2);
  clock_t end = clock();

  printf("Step multiplier: %d\n", step_multiplier);
  printf("Execution time: %.4f seconds\n", (double)(end - begin) / CLOCKS_PER_SEC);

  return 0;
}
