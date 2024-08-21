#ifndef WINDOWED_RUNNING_GRADIENT_H
#define WINDOWED_RUNNING_GRADIENT_H

#include <stddef.h>
#include <stdint.h>
#include <stdbool.h>

#define WINDOW_SIZE 10
#define GRADIENT_INCREASE_THRESHOLD 0.05
#define GRADIENT_DECREASE_THRESHOLD -0.05
#define DECREASE_TREND_LENGTH 3
#define NO_PEAK_FOUND UINT16_MAX

#define FORGETTING_FACTOR 0.5
#define MAX_SMALL_DECREASE_TREND 3
#define PEAK_IDENTIFIER_GRADIENT 0.8

typedef struct {
    size_t num_points;
    size_t max_points;
    double coefficients[2];
    double residual_sum_squares;
    double inverse_cov_matrix[2][2];
    double x[WINDOW_SIZE];
    double y[WINDOW_SIZE];
    double forgetting_factor;
} RunningGradient;

typedef struct {
    double prob_increase;
    size_t increase_count;
    bool potential_peak;
} ChunkAnalysisResult;

void init_running_gradient(RunningGradient *rg);
void add_data_point(RunningGradient * const rg, const double y);
double calculate_gradient(const RunningGradient * const rg);
double find_upper_quantile(const double * const data, const size_t size, const double quantile);

ChunkAnalysisResult process_chunk_left(const double *data, size_t chunk_size, RunningGradient *rg, size_t actual_start_idx);
ChunkAnalysisResult process_chunk_right(const double *data, size_t chunk_size, RunningGradient *rg, size_t actual_start_idx);

double calculate_probability(RunningGradient *rg);

#endif // WINDOWED_RUNNING_GRADIENT_H
