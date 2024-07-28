#ifndef RUNNING_GRADIENT_H
#define RUNNING_GRADIENT_H

#include <stddef.h>
#include <stdint.h>

// Structure to hold the running gradient calculation data
typedef struct {
    size_t num_points;
    size_t max_points;
    double *x;
    double *y;
    double inverse_cov_matrix[2][2];
    double coefficients[2];
    double residual_sum_squares;
} RunningGradient;

#define NO_PEAK_FOUND UINT16_MAX
#define WINDOW_SIZE 5
#define GRADIENT_INCREASE_THRESHOLD 0.05
#define GRADIENT_DECREASE_THRESHOLD -0.05
#define DECREASE_TREND_LENGTH 3

void init_running_gradient(RunningGradient *rg, size_t window_size);

void add_data_point(RunningGradient *rg, double y);

double calculate_gradient(const RunningGradient *rg);

void free_running_gradient(RunningGradient *rg);

/**
 * @brief Processes the data in chunks of 10 and calculates various statistics for each chunk.
 * 
 * @param data The data points.
 * @param size The number of data points.
 */
void process_windowed_data_in_chunks(const double *data, size_t size);
// I need to spend more time on this. 
uint16_t find_peak(const double *data, size_t start_idx, size_t end_idx, size_t  window_size);

uint16_t find_peak_every_fifth(const double *data, size_t start_idx, size_t end_idx, size_t window_size);

uint16_t find_overall_increasing_trend(const double *data, size_t size, size_t start_idx, size_t window_size);

int determine_trend_direction(const double *data, size_t size, size_t start_idx, size_t window_size, size_t trend_length, double gradient_increase_threshold, double gradient_decrease_threshold);

int determine_trend_direction_robust(const double *data, size_t size, size_t start_idx, size_t window_size, size_t trend_length, 
                              double gradient_increase_threshold, double gradient_decrease_threshold, double quantile_discard, double probability_of_increase);
                              

#endif // RUNNING_GRADIENT_H
