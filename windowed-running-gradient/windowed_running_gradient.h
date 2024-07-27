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
    double coefficients[2];
    double residual_sum_squares;
} RunningGradient;


/**
 * @brief Processes the data in chunks of 10 and calculates various statistics for each chunk.
 * 
 * @param data The data points.
 * @param size The number of data points.
 */
void process_windowed_data_in_chunks(const double *data, size_t size);
// I need to spend more time on this. 
//uint16_t find_peak(const double *data, size_t start_idx, size_t end_idx);


#endif // RUNNING_GRADIENT_H