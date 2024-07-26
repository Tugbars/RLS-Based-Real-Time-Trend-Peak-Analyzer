#ifndef RUNNING_GRADIENT_H
#define RUNNING_GRADIENT_H

#include <stddef.h>

typedef struct {
    size_t num_points;
    double inverse_cov_matrix[2][2];
    double coefficients[2];
    double residual_sum_squares;
} RunningGradient;

void init_running_gradient(RunningGradient *rg);
void add_data_point(RunningGradient *rg, double y);
double calculate_gradient(const RunningGradient *rg);
double calculate_standard_error(const RunningGradient *rg);
double probability_gradient_less_than(const RunningGradient *rg, double thresh);
double probability_gradient_greater_than(const RunningGradient *rg, double thresh);
double find_upper_quantile(const double* data, size_t size, double quantile);
void process_data_in_chunks(const double *data, size_t size);

#endif // RUNNING_GRADIENT_H
