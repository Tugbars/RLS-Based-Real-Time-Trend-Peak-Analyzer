#ifndef RUNNING_GRADIENT_PARAMETERS_H
#define RUNNING_GRADIENT_PARAMETERS_H

#include <stddef.h>

typedef struct {
    double gradient_increase_threshold;
    double gradient_decrease_threshold;
    size_t trend_length;
    double quantile_discard;
    size_t chunk_size;
} RunningGradientParameters;

extern RunningGradientParameters g_running_gradient_params;

void init_running_gradient_parameters(double gradient_increase_threshold, double gradient_decrease_threshold, size_t trend_length, double quantile_discard, size_t chunk_size);
void update_running_gradient_parameters(double gradient_increase_threshold, double gradient_decrease_threshold, size_t trend_length, double quantile_discard, size_t chunk_size);

#endif // RUNNING_GRADIENT_PARAMETERS_H
