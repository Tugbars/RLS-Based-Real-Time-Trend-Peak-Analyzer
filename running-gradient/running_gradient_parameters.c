#include "running_gradient_parameters.h"
#include <stdio.h>
#include <stdlib.h>

// Define the global parameters object
RunningGradientParameters g_running_gradient_params;

/**
 * @brief Initializes the global RunningGradientParameters structure with the given values.
 *
 * This function sets the initial values for the global RunningGradientParameters structure and performs basic validation
 * to ensure that the parameters are within acceptable ranges.
 *
 * @param gradient_increase_threshold The threshold above which a gradient is considered an increase.
 * @param gradient_decrease_threshold The threshold below which a gradient is considered a decrease.
 * @param trend_length The minimum number of consecutive steps with a noticeable increase required to determine a trend.
 * @param quantile_discard The upper quantile of data points to discard.
 * @param chunk_size The number of data points in each chunk.
 */
void init_running_gradient_parameters(double gradient_increase_threshold, double gradient_decrease_threshold, size_t trend_length, double quantile_discard, size_t chunk_size) {
    // Basic validation
    if (gradient_increase_threshold <= 0.0) {
        fprintf(stderr, "Error: Gradient thresholds must be positive.\n");
        exit(EXIT_FAILURE);
    }
    if (trend_length <= 0) {
        fprintf(stderr, "Error: Trend length must be greater than zero.\n");
        exit(EXIT_FAILURE);
    }
    if (quantile_discard < 0.0 || quantile_discard > 1.0) {
        fprintf(stderr, "Error: Quantile discard must be between 0 and 1.\n");
        exit(EXIT_FAILURE);
    }
    if (chunk_size <= 0) {
        fprintf(stderr, "Error: Chunk size must be greater than zero.\n");
        exit(EXIT_FAILURE);
    }

    // Initialize the global parameters
    g_running_gradient_params.gradient_increase_threshold = gradient_increase_threshold;
    g_running_gradient_params.gradient_decrease_threshold = gradient_decrease_threshold;
    g_running_gradient_params.trend_length = trend_length;
    g_running_gradient_params.quantile_discard = quantile_discard;
    g_running_gradient_params.chunk_size = chunk_size;
}

/**
 * @brief Updates the global RunningGradientParameters structure with new values.
 *
 * This function updates the values of the global RunningGradientParameters structure and performs basic validation
 * to ensure that the new parameters are within acceptable ranges.
 *
 * @param gradient_increase_threshold The new threshold above which a gradient is considered an increase.
 * @param gradient_decrease_threshold The new threshold below which a gradient is considered a decrease.
 * @param trend_length The new minimum number of consecutive steps with a noticeable increase required to determine a trend.
 * @param quantile_discard The new upper quantile of data points to discard.
 * @param chunk_size The new number of data points in each chunk.
 */
void update_running_gradient_parameters(double gradient_increase_threshold, double gradient_decrease_threshold, size_t trend_length, double quantile_discard, size_t chunk_size) {
    // Basic validation
    if (gradient_increase_threshold <= 0.0 || gradient_decrease_threshold <= 0.0) {
        fprintf(stderr, "Error: Gradient thresholds must be positive.\n");
        exit(EXIT_FAILURE);
    }
    if (trend_length <= 0) {
        fprintf(stderr, "Error: Trend length must be greater than zero.\n");
        exit(EXIT_FAILURE);
    }
    if (quantile_discard < 0.0 || quantile_discard > 1.0) {
        fprintf(stderr, "Error: Quantile discard must be between 0 and 1.\n");
        exit(EXIT_FAILURE);
    }
    if (chunk_size <= 0) {
        fprintf(stderr, "Error: Chunk size must be greater than zero.\n");
        exit(EXIT_FAILURE);
    }

    // Update the global parameters
    g_running_gradient_params.gradient_increase_threshold = gradient_increase_threshold;
    g_running_gradient_params.gradient_decrease_threshold = gradient_decrease_threshold;
    g_running_gradient_params.trend_length = trend_length;
    g_running_gradient_params.quantile_discard = quantile_discard;
    g_running_gradient_params.chunk_size = chunk_size;
}
