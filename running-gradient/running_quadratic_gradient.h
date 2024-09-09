#ifndef RUNNING_QUADRATIC_GRADIENT_H
#define RUNNING_QUADRATIC_GRADIENT_H

#include <stdint.h>  // For uint16_t
#include <math.h>
#include "mqs_def.h"
#include "running_cubic_gradient.h"

// Constants
#define RLS_WINDOW 30   // Define a constant for the window size

/**
 * @brief Structure to store the result of the concavity analysis.
 *
 * This structure contains the sum of second-order gradients for every 10 points within the RLS window.
 */
typedef struct {
    double sums[RLS_WINDOW / 10]; /**< Array containing the sum of second-order gradients for each 10-point segment. */
} ConcavityAnalysisResult;

/**
 * @brief Enum to represent the different patterns of concavity across the three segments.
 */
typedef struct {
    bool isNoisy;          /**< Indicates if the segment is noisy */
    bool moveLeft;         /**< Indicates if the analysis suggests moving left */
    bool moveRight;        /**< Indicates if the analysis suggests moving right */
    bool stay;             /**< Indicates if the analysis suggests staying in the current position */
    bool isPotentialPeak;  /**< Indicates if a potential peak is detected */
    bool isTruePeak;       /**< Indicates if a true peak is detected */
} ConcavityAnalysisOutput;

double* compute_second_order_gradients(const MqsRawDataPoint_t *values, uint16_t length, uint16_t start_index, double forgetting_factor);

ConcavityAnalysisResult initial_concavity_analysis(const MqsRawDataPoint_t *values, uint16_t length, uint16_t start_index, double forgetting_factor, bool reinitialize_after_each_segment);

ConcavityAnalysisOutput analyze_concavity_segments(const ConcavityAnalysisResult *concavity_result);

GradientTrendResult track_gradient_trends_with_quadratic_regression(const MqsRawDataPoint_t *values, uint16_t length, uint16_t start_index, uint16_t window_size, double forgetting_factor);

#endif // RUNNING_QUADRATIC_GRADIENT_H

