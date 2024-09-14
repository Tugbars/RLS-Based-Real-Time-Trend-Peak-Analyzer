#ifndef RUNNING_QUADRATIC_GRADIENT_H
#define RUNNING_QUADRATIC_GRADIENT_H

#include <stdint.h>  // For uint16_t
#include <math.h>
#include "mqs_def.h"
#include "running_cubic_gradient.h"

// Constants
#define RLS_WINDOW 30   // Define a constant for the window size
#define MINIMUM_REQUIRED_TREND_COUNT 5
#define ALLOWABLE_INCONSISTENCY_COUNT 2  // Allowable number of inconsistencies

typedef struct {
    bool peak_found;
    uint16_t peak_index;
} QuadraticPeakAnalysisResult;

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


GradientTrendResult track_gradient_trends_with_quadratic_regression(const MqsRawDataPoint_t *values, uint16_t length, uint16_t start_index, uint16_t window_size, double forgetting_factor);

double compute_total_second_order_gradient(const MqsRawDataPoint_t *values, uint16_t length, uint16_t start_index, double forgetting_factor);

/**
 * @brief Finds and verifies a peak in the data using quadratic RLS.
 */
QuadraticPeakAnalysisResult find_and_verify_quadratic_peak(const MqsRawDataPoint_t *values, uint16_t length, uint16_t start_index, double forgetting_factor);

#endif // RUNNING_QUADRATIC_GRADIENT_H

