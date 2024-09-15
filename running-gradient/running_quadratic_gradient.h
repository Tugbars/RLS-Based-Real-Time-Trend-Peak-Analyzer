#ifndef RUNNING_QUADRATIC_GRADIENT_H
#define RUNNING_QUADRATIC_GRADIENT_H

#include <stdint.h>  // For uint16_t
#include <math.h>
#include "mqs_def.h"
#include "running_cubic_gradient.h"

/** @brief Define a constant for the RLS window size. */
#define RLS_WINDOW 30

// #define DEBUG_GRADIENT_CALC
// #define DEBUG_QUADRATIC

/**
 * @brief Structure to store the result of quadratic peak analysis.
 */
typedef struct {
    bool peak_found;     /**< Flag indicating if a peak was found and verified. */
    uint16_t peak_index; /**< Index of the detected peak in the data array. */
} QuadraticPeakAnalysisResult;

/**
 * @brief Structure to represent the output of concavity analysis.
 */
typedef struct {
    bool isNoisy;          /**< Indicates if the segment is noisy. */
    bool moveLeft;         /**< Indicates if the analysis suggests moving left. */
    bool moveRight;        /**< Indicates if the analysis suggests moving right. */
    bool stay;             /**< Indicates if the analysis suggests staying in the current position. */
    bool isPotentialPeak;  /**< Indicates if a potential peak is detected. */
    bool isTruePeak;       /**< Indicates if a true peak is detected. */
} ConcavityAnalysisOutput;

/**
 * @brief Tracks the second-order gradient trends within a specified window using quadratic regression.
 *
 * This function calculates second-order gradients (curvatures) at each point within a specified window,
 * identifies regions of consistent increase and decrease, and returns the information about these trends.
 *
 * @param values           Array of data points.
 * @param length           Length of the data array.
 * @param start_index      Starting index in the data array from which to begin calculations.
 * @param window_size      The number of points to include in the gradient calculation.
 * @param forgetting_factor The forgetting factor for the RLS algorithm.
 * @return GradientTrendResult Structure containing information about increasing and decreasing trends.
 */
GradientTrendResult track_gradient_trends_with_quadratic_regression(
    const MqsRawDataPoint_t *values,
    uint16_t length,
    uint16_t start_index,
    uint16_t window_size,
    double forgetting_factor
);

/**
 * @brief Computes the total sum of second-order gradients within a specified window.
 *
 * This function calculates the second-order gradients (curvatures) at each point within a specified window,
 * and computes their total sum. It is useful for assessing the overall concavity of the data in the window.
 *
 * @param values           Array of data points.
 * @param length           Length of the data array.
 * @param start_index      Starting index in the data array from which to begin calculations.
 * @param forgetting_factor The forgetting factor for the RLS algorithm.
 * @return The total sum of second-order gradients within the window.
 */
double compute_total_second_order_gradient(
    const MqsRawDataPoint_t *values,
    uint16_t length,
    uint16_t start_index,
    double forgetting_factor
);

/**
 * @brief Finds and verifies a peak in the data using quadratic RLS.
 *
 * This function first detects a peak using second-order gradients calculated from quadratic regression,
 * and then verifies the peak using trend consistency criteria on both sides of the peak.
 *
 * @param values           Array of data points.
 * @param length           Length of the data array.
 * @param start_index      The start index in the data array from which to begin the gradient calculation.
 * @param forgetting_factor The forgetting factor used in the RLS algorithm.
 * @return QuadraticPeakAnalysisResult Structure containing the peak detection status and the peak index if found and verified.
 */
QuadraticPeakAnalysisResult find_and_verify_quadratic_peak(
    const MqsRawDataPoint_t *values,
    uint16_t length,
    uint16_t start_index,
    double forgetting_factor
);

#endif // RUNNING_QUADRATIC_GRADIENT_H