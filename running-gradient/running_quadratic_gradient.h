#ifndef RUNNING_QUADRATIC_GRADIENT_H
#define RUNNING_QUADRATIC_GRADIENT_H

#include <stddef.h>  // For size_t
#include <math.h>

// Constants
#define RLS_WINDOW 30   // Define a constant for the window size

#define GRADIENT_THRESHOLD 1e-6  // Define a threshold to determine if a gradient is close to zero

typedef enum {
    INDEX_LEFT_TO_PEAK,
    INDEX_RIGHT_TO_PEAK,
    INDEX_AT_PEAK,
    INDEX_POSITION_UNDETERMINED
} IndexRelativeToPeak;

/**
 * @brief Structure to store the result of the peak detection.
 */
typedef struct {
    bool peak_found; /**< Indicates whether a peak was found. */
    size_t peak_index; /**< The index of the detected peak. */
} PeakAnalysisResult;

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
typedef enum {
    INCREASE_INCREASE_DECREASE,
    INCREASE_DECREASE_INCREASE,
    INCREASE_DECREASE_DECREASE,
    DECREASE_INCREASE_INCREASE,
    DECREASE_INCREASE_DECREASE,
    DECREASE_DECREASE_INCREASE,
    DECREASE_DECREASE_DECREASE,
    UNDETERMINED_PATTERN
} ConcavityPattern;


double* compute_second_order_gradients(const double *values, size_t length, size_t start_index, double forgetting_factor);

double* compute_first_order_gradients(const double *values, size_t length, size_t start_index, double forgetting_factor);

IndexRelativeToPeak determine_position_relative_to_peak(const double *values, size_t length, size_t start_index, size_t target_index, double forgetting_factor);

PeakAnalysisResult find_and_verify_peak(const double *values, size_t length, size_t start_index, double forgetting_factor);

ConcavityAnalysisResult initial_concavity_analysis(const double *values, size_t length, size_t start_index, double forgetting_factor, bool reinitialize_after_each_segment);

ConcavityPattern analyze_concavity_segments(const ConcavityAnalysisResult *concavity_result, bool *isPotentialPeak, bool *isTruePeak);

bool find_true_peak_with_concavity_analysis(const double *values, size_t length, size_t start_index, double forgetting_factor, bool reinitialize_after_each_segment, size_t *start_idx, size_t *end_idx);


#endif // RUNNING_QUADRATIC_GRADIENT_H
