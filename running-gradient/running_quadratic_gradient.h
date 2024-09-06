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

/*
ConcavityPattern analyze_concavity_segments(const ConcavityAnalysisResult *concavity_result, bool *isPotentialPeak, bool *isTruePeak) {
    if (isnan(concavity_result->sums[0]) || isnan(concavity_result->sums[1]) || isnan(concavity_result->sums[2])) {
        printf("Insufficient data for determining concavity pattern.\n");
        *isPotentialPeak = false;
        *isTruePeak = false;
        return UNDETERMINED_PATTERN;
    }

    *isPotentialPeak = false;
    *isTruePeak = false;

    // Define thresholds for categorizing the gradient sums
    double increase_threshold = quadratic_analysis_params.minimum_second_order_gradient_sum;
    double small_increase_threshold = 0.0;

    // Map to 2 (increase), 1 (small increase), or 0 (decrease)
    int first_segment = concavity_result->sums[0] > increase_threshold ? 2 :
                        (concavity_result->sums[0] > small_increase_threshold ? 1 : 0);
    int second_segment = concavity_result->sums[1] > increase_threshold ? 2 :
                         (concavity_result->sums[1] > small_increase_threshold ? 1 : 0);
    int third_segment = concavity_result->sums[2] > increase_threshold ? 2 :
                        (concavity_result->sums[2] > small_increase_threshold ? 1 : 0);
                        
    printf("%d %d %d\n", first_segment, second_segment, third_segment);

    // Create a ternary pattern based on the segments (using 2 as increase, 1 as small increase, 0 as decrease)
    int pattern = (first_segment * 9) + (second_segment * 3) + third_segment;

    // Switch on the ternary pattern to determine the concavity pattern
    switch (pattern) {
        case 25: // 2, 2, 1: INCREASE, INCREASE, SMALL INCREASE //can be peak. go right. go to right to center it. 
        case 24: // 2, 2, 0: INCREASE, INCREASE, DECREASE   //can be peak. should be peak! go to right to center it. 
            *isPotentialPeak = true;
            return INCREASE_INCREASE_DECREASE;
        case 23: // 2, 1, 2: INCREASE, SMALL INCREASE, INCREASE   //very noisy left to curve. go to right. 
        case 22: // 2, 1, 1: INCREASE, SMALL INCREASE, SMALL INCREASE  // noisy curve. potentially can be a weird peak, go to right. 
        case 21: // 2, 1, 0: INCREASE, SMALL INCREASE, DECREASE // on the peak, potentially. go right to center it. 
        case 20: // 2, 0, 2: INCREASE, DECREASE, INCREASE //noisy and probably on the peak. centered. stay. 
            *isPotentialPeak = true;
            return INCREASE_DECREASE_INCREASE;
        case 18: // 2, 0, 0: INCREASE, DECREASE, DECREASE //on the peak. go to left to center it. no noise. 
            *isPotentialPeak = true;
            *isTruePeak = true;
            return INCREASE_DECREASE_DECREASE;   
        case 14: // 1, 2, 2: SMALL INCREASE, INCREASE, INCREASE //close to peak, go right side. no noise. 
            return DECREASE_INCREASE_INCREASE;
        case 11: // 1, 0, 2: SMALL INCREASE, DECREASE, INCREASE //noisy curve, go to left. not on the peak. 
        case 9:  // 0, 2, 2: DECREASE, INCREASE, INCREASE   //right to peak, go left. 
            return DECREASE_INCREASE_INCREASE;
        case 10: // 1, 0, 0: SMALL INCREASE, DECREASE, DECREASE  noisy but on the peak. go to left. 
        case 7:  // 0, 2, 0: DECREASE, INCREASE, DECREASE    //on the peak, noisy. go to left. 
            *isPotentialPeak = true;
            return DECREASE_INCREASE_DECREASE;
        case 4:  // 0, 0, 2: DECREASE, DECREASE, INCREASE  //right sidde of the peak or noise. go left. 
            return DECREASE_DECREASE_INCREASE;
        case 0:  // 0, 0, 0: DECREASE, DECREASE, DECREASE  //right sidee of the peak.  go to left. 
            *isTruePeak = true;
            return DECREASE_DECREASE_DECREASE;   
        case 26: // 2, 2, 2: INCREASE, INCREASE, INCREASE  //left sie of the peak. go right
            return INCREASE_INCREASE_INCREASE;
        default:
            return UNDETERMINED_PATTERN;
    }
}
*/