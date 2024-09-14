#ifndef SLIDING_WINDOW_ANALYSIS_H
#define SLIDING_WINDOW_ANALYSIS_H

#include <stdint.h>
#include <stdbool.h>

// Dependencies
#include "windowed_running_gradient.h"
#include "running_quadratic_gradient.h"
#include "running_peak_analysis.h"
#include "running_cubic_gradient.h"
#include "rls_analysis_parameters.h"

// Type Definitions
typedef void (*Callback_t)(void);  // Callback function type

typedef enum {
    SWP_INITIAL_ANALYSIS,
    SWP_SEGMENT_ANALYSIS,
    SWP_UPDATE_BUFFER_DIRECTION,
    SWP_UNDECIDED_TREND_CASE,
    SWP_PEAK_CENTERING,           // New state
    SWP_PEAK_FINDING_ANALYSIS,
    SWP_WAITING,
    SWP_STATE_LAST
} SwpState_t;

typedef struct {
    PeakPosition direction;       // Current direction (LEFT, RIGHT, ON_PEAK, UNDECIDED)
    Callback_t callback;          // Callback function
    const double* phaseAngles;    // Pointer to the phase angles array
    uint16_t phase_angle_size;    // Size of the phase angles array
} SlidingWindowAnalysisContext;

// Function Declarations
void startSlidingWindowAnalysis(const double* phaseAngles, uint16_t phase_angle_size, Callback_t callback);

#endif // SLIDING_WINDOW_ANALYSIS_H

//documentation

/*

******compare_gradient_parts************


//  if (gradient > gradient_analysis_params.gradient_threshold) {  // Use globally defined threshold 
 (*increase_count)++;
buna göre countlanıyor. 

double min_total = gradient_analysis_params.minimum_gradient_total;

 // Check if both gradients are insignificant
    if ((result->total_gradient_first_part > -min_total && result->total_gradient_first_part < min_total) &&
        (result->total_gradient_second_part > -min_total && result->total_gradient_second_part < min_total)) {
        return UNDECIDED;  // Both sides have insignificant gradients
    }


*************detect_significant_gradient_trends***************


*/