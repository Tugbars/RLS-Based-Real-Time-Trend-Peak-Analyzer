/**
 * @file sliding_window_analysis.h
 * @brief Header file for sliding window analysis of phase angle data.
 *
 * This header defines the structures, enums, and function prototypes used for performing
 * sliding window analysis on phase angle data collected from an impedance analyzer.
 * The analysis is used to detect trends, peaks, and to adjust the analysis window dynamically.
 *
 * Dependencies:
 * - windowed_running_gradient.h
 * - running_quadratic_gradient.h
 * - running_peak_analysis.h
 * - running_cubic_gradient.h
 * - rls_analysis_parameters.h
 */

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
#include "buffer_manager.h"
#include "mes_buffers.h"

/**
 * @typedef Callback_t
 * @brief Typedef for a callback function used in sliding window analysis.
 *
 * This callback is invoked upon completion or at specific stages of the sliding window analysis.
 */
typedef void (*Callback_t)(void);

extern bool boundaryErrorOccurred;

/**
 * @enum SwpState_t
 * @brief Enumeration of states in the sliding window analysis state machine.
 *
 * This enum defines the various states that the sliding window analysis can be in.
 */
typedef enum {
    SWP_INITIAL_ANALYSIS,
    SWP_SEGMENT_ANALYSIS,
    SWP_UPDATE_BUFFER_DIRECTION,
    SWP_UNDECIDED_TREND_CASE,
    SWP_PEAK_CENTERING,
    SWP_PEAK_FINDING_ANALYSIS,
    SWP_PEAK_TRUNCATION_HANDLING,  // New state
    SWP_WAITING,
    SWP_STATE_LAST  // This should always be last
} SwpState_t;

/**
 * @struct SlidingWindowAnalysisContext
 * @brief Context structure for the sliding window analysis.
 *
 * This structure holds the current state and data required for performing the sliding window analysis.
 */
typedef struct {
    PeakPosition direction;       /**< Current direction based on peak analysis (LEFT, RIGHT, ON_PEAK, UNDECIDED). */
    Callback_t callback;          /**< Callback function to be invoked during the analysis. */
    bool isTruncatedLeft;         // New field
    bool isTruncatedRight;         // New field
    const double* phaseAngles;    /**< Pointer to the array of phase angles to be analyzed. */
    uint16_t phase_angle_size;    /**< Size of the phase angles array. */
} SlidingWindowAnalysisContext;

/**
 * @brief Starts the sliding window analysis on the provided phase angle data.
 *
 * This function initializes and begins the sliding window analysis using the provided phase angle data.
 * It sets up the analysis context and invokes the analysis state machine, which progresses through
 * various states to detect trends, peaks, and adjust the analysis window accordingly.
 *
 * @param phaseAngles       Pointer to the array of phase angle data.
 * @param phase_angle_size  The size (number of elements) of the phase angles array.
 * @param callback          Callback function to be called during the analysis (optional).
 *
 * @note
 * - The phaseAngles array should contain the phase angle measurements collected from the impedance analyzer.
 * - The callback function can be used to perform actions upon certain events or completion of the analysis.
 */
void startSlidingWindowAnalysis(MesSweep_t *sweep, const double* phaseAngles, uint16_t phase_angle_size, Callback_t callback);

void set_boundary_error_flag(uint8_t flag);

#endif // SLIDING_WINDOW_ANALYSIS_H