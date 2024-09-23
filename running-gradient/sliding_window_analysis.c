#include "sliding_window_analysis.h"
#include <stdio.h>

/**
 * @file sliding_window_analysis.h
 * @brief Header file for sliding window analysis of phase angles.
 *
 * This file contains function declarations and data structures used for analyzing phase angle data
 * using sliding window techniques.
 */


/******************************************************************************/
/* Global Variables */
/******************************************************************************/

/**
 * @brief Union representing the status of the sliding window analysis process.
 *
 * This union contains flags that indicate the current status of the analysis,
 * such as whether a peak is found, if the analysis is undecided, or if the peak is centered.
 */
typedef union {
    uint8_t value; /**< The combined value of all status flags. */
    struct {
        uint8_t isPeakFound : 1;            /**< Flag indicating if a peak is found. */
        uint8_t isUndecided : 1;            /**< Flag indicating if the analysis is undecided. */
        uint8_t isCentered : 1;             /**< Flag indicating if the peak is centered. */
        uint8_t isSweepRequested : 1;       /**< Flag indicating if a sweep is requested. */
        uint8_t isSweepDone : 1;            /**< Flag indicating if the sweep is done. */
        uint8_t isNotCentered : 1;          /**< Flag indicating if the peak is still not centered. */
        uint8_t isVerificationFailed : 1;   /**< Flag indicating if peak verification failed. */
        uint8_t reserved : 1;               /**< Reserved bit for future use. */
    };
} SwpStatus_t;

/** @brief Context for sliding window analysis. */
SlidingWindowAnalysisContext ctx; // Sliding window analysis context

/** @brief Current status of the sliding window analysis. */
SwpStatus_t currentStatus = { .value = 0x0 };  // Initialize all flags to 0

/** @brief Current state of the sliding window analysis state machine. */
SwpState_t currentState = SWP_WAITING;  // Initialize to SWP_WAITING state

/******************************************************************************/
/* Function Prototypes (Internal Functions) */
/******************************************************************************/

/* OnEntry and OnExit function declarations for each state */
static void OnEntryInitialAnalysis(void);
static void OnEntrySegmentAnalysis(void);
static void OnEntryUpdateBufferDirection(void);
static void OnEntryUndecidedTrendCase(void);
static void OnEntryPeakFindingAnalysis(void);
static void OnEntryWaiting(void);
static void OnEntryPeakTruncationHandling(void);
static void OnEntryPeakCentering(void);

static void OnExitInitialAnalysis(void);
static void OnExitSegmentAnalysis(void);
static void OnExitUpdateBufferDirection(void);
static void OnExitUndecidedTrendCase(void);
static void OnExitPeakFindingAnalysis(void);
static void OnExitWaiting(void);
static void OnExitPeakCentering(void);
static void OnExitPeakTruncationHandling(void);

void SwpProcessStateChange(void);
static void initBufferManager(MqsRawDataPoint_t* dataBuffer);
static SwpState_t NextState(SwpState_t state);

/******************************************************************************/
/* State Function Definitions */
/******************************************************************************/

/**
 * @brief Structure containing function pointers and completion status for state functions.
 */
typedef struct {
    void (*onEntry)(void);
    void (*onExit)(void);
    bool isComplete;
} StateFuncs_t;

/**
 * @brief Array of state functions and their completion status.
 *
 * This array maps each state to its corresponding onEntry and onExit functions, along with a flag
 * indicating whether the state's processing is complete.
 */
static StateFuncs_t STATE_FUNCS[SWP_STATE_LAST] = {
    {OnEntryInitialAnalysis, OnExitInitialAnalysis, false},         // SWP_INITIAL_ANALYSIS
    {OnEntrySegmentAnalysis, OnExitSegmentAnalysis, false},         // SWP_SEGMENT_ANALYSIS
    {OnEntryUpdateBufferDirection, OnExitUpdateBufferDirection, false}, // SWP_UPDATE_BUFFER_DIRECTION
    {OnEntryUndecidedTrendCase, OnExitUndecidedTrendCase, false},   // SWP_UNDECIDED_TREND_CASE
    {OnEntryPeakCentering, OnExitPeakCentering, false},             // SWP_PEAK_CENTERING
    {OnEntryPeakFindingAnalysis, OnExitPeakFindingAnalysis, false}, // SWP_PEAK_FINDING_ANALYSIS
    {OnEntryPeakTruncationHandling, OnExitPeakTruncationHandling, false}, // SWP_PEAK_TRUNCATION_HANDLING
    {OnEntryWaiting, OnExitWaiting, false}                          // SWP_WAITING
};
/******************************************************************************/
/* Sliding Window Analysis Functions */
/******************************************************************************/


/**
 * @brief Initializes the buffer manager for sliding window analysis.
 *
 * This function sets up the buffer manager with the appropriate parameters, including the buffer size,
 * window size, starting index, and frequency parameters.
 */
static void initBufferManager(MqsRawDataPoint_t* dataBuffer) {
    // Initialize the buffer manager with appropriate parameters
    // Arguments:
    // buffer          -> The buffer that will store the phase and impedance values (array of MqsRawDataPoint_t)
    // BUFFER_SIZE     -> Total size of the buffer (maximum number of data points the buffer can hold)
    // 30              -> Sliding window size (number of data points to analyze in each window)
    // 0               -> Starting index in the phaseAngles array (this could be any index where you want to start the analysis)
    // 11300.0         -> Starting frequency for the frequency sweep (in Hz, e.g., starting at 11300 Hz)
    // 1.0             -> Frequency increment per step (in Hz, e.g., increment by 1 Hz for each data point)
    int start_index = 207;
    init_buffer_manager(dataBuffer, BUFFER_SIZE, WINDOW_SIZE, start_index, 11300.0, 1.0);

    // Debugging: Print initialization state
    printf("[DEBUG] Buffer Manager initialized:\n");
    printf("Buffer size: %d, Window size: %d, Start index: %d\n", BUFFER_SIZE, WINDOW_SIZE, start_index);
}


/**
 * @brief Starts the sliding window analysis process.
 *
 * This function initializes the context for the sliding window analysis, sets up the buffer manager,
 * and starts the state machine for processing.
 *
 * @param phaseAngles Pointer to the array of phase angles.
 * @param phase_angle_size Size of the phase angle array.
 * @param callback Function pointer to the callback function to be executed after analysis.
 */
void startSlidingWindowAnalysis(MesSweep_t *sweep, const double* phaseAngles, uint16_t phase_angle_size, Callback_t callback) {
    ctx.phaseAngles = phaseAngles;
    ctx.phase_angle_size = phase_angle_size;
    ctx.callback = callback;
    ctx.isTruncatedLeft = false;
    ctx.isTruncatedRight = false;

    // Initialize the buffer manager
    initBufferManager(sweep->data);
    
    	// detect_significant_gradient_trends, determine_trend_direction
	init_cubic_rls_analysis_parameters(                                                                                                        //CHECKS IF THERE IS A PEAK VIA COUNTING/SUMMING CONSISTENT GRADIENT INCREASES 
        12.0,   // significance_thresh: Threshold for determining significant cubic trends.                                                     
               // Trends with a sum of gradients above this threshold are considered significant.
        5,     // duration_thresh: Minimum number of consecutive points required for a trend to be considered significant.                     
               // Ensures that only sustained trends are analyzed.
        2,     // min_trend_count: Minimum number of consistent trends required for cubic analysis.
               // Helps filter out noise and minor fluctuations.
        2,     // max_third_order_trend_decrease_count: Maximum number of consecutive decreases allowed in a cubic trend.           
               // Allows for minor fluctuations without discarding a potentially valid trend.
        2      // max_third_order_trend_increase_count: Maximum number of consecutive increases allowed in a cubic trend.
               // Helps to filter out noise when identifying a cubic trend.
    );
    
    // Initialize the quadratic RLS analysis parameters with the necessary thresholds.                                                          
    init_quadratic_rls_analysis_parameters(
        1.0,  // centered_gradient_sum: If the total second-order gradient sum is less than or equal to this value,
              // the peak is considered centered based on the gradient analysis.
        2,    // max_decrease_count: Maximum number of consecutive negative gradients allowed
              // when tracking an increasing trend. Allows for minor fluctuations without discarding a valid trend.
        2,    // max_increase_count: Maximum number of consecutive positive gradients allowed
              // when tracking a decreasing trend. Helps to filter out noise when identifying a decreasing trend.
        5,    // min_trend_count: Minimum number of consistent trends required to consider a peak valid.
              // Ensures that only sustained trends are analyzed.
        2     // allowable_inconsistency_count: Allowable number of inconsistencies in trend detection.
              // Permits minor deviations without discarding the trend.
    );
    
    // Initialize the on-peak analysis parameters for average gradient thresholds and consistent trend count
   // // INCREASELERIN BŞALADIKLARI YERDEN YAP BUNU?!
    init_on_peak_analysis_parameters( 
        0.8,  // min_avg_increase: Minimum average increase required to consider a significant increasing trend. LEFT SIDE OF THE PEAK            //KULLANILIYOR MU?
               // Ensures that only trends with substantial average gradient are flagged as significant.
        -0.19, // min_avg_decrease: Minimum average decrease required to consider a significant decreasing trend. RIGHT SIDE OF THE PEAK          //KULLANILIYOR MU?
               // Ensures that only trends with substantial average gradient are flagged as significant.
        5      // min_consistent_trend_count: Minimum number of consecutive data points required for a trend to be considered consistent. 
               // Helps ensure that only sustained trends are analyzed.
    );
    
    // Initialize gradient analysis parameters
    init_gradient_analysis_params(
        0.1,  // gradient_thresh: Threshold for determining a significant gradient increase.                                                    //COMPARE GRADIENT PARTSIN UYESI. SADECE BUYUK GRADIENT ARTIŞLARINI SAYIYOR. 
        0.6   // min_gradient_total: Minimum total gradient to flag a significant trend.                                                        //BUNUN ALTINDAKİLERE UNDECIDED DİYORUZ. 
    );
    
        // Enable sweep request
    currentStatus.isSweepRequested = 1;  // Start sweep request

    // Set the initial state to waiting
    currentState = SWP_WAITING;

    // Set the initial state and start the state machine
    SwpProcessStateChange();  // Start the state machine
}

/******************************************************************************/
/* State Machine Entry Functions */
/******************************************************************************/

/**
 * @brief Entry function for the SWP_INITIAL_ANALYSIS state.
 *
 * This function loads the initial buffer with phase angle data and sets the sweep request flag.
 */
static void OnEntryInitialAnalysis(void) {
    //printf("→→→→→Entering SWP_INITIAL_ANALYSIS state.\n");

    load_initial_buffer(ctx.phaseAngles, ctx.phase_angle_size);
    
      // Debugging buffer state after load
    //printf("[DEBUG] Buffer Manager State after initial loading:\n");
    //printf("Current buffer index: %d\n", buffer_manager.current_buffer_index);
    //printf("Current phase index: %d\n", buffer_manager.current_phase_index);
    //printf("Window size: %d\n", buffer_manager.window_size);
    //printf("Buffer size: %d\n", buffer_manager.buffer_size);
    currentStatus.isSweepRequested = true;
    STATE_FUNCS[SWP_INITIAL_ANALYSIS].isComplete = true;
    SwpProcessStateChange();  // Move to next state
}


/**
 * @brief Entry function for the SWP_SEGMENT_ANALYSIS state.
 *
 * This function performs segment analysis on the current window of data to determine the direction of movement.
 */
static void OnEntrySegmentAnalysis(void) {
    //printf("→→→→→Entering SWP_SEGMENT_ANALYSIS state.\n");
    
    float forgetting_factor = 0.5f;

    SegmentAnalysisResult result = segment_trend_and_concavity_analysis(
        &buffer_manager.buffer[buffer_manager.current_buffer_index],
        buffer_manager.window_size,
        forgetting_factor
    );

    ctx.direction = result.nextDirection;

    if (ctx.direction == ON_PEAK) {
        currentStatus.isPeakFound = 1;
    } else if (ctx.direction == UNDECIDED) {
        currentStatus.isUndecided = 1;
    }

    STATE_FUNCS[SWP_SEGMENT_ANALYSIS].isComplete = true;
    SwpProcessStateChange();
}

/**
 * @brief Entry function for the SWP_UPDATE_BUFFER_DIRECTION state.
 *
 * This function updates the buffer based on the determined direction from the segment analysis.
 */
static void OnEntryUpdateBufferDirection(void) {
    // Pass the movement amount as an argument
    update_buffer_for_direction(ctx.phaseAngles, ctx.direction, buffer_manager.window_size / 2);

    STATE_FUNCS[SWP_UPDATE_BUFFER_DIRECTION].isComplete = true;
    SwpProcessStateChange();  // Move to next state
}

/**
 * @brief Entry function for the SWP_UNDECIDED_TREND_CASE state.
 *
 * This function handles the case when the segment analysis result is undecided by moving the window forward.
 */
static void OnEntryUndecidedTrendCase(void) {  //ÜZERİNE YAZIYOR.
    //printf("→→→→→→Entering SWP_UNDECIDED_TREND_CASE state.\n");

    handle_undecided_case(ctx.phaseAngles, ctx.phase_angle_size);

    STATE_FUNCS[SWP_UNDECIDED_TREND_CASE].isComplete = true;
    SwpProcessStateChange();
}

/**
 * @brief Entry function for the SWP_PEAK_CENTERING state.
 *
 * This function attempts to center the detected peak within the current window of data.
 * It performs the following steps:
 *
 * ### Function Workflow:
 * 1. **Reset State Completion Flag**:
 *    - Resets the `isComplete` flag for the current state to indicate processing has started.
 * 2. **Compute Total Sum of Second-Order Gradients**:
 *    - Calculates the total sum of second-order gradients over the current data window.
 *    - Uses the `compute_total_second_order_gradient` function with a forgetting factor.
 *    - The second-order gradient provides information about the curvature of the data.
 * 3. **Check if Peak is Centered**:
 *    - Compares the total gradient sum with a predefined threshold (`centered_gradient_sum`).
 *    - If the sum is less than or equal to this threshold, the peak is considered centered.
 *    - Sets `currentStatus.isCentered` to `1` and marks the state as complete.
 * 4. **Analyze Gradient Trends**:
 *    - If the peak is not centered, it analyzes the gradient trends using quadratic regression.
 *    - Calls `track_gradient_trends_with_quadratic_regression` to get increasing and decreasing trend information.
 * 5. **Validate Trend Data**:
 *    - Checks if both increasing and decreasing trends are valid.
 *    - If not, considers the peak centered to prevent infinite loops and marks the state as complete.
 * 6. **Calculate Trend Durations**:
 *    - Calculates the durations of the increasing and decreasing trends.
 *    - Adjusts for buffer wrapping using modulo arithmetic.
 * 7. **Determine Shift Direction and Amount**:
 *    - Compares the durations of the trends to decide the direction to shift.
 *    - If the increasing duration is longer, the peak is to the left (shift right).
 *    - If the decreasing duration is longer, the peak is to the right (shift left).
 *    - Calculates the shift amount as the difference between the durations.
 * 8. **Adjust Data Window**:
 *    - Calls `update_buffer_for_direction` to shift the data window accordingly.
 *    - Does not mark the state as complete to allow re-entry and re-evaluation.
 * 9. **Handle Peak Centered Cases**:
 *    - If the durations are equal or the shift amount is zero, considers the peak centered.
 *    - Marks the state as complete.
 * 10. **State Transition**:
 *     - If the state is complete, calls `SwpProcessStateChange` to transition to the next state.
 *
 * ### Mathematical Background:
 * - **Second-Order Gradients**:
 *   - Provide information about the curvature (concavity/convexity) of the data.
 *   - Calculated using quadratic regression over the data window.
 *   - A total sum close to zero indicates symmetry around the peak.
 * - **Quadratic Regression**:
 *   - Fits a quadratic model to the data: \( y = ax^2 + bx + c \).
 *   - The second derivative (second-order gradient) is \( 2a \), representing curvature.
 * - **Trend Analysis**:
 *   - Identifies consistent increasing or decreasing trends in the second-order gradients.
 *   - Uses durations of trends to determine the asymmetry of the peak.
 * - **Buffer Shifting**:
 *   - Adjusts the data window by shifting it left or right to center the peak.
 *   - The shift amount is based on the difference in trend durations.
 *
 * @note
 * - The function uses global variables like `buffer_manager`, `quadratic_analysis_params`, and `currentStatus`.
 * - It is part of a state machine managing the sliding window analysis process.
 * - The state machine will re-enter this state if the peak is not yet centered after shifting.
 *
 * @see compute_total_second_order_gradient
 * @see track_gradient_trends_with_quadratic_regression
 * @see update_buffer_for_direction
 * @see SwpProcessStateChange
 */
static void OnEntryPeakCentering(void) {
    // Reset the isComplete flag
    STATE_FUNCS[SWP_PEAK_CENTERING].isComplete = false;
    
    // Reset the isNotCentered flag
    currentStatus.isNotCentered = 0;

    // Get the start index in the buffer
    uint16_t start_index = buffer_manager.current_buffer_index;
    uint16_t window_size = buffer_manager.window_size;

    // Compute the total sum of second-order gradients using the buffer directly
    double total_gradient_sum = compute_total_second_order_gradient(
        buffer_manager.buffer,
        buffer_manager.buffer_size,
        start_index,
        0.2 // Forgetting factor
    );

    printf("Total sum of second-order gradients: %.6f\n", total_gradient_sum);

    // If total_gradient_sum <= centered_gradient_sum, consider the peak centered
    if (total_gradient_sum <= quadratic_analysis_params.centered_gradient_sum &&
        total_gradient_sum > -quadratic_analysis_params.centered_gradient_sum) {
        currentStatus.isCentered = 1;
        printf("Peak is centered based on total_gradient_sum.\n");
        STATE_FUNCS[SWP_PEAK_CENTERING].isComplete = true;

    } else {
        // Else, we need to adjust the buffer to center the peak
        // Use track_gradient_trends_with_quadratic_regression to get trend info
        GradientTrendResult gradient_trends = track_gradient_trends_with_quadratic_regression(
            buffer_manager.buffer,
            buffer_manager.buffer_size,
            start_index,
            window_size,
            0.5 // Forgetting factor
        );

        // Check if both trends are valid
        if (!gradient_trends.increase_info.valid || !gradient_trends.decrease_info.valid) {
            printf("Invalid trend data. Cannot proceed with centering.\n");
            currentStatus.isCentered = 1; // Consider it centered for now
            STATE_FUNCS[SWP_PEAK_CENTERING].isComplete = true;
        } else {
            // Calculate durations of increase and decrease trends
            uint16_t increase_start = gradient_trends.increase_info.start_index;
            uint16_t increase_end = gradient_trends.increase_info.end_index;
            uint16_t decrease_start = gradient_trends.decrease_info.start_index;
            uint16_t decrease_end = gradient_trends.decrease_info.end_index;

            // Adjust indices relative to the buffer
            uint16_t increase_duration = (increase_end + buffer_manager.buffer_size - increase_start) % buffer_manager.buffer_size;
            uint16_t decrease_duration = (decrease_end + buffer_manager.buffer_size - decrease_start) % buffer_manager.buffer_size;

            int shift_amount = 0;
            int direction = UNDECIDED; // LEFT_SIDE or RIGHT_SIDE

            if (increase_duration > decrease_duration) {
                // Peak is to the left; we need to move right
                shift_amount = (increase_duration - decrease_duration) + 1;
                direction = RIGHT_SIDE;
                printf("Increase duration (%u) > decrease duration (%u). Moving right by %d.\n",
                       increase_duration, decrease_duration, shift_amount);
                       
            } else if (decrease_duration > increase_duration) {
                // Peak is to the right; we need to move left
                shift_amount = (decrease_duration - increase_duration) + 1;
                direction = LEFT_SIDE;
                printf("Decrease duration (%u) > increase duration (%u). Moving left by %d.\n",
                       decrease_duration, increase_duration, shift_amount);
            } else {
                // Durations are equal; consider the peak centered
                currentStatus.isCentered = 1;
                printf("Increase and decrease durations are equal. Peak is centered.\n");
                STATE_FUNCS[SWP_PEAK_CENTERING].isComplete = true;
            }

            if (shift_amount == 0) {
                // No need to shift
                currentStatus.isCentered = 1;
                printf("No shift needed. Peak is centered.\n");
                STATE_FUNCS[SWP_PEAK_CENTERING].isComplete = true;
            } else {
                // Use move_window_and_update_if_needed to shift the window and update if necessary
                move_window_and_update_if_needed(ctx.phaseAngles, direction, shift_amount);  
                currentStatus.isCentered = 1;  // Peak is considered centered after window move

                // The state machine will re-enter this state to check again
                STATE_FUNCS[SWP_PEAK_CENTERING].isComplete = true;
            }
        }
    }

    // The state machine will handle the transition based on the isComplete flag
    if (STATE_FUNCS[SWP_PEAK_CENTERING].isComplete) {
        SwpProcessStateChange();
    }
}

/**
 * @brief Entry function for the SWP_PEAK_FINDING_ANALYSIS state.
 *
 * This function performs several key tasks to verify the centering and validity of the detected peak:
 * 
 * ### Step 1: Re-check the Centration of the Peak
 * - The function first computes the total sum of the second-order gradients across the current data window
 *   using a quadratic regression analysis. This sum reflects the curvature of the data and helps determine
 *   whether the peak is centered in the window.
 * - If the sum is outside the predefined threshold range (`centered_gradient_sum`), it indicates that the peak is not yet
 *   centered, and the function returns the state machine to the peak-centering state.
 *
 * ### Step 2: Perform Peak Verification
 * - If the sum of second-order gradients is within the acceptable threshold, the function proceeds to verify the peak.
 * - The verification process, done via `find_and_verify_quadratic_peak`, checks if the data exhibits a valid peak shape:
 *   - Increasing trends on the left side of the peak.
 *   - Decreasing trends on the right side of the peak.
 * - It also accounts for truncation scenarios. If the verification process encounters the boundaries of the window before
 *   completing the necessary trend checks, truncation flags (`is_truncated_left` and `is_truncated_right`) are set.
 * 
 * ### Step 3: Handle Truncation Scenarios
 * - If the peak verification is truncated on either the left or right side, this information is logged.
 * - The function prints specific messages indicating whether truncation occurred on the left, right, or both sides
 *   of the data window.
 *
 * ### Step 4: Final Peak Verification and State Transition
 * - If the peak is successfully verified based on the trend analysis, the function prints the sliding window analysis
 *   interval and proceeds to mark the sweep as complete (`isSweepDone` flag).
 * - If the peak verification fails (e.g., the data does not exhibit the required trend patterns), the function returns
 *   to the peak-centering state by setting the `isNotCentered` flag.
 *
 * ### Summary of Key Functionality:
 * - Re-checks peak centration based on the sum of second-order gradients.
 * - Verifies the peak using quadratic regression and trend analysis.
 * - Handles boundary truncations during verification.
 * - Directs the state machine to either the peak-centering or peak-verification complete state.
 *
 * @see compute_total_second_order_gradient
 * @see find_and_verify_quadratic_peak
 * @see print_analysis_interval
 * @see SwpProcessStateChange
 */
/**
 * @brief Entry function for the SWP_PEAK_FINDING_ANALYSIS state.
 *
 * This function verifies the peak at the current buffer index and handles truncation cases.
 */
static void OnEntryPeakFindingAnalysis(void) {
    // Print the total second-order gradient sum again to verify the peak
    double total_gradient_sum = compute_total_second_order_gradient(
        buffer_manager.buffer,
        buffer_manager.buffer_size,
        buffer_manager.current_buffer_index,
        0.2 // Forgetting factor
    );
    
    printf("Total sum of second-order gradients during peak verification: %.6f\n", total_gradient_sum);

    // Check if the sum of second-order gradients is still outside the threshold
    if (total_gradient_sum > quadratic_analysis_params.centered_gradient_sum || 
        total_gradient_sum < -quadratic_analysis_params.centered_gradient_sum) {
        // Peak is not centered yet, flag it and return to peak centering
        printf("Peak not centered, returning to peak centering state.\n");
        currentStatus.isNotCentered = 1;  // Set flag to indicate it's not centered
    } else {
        // Now handle the truncation flags in case the verification was truncated
        QuadraticPeakAnalysisResult verification_result = find_and_verify_quadratic_peak(
            buffer_manager.buffer,
            buffer_manager.buffer_size,
            buffer_manager.current_buffer_index,  // Peak index
            0.5 // Forgetting factor
        );

        // Set truncation flags based on the verification result
      if (verification_result.is_truncated_left || verification_result.is_truncated_right) {
            if (verification_result.is_truncated_left) {
            printf("Peak verification truncated on the left side.\n");
            ctx.isTruncatedLeft = true;
            }
            if (verification_result.is_truncated_right) {
            printf("Peak verification truncated on the right side.\n");
            ctx.isTruncatedRight = true;
            }
        } else {
            // Peak verification successful
            if (verification_result.peak_found) {
                printf("Peak verification successful, peak is centered.\n");
                print_analysis_interval(ctx.phaseAngles, ctx.phase_angle_size);  // Print the buffer interval
                currentStatus.isSweepDone = 1;  // Mark the sweep as done
            } else {
                printf("Peak verification failed, returning to peak centering.\n");
                currentStatus.isNotCentered = 1;  // Set flag to indicate it's not centered
            }
        }
    }

    // Mark the state as complete
    STATE_FUNCS[SWP_PEAK_FINDING_ANALYSIS].isComplete = true;
    SwpProcessStateChange();  // Process the next state
}

/**
 * @brief Entry function for the SWP_PEAK_TRUNCATION_HANDLING state.
 *
 * This function handles the scenario where peak verification during the peak finding analysis
 * is incomplete due to data truncation on the left or right side of the data window.
 * It adjusts the data window by moving it in the appropriate direction and attempts to
 * re-verify the peak after the adjustment.
 *
 * ### Function Workflow:
 * 1. **Initialization**:
 *    - Initializes a flag `peak_verified_after_truncation` to track if the peak is verified
 *      after handling truncation.
 * 2. **Left Side Truncation Handling**:
 *    - Checks if the peak verification was truncated on the left side (`ctx.isTruncatedLeft`).
 *    - If so, it moves the data window to the left by a predefined amount (e.g., 5 positions)
 *      using `move_window_and_update_if_needed`.
 *    - Re-verifies the peak at the current buffer index by calling `verify_peak_at_index`.
 *    - Resets the `isTruncatedLeft` flag in the context.
 * 3. **Right Side Truncation Handling**:
 *    - Checks if the peak verification was truncated on the right side (`ctx.isTruncatedRight`).
 *    - If so, it moves the data window to the right by a predefined amount (e.g., 5 positions)
 *      using `move_window_and_update_if_needed`.
 *    - Re-verifies the peak at the current buffer index by calling `verify_peak_at_index`.
 *    - Resets the `isTruncatedRight` flag in the context.
 * 4. **Peak Verification Post-Truncation**:
 *    - If the peak is successfully verified after adjustments, it marks the sweep as done
 *      by setting `currentStatus.isSweepDone`.
 *    - If the peak verification fails after adjustments, it sets `currentStatus.isVerificationFailed`
 *      to indicate that the peak verification failed, prompting a return to peak finding analysis.
 * 5. **State Transition**:
 *    - Marks the current state as complete (`isComplete = true`).
 *    - Calls `SwpProcessStateChange()` to transition to the next state based on the updated status flags.
 *
 * ### Related Functions:
 * - **`find_and_verify_quadratic_peak`**:
 *   - Detects and verifies peaks using second-order gradients from quadratic regression.
 *   - Returns a `QuadraticPeakAnalysisResult` containing peak verification status and truncation flags.
 * - **`verify_quadratic_peak`**:
 *   - Verifies a detected peak by checking for consistent increasing trends on the left side and
 *     decreasing trends on the right side of the peak.
 *   - Considers truncation scenarios if the data window does not fully encompass the required trends.
 *
 * ### Notes:
 * - The function relies on the context (`ctx`) to check for truncation flags (`isTruncatedLeft`, `isTruncatedRight`)
 *   and to reset them after handling.
 * - The function modifies `currentStatus` to communicate the result of the truncation handling to the state machine.
 * - The data window adjustments are crucial for ensuring that the peak is fully captured within the analysis window,
 *   which is essential for accurate peak verification.
 *
 * @see find_and_verify_quadratic_peak
 * @see verify_quadratic_peak
 * @see move_window_and_update_if_needed
 * @see SwpProcessStateChange
 */
static void OnEntryPeakTruncationHandling(void) {
    bool peak_verified_after_truncation = false;

    if (ctx.isTruncatedLeft) {
        printf("Handling truncation on the left side.\n");
        move_window_and_update_if_needed(ctx.phaseAngles, LEFT_SIDE, 5);
        peak_verified_after_truncation = verify_peak_at_index(buffer_manager.current_buffer_index);
        ctx.isTruncatedLeft = false;  // Reset flag
    }

    if (ctx.isTruncatedRight) {
        printf("Handling truncation on the right side.\n");
        move_window_and_update_if_needed(ctx.phaseAngles, RIGHT_SIDE, 5);
        peak_verified_after_truncation = verify_peak_at_index(buffer_manager.current_buffer_index);
        ctx.isTruncatedRight = false;  // Reset flag
    }

    // Check if the peak is verified after adjustments (post-truncation)
    if (peak_verified_after_truncation) {
        printf("Peak verification successful after truncation handling, peak is centered.\n");
        print_analysis_interval(ctx.phaseAngles, ctx.phase_angle_size);  // Print the buffer interval
        currentStatus.isSweepDone = 1;  // Mark the sweep as done
    } else {
        // If peak is not verified, go back to peak finding analysis
        printf("Peak verification failed after truncation handling, returning to peak finding analysis.\n");
        currentStatus.isVerificationFailed = 1;  // Set flag to indicate verification failed
    }

    // Mark the state as complete
    STATE_FUNCS[SWP_PEAK_TRUNCATION_HANDLING].isComplete = true;
    SwpProcessStateChange();  // Process the next state
}


/**
 * @brief Entry function for the SWP_WAITING state.
 *
 * This function waits for a sweep request and executes the callback if provided.
 */
static void OnEntryWaiting(void) {

    if (currentStatus.isSweepRequested) {
        STATE_FUNCS[SWP_WAITING].isComplete = true;
    }
}

/******************************************************************************/
/* Exit Functions */
/******************************************************************************/
static void OnExitInitialAnalysis(void) {
    //printf("Exiting SWP_INITIAL_ANALYSIS state.\n");
}

static void OnExitSegmentAnalysis(void) {
    //printf("Exiting SWP_SEGMENT_ANALYSIS state.\n");
}

static void OnExitUpdateBufferDirection(void) {
    //printf("Exiting SWP_UPDATE_BUFFER_DIRECTION state.\n");
}

static void OnExitUndecidedTrendCase(void) {
    //printf("Exiting SWP_UNDECIDED_TREND_CASE state.\n");
}

static void OnExitPeakFindingAnalysis(void) {
    //printf("Exiting SWP_PEAK_FINDING_ANALYSIS state.\n");
}

static void OnExitWaiting(void) {
    //printf("Exiting SWP_WAITING state.\n");
}

static void OnExitPeakCentering(void) {
    //printf("Exiting SWP_PEAK_CENTERING state.\n");
}

static void OnExitPeakTruncationHandling(void) {
    // Any cleanup code can go here
    // For now, we'll leave it empty
}

/******************************************************************************/
/* State Machine Process and Transitions */
/******************************************************************************/

/**
 * @brief Processes state changes in the sliding window analysis state machine.
 *
 * This function manages transitions between states based on the current status flags and ensures that
 * the appropriate onEntry and onExit functions are called for each state.
 */
void SwpProcessStateChange(void) {
    bool stateChanged;
    
    // The state machine should process until no further state changes occur
    do {
        // Determine the next state from the current state
        SwpState_t nextState = NextState(currentState);  
        stateChanged = (nextState != currentState);

        // Reset the current status flags before moving to the next state
        currentStatus.value = 0x0;

        if (stateChanged) {
            //printf("[DEBUG] Transitioning from state %d to state %d\n", currentState, nextState);

            // Exit the current state (if applicable)
            if (STATE_FUNCS[currentState].onExit != NULL) {
                STATE_FUNCS[currentState].onExit();
            }

            // Move to the new state
            currentState = nextState;

            // Reset the completion flag for the new state
            STATE_FUNCS[currentState].isComplete = false;

            // Enter the new state
            if (STATE_FUNCS[currentState].onEntry != NULL) {
                STATE_FUNCS[currentState].onEntry();
            }
        }

        // If the new state is complete, proceed to the next state
        if (STATE_FUNCS[currentState].isComplete) {
            printf("[DEBUG] State %d is complete. Proceeding to next state...\n", currentState);
        }

    } while (stateChanged && STATE_FUNCS[currentState].isComplete);

    // Handle the callback execution when we reach the waiting state
    if (currentState == SWP_WAITING && ctx.callback) {
        printf("[DEBUG] Executing callback as we have reached SWP_WAITING.\n");
        ctx.callback();  // Execute the callback function
        ctx.callback = NULL;  // Clear the callback to avoid re-execution
    }
}

/******************************************************************************/
/* Helper Functions for State Transitioning */
/******************************************************************************/

/**
 * @brief Determines the next state in the state machine based on the current state and status flags.
 *
 * @param state The current state.
 * @return SwpState_t The next state to transition to.
 */
static SwpState_t NextState(SwpState_t state) {
    switch (state) {
        case SWP_INITIAL_ANALYSIS:
            if (currentStatus.isPeakFound) {
                return SWP_PEAK_CENTERING;
            }
            if (currentStatus.isUndecided) {
                return SWP_UNDECIDED_TREND_CASE;
            }
            return SWP_SEGMENT_ANALYSIS;

        case SWP_SEGMENT_ANALYSIS:
            if (currentStatus.isUndecided) {
                return SWP_UNDECIDED_TREND_CASE;
            }
            if (currentStatus.isPeakFound) {
                return SWP_PEAK_CENTERING;
            }
            return SWP_UPDATE_BUFFER_DIRECTION;

        case SWP_PEAK_CENTERING:
            if (currentStatus.isCentered) {
                return SWP_PEAK_FINDING_ANALYSIS;
            } else {
                return SWP_PEAK_CENTERING;
            }

        case SWP_UPDATE_BUFFER_DIRECTION:
            return SWP_SEGMENT_ANALYSIS;

        case SWP_UNDECIDED_TREND_CASE:
            return SWP_SEGMENT_ANALYSIS;

        case SWP_PEAK_FINDING_ANALYSIS:
            if (ctx.isTruncatedLeft || ctx.isTruncatedRight) {
                return SWP_PEAK_TRUNCATION_HANDLING;
            }
            if (currentStatus.isNotCentered) {
                return SWP_PEAK_CENTERING;
            }
            if (currentStatus.isSweepDone) {
                return SWP_WAITING;
            }
            return SWP_PEAK_FINDING_ANALYSIS;

        case SWP_PEAK_TRUNCATION_HANDLING:
            if (currentStatus.isSweepDone) {
                return SWP_WAITING;
            }
            if (currentStatus.isVerificationFailed) {
                return SWP_PEAK_FINDING_ANALYSIS;
            }
            // Remain in truncation handling if no flags are set
            return SWP_PEAK_TRUNCATION_HANDLING;

        case SWP_WAITING:
            if (currentStatus.isSweepRequested) {
                return SWP_INITIAL_ANALYSIS;
            }
            return SWP_WAITING;

        default:
            return SWP_WAITING;
    }
}


