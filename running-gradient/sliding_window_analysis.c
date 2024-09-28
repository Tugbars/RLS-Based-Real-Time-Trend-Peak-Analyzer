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
        uint8_t isBoundaryError : 1;        /**< Flag indicating if a boundary error occurred. */
    };
} SwpStatus_t;

/** @brief Context for sliding window analysis. */
SlidingWindowAnalysisContext ctx; // Sliding window analysis context

/** @brief Current status of the sliding window analysis. */
SwpStatus_t currentStatus = { .value = 0x0 };  // Initialize all flags to 0

/** @brief Current state of the sliding window analysis state machine. */
SwpState_t currentState = SWP_WAITING;  // Initialize to SWP_WAITING state

// Add a global flag to track boundary error
bool boundaryErrorOccurred = false;

/*******************************************************************************
 * Undecided deadlock error control
 ******************************************************************************/

/**
 * @brief Tracks the number of consecutive times the state machine enters the undecided case.
 *un
 * This variable is incremented each time the state machine enters the `SWP_UNDECIDED_TREND_CASE` state.
 * It limits the number of times the system can enter the undecided state within a session. If the counter
 * exceeds a predefined limit (e.g., 7), the system will transition to the `SWP_WAITING` state to prevent
 * the state machine from getting stuck in an undecided loop. This helps ensure that the analysis progresses
 * or halts gracefully if a decision cannot be reached.
 */
static uint8_t undecided_case_counter = 0;

/**
 * @brief Flag indicating that an error occurred due to exceeding the undecided case limit.
 *
 * This flag is set to `true` when the `undecided_case_counter` exceeds the maximum allowed value
 * (e.g., 7 undecided states). When this flag is set, the system will transition to the `SWP_WAITING` state,
 * and a warning will be displayed in the callback function executed by the state machine. The flag
 * will be reset when the `SWP_WAITING` state is entered, indicating that the system is ready for
 * a new session or sweep.
 */
static bool undecided_error_flag = false;

/** @brief Maximum number of times the state machine can enter the undecided case within a session.
 *
 * This variable defines the limit on how many times the state machine can enter the `SWP_UNDECIDED_TREND_CASE`
 * state before it triggers an error. If the `undecided_case_counter` exceeds this value, the system will
 * transition to the `SWP_WAITING` state and raise an error flag, indicating that a resolution could not be reached.
 */
static const uint8_t MAX_UNDECIDED_CASE_ATTEMPTS = 7;

/*******************************************************************************
 * Peak centering deadlock error control
 ******************************************************************************/

/**
 * @brief Defines the maximum number of peak centering attempts.
 * 
 * This constant sets a limit on the number of times the algorithm will attempt to center the peak.
 * After reaching this limit, the system will stop further attempts to prevent an infinite loop in cases
 * where the peak cannot be accurately centered. This ensures that the process either succeeds within 
 * a reasonable number of attempts or gracefully exits if centering is not possible.
 */
#define MAX_CENTERING_ATTEMPTS 4

/**
 * @brief Tracks the number of attempts made to center the peak.
 * 
 * This variable is used to limit the number of centering attempts in the peak centering state. 
 * Each time the peak centering process fails to accurately center the peak, this counter is 
 * incremented. Once it reaches a predefined maximum number of attempts (MAX_CENTERING_ATTEMPTS),
 * the system will exit the centering state and transition to a waiting or error state to prevent 
 * an infinite loop in cases where the peak cannot be centered properly.
 */
static uint8_t centering_attempts = 0;

/**
 * @brief Adjusts the forgetting factor used in peak centering analysis.
 * 
 * The forgetting factor controls how much weight is given to more recent data points during the 
 * gradient analysis for peak centering. This variable is incremented by 0.1 after each centering 
 * attempt to give more weight to the newer data points, allowing the algorithm to make progressively 
 * larger adjustments in subsequent attempts. This helps improve the accuracy of peak centering 
 * over multiple attempts.
 */
static double centering_forgetting_factor = 0.7;

/*******************************************************************************
 * Repetitive shift deadlock error control
 ******************************************************************************/
 
 /**
 * @brief Size of the shift tracker array used to detect alternating left and right shifts.
 *
 * This constant defines the size of the `shift_tracker` array, which stores the last two
 * shift directions (LEFT, RIGHT) made during the sliding window analysis. The purpose is to
 * detect whether the analysis is getting stuck in an alternating pattern of shifts.
 */
#define SHIFT_TRACKER_SIZE 2

/**
 * @brief Array to track the last two shift directions.
 *
 * This array stores the most recent shift directions (LEFT or RIGHT) that occurred during the
 * sliding window analysis. If the shift directions alternate between LEFT and RIGHT over the
 * last two shifts, the system concludes that the analysis is stuck and raises an error.
 */
static PeakPosition shift_tracker[SHIFT_TRACKER_SIZE] = {UNDECIDED, UNDECIDED};

/**
 * @brief Index for tracking the current position in the `shift_tracker` array.
 *
 * This variable tracks the current index in the `shift_tracker` array, which stores the last two
 * shift directions. The array functions as a circular buffer, where the index wraps around after
 * reaching the size of the array (`SHIFT_TRACKER_SIZE`).
 */
static uint8_t shift_tracker_index = 0;

/**
 * @brief Updates the shift tracker with the latest direction.
 *
 * This function updates the `shift_tracker` array with the most recent shift direction. The array
 * functions as a circular buffer, meaning that once it reaches its size limit, it starts overwriting
 * the oldest direction. The function also checks whether the last two shift directions alternated
 * between LEFT and RIGHT. If this condition is met, an error is raised, and the analysis transitions
 * to the `SWP_WAITING` state.
 *
 * @param new_direction The new shift direction (LEFT, RIGHT, or UNDECIDED).
 */
static void update_shift_tracker(PeakPosition new_direction) {
    // Update the shift tracker with the latest direction
    shift_tracker[shift_tracker_index] = new_direction;
    shift_tracker_index = (shift_tracker_index + 1) % SHIFT_TRACKER_SIZE;

    // Check if the analysis is stuck (LEFT → RIGHT → LEFT or RIGHT → LEFT → RIGHT)
    if ((shift_tracker[0] == LEFT_SIDE && shift_tracker[1] == RIGHT_SIDE) ||
        (shift_tracker[0] == RIGHT_SIDE && shift_tracker[1] == LEFT_SIDE)) {
        printf("Error: Sliding window analysis is stuck alternating between left and right shifts.\n");
        boundaryErrorOccurred = true;  // Raise an error flag
        currentStatus.isSweepDone = 1; // End the analysis and go to SWP_WAITING
    }
}

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
 * @brief Performs a buffer update if one is required.
 *
 * This function checks whether a buffer update is necessary based on the `buffer_update_info` structure,
 * which contains information about pending buffer updates. If a buffer update is required (indicated by the
 * `needs_update` flag), the function proceeds to update the phase angles in the buffer by calling 
 * `update_phaseAngle_to_buffer`. Once the update is completed, the function triggers a state transition
 * by calling `SwpProcessStateChange()`. If no buffer update is needed, it directly calls `SwpProcessStateChange()`.
 *
 * @param phaseAngles Pointer to the array of phase angles that are used in the analysis.
 */
void SweepSampleCb(const double* phaseAngles) {
    if (buffer_update_info.needs_update) {
        AdptSweepAddDataPoint(                                                    //MESSWEEPADDDATAPOINT YAPISINA UYMALI
            phaseAngles,
            buffer_update_info.phase_index_start,
            buffer_update_info.buffer_start_index,
            buffer_update_info.move_amount
        );

        // Reset the update flag
        buffer_update_info.needs_update = false;

        // Now proceed to the next state
        SwpProcessStateChange();
    } else {
        // No update needed, proceed to next state
        SwpProcessStateChange();
    }
}

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
    int start_index = 162;                                                                                                                    //START FREQUENCY  - MY OWN THRESHOLD. 
    init_buffer_manager(dataBuffer, BUFFER_SIZE, WINDOW_SIZE, start_index, 11300.0, 1.0);                                                    //START INDEX. START FREQUENCY - MY OWN THRESHOLD OLMALI. PEAK BURADAN HESAPLANACAK MESELA 11300 ISE PEAK 110'da CIKTIYSA 11410 olacak.

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
void startSlidingWindowAnalysis(MesSweep_t *sweep, const double* phaseAngles, uint16_t phase_angle_size, Callback_t callback) { //ILK MES START OLARAK GOREV ALABILIR. 
    ctx.phaseAngles = phaseAngles;      //GEREK YOK.
    ctx.phase_angle_size = phase_angle_size; //GEREK YOK. 
    ctx.callback = callback;
    ctx.isTruncatedLeft = false;
    ctx.isTruncatedRight = false;

    // Initialize the buffer manager
    initBufferManager(sweep->data);
    
    	// detect_significant_gradient_trends, determine_trend_direction
	init_cubic_rls_analysis_parameters(                                                                                                      
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
        1,    // max_decrease_count: Maximum number of consecutive negative gradients allowed
              // when tracking an increasing trend. Allows for minor fluctuations without discarding a valid trend.
        1,    // max_increase_count: Maximum number of consecutive positive gradients allowed
              // when tracking a decreasing trend. Helps to filter out noise when identifying a decreasing trend.
        5,    // min_trend_count: Minimum number of consistent trends required to consider a peak valid.
              // Ensures that only sustained trends are analyzed.
        1     // allowable_inconsistency_count: Allowable number of inconsistencies in trend detection.
              // Permits minor deviations without discarding the trend.
    );
    
    // Initialize the on-peak analysis parameters for average gradient thresholds and consistent trend count
   // // INCREASELERIN BŞALADIKLARI YERDEN YAP BUNU?!
    init_on_peak_analysis_parameters( 
        0.8,  // min_avg_increase: Minimum average increase required to consider a significant increasing trend. LEFT SIDE OF THE PEAK            
               // Ensures that only trends with substantial average gradient are flagged as significant.
        -0.19, // min_avg_decrease: Minimum average decrease required to consider a significant decreasing trend. RIGHT SIDE OF THE PEAK          
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
    // Set up buffer update info instead of directly loading the buffer
    buffer_update_info.needs_update = true;
    buffer_update_info.phase_index_start = buffer_manager.current_phase_index;
    buffer_update_info.buffer_start_index = buffer_manager.current_buffer_index;
    buffer_update_info.move_amount = buffer_manager.window_size;

    currentStatus.isSweepRequested = true;
    STATE_FUNCS[SWP_INITIAL_ANALYSIS].isComplete = true;

    // After setting up, we can perform the buffer update
    SweepSampleCb(ctx.phaseAngles);
}

/**
 * @brief Entry function for the SWP_SEGMENT_ANALYSIS state.
 *
 * This function performs segment analysis on the current window of data to determine the direction of movement.
 * If the analysis is stuck alternating between left and right shifts, it raises an error and transitions to
 * the `SWP_WAITING` state to prevent infinite looping.
 */
static void OnEntrySegmentAnalysis(void) {
    float forgetting_factor = 0.5f;

    SegmentAnalysisResult result = segment_trend_and_concavity_analysis(
        &buffer_manager.buffer[buffer_manager.current_buffer_index],
        buffer_manager.window_size,
        forgetting_factor
    );

    ctx.direction = result.nextDirection;

    // Update the shift tracker with the new direction
    update_shift_tracker(ctx.direction);

    // Check if we are on the peak
    if (ctx.direction == ON_PEAK) {
        currentStatus.isPeakFound = 1;
    } else if (ctx.direction == UNDECIDED || ctx.direction == NEGATIVE_UNDECIDED) {
        currentStatus.isUndecided = 1;  // Raise undecided flag for both cases
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
    // Set up buffer update info
    update_buffer_for_direction(ctx.phaseAngles, ctx.direction, buffer_manager.window_size / 2);

    // After setting up, perform the buffer update if needed
    SweepSampleCb(ctx.phaseAngles);                                                                                              

    STATE_FUNCS[SWP_UPDATE_BUFFER_DIRECTION].isComplete = true;
    SwpProcessStateChange();  // Move to next state
}

/**
 * @brief Entry function for the SWP_UNDECIDED_TREND_CASE state.
 *
 * This function handles the case when the segment analysis result is undecided by moving the window forward.
 * If it's a negative undecided case, it moves the window backward instead.
 */
static void OnEntryUndecidedTrendCase(void) {
    
    if (undecided_case_counter >= MAX_UNDECIDED_CASE_ATTEMPTS) {
        // Set the error flag and transition to SWP_WAITING
        printf("Undecided case limit exceeded. Setting error flag and returning to SWP_WAITING.\n");
        undecided_error_flag = true;
        currentState = SWP_WAITING;
        SwpProcessStateChange();
        return;
    }
    
    if (ctx.direction == NEGATIVE_UNDECIDED) {
        handle_negative_undecided_case(ctx.phaseAngles, ctx.phase_angle_size);  // Move the window backward
    } else {
        handle_undecided_case(ctx.phaseAngles, ctx.phase_angle_size);  // Move the window forward
    }

    // After setting up, perform the buffer update if needed
    SweepSampleCb(ctx.phaseAngles);                                                                                                  

    STATE_FUNCS[SWP_UNDECIDED_TREND_CASE].isComplete = true;
    SwpProcessStateChange();
}


/**
 * @brief Adjusts the global forgetting factor based on gradient sums to improve peak centering.
 *
 * This function calculates the total second-order gradient sum within the buffer window and 
 * adjusts the global centering forgetting factor to better center the peak. It first checks 
 * whether increasing or decreasing the forgetting factor results in a more favorable gradient 
 * sum. The function then modifies the forgetting factor iteratively to emphasize newer or 
 * older data points, improving the accuracy of the peak centering process.
 *
 * ### Parameters:
 * @param buffer Pointer to the data buffer containing the phase angle values.
 * @param buffer_size Size of the data buffer.
 * @param start_index Starting index in the buffer for the gradient calculation.
 * @param previous_gradient_sum Pointer to the previous total gradient sum for comparison.
 * @param centering_attempts Current number of centering attempts.
 *
 * @return double The updated global forgetting factor after the adjustment.
 *
 * ### Intention:
 * - The goal of this function is to adjust the forgetting factor used during peak centering
 *   to achieve better alignment of the data window with the actual peak in the signal.
 * - It does this by comparing the total sum of second-order gradients (a measure of the curvature
 *   of the data) across different forgetting factor values.
 * - If increasing or decreasing the forgetting factor improves the total gradient sum (indicating 
 *   better symmetry around the peak), the function adjusts the global forgetting factor accordingly.
 * - The function is designed to handle both symmetric and asymmetric peaks, refining the 
 *   centering process over multiple attempts.
 *
 * ### Handling Peak Asymmetry:
 * - In real-world signals, peaks can often be **asymmetric** due to noise or inherent system behavior.
 * - **Increasing** the forgetting factor helps correct asymmetry by giving more weight to recent data points.
 * - **Decreasing** the forgetting factor allows the algorithm to smooth the influence of newer points
 *   and gives more emphasis to older data, which helps when the peak asymmetry lies in earlier data.
 * - The function dynamically adjusts the forgetting factor to **adapt** to the specific peak characteristics
 *   and asymmetries, refining the centering process until the second-order gradients are minimized.
 * 
 * - The iterative adjustment helps ensure that the peak is centered by minimizing the total sum 
 *   of the second-order gradients, effectively handling any peak asymmetry.
 */
double adjust_forgetting_factor(const double* buffer, uint16_t buffer_size, uint16_t start_index, 
                                double* previous_gradient_sum, 
                                uint8_t centering_attempts) 
{
    // Get the total gradient sum with the global centering_forgetting_factor
    double total_gradient_sum = compute_total_second_order_gradient(buffer, buffer_size, start_index, centering_forgetting_factor);

    // Check if this is the first centering attempt
    if (*previous_gradient_sum == -1.0) {
        // First attempt, no comparison to previous value
        *previous_gradient_sum = total_gradient_sum;
    } else if (centering_attempts == 1) {
        // Second attempt, try both increasing and decreasing the forgetting factor

        // Store the current total_gradient_sum and forgetting factor
        double increased_forgetting_factor = centering_forgetting_factor + 0.2;
        double decreased_forgetting_factor = centering_forgetting_factor - 0.2;

        // Compute the total gradient sum for increased forgetting factor
        double increased_gradient_sum = compute_total_second_order_gradient(
            buffer, buffer_size, start_index, increased_forgetting_factor
        );

        // Compute the total gradient sum for decreased forgetting factor
        double decreased_gradient_sum = compute_total_second_order_gradient(
            buffer, buffer_size, start_index, decreased_forgetting_factor
        );

        // Compare which one is better
        if (fabs(increased_gradient_sum) < fabs(*previous_gradient_sum)) {
            // Use increased forgetting factor if it's better
            centering_forgetting_factor = increased_forgetting_factor;
            total_gradient_sum = increased_gradient_sum;
            printf("Increased forgetting factor improved centering: %.2f\n", centering_forgetting_factor);
        } else if (fabs(decreased_gradient_sum) < fabs(*previous_gradient_sum)) {
            // Use decreased forgetting factor if it's better
            centering_forgetting_factor = decreased_forgetting_factor;
            total_gradient_sum = decreased_gradient_sum;
            printf("Decreased forgetting factor improved centering: %.2f\n", centering_forgetting_factor);
        } else {
            // If neither is better, stick with the current one
            printf("Neither increased nor decreased forgetting factor improved centering, keeping it: %.2f\n", centering_forgetting_factor);
        }
    } else {
        // For subsequent attempts, compare the current and previous gradient sums
        if (fabs(total_gradient_sum) > fabs(*previous_gradient_sum)) {
            // If the new gradient sum is worse, decrease the forgetting factor
            centering_forgetting_factor -= 0.2;
            if (centering_forgetting_factor < 0.1) {
                // Ensure the factor does not go below a certain minimum (e.g., 0.1)
                centering_forgetting_factor = 0.1;
            }
            printf("Centering worsened, reducing forgetting factor to: %.2f\n", centering_forgetting_factor);
        } else {
            // If the new gradient sum is better, increase the forgetting factor
            centering_forgetting_factor += 0.2;
            printf("Centering improved, increasing forgetting factor to: %.2f\n", centering_forgetting_factor);
        }
    }

    // Update the previous gradient sum with the current one for the next comparison
    *previous_gradient_sum = total_gradient_sum;

    return centering_forgetting_factor;  // Return the updated forgetting factor
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
    
     // Add this at the top of the function (initialize a static variable to hold the previous sum)
    static double previous_gradient_sum = -1.0;
    
    // Reset the isComplete flag for this state
    STATE_FUNCS[SWP_PEAK_CENTERING].isComplete = false;

    // Reset the isNotCentered flag
    currentStatus.isNotCentered = 0;
    
    // Check if we've exceeded the max number of centering attempts
    if (centering_attempts >= MAX_CENTERING_ATTEMPTS) {
        printf("Max centering attempts reached. Transitioning to WAITING state.\n");
        currentStatus.isSweepDone = 1; // Mark the sweep as done
        STATE_FUNCS[SWP_PEAK_CENTERING].isComplete = true;
        SwpProcessStateChange();
        return;
    }

    // Perform the rest of the peak centering logic as before
    double total_gradient_sum = compute_total_second_order_gradient(
        buffer_manager.buffer,
        buffer_manager.buffer_size,
        buffer_manager.current_buffer_index,
        centering_forgetting_factor // Use the updated global forgetting factor
    );
    
        // Adjust the global forgetting factor based on the total gradient sum
    adjust_forgetting_factor(
        buffer_manager.buffer,
        buffer_manager.buffer_size,
        buffer_manager.current_buffer_index,
        &previous_gradient_sum,
        centering_attempts
    );
    
    // Increment the centering attempts counter
    centering_attempts++;

    // Check if the peak is centered based on the total gradient sum
    if (total_gradient_sum <= quadratic_analysis_params.centered_gradient_sum &&
        total_gradient_sum >= -quadratic_analysis_params.centered_gradient_sum) {

        currentStatus.isCentered = 1;
        printf("Peak is centered based on total_gradient_sum.\n");
        STATE_FUNCS[SWP_PEAK_CENTERING].isComplete = true;

    } else {
        // Need to adjust the buffer to center the peak
        GradientTrendResult gradient_trends = track_gradient_trends_with_quadratic_regression(
            buffer_manager.buffer,
            buffer_manager.buffer_size,
            buffer_manager.current_buffer_index,
            buffer_manager.window_size,
            centering_forgetting_factor  // Use the global forgetting factor
        );

        // Check if both increasing and decreasing trends are valid
        if (!gradient_trends.increase_info.valid || !gradient_trends.decrease_info.valid) {
            printf("Invalid trend data. Cannot proceed with centering.\n");
            currentStatus.isCentered = 1; // Consider it centered for now to prevent infinite loops
            STATE_FUNCS[SWP_PEAK_CENTERING].isComplete = true;
        } else {
            // Calculate durations of increasing and decreasing trends
            uint16_t increase_start = gradient_trends.increase_info.start_index;
            uint16_t increase_end = gradient_trends.increase_info.end_index;
            uint16_t decrease_start = gradient_trends.decrease_info.start_index;
            uint16_t decrease_end = gradient_trends.decrease_info.end_index;

            // Adjust indices relative to the buffer size to handle wrap-around
            uint16_t increase_duration = (increase_end + buffer_manager.buffer_size - increase_start) % buffer_manager.buffer_size;
            uint16_t decrease_duration = (decrease_end + buffer_manager.buffer_size - decrease_start) % buffer_manager.buffer_size;

            int shift_amount = 0;
            int direction = UNDECIDED; // Initialize direction

            if (increase_duration > decrease_duration) {
                // Peak is to the left; shift window to the right
                shift_amount = (increase_duration - decrease_duration) / 2;
                direction = RIGHT_SIDE;
                printf("Increase duration (%u) > decrease duration (%u). Moving right by %d.\n",
                       increase_duration, decrease_duration, shift_amount);

            } else if (decrease_duration > increase_duration) {
                // Peak is to the right; shift window to the left
                shift_amount = (decrease_duration - increase_duration) / 2;
                direction = LEFT_SIDE;
                printf("Decrease duration (%u) > increase duration (%u). Moving left by %d.\n",
                       decrease_duration, increase_duration, shift_amount);
            } else {
                // Durations are equal; peak is centered
                currentStatus.isCentered = 1;
                printf("Increase and decrease durations are equal. Peak is centered.\n");
                STATE_FUNCS[SWP_PEAK_CENTERING].isComplete = true;
            }

            if (shift_amount == 0) {
                // No need to shift; peak is centered
                currentStatus.isCentered = 1;
                printf("No shift needed. Peak is centered.\n");
                STATE_FUNCS[SWP_PEAK_CENTERING].isComplete = true;
        
                // Proceed to next state
                SwpProcessStateChange();

            } else {
                
                printf("Peak not centered. Attempt %d of %d\n", centering_attempts + 1, MAX_CENTERING_ATTEMPTS);
    
                // Move the window and set up buffer update if necessary
                move_window_and_update_if_needed(ctx.phaseAngles, direction, shift_amount);
        
                // After moving the window, perform buffer update if needed
                SweepSampleCb(ctx.phaseAngles);                                                                                                    

                currentStatus.isCentered = 1;

                // The state will re-enter to check centering again
                STATE_FUNCS[SWP_PEAK_CENTERING].isComplete = true;
                
                SwpProcessStateChange();
            }
        }
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
static void OnEntryPeakFindingAnalysis(void) {                                                                                                 // yok 
    // Print the total second-order gradient sum again to verify the peak
    double total_gradient_sum = compute_total_second_order_gradient(
        buffer_manager.buffer,
        buffer_manager.buffer_size,
        buffer_manager.current_buffer_index,
        0.5 // Forgetting factor
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
                // If the peak is successfully centered, reset the counter and proceed
                centering_attempts = 0;
                centering_forgetting_factor = 0.7; // Reset the forgetting factor
                
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
    // Since this state may involve multiple transitions based on truncation handling, 
    // we keep isComplete = false initially
    STATE_FUNCS[SWP_PEAK_TRUNCATION_HANDLING].isComplete = false; 

    if (ctx.isTruncatedLeft) {
        printf("Handling truncation on the left side.\n");
        move_window_and_update_if_needed(ctx.phaseAngles, LEFT_SIDE, 5);

        // Perform buffer update if needed
        SweepSampleCb(ctx.phaseAngles);
    }

    if (ctx.isTruncatedRight) {
        printf("Handling truncation on the right side.\n");
        move_window_and_update_if_needed(ctx.phaseAngles, RIGHT_SIDE, 5);

        // Perform buffer update if needed
        SweepSampleCb(ctx.phaseAngles);
    }

    // Don't mark the state complete here since we will evaluate the outcome in OnExit
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

/**
 * @brief Exit function for the SWP_PEAK_TRUNCATION_HANDLING state.
 *
 * This function handles the peak verification after truncation and buffer update are done.
 */
static void OnExitPeakTruncationHandling(void) {
    // Verify if the peak was successfully centered after truncation
    bool peak_verified_after_truncation = verify_peak_at_index(buffer_manager.current_buffer_index);
    if (peak_verified_after_truncation) {
        printf("Peak verification successful after truncation handling, peak is centered.\n");
        print_analysis_interval(ctx.phaseAngles, ctx.phase_angle_size);  // Print the buffer interval
        currentStatus.isSweepDone = 1;  // Mark the sweep as done
        // Now the state is ready to transition, mark it complete

    } else {
        // If peak is not verified, return to peak finding analysis
        printf("Peak verification failed after truncation handling, returning to peak finding analysis.\n");
        currentStatus.isVerificationFailed = 1;  // Set flag to indicate verification failed
       
    }
    
    STATE_FUNCS[SWP_PEAK_TRUNCATION_HANDLING].isComplete = true;
}

/******************************************************************************/
/* State Machine Process and Transitions */
/******************************************************************************/

/**
 * @brief Processes state changes in the sliding window analysis state machine.
 *
 * This function manages the state transitions in the sliding window analysis state machine.
 * It ensures that the appropriate OnExit and OnEntry functions are called during state transitions.
 * The state machine continues processing until no further state changes occur.
 * 
 * ### Structure and Purpose:
 * - The function uses a loop to repeatedly check if a state transition is needed.
 * - Each state is represented by a set of functions: `OnEntry`, `OnExit`, and a `isComplete` flag that indicates
 *   whether the state has finished its processing.
 * - State transitions are triggered based on the `NextState` function, which determines the appropriate state to move to
 *   after completing the current one.
 *
 * ### Key Concepts:
 * - **State Completion (`isComplete`)**: Each state has an `isComplete` flag to indicate whether the state's processing
 *   is done. The state machine will only transition to the next state if this flag is set. If the state is not marked
 *   as complete, the next state is often determined by the logic inside the `OnExit()` function of the current state,
 *   allowing flexibility for more complex state transitions that are dependent on the exit logic of a state.
 * 
 * - **OnExit/OnEntry Functions**: 
 *   - When a state is marked as incomplete (`!isComplete`), the `OnExit()` function (if defined) is called to handle
 *     cleanup, and potentially dictate what the next state will be.
 *   - Once a state transition is made, the `OnEntry()` function for the new state is invoked to initialize the new state's processing.
 *
 * - **State Transition (`NextState`)**: The `NextState()` function determines the next state to transition into
 *   based on the current state and flags. This function, in combination with the `OnExit()` function, ensures the
 *   correct state is selected, especially for states that require further processing or specific conditions before moving forward.
 * 
 * ### Intention:
 * - The function aims to provide flexibility in state transitions by allowing states to either fully complete and move on
 *   to the next state, or re-enter themselves if they require further processing, such as in asynchronous operations.
 * - The `OnExit()` logic allows states to determine the appropriate next state during exit, particularly useful in scenarios
 *   where a state's exit depends on external factors or complex logic.
 * - The loop ensures the state machine processes transitions in sequence without breaking the flow, making it adaptable for real-time or multi-step operations.
 * 
 * ### Working Mechanism:
 * 1. **State Completion Check**:
 *    - If the current state is incomplete (`!STATE_FUNCS[currentState].isComplete`), the function calls the `OnExit()` function (if defined)
 *      to handle cleanup and determine any next steps.
 *    - The `NextState()` function is called to determine the next state, and `stateChanged` tracks if the state has changed.
 *    - This allows a state to finalize any pending tasks in its `OnExit()` before transitioning to the next one.
 * 
 * 2. **State Transition**:
 *    - After calling `OnExit()`, the `isComplete` flag for the new state is reset, and the `OnEntry()` function for the new state is invoked.
 *    - The `currentStatus.value` is reset to ensure no status flags carry over to the next state.
 * 
 * 3. **Complete State Handling**:
 *    - If a state is marked as complete (`STATE_FUNCS[currentState].isComplete`), the function checks for the next state and processes it.
 *    - The loop continues processing state transitions until no more state changes occur.
 * 
 * 4. **Callback Handling**:
 *    - If the state machine reaches the `SWP_WAITING` state and a callback is set (`ctx.callback`), the callback is executed and cleared.
 * 
 * ### Example Use Case:
 * The function is suitable for scenarios where:
 * - Some states need multiple passes to complete their processing (e.g., waiting for an asynchronous event).
 * - Transitions between states depend on conditions evaluated during the `OnExit()` function.
 * - A callback function is executed after the state machine reaches a specific terminal state (e.g., `SWP_WAITING`).
 * 
 * ### Flexibility:
 * The function handles both synchronous and asynchronous transitions. States that require multiple iterations
 * to complete can use the `OnExit()` logic to determine the next state dynamically. Meanwhile, states that do not
 * require re-entry can transition immediately.
 * 
 * @note The state machine continues processing until no further state changes occur, ensuring all transitions are properly handled.
 *
 */
void SwpProcessStateChange(void) {
    bool stateChanged;

    // The state machine should process until no further state changes occur
    do {
        // If the current state is incomplete, handle its exit
        if (!STATE_FUNCS[currentState].isComplete) {
            // Ensure this debug message is only printed during truncation handling or incomplete states
         
            // Execute OnExit if it's not complete
            if (STATE_FUNCS[currentState].onExit != NULL) {
                STATE_FUNCS[currentState].onExit();
            }

            // Determine the next state after OnExit
            SwpState_t nextState = NextState(currentState);
            stateChanged = (nextState != currentState);

            // Reset the currentStatus before transitioning
            currentStatus.value = 0x0;

            if (stateChanged) {
                
                // Ensure the old state's completion flag is set to false
                STATE_FUNCS[currentState].isComplete = false;
                // Move to the new state
                currentState = nextState;

                STATE_FUNCS[currentState].isComplete = false;

                // Call the new state's OnEntry function
                if (STATE_FUNCS[currentState].onEntry != NULL) {
                    STATE_FUNCS[currentState].onEntry();
                }
            }
        }

        // Process the state if it is marked complete
        if (STATE_FUNCS[currentState].isComplete) {
            //printf("[DEBUG] State %d is complete. Proceeding to next state...\n", currentState);

            // Determine the next state
            SwpState_t nextState = NextState(currentState);
            stateChanged = (nextState != currentState);
            
            // Clear currentStatus flags
            currentStatus.value = 0x0;

            if (stateChanged) {
                // Transition to the new state
                currentState = nextState;

                // Ensure the new state's completion flag is set to false
                STATE_FUNCS[currentState].isComplete = false;

                // Call the OnEntry function for the new state
                if (STATE_FUNCS[currentState].onEntry != NULL) {
                    STATE_FUNCS[currentState].onEntry();
                }
            }
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
 * This function is used to manage state transitions in the sliding window analysis process. It checks for various flags
 * such as boundary errors, sweep completion, and centering conditions to determine the appropriate next state.
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
            if (currentStatus.isSweepDone) {
                return SWP_WAITING;  // Transition to WAITING state if sweep is done
            }
        return SWP_UPDATE_BUFFER_DIRECTION;

        case SWP_UPDATE_BUFFER_DIRECTION:
            if (currentStatus.isBoundaryError) {
                return SWP_WAITING;  // Transition to waiting if buffer boundary error occurred
            }
            return SWP_SEGMENT_ANALYSIS;

        case SWP_PEAK_CENTERING:
            if (currentStatus.isBoundaryError) {
                return SWP_WAITING;  // Transition to waiting if buffer boundary error occurred
            }
            if (currentStatus.isCentered) {
                return SWP_PEAK_FINDING_ANALYSIS;
            }
            if (currentStatus.isSweepDone) {
                return SWP_WAITING;
            }
            return SWP_PEAK_CENTERING;

        case SWP_UNDECIDED_TREND_CASE:
            // Check if undecided case counter has reached the limit
            if (undecided_case_counter >= MAX_UNDECIDED_CASE_ATTEMPTS) {
                return SWP_WAITING;  // Transition to SWP_WAITING if limit exceeded
            }
            if (currentStatus.isBoundaryError) {
                return SWP_WAITING;  // Transition to waiting if buffer boundary error occurred
            }
            return SWP_SEGMENT_ANALYSIS;

        case SWP_PEAK_FINDING_ANALYSIS:
            if (currentStatus.isBoundaryError) {
                return SWP_WAITING;  // Transition to waiting if buffer boundary error occurred
            }
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
            if (currentStatus.isBoundaryError) {
                return SWP_WAITING;  // Transition to waiting if buffer boundary error occurred
            }
            if (currentStatus.isSweepDone) {
                return SWP_WAITING;
            }
            if (currentStatus.isVerificationFailed) {
                return SWP_PEAK_FINDING_ANALYSIS;
            }
            return SWP_PEAK_FINDING_ANALYSIS;

        case SWP_WAITING:
            if (currentStatus.isSweepRequested) {
                return SWP_INITIAL_ANALYSIS;
            }
            return SWP_WAITING;

        default:
            return SWP_WAITING;
    }
}

/**
 * @brief Sets the boundary error flag in the currentStatus.
 *
 * This function allows external modules (like buffer_manager) to set the boundary error flag in the sliding window
 * analysis state machine, which will trigger a state transition to SWP_WAITING.
 *
 * @param flag The value to set for the boundary error flag (1 for error, 0 to clear).
 */
// Function to set the boundary error flag
void set_boundary_error_flag(uint8_t flag) {
    currentStatus.isBoundaryError = flag;
    if (flag) {
        boundaryErrorOccurred = true; // Set the persistent boundary error flag
    }
}


//9.1

//139, 180 hata veriyor. 