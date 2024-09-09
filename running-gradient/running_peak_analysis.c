#include <stdio.h>
#include <time.h>
#include <stdlib.h>

#include "windowed_running_gradient.h"
#include "running_quadratic_gradient.h"
#include "running_peak_analysis.h"
#include "running_cubic_gradient.h"
#include "rls_analysis_parameters.h"


	 // Prepare a buffer for analysis (using 200 data points as an example)
    MqsRawDataPoint_t buffer[BUFFER_SIZE];

/**
 * @brief Determines the directional trend within a given window based on the analysis of peak trends.
 *
 * This function is designed to analyze trends within a specific window of data and determine the direction
 * in which the data is moving. The analysis is focused on identifying significant increasing or decreasing trends,
 * as well as detecting consistent trends that, while not significant, suggest a directional movement over a period.
 *
 * Based on these criteria, the function sets flags that indicate the direction of movement (left or right)
 * and the proximity to a peak. The flags help in determining whether the analysis is close to identifying
 * a peak or still far from one.
 *
 * ### Function Workflow:
 * 1. **Initialization**: The function starts by initializing the `TrendDirectionFlags` structure with default values.
 *    The `far_to_peak` flag is set to `true` by default, indicating an assumption that the data is far from a peak.
 * 2. **Significant Trend Detection**: The function first checks for significant trends. If a significant increase is detected,
 *    particularly towards the right side of the window, it sets the `move_to_right` flag, indicating the need to shift the
 *    analysis to the right. Similarly, if a significant decrease is detected towards the left, it sets the `go_to_left` flag.
 *    In both cases, the `close_to_peak` flag is set to `true`, and `far_to_peak` is set to `false`.
 * 3. **Low Consistent Trend Detection**: If no significant trends are detected, the function then checks for low consistent trends.
 *    If the data shows a consistent increase or decrease over more than 10 data points, it sets the corresponding
 *    `low_consistent_right` or `low_consistent_left` flags. This step helps in identifying gradual but steady trends
 *    that might not be as pronounced as significant trends but still indicate a directional movement.
 * 4. **Finalization**: After analyzing the trends, the function returns the `TrendDirectionFlags` structure containing
 *    the flags that indicate the detected trends and the recommended direction for further analysis.
 *
 * ### Function Purpose:
 * The primary goal of the `determine_trend_direction` function is to guide further analysis by determining the
 * direction in which the data is trending within a given window. By identifying significant and low consistent trends,
 * the function helps in deciding whether to move the analysis window to the left or right or to stay within the current window.
 * This is particularly useful in real-time data analysis, where the identification of trends and peaks is crucial
 * for making informed decisions based on the data's behavior over time.
 *
 * ### Parameters:
 * @param trends A pointer to a `PeakTrendAnalysisResult` structure that contains detailed information about
 *        the detected significant increase and decrease trends within the data.
 * @param window_cubic The size of the window used for cubic regression analysis, representing the number of data points
 *        considered for trend detection within the window.
 * @param start_cubic_index The starting index of the window in the data array, indicating the position within the
 *        overall data set where the current window begins.
 *
 */
static TrendDirectionFlags determine_trend_direction(const PeakTrendAnalysisResult *trends, uint16_t window_cubic, uint16_t start_cubic_index) {
    TrendDirectionFlags flags = {false, false, false, true, false, false, false}; // Initialize with on_the_peak

    printf("Window start index: %u, Window size: %u\n", start_cubic_index, window_cubic);

    double average_increase = 0.0;
    double average_decrease = 0.0;
    uint16_t increase_duration = 0;
    uint16_t decrease_duration = 0;

    // Check for significant increasing trend
    if (trends->significant_increase) { //SIGNIFICANT INCREASE THRESHOLDS SHOULD BE CHECKED 
        increase_duration = trends->increase_info.end_index - trends->increase_info.start_index;
        average_increase = trends->increase_info.max_sum / (double)increase_duration;

        printf("Significant increase detected. Increase end index: %u\n", trends->increase_info.end_index);
        printf("Average Increase: %.6f over interval [%u - %u]\n", average_increase, trends->increase_info.start_index, trends->increase_info.end_index);

        // Use the parameters from peak_analysis_params instead of hardcoded values
        if (increase_duration > peak_analysis_params.min_consistent_trend_count && average_increase > peak_analysis_params.min_average_increase) {
            flags.move_to_right = true;
            flags.close_to_peak = true;
            flags.far_to_peak = false;
            printf("Average increase over threshold. Flags set: move_to_right, close_to_peak.\n");
        }

        // Also check if the significant increase is towards the right side of the window
        if (trends->increase_info.end_index > start_cubic_index + (window_cubic / 2)) {
            if (trends->increase_info.end_index >= start_cubic_index + window_cubic - 1) {
                flags.move_to_right = true;
                flags.close_to_peak = true;
                flags.far_to_peak = false;
                printf("Significant increasing trend extends beyond the window. Flags set: move_to_right, close_to_peak.\n");
            }
        }
    }

    // Check for significant decreasing trend
    if (trends->significant_decrease) {
        decrease_duration = trends->decrease_info.end_index - trends->decrease_info.start_index;
        average_decrease = trends->decrease_info.max_sum / (double)decrease_duration;

        printf("Significant decrease detected. Decrease start index: %u\n", trends->decrease_info.start_index);
        printf("Average Decrease: %.6f over interval [%u - %u]\n", average_decrease, trends->decrease_info.start_index, trends->decrease_info.end_index);

        // Use the parameters from peak_analysis_params instead of hardcoded values
        if (decrease_duration > peak_analysis_params.min_consistent_trend_count && average_decrease < peak_analysis_params.min_average_decrease) {
            flags.go_to_left = true;
            flags.close_to_peak = true;
            flags.far_to_peak = false;
            printf("Average decrease below threshold. Flags set: go_to_left, close_to_peak.\n");
        }

        // Also check if the significant decrease is towards the left side of the window
        if (trends->decrease_info.start_index < start_cubic_index + (window_cubic / 2)) {
            flags.go_to_left = true;
            flags.close_to_peak = true;
            flags.far_to_peak = false;
            printf("Significant decreasing trend in the first half of the window. Flags set: go_to_left, close_to_peak.\n");
        }
    }

    // Set on_the_peak flag if both significant increase and decrease are detected
    if (trends->significant_increase && trends->significant_decrease) {
        flags.on_the_peak = true;
        printf("On the peak detected. Flags set: on_the_peak.\n");
    }
    //HATA BURADAN GELIOR. 
    // If no significant trend is detected, check for low consistent trends
    
    /*
    if (!trends->significant_increase && !trends->significant_decrease) {
        increase_duration = trends->increase_info.end_index - trends->increase_info.start_index;
        decrease_duration = trends->decrease_info.end_index - trends->decrease_info.start_index;

        // Use a parameterized value for low consistent trends (replace with appropriate value)
        if (increase_duration > (WINDOW_SIZE/2)) { //SHOULD NOT BE HARDCODED. 
            flags.low_consistent_right = true;
            printf("Low consistent right trend detected over %u indices.\n", increase_duration);
        }

        if (decrease_duration > (WINDOW_SIZE/2)) {
            flags.low_consistent_left = true;
            printf("Low consistent left trend detected over %u indices.\n", decrease_duration);
        }
    }
    */

    printf("Trend direction determination complete.\n");
    return flags;
}

void log_peak_detection(bool on_peak) {
    if (on_peak) {
        printf("--> On the peak detected. Skipping direction analysis and proceeding with concavity analysis.\n");
    }
}

void log_direction_determined(int direction) {
    switch (direction) {
        case RIGHT_SIDE:
            printf("--> Direction determined: RIGHT_SIDE.\n");
            break;
        case LEFT_SIDE:
            printf("--> Direction determined: LEFT_SIDE.\n");
            break;
        case ON_PEAK:
            printf("--> Direction determined: ON_PEAK.\n");
            break;
        case UNDECIDED:
        default:
            printf("--> Direction determined: UNDECIDED.\n");
            break;
    }
}

void log_gradient_comparison(GradientComparisonResult gradient_result) {
    printf("--> Total Gradient First Part: %.6f\n", gradient_result.total_gradient_first_part);
    printf("--> Total Gradient Second Part: %.6f\n", gradient_result.total_gradient_second_part);
    printf("ooo-->Dominant Side: %d\n", gradient_result.dominant_side);

    if (gradient_result.dominant_side == LEFT_SIDE) {
        printf("--> Direction overridden by gradient comparison: LEFT_SIDE.\n");
    } else if (gradient_result.dominant_side == RIGHT_SIDE) {
        printf("--> Direction overridden by gradient comparison: RIGHT_SIDE.\n");
    } else {
        printf("--> Gradient comparison result was UNDECIDED, not overriding the existing direction.\n");
    }
}

void log_significant_trend(const char* type, double max_sum) {
    printf("--> Significant %s trend with max sum %.6f detected, setting ON_PEAK.\n", type, max_sum);
}

void log_concavity_analysis(ConcavityAnalysisOutput concavity_output) {
    if (concavity_output.isNoisy) {
        printf("--> The segment is noisy.\n");
    }
    if (concavity_output.moveLeft) {
        printf("--> Suggested Action: Move LEFT\n");
    }
    if (concavity_output.moveRight) {
        printf("--> Suggested Action: Move RIGHT\n");
    }
    if (concavity_output.stay) {
        printf("--> Suggested Action: Stay in the current position\n");
    }
    if (concavity_output.isPotentialPeak) {
        printf("--> A potential peak is detected.\n");
    }
    if (concavity_output.isTruePeak) {
        printf("--> A true peak is detected.\n");
    }
}

void log_final_direction(int nextDirection) {
    printf("Final Direction: ");
    switch (nextDirection) {
        case RIGHT_SIDE:
            printf("RIGHT_SIDE\n");
            break;
        case LEFT_SIDE:
            printf("LEFT_SIDE\n");
            break;
        case ON_PEAK:
            printf("ON_PEAK\n");
            break;
        case UNDECIDED:
        default:
            printf("UNDECIDED\n");
            break;
    }
}

/**
 * @brief Analyzes the trend and concavity within a window of data to determine the direction of movement.
 *
 * This function will now work on a subset of the buffer, starting from the provided window start index.
 *
 * @param data The pointer to the start of the window (within the larger buffer).
 * @param window_size The number of data points in the window (subset of the buffer).
 * @param forgetting_factor A double value typically close to 1 that determines the weight given to newer data points in recursive analysis.
 * @return SegmentAnalysisResult A structure containing information about the direction to move (left, right, undecided),
 *         whether a potential or true peak was found, and the detected concavity pattern.
 */
SegmentAnalysisResult segment_trend_and_concavity_analysis(const MqsRawDataPoint_t *data, uint16_t window_size, double forgetting_factor) {
    SegmentAnalysisResult result = {
        .isPotentialPeak = false,
        .isTruePeak = false,
        .nextDirection = UNDECIDED,
        .concavityOutput = {
            .isNoisy = false,
            .moveLeft = false,
            .moveRight = false,
            .stay = false,
            .isPotentialPeak = false,
            .isTruePeak = false
        }
    };

    printf("****************************************** Cubic regression increase/decrease interval analysis...\n");
    PeakTrendAnalysisResult significant_trends = detect_significant_gradient_trends(data, window_size, 0, 30, forgetting_factor);

    printf("****************************************** Determining cubic regression linear gradient sums...\n");
    TrendDirectionFlags direction_flags = determine_trend_direction(&significant_trends, window_size, 0);

    if (direction_flags.on_the_peak) {
        result.nextDirection = ON_PEAK;
        log_peak_detection(true);
        goto concavity_analysis;
    }

    if (direction_flags.move_to_right) {
        result.nextDirection = RIGHT_SIDE;
        log_direction_determined(RIGHT_SIDE);
    } else if (direction_flags.go_to_left) {
        result.nextDirection = LEFT_SIDE;
        log_direction_determined(LEFT_SIDE);
    } else if (direction_flags.low_consistent_right) {
        result.nextDirection = RIGHT_SIDE;
        log_direction_determined(RIGHT_SIDE);
    } else if (direction_flags.low_consistent_left) {
        result.nextDirection = LEFT_SIDE;
        log_direction_determined(LEFT_SIDE);
    }

    if (direction_flags.far_to_peak) {
        printf("|||--> Far from peak. Determining direction by comparing gradient parts...\n");
        GradientComparisonResult gradient_result = compare_gradient_parts(data, 0, forgetting_factor); //BURADAKI THRESHOLDLARA DA BAKMAMIZ LAZIM.
        log_gradient_comparison(gradient_result);

        if (gradient_result.dominant_side != UNDECIDED) {
            result.nextDirection = gradient_result.dominant_side;
        }
    } else {
concavity_analysis: 
        printf("****************************************** Performing concavity analysis...\n");

        GradientTrendResult gradient_trends = track_gradient_trends_with_quadratic_regression(data, window_size, 0, 30, forgetting_factor);

        if (gradient_trends.increase_info.valid && gradient_trends.increase_info.max_sum > 2.5) {
            log_significant_trend("increasing", gradient_trends.increase_info.max_sum);
            result.nextDirection = ON_PEAK;
            goto end_analysis;
        }
        
        if (gradient_trends.decrease_info.valid && gradient_trends.decrease_info.max_sum < -2.5) {
            log_significant_trend("decreasing", gradient_trends.decrease_info.max_sum);
            result.nextDirection = ON_PEAK;
            goto end_analysis;
        }

        printf("*********************** CONCAVITY SEGMENT ANALYSIS ************************\n");
        //ConcavityAnalysisResult concavity_result = initial_concavity_analysis(data, window_size, 0, forgetting_factor, true); //NOT NECESSARY ANYMORE
        //ConcavityAnalysisOutput concavity_output = analyze_concavity_segments(&concavity_result);
        //log_concavity_analysis(concavity_output);
    }

end_analysis:
    log_final_direction(result.nextDirection);
    return result;
}

BufferManager buffer_manager;

/**
 * @brief Initializes the BufferManager struct with starting values.
 *
 * This function sets up the buffer manager by initializing its buffer, buffer size, window size, 
 * current phase index, and frequency parameters. The buffer index is set to the middle of the buffer.
 *
 * @param buffer A pointer to the phase data buffer to be managed.
 * @param buffer_size Size of the buffer.
 * @param window_size Size of the sliding window.
 * @param start_index Starting index for phase data analysis.
 * @param start_frequency The starting frequency for the phase index.
 * @param frequency_increment Increment value to calculate frequency for each phase index.
 *
 * @note This function prints out debugging information, including the initial state of the buffer manager.
 */
void init_buffer_manager(MqsRawDataPoint_t* buffer, uint16_t buffer_size, uint16_t window_size, int16_t start_index, double start_frequency, double frequency_increment) {
    buffer_manager.buffer = buffer;
    buffer_manager.buffer_size = buffer_size;
    buffer_manager.window_size = window_size;
    buffer_manager.current_phase_index = start_index;
    buffer_manager.current_buffer_index = buffer_size / 2; // Start from the middle of the buffer
    
    buffer_manager.mapping_start_index = start_index;
    buffer_manager.mapping_end_index = start_index + window_size - 1;
    
    buffer_manager.start_frequency = start_frequency;
    buffer_manager.frequency_increment = frequency_increment;

    // Debugging: Print initialization state
    printf("Buffer Manager initialized:\n");
    printf("Buffer size: %d, Window size: %d, Start index: %d\n", buffer_size, window_size, start_index);
    printf("Start frequency: %.2f Hz, Frequency increment: %.2f Hz\n", start_frequency, frequency_increment);
}


/**
 * @brief Loads the initial window of phase angle data into the buffer.
 *
 * @param phaseAngles      Pointer to the array of phase angles.
 * @param phase_angle_size Size of the phase angle array.
 */
void load_initial_buffer(const double* phaseAngles, uint16_t phase_angle_size) {
    printf("Loading initial window of data into the buffer.\n");
    
    for (uint16_t i = 0; i < buffer_manager.window_size; ++i) {
        int16_t buffer_index = (buffer_manager.current_buffer_index + i) % buffer_manager.buffer_size;
        int16_t phase_index = buffer_manager.current_phase_index + i;

        // Ensure phase_index is within bounds
        if (phase_index >= phase_angle_size) {
            //printf("[load_initial_buffer] Phase index %d out of bounds. Setting phaseAngle to 0.0.\n", phase_index);
            buffer_manager.buffer[buffer_index].phaseAngle = 0.0;
        } else {
            buffer_manager.buffer[buffer_index].phaseAngle = phaseAngles[phase_index];
            buffer_manager.buffer[buffer_index].impedance = 0.0; // Initialize impedance as needed
            //printf("[load_initial_buffer] Buffer[%d] set to phaseAngle %.6f\n", buffer_index, phaseAngles[phase_index]);
        }
    }

    printf("Initial buffer loading complete.\n");
}

 
/**
 * @brief Updates the buffer with phaseAngle values and optionally sets impedance.
 *
 * This function maps phase index values to buffer indices and updates the buffer with the corresponding
 * phaseAngle values. It also logs the updates made to the buffer.
 *
 * @param phaseAngles A pointer to the array of phase angles.
 * @param phase_index_start The starting index of the phase angles to update from.
 * @param buffer_start_index The starting index of the buffer where the phase angles will be stored.
 */
void update_phaseAngle_to_buffer(const double* phaseAngles, int16_t phase_index_start, int16_t buffer_start_index) {
    // Loop to update the buffer with new phaseAngle values
    for (uint16_t i = 0; i < buffer_manager.window_size; ++i) {
        int16_t phase_index = phase_index_start + i;
        int16_t buffer_index = (buffer_start_index + i) % buffer_manager.buffer_size;

        // Debugging: Log the details of buffer update before actual update
        printf("[update_phaseAngle_to_buffer] [DEBUG] Before update -> Buffer index: %d, Phase index: %d, PhaseAngle: %.6f\n", buffer_index, phase_index, phaseAngles[phase_index]);

        // Copy phaseAngle and optionally set impedance
        buffer_manager.buffer[buffer_index].phaseAngle = phaseAngles[phase_index];
        buffer_manager.buffer[buffer_index].impedance = 0.0; // Mock impedance, update if needed

        // Debugging: Log the update in the buffer
        //printf("[update_phaseAngle_to_buffer] Buffer updated at phase index %d with phaseAngle %.6f (buffer index: %d)\n", phase_index, phaseAngles[phase_index], buffer_index);
    }
}

/**
 * @brief Updates the buffer with phase angle data based on the movement direction.
 *
 * This function moves the sliding window left or right based on the given direction and updates the buffer
 * with new phase angle values. The phase index is adjusted accordingly. If the direction is left, 
 * the window slides left, and if the direction is right, the window slides right. The updated phase angle 
 * values are mapped to their corresponding buffer index.
 *
 * @param phaseAngles A pointer to the array of phase angles.
 * @param direction An integer representing the direction of movement. Use LEFT_SIDE for left, 
 * RIGHT_SIDE for right, and UNDECIDED for no movement.
 *
 * @note The function performs bounds checking to ensure that the next phase index stays within valid limits.
 */
void update_buffer_for_direction(const double* phaseAngles, int direction) {
    int16_t phase_index_start = buffer_manager.current_phase_index;
    int16_t phase_index_end = phase_index_start + buffer_manager.window_size;
    int16_t buffer_start_index = buffer_manager.current_buffer_index;

    // update_buffer: Log initial indices before movement
    printf("[update_buffer] Initial phase index start: %d, end: %d, buffer start index: %d\n", phase_index_start, phase_index_end, buffer_start_index);

    // Calculate the shift direction (LEFT or RIGHT)
    if (direction == LEFT_SIDE) {
        phase_index_start -= buffer_manager.window_size / 2;
        phase_index_end = phase_index_start + buffer_manager.window_size;
        buffer_start_index -= buffer_manager.window_size / 2;

        if (buffer_start_index < 0) buffer_start_index += buffer_manager.buffer_size;
        printf("[update_buffer] Moving left. New phase index range: [%d - %d], buffer start index: %d\n", phase_index_start, phase_index_end, buffer_start_index);
    } else if (direction == RIGHT_SIDE) {
        phase_index_start += buffer_manager.window_size / 2;
        phase_index_end = phase_index_start + buffer_manager.window_size;
        buffer_start_index += buffer_manager.window_size / 2;

        if (buffer_start_index >= buffer_manager.buffer_size) buffer_start_index -= buffer_manager.buffer_size;
        printf("[update_buffer] Moving right. New phase index range: [%d - %d], buffer start index: %d\n", phase_index_start, phase_index_end, buffer_start_index);
    } else {
        printf("[update_buffer] Direction undecided. No movement.\n");
        return;
    }

    // Ensure phase index stays within bounds of phaseAngles
    if (phase_index_start < 0 || phase_index_end > buffer_manager.buffer_size) {
        printf("[update_buffer] Out of bounds of phase angles (index: [%d - %d]). Operation aborted.\n", phase_index_start, phase_index_end);
        return;
    }

    // Call external function to update buffer with new phase angles
    update_phaseAngle_to_buffer(phaseAngles, phase_index_start, buffer_start_index);

    // Update the current phase and buffer indices
    buffer_manager.current_phase_index = phase_index_start;
    buffer_manager.current_buffer_index = buffer_start_index;

    // Log final indices after update
    printf("[update_buffer] Updated buffer indices. Current phase index: %d, Buffer start index: %d\n", buffer_manager.current_phase_index, buffer_manager.current_buffer_index);
}

/**
 * @brief Handles the case when the segment analysis result is undecided.
 *
 * This function jumps forward in the phaseAngle array by the window size, discards the previous analysis,
 * and reloads the buffer from the new phase index. The new portion of the phaseAngle array is loaded into the buffer,
 * overwriting the previous values in the middle section of the buffer.
 *
 * @param phaseAngles A pointer to the array of phase angles.
 * @param phase_angle_size The total size of the phase angle array.
 */
void handle_undecided_case(const double* phaseAngles, uint16_t phase_angle_size) {
    // Flag to indicate entry into the function
    printf("[handle_undecided_case] Entering undecided case handler.\n");

    // Jump forward by the window size in the phaseAngle array
    buffer_manager.current_phase_index += buffer_manager.window_size;

    // Ensure we don't exceed the phase angle size
    if (buffer_manager.current_phase_index + buffer_manager.window_size > phase_angle_size) {
        printf("[handle_undecided_case] Phase index exceeds phase angle size. Stopping analysis.\n");
        return;
    }

    // Load new phase angle values from the phaseAngle array into the middle section of the buffer
    printf("[handle_undecided_case] Jumping to phase index: %d\n", buffer_manager.current_phase_index);
    update_phaseAngle_to_buffer(phaseAngles, buffer_manager.current_phase_index, buffer_manager.current_buffer_index);

    // Log the buffer after loading new data
    printf("[handle_undecided_case] Loaded new phase angle values into the buffer. Starting analysis from buffer index %d.\n", buffer_manager.current_buffer_index);
}

/**
 * @brief Verifies the detected peak using the find_and_verify_cubic_peak function.
 *
 * This function is triggered when the segment analysis detects a peak (ON_PEAK). It uses the current phase index
 * as the starting point for verifying the peak by calling the find_and_verify_cubic_peak function. If the peak
 * is verified successfully, it logs the result.
 *
 * @param phaseAngles A pointer to the array of phase angles.
 * @param phase_angle_size The size of the phase angle array.
 * @param start_index The starting index for peak verification.
 */
void verify_peak_at_index(const MqsRawDataPoint_t *phaseAngles, uint16_t phase_angle_size, uint16_t start_index) {
    printf("[verify_peak_at_index] Verifying peak at phase index: %d\n", start_index);
    
    // Call the find_and_verify_cubic_peak function
    CubicPeakAnalysisResult peak_result = find_and_verify_cubic_peak(
        phaseAngles,           // Pointer to phase angle data
        phase_angle_size,      // Total size of the phase angle data
        start_index,           // Start index for peak verification
        0.5                    // Forgetting factor for RLS analysis
    );

    if (peak_result.peak_found) {
        printf("[verify_peak_at_index] Verified peak found at index: %d\n", peak_result.peak_index);
    } else {
        printf("[verify_peak_at_index] Peak verification failed at index: %d\n", start_index);
    }
}


/** calculate_frequency
 * @brief Performs sliding window analysis with buffer management and trend detection.
 *
 * This function controls the sliding window logic. It first initializes the buffer with the initial phase data 
 * and then enters a loop that shifts the window based on the direction determined by the segment analysis. 
 * The window continues to slide left or right until a peak is detected or the data bounds are exceeded.
 *
 * @param phaseAngles A pointer to the array of phase angles.
 * @param phase_angle_size The size of the phase angle array.
 *
 * @note Debugging information is printed at every major step, including buffer contents, window updates, 
 * direction determinations, and peak detection.
 */
void perform_sliding_window_analysis(const double* phaseAngles, uint16_t phase_angle_size) {
    int16_t direction = UNDECIDED;

    // Step 1: Load the initial window of values into the buffer
    load_initial_buffer(phaseAngles, phase_angle_size);

    // Step 2: Debugging: Log the buffer manager state
    printf("[DEBUG] Buffer Manager State after initial loading:\n");
    printf("Current buffer index: %d\n", buffer_manager.current_buffer_index);
    printf("Current phase index: %d\n", buffer_manager.current_phase_index);
    printf("Window size: %d\n", buffer_manager.window_size);
    printf("Buffer size: %d\n", buffer_manager.buffer_size);

    // Step 3: Pass the starting point of the buffer (window) to the segment analysis
    SegmentAnalysisResult result = segment_trend_and_concavity_analysis(
        &buffer_manager.buffer[buffer_manager.current_buffer_index],
        buffer_manager.window_size,
        0.5
    );

    // Step 4: Debugging: Log the result of the initial analysis
    printf("[DEBUG] Initial segment analysis result: Direction: %d, Potential peak: %s, True peak: %s\n",
           result.nextDirection,
           result.isPotentialPeak ? "Yes" : "No",
           result.isTruePeak ? "Yes" : "No");

    // Step 5: Determine initial direction based on analysis
    direction = result.nextDirection;

    // Step 6: If a peak is detected right away, log and exit
    if (direction == ON_PEAK) {
        printf("Peak detected immediately at phase index: %d\n", buffer_manager.current_phase_index);
        return; // No need to enter the while loop if peak is already found
    }

    printf("Initial direction determined: %d\n", direction);

    // Step 7: Enter the sliding window logic if peak wasn't found immediately
    bool run_once = false;  // Set to true for single iteration, false for indefinite execution

    while ((direction != ON_PEAK && buffer_manager.current_phase_index >= 0 && buffer_manager.current_phase_index < phase_angle_size)) {

        // Step 8: If the result is undecided, handle it by jumping forward in the phaseAngle array
        if (direction == UNDECIDED) {
            printf("[perform_sliding_window_analysis] Result was undecided, handling undecided case.\n");
            handle_undecided_case(phaseAngles, phase_angle_size);
            result = segment_trend_and_concavity_analysis(
                &buffer_manager.buffer[buffer_manager.current_buffer_index],
                buffer_manager.window_size,
                0.5
            );
            direction = result.nextDirection; // Get new direction after handling undecided case
            continue; // Restart the loop with new data
        }

        // Debugging: Log the current window range before update
        printf("[DEBUG] Current window range in buffer: [%d - %d], Buffer manager current phase index: %d\n",
               buffer_manager.current_buffer_index,
               (buffer_manager.current_buffer_index + buffer_manager.window_size - 1) % buffer_manager.buffer_size,
               buffer_manager.current_phase_index);

        // Step 9: Load the next window of values into the buffer based on the direction
        printf("Updating buffer for direction: %d\n", direction);
        update_buffer_for_direction(phaseAngles, direction);

        // Debugging: Log the buffer contents after the update
        /*
        printf("Buffer contents after updating:\n");
        for (uint16_t i = 0; i < buffer_manager.buffer_size; ++i) {
            printf("Buffer[%d]: phaseAngle = %.6f, impedance = %.6f\n", i, buffer_manager.buffer[i].phaseAngle, buffer_manager.buffer[i].impedance);
        }
        */

        // Debugging: Log updated buffer manager state
        printf("[DEBUG] Updated Buffer Manager State:\n");
        printf("Current buffer index: %d\n", buffer_manager.current_buffer_index);
        printf("Current phase index: %d\n", buffer_manager.current_phase_index);
        printf("Window size: %d\n", buffer_manager.window_size);
        printf("Buffer size: %d\n", buffer_manager.buffer_size);

        // Step 10: Pass the correct window into the analysis (sliding window)
        result = segment_trend_and_concavity_analysis(
            &buffer_manager.buffer[buffer_manager.current_buffer_index],
            buffer_manager.window_size,
            0.5
        );

        // Debugging: Log the analysis result
        printf("[DEBUG] Segment analysis result: Next direction: %d, Potential peak: %s, True peak: %s\n",
               result.nextDirection,
               result.isPotentialPeak ? "Yes" : "No",
               result.isTruePeak ? "Yes" : "No");

        // Step 11: Update direction based on the result of analysis
        direction = result.nextDirection;

        // If a peak is detected, log and break the loop
        if (direction == ON_PEAK) {
            // Call verify_peak_at_index to verify the detected peak
            verify_peak_at_index(phaseAngles, phase_angle_size, buffer_manager.current_phase_index);
            break; // Exit the loop after peak detection
            
        }

        printf("Moving to the %s based on analysis.\n", (direction == LEFT_SIDE) ? "left" : "right");

        // Step 12: If the flag is set to run once, break after one iteration
        if (run_once) {
            printf("Debug mode: Loop will exit after one iteration.\n");
            break; // Exit after one loop iteration when in "run once" mode
        }
    }

    // Final Debugging: Log the final state after sliding window process
    printf("[DEBUG] Final Buffer Manager State:\n");
    printf("Current buffer index: %d\n", buffer_manager.current_buffer_index);
    printf("Current phase index: %d\n", buffer_manager.current_phase_index);
}



//indexleri düzgün ayarla
//undecided için bir case ayarla. 