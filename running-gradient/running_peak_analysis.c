#include <stdio.h>
#include <time.h>
#include <stdlib.h>

#include "windowed_running_gradient.h"
#include "running_quadratic_gradient.h"
#include "running_peak_analysis.h"
#include "running_cubic_gradient.h"
#include "rls_analysis_parameters.h"

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
    if (trends->significant_increase) {
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

    // If no significant trend is detected, check for low consistent trends
    if (!trends->significant_increase && !trends->significant_decrease) {
        increase_duration = trends->increase_info.end_index - trends->increase_info.start_index;
        decrease_duration = trends->decrease_info.end_index - trends->decrease_info.start_index;

        // Use a parameterized value for low consistent trends (replace with appropriate value)
        if (increase_duration > (WINDOW_SIZE/3)) { //SHOULD NOT BE HARDCODED. 
            flags.low_consistent_right = true;
            printf("Low consistent right trend detected over %u indices.\n", increase_duration);
        }

        if (decrease_duration > (WINDOW_SIZE/3)) {
            flags.low_consistent_left = true;
            printf("Low consistent left trend detected over %u indices.\n", decrease_duration);
        }
    }

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
 * @brief Analyzes the trend and concavity within a data segment to determine the direction of movement.
 *
 * This function performs a comprehensive analysis of a given data segment to determine the appropriate direction to move.
 * It starts by analyzing peak trends to establish an initial direction. Depending on whether the data is close to a peak
 * or far from it, it either compares gradient parts or performs concavity analysis to refine the direction. The function
 * considers various patterns in the concavity and flags indicating proximity to a peak to decide whether to move left, right,
 * or to remain undecided.
 *
 * @param data The array of `MqsRawDataPoint_t` values representing the data points.
 * @param size The total number of data points in the array.
 * @param start_index The starting index within the data array where the analysis should begin.
 * @param forgetting_factor A double value typically close to 1 that determines the weight given to newer data points in recursive analysis.
 * @return SegmentAnalysisResult A structure containing information about the direction to move (left, right, undecided),
 *         whether a potential or true peak was found, and the detected concavity pattern.
 */
SegmentAnalysisResult segment_trend_and_concavity_analysis(const MqsRawDataPoint_t *data, uint16_t size, uint16_t start_index, double forgetting_factor) {
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
    PeakTrendAnalysisResult significant_trends = detect_significant_gradient_trends(data, size, start_index, 30, forgetting_factor);

    printf("****************************************** Determining cubic regression linear gradient sums...\n");
    TrendDirectionFlags direction_flags = determine_trend_direction(&significant_trends, WINDOW_SIZE, start_index);

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
        GradientComparisonResult gradient_result = compare_gradient_parts(data, start_index, forgetting_factor);
        log_gradient_comparison(gradient_result);

        if (gradient_result.dominant_side != UNDECIDED) {
            result.nextDirection = gradient_result.dominant_side;
        }
    } else {
concavity_analysis: 
        printf("****************************************** Performing concavity analysis...\n");

        GradientTrendResult gradient_trends = track_gradient_trends_with_quadratic_regression(data, size, start_index, 30, forgetting_factor);

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
        ConcavityAnalysisResult concavity_result = initial_concavity_analysis(data, size, start_index, forgetting_factor, true);
        ConcavityAnalysisOutput concavity_output = analyze_concavity_segments(&concavity_result);
        log_concavity_analysis(concavity_output);
    }

end_analysis:
    log_final_direction(result.nextDirection);
    return result;
}

BufferManager buffer_manager;

void init_buffer_manager(MqsRawDataPoint_t* buffer, uint16_t buffer_size, uint16_t window_size, int16_t start_index, double start_frequency, double frequency_increment) {
    buffer_manager.buffer = buffer;
    buffer_manager.buffer_size = buffer_size;
    buffer_manager.window_size = window_size;
    buffer_manager.current_phase_index = start_index;
    buffer_manager.current_buffer_index = buffer_size / 2; // Start from the middle
    buffer_manager.mapping_start_index = start_index;
    buffer_manager.mapping_end_index = start_index + window_size - 1;
    buffer_manager.start_frequency = start_frequency;
    buffer_manager.frequency_increment = frequency_increment;
}

/**
 * @brief Updates the buffer with phase angle data based on direction and current buffer index.
 */
void update_buffer_for_direction(const double* phaseAngles, int direction) {
    // Determine the next range of values to load based on direction
    int16_t next_phase_index = buffer_manager.current_phase_index;

    // Move left or right based on direction
    if (direction == LEFT_SIDE) {
        next_phase_index -= buffer_manager.window_size / 2; // Move left
    } else if (direction == RIGHT_SIDE) {
        next_phase_index += buffer_manager.window_size / 2; // Move right
    }

    // Ensure the phase index stays within bounds of phaseAngles
    if (next_phase_index < 0 || next_phase_index >= buffer_manager.buffer_size) {
        printf("Out of bounds of phase angles.\n");
        return;
    }

    // Update the buffer with new phaseAngle values
    for (uint16_t i = 0; i < buffer_manager.window_size; ++i) {
        int16_t phase_index = next_phase_index + i;  // New phase index based on window size
        int16_t buffer_index = buffer_manager.current_buffer_index + (direction == LEFT_SIDE ? -i : i); // Adjust buffer index for left or right

        // Ensure buffer index stays within bounds
        if (buffer_index < 0 || buffer_index >= buffer_manager.buffer_size) {
            printf("Buffer out of bounds!\n");
            break;
        }

        // Copy phaseAngle and set impedance to 0.0 (mock value)
        buffer_manager.buffer[buffer_index].phaseAngle = phaseAngles[phase_index];
        buffer_manager.buffer[buffer_index].impedance = 0.0; // Mock impedance, can be changed
    }

    // Update the current phase index and buffer index based on direction
    buffer_manager.current_phase_index = next_phase_index;
    buffer_manager.current_buffer_index = (direction == LEFT_SIDE) 
        ? buffer_manager.current_buffer_index - buffer_manager.window_size / 2
        : buffer_manager.current_buffer_index + buffer_manager.window_size / 2;

    // Update mapping start and end indices
    buffer_manager.mapping_start_index = next_phase_index;
    buffer_manager.mapping_end_index = next_phase_index + buffer_manager.window_size - 1;
}

/**
 * @brief Calculate the frequency for a given index.
 */
double calculate_frequency(int16_t index) {
    return buffer_manager.start_frequency + index * buffer_manager.frequency_increment;
}

/**
 * @brief Performs sliding window analysis with buffer management.
 */
void perform_sliding_window_analysis(const double* phaseAngles, uint16_t phase_angle_size) {
    int16_t direction = UNDECIDED;

    // Continue until ON_PEAK is found or we run out of data
    while (direction != ON_PEAK && buffer_manager.current_phase_index >= 0 && buffer_manager.current_phase_index < phase_angle_size) {
        // Load the next window of values into the buffer
        update_buffer_for_direction(phaseAngles, direction);

        // Analyze the loaded values (e.g., using segment_trend_and_concavity_analysis)
        SegmentAnalysisResult result = segment_trend_and_concavity_analysis(buffer_manager.buffer, buffer_manager.buffer_size, 
                                                                            buffer_manager.current_buffer_index, 0.5);

        // Update direction based on the result of analysis
        direction = result.nextDirection;

        if (direction == ON_PEAK) {
            printf("Peak detected at phase index: %d\n", buffer_manager.current_phase_index);
            printf("Mapped Start Index: %d, Frequency: %.2f Hz\n", buffer_manager.mapping_start_index, calculate_frequency(buffer_manager.mapping_start_index));
            printf("Mapped End Index: %d, Frequency: %.2f Hz\n", buffer_manager.mapping_end_index, calculate_frequency(buffer_manager.mapping_end_index));
            break;
        }

        printf("Moving to the %s based on analysis.\n", (direction == LEFT_SIDE) ? "left" : "right");
    }
}