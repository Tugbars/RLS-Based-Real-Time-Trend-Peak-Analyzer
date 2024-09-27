#include <stdio.h>
#include <time.h>
#include <stdlib.h>

#include "windowed_running_gradient.h"
#include "running_quadratic_gradient.h"
#include "running_peak_analysis.h"
#include "running_cubic_gradient.h"
#include "rls_analysis_parameters.h"
#include "buffer_manager.h"

/**a
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

    //printf("Window start index: %u, Window size: %u\n", start_cubic_index, window_cubic);

    double average_increase = 0.0;
    double average_decrease = 0.0;
    uint16_t increase_duration = 0;
    uint16_t decrease_duration = 0;

    // Check for significant increasing trend
    if (trends->significant_increase) { //SIGNIFICANT INCREASE THRESHOLDS SHOULD BE CHECKED 
        increase_duration = trends->increase_info.end_index - trends->increase_info.start_index;
        average_increase = trends->increase_info.max_sum / (double)increase_duration;

        //printf("Significant increase detected. Increase end index: %u\n", trends->increase_info.end_index);
        //printf("Average Increase: %.6f over interval [%u - %u]\n", average_increase, trends->increase_info.start_index, trends->increase_info.end_index);

        // Use the parameters from peak_analysis_params instead of hardcoded values
        if (increase_duration > peak_analysis_params.min_consistent_trend_count || average_increase > peak_analysis_params.min_average_increase) {
            flags.move_to_right = true;
            flags.close_to_peak = true;
            flags.far_to_peak = false;
            printf("Flags set: move_to_right, close_to_peak.\n");
        }

        // Also check if the significant increase is towards the right side of the window
        if (trends->increase_info.end_index > start_cubic_index + (window_cubic / 2)) {
            if (trends->increase_info.end_index >= start_cubic_index + window_cubic - 1) {
                flags.move_to_right = true;
                flags.close_to_peak = true;
                flags.far_to_peak = false;
                printf("Flags set: move_to_right, close_to_peak.\n");
            }
        }
    }

    // Check for significant decreasing trend
    if (trends->significant_decrease) {
        decrease_duration = trends->decrease_info.end_index - trends->decrease_info.start_index;
        average_decrease = trends->decrease_info.max_sum / (double)decrease_duration;

        //printf("Significant decrease detected. Decrease start index: %u\n", trends->decrease_info.start_index);
        //printf("Average Decrease: %.6f over interval [%u - %u]\n", average_decrease, trends->decrease_info.start_index, trends->decrease_info.end_index);

        // Use the parameters from peak_analysis_params instead of hardcoded values
        if (decrease_duration > peak_analysis_params.min_consistent_trend_count || average_decrease < peak_analysis_params.min_average_decrease) {
            flags.go_to_left = true;
            flags.close_to_peak = true;
            flags.far_to_peak = false;
            printf("Flags set: go_to_left, close_to_peak.\n");
        }

        // Also check if the significant decrease is towards the left side of the window
        if (trends->decrease_info.start_index < start_cubic_index + (window_cubic / 2)) {
            flags.go_to_left = true;
            flags.close_to_peak = true;
            flags.far_to_peak = false;
            printf("Flags set: go_to_left, close_to_peak.\n");
        }
    }

    // Set on_the_peak flag if both significant increase and decrease are detected
    if (trends->significant_increase && trends->significant_decrease) {
        flags.on_the_peak = true;
        printf("On the peak detected. Flags set: on_the_peak.\n");
    }

    //printf("Trend direction determination complete.\n");
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
            //printf("--> Direction determined: RIGHT_SIDE.\n");
            break;
        case LEFT_SIDE:
            //printf("--> Direction determined: LEFT_SIDE.\n");
            break;
        case ON_PEAK:
            //printf("--> Direction determined: ON_PEAK.\n");
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
    //printf("ooo-->Dominant Side: %d\n", gradient_result.dominant_side);

    if (gradient_result.dominant_side == LEFT_SIDE) {
        printf("--> Direction overridden by gradient comparison: LEFT_SIDE.\n");
    } else if (gradient_result.dominant_side == RIGHT_SIDE) {
        printf("--> Direction overridden by gradient comparison: RIGHT_SIDE.\n");
    } else {
        printf("--> Gradient comparison result was UNDECIDED, not overriding the existing direction.\n");
    }
}

void log_significant_trend(const char* type, double max_sum) {
    printf("--> Concavity analysis. Significant %s trend with max sum %.6f detected, setting ON_PEAK.\n", type, max_sum);
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
    
    printf("\n");
}

/**
 * @brief Analyzes the trend and concavity within a window of data to determine the direction of movement towards a peak.
 *
 * This function is designed to evaluate the trend and concavity within a windowed dataset using Recursive Least Squares (RLS)
 * regression techniques. Specifically, it fits the values within the window (size defined by `window_size`) to a third-order
 * polynomial model using RLS to detect rapid increases in the first-order gradients, which can indicate proximity to a peak.
 *
 * ### Function Workflow:
 * 
 * 1. **Gradient Trend Detection via Cubic RLS**:
 *    - The function initially applies a windowed Recursive Least Squares (RLS) fitting to a third-order polynomial.
 *    - The goal is to detect rapid increases in the first-order gradients, which suggest proximity to a peak.
 *    - Using the RLS method ensures that the fitting accounts for all the data within the window while continuously updating 
 *      as new data points are added. The window size is fixed and equals the defined `window_size`.
 *    - If significant increasing gradients are detected, it indicates we are nearing a peak.
 *
 * 2. **Evaluating Trends (Left vs. Right)**:
 *    - The function evaluates the trends of increasing and decreasing intervals by analyzing the gradients.
 *    - After analyzing the trend intervals, the function determines whether the ongoing trend suggests the peak is 
 *      to the left or right of the current window.
 *    - If the evaluation suggests that the peak is still far away (i.e., no rapid changes in the gradients), 
 *      the window is split in half, and the sum of gradients for each half is compared.
 *        - The first half corresponds to the earlier part of the window, and the second half to the latter.
 *        - The function checks whether the total sum of the first-order gradients in the first half is higher 
 *          or lower than that of the second half. This helps infer the direction in which the peak might be located.
 *    - If this test yields no conclusive result, the function defaults to using the significant gradient trends 
 *      derived from the `detect_significant_gradient_trends` function.
 *
 * 3. **Peak Detection and Validation**:
 *    - If the analysis concludes that the window is centered around a peak (indicated by `on_the_peak`), the function
 *      attempts to validate whether this is a true peak or a false one.
 *    - To validate the peak, the function fits the data in the window to a second-order polynomial using RLS.
 *    - The second-order gradient is then analyzed to determine the concavity of the curve.
 *    - If the second-order gradient is sufficiently high (indicating a sharp concave curve), the function concludes that
 *      the window is centered around a valid peak.
 *
 * ### Function Details:
 * 
 * - **Input Data**: The function accepts a pointer to the data within a window (`data`), the size of the window (`window_size`), 
 *   and a `forgetting_factor` which influences how much weight is given to newer data points during the RLS fitting.
 * - **Window Splitting**: If the peak is not close, the window is split into two halves. The function then compares the 
 *   gradients of the first and second halves to infer the direction of the peak.
 * - **Trend Direction**: Based on the analysis, the function assigns a direction indicating whether the peak is likely 
 *   to be to the left or right of the current window. If the results are inconclusive, it returns `UNDECIDED`.
 * - **Concavity Analysis**: When a peak is suspected, the function verifies whether the second-order gradient 
 *   supports the presence of a true peak, marking the peak as valid if so.
 *
 * @param data The pointer to the start of the window (within the larger buffer).
 * @param window_size The number of data points in the window (subset of the buffer).
 * @param forgetting_factor A double value, typically close to 1, that determines the weight given to newer data points in the RLS analysis.
 * @return SegmentAnalysisResult A structure containing information about the direction to move (left, right, undecided),
 *         whether a potential or true peak was found, and the detected concavity pattern.
 *
 * @see detect_significant_gradient_trends
 * @see determine_trend_direction
 * @see compute_total_second_order_gradient
 * @see track_gradient_trends_with_quadratic_regression
 */
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

    printf("[analysis]-->Cubic regression increase/decrease...\n");
    PeakTrendAnalysisResult significant_trends = detect_significant_gradient_trends(data, window_size, 0, CUBIC_RLS_WINDOW, forgetting_factor);

    TrendDirectionFlags direction_flags = determine_trend_direction(&significant_trends, window_size, 0);

    // Check if we are on the peak
    if (direction_flags.on_the_peak) {
        result.nextDirection = ON_PEAK;
        log_peak_detection(true);
        goto concavity_analysis;
    }

    // Check for directional movements based on the trends
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
    
    printf("\n");

    // Check if far from peak, compare the left and right parts of the window
    if (direction_flags.far_to_peak) {
        
        printf("[analysis]--> compare gradient parts.\n");
        GradientComparisonResult gradient_result = compare_gradient_parts(data, 0, forgetting_factor); 
        log_gradient_comparison(gradient_result);
        
        printf("\n");

        // Handle the NEGATIVE_UNDECIDED case
        if (gradient_result.dominant_side == NEGATIVE_UNDECIDED) {
            //printf("[analysis]--> Both sides have negative gradients with no significant difference. Negative undecided.\n");
            result.nextDirection = NEGATIVE_UNDECIDED;
        } else if (gradient_result.dominant_side != UNDECIDED) {
            result.nextDirection = gradient_result.dominant_side;
        }
    } else {
concavity_analysis: 
        printf("[analysis]--> Performing concavity analysis...\n");

        GradientTrendResult gradient_trends = track_gradient_trends_with_quadratic_regression(data, window_size, 0, RLS_WINDOW, forgetting_factor); 
        
        // If the increase trend is strong enough, mark it as ON_PEAK
        if (gradient_trends.increase_info.valid && gradient_trends.increase_info.max_sum > 2.5) { 
            log_significant_trend("increasing", gradient_trends.increase_info.max_sum);
            result.nextDirection = ON_PEAK;
            
            printf("\n");
            goto end_analysis;
        }
        
        // If the decrease trend is strong enough, mark it as ON_PEAK
        if (gradient_trends.decrease_info.valid && gradient_trends.decrease_info.max_sum < -2.5) {
            log_significant_trend("decreasing", gradient_trends.decrease_info.max_sum);
            result.nextDirection = ON_PEAK;
            
            printf("\n");
            goto end_analysis;
        }
    }

end_analysis:
    log_final_direction(result.nextDirection);
    return result;
}



/**
 * @brief Verifies the detected peak using the find_and_verify_cubic_peak function.
 *
 * This function is triggered when the segment analysis detects a peak (ON_PEAK). It uses the current buffer index
 * as the starting point for verifying the peak by calling the find_and_verify_cubic_peak function. If the peak
 * is verified successfully, it logs the result.
 *
 * @param phaseAngles A pointer to the array of phase angles.
 * @param phase_angle_size The size of the phase angle array.
 * @param start_index The starting index for peak verification.
 */
bool verify_peak_at_index(uint16_t buffer_start_index) {
    // Print the buffer start index being verified
    printf("[verify_peak_at_index] Verifying peak at buffer start index: %d\n", buffer_start_index);

    // Call the find_and_verify_quadratic_peak function using the current buffer index
    QuadraticPeakAnalysisResult peak_result = find_and_verify_quadratic_peak(
        buffer_manager.buffer,  // Pointer to the buffer
        buffer_manager.buffer_size,  // Total buffer size
        buffer_start_index,      // Current buffer index for peak verification
        0.5                      // Forgetting factor for RLS analysis
    );

    // Debugging: Print out the results from the find_and_verify_quadratic_peak function
    //printf("[verify_peak_at_index] Peak result:\n");
    //printf("  peak_found: %d\n", peak_result.peak_found);
    //printf("  peak_index: %d\n", peak_result.peak_index);
    //printf("  is_truncated_left: %d\n", peak_result.is_truncated_left);
    //printf("  is_truncated_right: %d\n", peak_result.is_truncated_right);

    // If a peak is found, adjust the peak index to match the real phaseAngle array
    if (peak_result.peak_found) {
        // Debugging: Print current buffer and phase information
        //printf("[verify_peak_at_index] Debugging current buffer state:\n");
        //printf("  buffer_manager.current_phase_index: %d\n", buffer_manager.current_phase_index);
        //printf("  buffer_start_index: %d\n", buffer_start_index);
        //printf("  peak_result.peak_index: %d\n", peak_result.peak_index);

        // Adjust the peak index to match the real phaseAngle array
        uint16_t real_peak_index = buffer_manager.current_phase_index + (peak_result.peak_index);

        // Print the real peak index in the phaseAngle array
        printf("[verify_peak_at_index] Verified peak found at real phaseAngle index: %d\n", real_peak_index);

        return true;
    } else {
        printf("[verify_peak_at_index] Peak verification failed at buffer index: %d\n", buffer_start_index);
        return false;
    }
}

