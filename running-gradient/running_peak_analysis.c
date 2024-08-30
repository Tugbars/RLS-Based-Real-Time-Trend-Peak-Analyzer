#include <stdio.h>
#include <time.h>
#include <stdlib.h>

#include "windowed_running_gradient.h"
#include "running_quadratic_gradient.h"
#include "running_peak_analysis.h"
#include "running_cubic_gradient.h"

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
static TrendDirectionFlags determine_trend_direction(const PeakTrendAnalysisResult *trends, size_t window_cubic, size_t start_cubic_index) {
    TrendDirectionFlags flags = {false, false, false, true, false, false, false}; // Initialize with on_the_peak

    printf("Window start index: %zu, Window size: %zu\n", start_cubic_index, window_cubic);

    // Check for significant increasing trend
    if (trends->significant_increase) {
        size_t increase_end_index = trends->increase_info.end_index;

        printf("Significant increase detected. Increase end index: %zu\n", increase_end_index);

        // Determine if the significant increase is towards the right side of the window
        if (increase_end_index > start_cubic_index + (window_cubic / 2)) {
            if (increase_end_index >= start_cubic_index + window_cubic - 1) {
                flags.move_to_right = true;
                flags.close_to_peak = true;
                flags.far_to_peak = false;
                printf("Significant increasing trend extends beyond the window. Flags set: move_to_right, close_to_peak.\n");
            }
        }
    }

    // Check for significant decreasing trend
    if (trends->significant_decrease) {
        size_t decrease_start_index = trends->decrease_info.start_index;

        printf("Significant decrease detected. Decrease start index: %zu\n", decrease_start_index);

        // Determine if the significant decrease is towards the left side of the window
        if (decrease_start_index < start_cubic_index + (window_cubic / 2)) {
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
        size_t increase_duration = trends->increase_info.end_index - trends->increase_info.start_index;
        size_t decrease_duration = trends->decrease_info.end_index - trends->decrease_info.start_index;

        if (increase_duration > 15) {  // Replace with parameterized value later
            flags.low_consistent_right = true;
            printf("Low consistent right trend detected over %zu indices.\n", increase_duration);
        }

        if (decrease_duration > 15) {  // Replace with parameterized value later
            flags.low_consistent_left = true;
            printf("Low consistent left trend detected over %zu indices.\n", decrease_duration);
        }
    }

    printf("Trend direction determination complete.\n");
    return flags;
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
 * @param data The array of double values representing the data points.
 * @param size The total number of data points in the array.
 * @param start_index The starting index within the data array where the analysis should begin.
 * @param forgetting_factor A double value typically close to 1 that determines the weight given to newer data points in recursive analysis.
 * @return SegmentAnalysisResult A structure containing information about the direction to move (left, right, undecided), 
 *         whether a potential or true peak was found, and the detected concavity pattern.
 */
SegmentAnalysisResult segment_trend_and_concavity_analysis(const double *data, size_t size, size_t start_index, double forgetting_factor) {
    SegmentAnalysisResult result = {false, false, UNDECIDED, UNDETERMINED_PATTERN};
    
    printf("****************************************** Cubic regression increase/decrease interval analysis...\n");
    PeakTrendAnalysisResult significant_trends = detect_significant_gradient_trends(data, size, start_index, 30, forgetting_factor);
    
    printf("****************************************** Determining cubic regression linear gradient sums...\n");
    TrendDirectionFlags direction_flags = determine_trend_direction(&significant_trends, 30, start_index);

    if (direction_flags.on_the_peak) {
        result.nextDirection = ON_PEAK;
        printf("--> On the peak detected. Skipping direction analysis and proceeding with concavity analysis.\n");
        goto concavity_analysis;
    }

    if (direction_flags.move_to_right) {
        result.nextDirection = RIGHT_SIDE;
        printf("--> Direction determined: RIGHT_SIDE (move_to_right flag is set).\n");
    } else if (direction_flags.go_to_left) {
        result.nextDirection = LEFT_SIDE;
        printf("--> Direction determined: LEFT_SIDE (go_to_left flag is set).\n");
    } else if (direction_flags.low_consistent_right) {
        result.nextDirection = RIGHT_SIDE;
        printf("--> Direction determined: RIGHT_SIDE (low consistent right trend detected).\n");
    } else if (direction_flags.low_consistent_left) {
        result.nextDirection = LEFT_SIDE;
        printf("--> Direction determined: LEFT_SIDE (low consistent left trend detected).\n");
    }

    if (direction_flags.far_to_peak) {
        printf("--> Far from peak. Determining direction by comparing gradient parts...\n");
        GradientComparisonResult gradient_result = compare_gradient_parts(data, start_index, forgetting_factor);

        printf("--> Total Gradient First Part: %.6f\n", gradient_result.total_gradient_first_part);
        printf("--> Total Gradient Second Part: %.6f\n", gradient_result.total_gradient_second_part);
        printf("--> Increase Count First Part: %zu\n", gradient_result.increase_count_first_part);
        printf("--> Increase Count Second Part: %zu\n", gradient_result.increase_count_second_part);
        printf("Dominant Side: %d\n", gradient_result.dominant_side);

        if (gradient_result.dominant_side != UNDECIDED) {
            result.nextDirection = gradient_result.dominant_side;
            if (gradient_result.dominant_side == LEFT_SIDE) {
                printf("--> Direction overridden by gradient comparison: LEFT_SIDE.\n");
            } else if (gradient_result.dominant_side == RIGHT_SIDE) {
                printf("--> Direction overridden by gradient comparison: RIGHT_SIDE.\n");
            }
        } else {
            printf("--> Gradient comparison result was UNDECIDED, not overriding the existing direction.\n");
        }
    } else {
        concavity_analysis:
        printf("--> Close to peak. Performing concavity analysis...\n");
        
        if (significant_trends.significant_increase) {
            size_t start_index = significant_trends.increase_info.start_index;
            size_t end_index = significant_trends.increase_info.end_index;
            size_t window_size = end_index - start_index;
    
            printf("--> Significant increase detected from index %zu to %zu. Calling quadratic regression...\n", start_index, end_index);
            track_gradient_trends_with_quadratic_regression(data, size, start_index, window_size, forgetting_factor);
        }
        
        if (significant_trends.significant_decrease) {
            size_t start_index = significant_trends.decrease_info.start_index;
            size_t end_index = significant_trends.decrease_info.end_index;
            size_t window_size = end_index - start_index;
    
            printf("--> Significant decrease detected from index %zu to %zu. Calling quadratic regression...\n", start_index, end_index);
            track_gradient_trends_with_quadratic_regression(data, size, start_index, window_size, forgetting_factor);
        }
        
        printf("*********************** CONCAVITY SEGMENT ANALYSIS ************************\n");
        ConcavityAnalysisResult concavity_result = initial_concavity_analysis(data, size, start_index, forgetting_factor, true);
        result.pattern = analyze_concavity_segments(&concavity_result, &result.isPotentialPeak, &result.isTruePeak);

        switch (result.pattern) {
            case INCREASE_INCREASE_DECREASE:
                printf("--> Pattern: INCREASE -> INCREASE -> DECREASE\n");
                break;
            case INCREASE_DECREASE_INCREASE:
                printf("--> Pattern: INCREASE -> DECREASE -> INCREASE\n");
                break;
            case INCREASE_DECREASE_DECREASE:
                printf("--> Pattern: INCREASE -> DECREASE -> DECREASE\n");
                break;
            case DECREASE_INCREASE_INCREASE:
                printf("--> Pattern: DECREASE -> INCREASE -> INCREASE\n");
                break;
            case DECREASE_INCREASE_DECREASE:
                printf("--> Pattern: DECREASE -> INCREASE -> DECREASE\n");
                break;
            case DECREASE_DECREASE_INCREASE:
                printf("--> Pattern: DECREASE -> DECREASE -> INCREASE\n");
                break;
            case DECREASE_DECREASE_DECREASE:
                printf("--> Pattern: DECREASE -> DECREASE -> DECREASE\n");
                break;
            case INCREASE_INCREASE_INCREASE:
                printf("--> Pattern: INCREASE -> INCREASE -> INCREASE\n");
                break;
            default:
                printf("--> Pattern: UNDETERMINED\n");
                break;
        }

        if (result.isTruePeak || result.isPotentialPeak) { 
            printf("--> Peak detected. Waiting for further analysis.\n");
        }
    }

    printf("Final Direction: ");
    switch (result.nextDirection) {
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

    return result;
}


