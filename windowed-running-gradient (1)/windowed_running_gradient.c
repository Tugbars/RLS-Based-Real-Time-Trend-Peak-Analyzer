#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <assert.h>
#include "windowed_running_gradient.h"

#define NO_PEAK_FOUND UINT16_MAX
#define GRADIENT_INCREASE_THRESHOLD 1.0
#define GRADIENT_DECREASE_THRESHOLD -1.0
#define WINDOW_SIZE 5
#define DECREASE_TREND_LENGTH 5


/**
 * @brief Initializes the RunningGradient structure.
 * 
 * @param rg Pointer to the RunningGradient structure.
 */
void init_running_gradient(RunningGradient *rg, size_t window_size) {
    rg->num_points = 0;
    rg->max_points = window_size;
    rg->x = (double *)malloc(window_size * sizeof(double));
    rg->y = (double *)malloc(window_size * sizeof(double));
    rg->coefficients[0] = 0;
    rg->coefficients[1] = 0;
    rg->residual_sum_squares = 0.0;

    for (size_t i = 0; i < window_size; ++i) {
        rg->x[i] = (double)i; // Pre-fill x values for the window
    }
}

/**
 * @brief Frees the allocated memory for the RunningGradient structure.
 * 
 * @param rg Pointer to the RunningGradient structure.
 */
void free_running_gradient(RunningGradient *rg) {
    free(rg->x);
    free(rg->y);
}

/**
 * @brief Partitions the array for quickselect.
 * 
 * @param arr The array.
 * @param left The left index.
 * @param right The right index.
 * @param pivotIndex The pivot index.
 * @return The index of the pivot after partition.
 */
size_t partition(double* arr, size_t left, size_t right, size_t pivotIndex) {
    double pivotValue = arr[pivotIndex];
    double temp = arr[pivotIndex];
    arr[pivotIndex] = arr[right];
    arr[right] = temp;  // Move pivot to end
    size_t storeIndex = left;
    for (size_t i = left; i < right; i++) {
        if (arr[i] < pivotValue) {
            temp = arr[i];
            arr[i] = arr[storeIndex];
            arr[storeIndex] = temp;
            storeIndex++;
        }
    }
    temp = arr[storeIndex];
    arr[storeIndex] = arr[right];
    arr[right] = temp;  // Move pivot to its final place
    return storeIndex;
}

/**
 * @brief Quickselect algorithm to find the k-th smallest element.
 * 
 * @param arr The array.
 * @param left The left index.
 * @param right The right index.
 * @param k The index to find.
 * @return The k-th smallest element.
 */
double quickselect(double* arr, size_t left, size_t right, size_t k) {
    if (left == right) {
        return arr[left];  // If the list contains only one element, return that element
    }
    size_t pivotIndex = left + (right - left) / 2;
    pivotIndex = partition(arr, left, right, pivotIndex);
    if (k == pivotIndex) {
        return arr[k];
    } else if (k < pivotIndex) {
        return quickselect(arr, left, pivotIndex - 1, k);
    } else {
        return quickselect(arr, pivotIndex + 1, right, k);
    }
}

/**
 * @brief Adds a new data point to the RunningGradient structure.
 * 
 * @param rg Pointer to the RunningGradient structure.
 * @param y The new data point value.
 */
void add_data_point(RunningGradient *rg, double y) {
    if (rg->num_points < rg->max_points) {
        rg->y[rg->num_points++] = y;
    } else {
        // Shift data to the left to make room for the new point
        for (size_t i = 0; i < rg->max_points - 1; ++i) {
            rg->y[i] = rg->y[i + 1];
        }
        rg->y[rg->max_points - 1] = y;
    }

    // Compute the least squares coefficients using the data in the window
    double sum_x = 0.0, sum_y = 0.0, sum_xx = 0.0, sum_xy = 0.0;
    for (size_t i = 0; i < rg->num_points; ++i) {
        sum_x += rg->x[i];
        sum_y += rg->y[i];
        sum_xx += rg->x[i] * rg->x[i];
        sum_xy += rg->x[i] * rg->y[i];
    }

    double denominator = rg->num_points * sum_xx - sum_x * sum_x;
    if (denominator != 0.0) {
        rg->coefficients[0] = (rg->num_points * sum_xy - sum_x * sum_y) / denominator;
        rg->coefficients[1] = (sum_y * sum_xx - sum_x * sum_xy) / denominator;
    } else {
        rg->coefficients[0] = 0.0;
        rg->coefficients[1] = 0.0;
    }

    // Update the residual sum of squares
    rg->residual_sum_squares = 0.0;
    for (size_t i = 0; i < rg->num_points; ++i) {
        double error = rg->y[i] - (rg->coefficients[0] * rg->x[i] + rg->coefficients[1]);
        rg->residual_sum_squares += error * error;
    }
}

/**
 * @brief Calculates the current gradient (slope) of the fitted line.
 * 
 * @param rg Pointer to the RunningGradient structure.
 * @return The gradient (slope) of the fitted line.
 */
double calculate_gradient(const RunningGradient *rg) {
    assert(rg->num_points > 1 && "You must add more values into this object before calling this function.");
    return rg->coefficients[0];
}

/**
 * @brief Calculates the standard error of the gradient estimate.
 * 
 * @param rg Pointer to the RunningGradient structure.
 * @return The standard error of the gradient estimate.
 */
double calculate_standard_error(const RunningGradient *rg) {
    assert(rg->num_points > 2 && "You must add more values into this object before calling this function.");
    double s = rg->residual_sum_squares / (rg->num_points - 2);
    double adjust = 12.0 / (pow(rg->num_points, 3.0) - rg->num_points);
    return sqrt(s * adjust);
}

/**
 * @brief Calculates the normal cumulative distribution function.
 * 
 * @param value The value to calculate the CDF for.
 * @param mean The mean of the normal distribution.
 * @param stddev The standard deviation of the normal distribution.
 * @return The cumulative probability for the given value.
 */
double normal_cdf(double value, double mean, double stddev) {
    if (stddev == 0) {
        if (value < mean)
            return 0;
        else if (value > mean)
            return 1;
        else
            return 0.5;
    }
    value = (value - mean) / stddev;
    return 0.5 * erfc(-value / sqrt(2.0));
}

/**
 * @brief Calculates the probability that the gradient is less than a threshold.
 * 
 * @param rg Pointer to the RunningGradient structure.
 * @param thresh The threshold for the gradient.
 * @return The probability that the gradient is less than the threshold.
 */
double probability_gradient_less_than(const RunningGradient *rg, double thresh) {
    return normal_cdf(thresh, calculate_gradient(rg), calculate_standard_error(rg));
}

/**
 * @brief Calculates the probability that the gradient is greater than a threshold.
 * 
 * @param rg Pointer to the RunningGradient structure.
 * @param thresh The threshold for the gradient.
 * @return The probability that the gradient is greater than the threshold.
 */
double probability_gradient_greater_than(const RunningGradient *rg, double thresh) {
    return 1.0 - probability_gradient_less_than(rg, thresh);
}

/**
 * @brief Finds the upper quantile in the data.
 * 
 * This function calculates the quantile value below which a specified percentage
 * of data points fall. It uses the quickselect algorithm to efficiently find the
 * quantile value without fully sorting the data.
 * 
 * The quantile value is determined by the formula:
 * idx_upper = round((size - 1) * (1 - quantile))
 * 
 * For example, a quantile of 0.9 will find the value below which 90% of the data points fall.
 * 
 * @param data The data points.
 * @param size The number of data points.
 * @param quantile The quantile to find (must be in the range [0, 1]).
 * @return The value such that quantile percent of the values are greater than it.
 */
double find_upper_quantile(const double* data, size_t size, double quantile) {
    assert(0 <= quantile && quantile <= 1.0 && "quantile must be in the range [0, 1]");
    assert(size > 0 && "The container must contain more than 0 elements.");

    double* sorted_data = (double*)malloc(size * sizeof(double));
    if (!sorted_data) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }
    for (size_t i = 0; i < size; ++i) {
        sorted_data[i] = data[i];
    }

    size_t idx_upper = (size_t)round((size - 1) * (1 - quantile));
    double upper_q = quickselect(sorted_data, 0, size - 1, idx_upper);

    free(sorted_data);
    return upper_q;
}

/**
 * @brief Counts the number of steps without a noticeable decrease, robust to outliers.
 * 
 * This function counts the number of consecutive data points, starting from the end of the array,
 * that do not show a significant decreasing trend. It discards the upper quantile of data points 
 * to mitigate the impact of outliers. The function returns the count of such steps.
 * 
 * @param data The array of data points.
 * @param size The number of data points.
 * @param probability_of_decrease The probability threshold to consider a decrease.
 * @param quantile_discard The upper quantile of data points to discard (e.g., 0.10 for the top 10%).
 * @return The number of steps without a noticeable decrease.
 * 
 * @note This function is useful for analyzing trends in noisy data sets, where outliers can distort the overall trend.
 */
size_t count_steps_without_decrease_robust(const double* data, size_t size, double probability_of_decrease, double quantile_discard) {
    assert(0 <= quantile_discard && quantile_discard <= 1 &&
           "quantile_discard must be in the range [0, 1]");
    assert(0.5 < probability_of_decrease && probability_of_decrease < 1 && 
           "probability_of_decrease must be in the range (0.5, 1)");

    if (size == 0)
        return 0;

    double quantile_thresh = find_upper_quantile(data, size, quantile_discard);

    RunningGradient rg;
    init_running_gradient(&rg, WINDOW_SIZE);  // Use the windowed approach

    size_t count = 0;
    for (size_t i = size; i-- > 0;) {
        if (data[i] <= quantile_thresh) {
            add_data_point(&rg, data[i]);
            if (rg.num_points > 1 && calculate_gradient(&rg) < 0) {
                free_running_gradient(&rg);
                return count;
            }
        }
        count++;
    }

    free_running_gradient(&rg);
    return size;
}


/**
 * @brief Calculates the probability that the values in the data set are increasing, robust to outliers.
 * 
 * This function estimates the probability that the data trend is increasing, while ignoring the upper quantile of outliers.
 * It fits a linear model to the data points, excluding the top quantile of data points to mitigate the impact of outliers.
 * The function returns a probability that indicates the likelihood of the data trend being positive.
 * 
 * @param data The array of data points.
 * @param size The number of data points.
 * @param quantile_discard The upper quantile of data points to discard (e.g., 0.10 for the top 10%).
 * @return The probability that the values are increasing.
 * 
 * @note This function is useful for analyzing trends in noisy data sets, where outliers can distort the overall trend.
 */
double probability_values_are_increasing_robust(const double* data, size_t size, double quantile_discard) {
    assert(size > 0 && "The container must contain more than 0 elements.");
    assert(0 <= quantile_discard && quantile_discard <= 1 &&
           "quantile_discard must be in the range [0, 1]");

    double quantile_thresh = find_upper_quantile(data, size, quantile_discard);

    RunningGradient rg;
    init_running_gradient(&rg, WINDOW_SIZE);  // Use the windowed approach

    for (size_t i = 0; i < size; ++i) {
        if (data[i] <= quantile_thresh) {
            add_data_point(&rg, data[i]);
        }
    }

    double probability = 0.5; // Default probability if not enough points
    if (rg.num_points > 2) {
        probability = probability_gradient_greater_than(&rg, 0);
    }

    free_running_gradient(&rg);
    return probability;
}


/**
 * @brief Processes the data in chunks of a specified size and calculates various statistics for each chunk.
 * 
 * @param data The data points.
 * @param size The number of data points.
 */
void process_windowed_data_in_chunks(const double *data, size_t size) {
    size_t chunk_size = 10;          // Define the chunk size
    size_t window_size = 5;          // Define the window size for RunningGradient
    double quantile_discard = 0.1;   // Adjust as needed
    double probability_of_decrease = 0.51; // Adjust as needed

    for (size_t i = 0; i < size; i += chunk_size) {
        size_t end = (i + chunk_size < size) ? i + chunk_size : size;

        RunningGradient rg;
        init_running_gradient(&rg, window_size);  // Initialize with window size

        for (size_t j = i; j < end; ++j) {
            add_data_point(&rg, data[j]);
        }

        if (rg.num_points > 1) {
            double gradient = calculate_gradient(&rg);
            double std_error = calculate_standard_error(&rg);
            double prob_less_than_0 = probability_gradient_less_than(&rg, 0);
            double prob_greater_than_0 = probability_gradient_greater_than(&rg, 0);
            size_t steps_without_decrease = count_steps_without_decrease_robust(data + i, end - i, probability_of_decrease, quantile_discard);
            double prob_increasing_robust = probability_values_are_increasing_robust(data + i, end - i, quantile_discard);

            printf("Chunk %zu-%zu: Gradient = %f, Std Error = %f, P(Grad < 0) = %f, P(Grad > 0) = %f, Steps without decrease = %zu, P(Values increasing robust) = %f\n", 
                    i, end - 1, gradient, std_error, prob_less_than_0, prob_greater_than_0, steps_without_decrease, prob_increasing_robust);
        } else {
            printf("Chunk %zu-%zu: Not enough data points to calculate gradient.\n", i, end - 1);
        }

        free_running_gradient(&rg);  // Free allocated memory
    }
}

/**
 * @brief Searches for a peak in the given array interval using a windowed running gradient approach.
 * 
 * @param data The array of data points.
 * @param start_idx The starting index of the interval.
 * @param end_idx The ending index of the interval.
 * @return The index of the peak if found, otherwise NO_PEAK_FOUND.
 */
uint16_t find_peak(const double *data, size_t start_idx, size_t end_idx) {
    assert(start_idx < end_idx && "Start index must be less than end index.");

    RunningGradient rg;
    init_running_gradient(&rg, WINDOW_SIZE);

    bool increasing = false;
    uint16_t peak_index = NO_PEAK_FOUND;
    size_t decrease_trend_count = 0;

    for (size_t i = start_idx; i <= end_idx; ++i) {
        add_data_point(&rg, data[i]);

        // printf("Added point at index %zu: %f\n", i, data[i]);

        if (rg.num_points > 1) {
            double gradient = calculate_gradient(&rg);
            printf("Gradient at index %zu: %f\n", i, gradient);

            if (gradient > GRADIENT_INCREASE_THRESHOLD) {
                increasing = true;
                peak_index = (uint16_t)i; // Update peak index to the current point
                decrease_trend_count = 0; // Reset the decrease trend count
            } else if (increasing) {
                // Check for a decreasing trend
                if (gradient < GRADIENT_DECREASE_THRESHOLD) {
                    decrease_trend_count++;
                } else {
                    decrease_trend_count = 0; // Reset if trend is broken
                }

                // If we've seen a sufficient decrease trend, declare the peak
                if (decrease_trend_count >= DECREASE_TREND_LENGTH) {
                    printf("Peak found at index: %u, value: %f\n", peak_index, data[peak_index]);
                    free_running_gradient(&rg);
                    return peak_index;
                }
            }
        }
    }

    printf("No peak found in the given interval.\n");
    free_running_gradient(&rg);
    return NO_PEAK_FOUND; // No peak found
}

