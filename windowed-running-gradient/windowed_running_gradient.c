#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <assert.h>
#include "windowed_running_gradient.h"

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
    rg->inverse_cov_matrix[0][0] = 1e6; rg->inverse_cov_matrix[0][1] = 0;
    rg->inverse_cov_matrix[1][0] = 0; rg->inverse_cov_matrix[1][1] = 1e6;

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
 * @brief Adds a new data point to the RunningGradient structure and updates the regression coefficients.
 * 
 * This function updates the regression coefficients using a numerically stable Recursive Least Squares (RLS) algorithm.
 * It maintains a window of recent data points and adjusts the coefficients incrementally as new data points are added.
 * The inverse covariance matrix is updated to ensure efficient and stable coefficient updates.
 * 
 * @param rg Pointer to the RunningGradient structure.
 * @param y The new data point value.
 */
void add_data_point(RunningGradient *rg, double y) {
    // Add the new data point to the windowed data set
    if (rg->num_points < rg->max_points) {
        // Add the new data point to the array and increment the number of data points
        rg->y[rg->num_points++] = y;
    } else {
        // Shift all existing data points one position to the left
        for (size_t i = 0; i < rg->max_points - 1; ++i) {
            rg->y[i] = rg->y[i + 1];
        }
        // Add the new data point at the end of the array
        rg->y[rg->max_points - 1] = y;
    }

    // Update the inverse covariance matrix and coefficients using RLS
    // Setup the data point vector x
    double x[2] = {rg->num_points, 1};

    // Calculate a scalar value 'temp' to ensure numerical stability
    double temp = 1.0 + (x[0] * rg->inverse_cov_matrix[0][0] + x[1] * rg->inverse_cov_matrix[1][0]) * x[0] + 
                       (x[0] * rg->inverse_cov_matrix[0][1] + x[1] * rg->inverse_cov_matrix[1][1]) * x[1];

    // Check if 'temp' is too close to zero to avoid numerical instability issues
    if (fabs(temp) < 1e-10) {
        fprintf(stderr, "Numerical stability issue: temp is too close to zero.\n");
        exit(EXIT_FAILURE);
    }

    // Compute temporary vector 'tmp' to adjust the inverse covariance matrix
    double tmp[2] = {rg->inverse_cov_matrix[0][0] * x[0] + rg->inverse_cov_matrix[0][1] * x[1], 
                     rg->inverse_cov_matrix[1][0] * x[0] + rg->inverse_cov_matrix[1][1] * x[1]};

    // Update the inverse covariance matrix by removing old data influence and adding new data point influence
    rg->inverse_cov_matrix[0][0] -= (tmp[0] * tmp[0]) / temp;
    rg->inverse_cov_matrix[0][1] -= (tmp[0] * tmp[1]) / temp;
    rg->inverse_cov_matrix[1][0] = rg->inverse_cov_matrix[0][1]; // Ensure symmetry
    rg->inverse_cov_matrix[1][1] -= (tmp[1] * tmp[1]) / temp;

    // Calculate the prediction error
    double prediction_error = y - (x[0] * rg->coefficients[0] + x[1] * rg->coefficients[1]);

    // Update the coefficients (slope and intercept) of the linear model using the new data point
    // The inverse covariance matrix helps determine the influence of the new data point on the coefficients
    rg->coefficients[0] += (rg->inverse_cov_matrix[0][0] * x[0] + rg->inverse_cov_matrix[0][1] * x[1]) * prediction_error;
    rg->coefficients[1] += (rg->inverse_cov_matrix[1][0] * x[0] + rg->inverse_cov_matrix[1][1]) * prediction_error;

    // Update the residual sum of squares (RSS)
    // RSS measures the discrepancy between the data and the regression model
    rg->residual_sum_squares = 0.0;
    for (size_t i = 0; i < rg->num_points; ++i) {
        // Calculate the prediction error for each data point within the window
        double error = rg->y[i] - (rg->coefficients[0] * rg->x[i] + rg->coefficients[1]);
        // Sum the squared prediction errors to obtain the residual sum of squares
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
 * @brief Counts the number of steps with a noticeable increase, robust to outliers.
 * 
 * This function counts the number of consecutive data points, starting from the start of the array,
 * that show a significant increasing trend. It discards the upper quantile of data points 
 * to mitigate the impact of outliers. The function returns the count of such steps.
 * 
 * @param data The array of data points.
 * @param size The number of data points.
 * @param probability_of_increase The probability threshold to consider an increase.
 * @param quantile_discard The upper quantile of data points to discard (e.g., 0.10 for the top 10%).
 * @return The number of steps with a noticeable increase.
 * 
 * @note This function is useful for analyzing trends in noisy data sets, where outliers can distort the overall trend.
 */
size_t count_steps_with_increase_robust(const double* data, size_t size, double probability_of_increase, double quantile_discard) {
    assert(0 <= quantile_discard && quantile_discard <= 1 &&
           "quantile_discard must be in the range [0, 1]");
    assert(0.5 < probability_of_increase && probability_of_increase < 1 && 
           "probability_of_increase must be in the range (0.5, 1)");

    if (size == 0)
        return 0;

    double quantile_thresh = find_upper_quantile(data, size, quantile_discard);

    RunningGradient rg;
    init_running_gradient(&rg, WINDOW_SIZE);  // Use the windowed approach

    size_t count = 0;
    for (size_t i = 0; i < size; ++i) {
        if (data[i] <= quantile_thresh) {
            add_data_point(&rg, data[i]);
            if (rg.num_points > 1 && calculate_gradient(&rg) > 0) {
                count++;
            } else {
                free_running_gradient(&rg);
                return count;
            }
        }
    }

    free_running_gradient(&rg);
    return count;
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
            //size_t steps_without_decrease = count_steps_without_decrease_robust(data + i, end - i, probability_of_decrease, quantile_discard);
            double prob_increasing_robust = probability_values_are_increasing_robust(data + i, end - i, quantile_discard);

            printf("Chunk %zu-%zu: Gradient = %f, Std Error = %f, P(Grad < 0) = %f, P(Grad > 0) = %f, P(Values increasing robust) = %f\n", 
                    i, end - 1, gradient, std_error, prob_less_than_0, prob_greater_than_0, prob_increasing_robust);
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
uint16_t find_peak(const double *data, size_t start_idx, size_t end_idx, size_t  window_size) {
    assert(start_idx < end_idx && "Start index must be less than end index.");

    RunningGradient rg;
    init_running_gradient(&rg, window_size);

    bool increasing = false;
    uint16_t peak_index = NO_PEAK_FOUND;
    size_t decrease_trend_count = 0;

    for (size_t i = start_idx; i <= end_idx; ++i) {
        add_data_point(&rg, data[i]);

         //printf("Added point at index %zu: %f\n", i, data[i]);

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
                    //printf("Peak found at index: %u, value: %f\n", peak_index, data[peak_index]);
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

/**
 * @brief Wraps the find_peak function to take every fifth item in the specified interval.
 * 
 * @param data The array of data points.
 * @param start_idx The starting index of the interval.
 * @param end_idx The ending index of the interval.
 * @return The index of the peak if found, otherwise NO_PEAK_FOUND.
 */
uint16_t find_peak_every_fifth(const double *data, size_t start_idx, size_t end_idx, size_t  window_size) {
    assert(start_idx < end_idx && "Start index must be less than end index.");
    
    size_t reduced_size = (end_idx - start_idx) / 5 + 1;
    double *reduced_data = (double *)malloc(reduced_size * sizeof(double));
    if (reduced_data == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }
    
    size_t idx = 0;
    for (size_t i = start_idx; i <= end_idx; i += 5) {
        reduced_data[idx++] = data[i];
    }

    uint16_t peak_index = find_peak(reduced_data, 0, idx - 1, window_size);
    
    if (peak_index != NO_PEAK_FOUND) {
        peak_index = start_idx + peak_index * 5;
    }

    free(reduced_data);
    return peak_index;
}


/**
 * @brief Counts the number of steps with a noticeable increase to the left, robust to outliers.
 * 
 * This function counts the number of consecutive data points, starting from the start of the array,
 * that show a significant increasing trend. It discards the upper quantile of data points 
 * to mitigate the impact of outliers. The function returns the count of such steps.
 * 
 * @param data The array of data points.
 * @param start_idx The starting index.
 * @param probability_of_increase The probability threshold to consider an increase.
 * @param quantile_discard The upper quantile of data points to discard (e.g., 0.10 for the top 10%).
 * @param gradient_increase_threshold The threshold above which a gradient is considered an increase.
 * @return The number of steps with a noticeable increase.
 * 
 * @note This function is useful for analyzing trends in noisy data sets, where outliers can distort the overall trend.
 */
size_t count_steps_with_increase_robust_left(const double* data, size_t start_idx, double quantile_discard, double gradient_increase_threshold, double gradient_decrease_threshold) {
    assert(0 <= quantile_discard && quantile_discard <= 1 &&
           "quantile_discard must be in the range [0, 1]");

    if (start_idx == 0)
        return 0;

    double quantile_thresh = find_upper_quantile(data, start_idx + 1, quantile_discard);

    RunningGradient rg;
    init_running_gradient(&rg, WINDOW_SIZE);  // Use the windowed approach

    size_t count = 0;
    for (size_t i = start_idx; i != SIZE_MAX; --i) {
        if (data[i] <= quantile_thresh) {
            add_data_point(&rg, data[i]);
            if (rg.num_points > 1) {
                double gradient = calculate_gradient(&rg);
                if (gradient > gradient_increase_threshold) {
                    count++;
                } else if (gradient > gradient_decrease_threshold) {
                    // Ignore small decreases
                    continue;
                } else {
                    break;
                }
            }
        }
        if (i == 0) break;  // Prevent underflow
    }

    free_running_gradient(&rg);
    return count;
}

/**
 * @brief Counts the number of steps with a noticeable increase to the right, robust to outliers.
 * 
 * This function counts the number of consecutive data points, starting from the start of the array,
 * that show a significant increasing trend. It discards the upper quantile of data points 
 * to mitigate the impact of outliers. The function returns the count of such steps.
 * 
 * @param data The array of data points.
 * @param start_idx The starting index.
 * @param probability_of_increase The probability threshold to consider an increase.
 * @param quantile_discard The upper quantile of data points to discard (e.g., 0.10 for the top 10%).
 * @param gradient_increase_threshold The threshold above which a gradient is considered an increase.
 * @return The number of steps with a noticeable increase.
 * 
 * @note This function is useful for analyzing trends in noisy data sets, where outliers can distort the overall trend.
 */
size_t count_steps_with_increase_robust_right(const double* data, size_t start_idx, size_t size, double quantile_discard, double gradient_increase_threshold, double gradient_decrease_threshold) {
    assert(0 <= quantile_discard && quantile_discard <= 1 &&
           "quantile_discard must be in the range [0, 1]");
    assert(start_idx < size && "start_idx must be less than size");

    if (start_idx >= size - 1)
        return 0;

    double quantile_thresh = find_upper_quantile(data + start_idx, size - start_idx, quantile_discard);

    RunningGradient rg;
    init_running_gradient(&rg, WINDOW_SIZE);  // Use the windowed approach

    size_t count = 0;
    for (size_t i = start_idx + 1; i < size; ++i) {
        if (data[i] <= quantile_thresh) {
            add_data_point(&rg, data[i]);
            if (rg.num_points > 1) {
                double gradient = calculate_gradient(&rg);
                if (gradient > gradient_increase_threshold) {
                    count++;
                } else if (gradient > gradient_decrease_threshold) {
                    // Ignore small decreases
                    continue;
                } else {
                    break;
                }
            }
        }
    }

    free_running_gradient(&rg);
    return count;
}

int determine_trend_direction_robust(const double *data, size_t size, size_t start_idx, size_t window_size, size_t trend_length, 
                              double gradient_increase_threshold, double gradient_decrease_threshold, double quantile_discard, double probability_of_increase) {
    if (start_idx >= size) {
        printf("start_idx is out of bounds\n");
        return 0;
    }

    printf("Starting trend analysis at index %zu with window size %zu and trend length %zu\n", start_idx, window_size, trend_length);

    size_t left_increase_count = 0, right_increase_count = 0;
    double left_prob_increase = 0.0, right_prob_increase = 0.0;

    // Analyze the left side
    if (start_idx > 0) {
        left_increase_count = count_steps_with_increase_robust_left(data, start_idx, quantile_discard, gradient_increase_threshold, gradient_decrease_threshold);
        left_prob_increase = probability_values_are_increasing_robust(data, start_idx, quantile_discard);
    }

    // Analyze the right side
    if (start_idx < size - 1) {
        right_increase_count = count_steps_with_increase_robust_right(data, start_idx, size, quantile_discard, gradient_increase_threshold, gradient_decrease_threshold);
        right_prob_increase = probability_values_are_increasing_robust(data + start_idx + 1, size - start_idx - 1, quantile_discard);
    }

    // Determine the direction based on the increase counts and probabilities
    if (left_increase_count >= trend_length) {
        printf("Trend increasing to the left detected.\n");
        return -1; // Increasing trend to the left
    } else if (right_increase_count >= trend_length) {
        printf("Trend increasing to the right detected.\n");
        return 1; // Increasing trend to the right
    } else if (left_prob_increase > right_prob_increase) {
        printf("Higher probability of increase to the left detected.\n");
        return -1; // Higher probability of increase to the left
    } else if (right_prob_increase > left_prob_increase) {
        printf("Higher probability of increase to the right detected.\n");
        return 1; // Higher probability of increase to the right
    }

    printf("Unclear trend detected.\n");
    return 0; // Unclear trend
}


