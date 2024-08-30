#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <assert.h>
#include "windowed_running_gradient.h"


/**
 * @brief Initializes the RunningGradient structure.
 *
 * This function initializes the RunningGradient structure, which includes setting initial values for 
 * coefficients, residual sum of squares, inverse covariance matrix, and pre-filling the x values for the window.
 *
 * @param rg Pointer to the RunningGradient structure.
 */
void init_running_gradient(RunningGradient *rg) {
    rg->num_points = 0;
    rg->max_points = WINDOW_SIZE;
    rg->coefficients[0] = 0;
    rg->coefficients[1] = 0;
    rg->residual_sum_squares = 0.0;
    rg->inverse_cov_matrix[0][0] = 1e9;
    rg->inverse_cov_matrix[0][1] = 0;
    rg->inverse_cov_matrix[1][0] = 0;
    rg->inverse_cov_matrix[1][1] = 1e9;
    rg->forgetting_factor = FORGETTING_FACTOR;

    for (size_t i = 0; i < WINDOW_SIZE; ++i) {
        rg->x[i] = (double)i; // Pre-fill x values for the window
        rg->y[i] = 0.0; // Initialize y values to 0
    }
}

/**
 * @brief Partitions the array for quickselect.
 *
 * This function partitions the array around a pivot element chosen during quickselect. 
 * Elements less than the pivot are moved to the left of the pivot, and elements greater 
 * than the pivot are moved to the right.
 *
 * @param arr The array to be partitioned.
 * @param left The left index for the partitioning.
 * @param right The right index for the partitioning.
 * @param pivotIndex The pivot index around which partitioning is done.
 * @return The index of the pivot after partitioning.
 */
static inline size_t partition(double* arr, size_t left, size_t right, size_t pivotIndex) {
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
 * This function implements the Quickselect algorithm to find the k-th smallest element in an array.
 * Quickselect is similar to Quicksort and works by repeatedly partitioning the array.
 *
 * @param arr The array in which to find the k-th smallest element.
 * @param left The left index for the partitioning.
 * @param right The right index for the partitioning.
 * @param k The index of the k-th smallest element to find.
 * @return The k-th smallest element in the array.
 */
static inline double quickselect(double* arr, size_t left, size_t right, size_t k) {
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
 * The Recursive Least Squares (RLS) algorithm updates the coefficients of the linear regression model by incorporating
 * new data points one at a time. It is a form of adaptive filter that adjusts the model parameters to minimize the 
 * prediction error. The inverse covariance matrix is used to calculate the adjustment to the coefficients.
 *
 * RLS is particularly suitable for running gradient calculations because it efficiently updates the model parameters
 * without the need to reprocess all previous data points. This is essential for real-time applications where data is
 * continuously arriving, and the gradient needs to be updated dynamically.
 *
 * @param rg Pointer to the RunningGradient structure.
 * @param y The new data point value.
 */
void add_data_point(RunningGradient *const rg, const double y) {
    // Check if the current window is not full
    if (rg->num_points < rg->max_points) {
        // Add new data point to the window
        rg->y[rg->num_points] = y;
        rg->x[rg->num_points] = (double)rg->num_points;
        rg->num_points++;
    } else {
        // Shift the window to the left and add new data point
        for (size_t i = 0; i < rg->max_points - 1; ++i) {
            rg->y[i] = rg->y[i + 1];
            rg->x[i] = rg->x[i + 1];
        }
        rg->y[rg->max_points - 1] = y;
        rg->x[rg->max_points - 1] = (double)(rg->num_points);
        rg->num_points++;
    }

    // Define the input vector for the new data point
    double x[2] = { (double)(rg->num_points - 1), 1.0 }; // Adjusting index for the new point

    // Calculate the denominator term in the RLS update equations
    const double temp = rg->forgetting_factor + 
                        (x[0] * rg->inverse_cov_matrix[0][0] + x[1] * rg->inverse_cov_matrix[1][0]) * x[0] +
                        (x[0] * rg->inverse_cov_matrix[0][1] + x[1] * rg->inverse_cov_matrix[1][1]) * x[1];

    // Check for numerical stability
    if (fabs(temp) < 1e-10) {
        fprintf(stderr, "Numerical stability issue: temp is too close to zero.\n");
        exit(EXIT_FAILURE);
    }

    // Calculate the temporary vector used for updating the inverse covariance matrix
    const double tmp[2] = { rg->inverse_cov_matrix[0][0] * x[0] + rg->inverse_cov_matrix[0][1] * x[1],
                            rg->inverse_cov_matrix[1][0] * x[0] + rg->inverse_cov_matrix[1][1] * x[1] };

    // Update the inverse covariance matrix using the Sherman-Morrison formula
    rg->inverse_cov_matrix[0][0] = (rg->inverse_cov_matrix[0][0] - (tmp[0] * tmp[0]) / temp) / rg->forgetting_factor;
    rg->inverse_cov_matrix[0][1] = (rg->inverse_cov_matrix[0][1] - (tmp[0] * tmp[1]) / temp) / rg->forgetting_factor;
    rg->inverse_cov_matrix[1][0] = rg->inverse_cov_matrix[0][1]; // Ensure symmetry
    rg->inverse_cov_matrix[1][1] = (rg->inverse_cov_matrix[1][1] - (tmp[1] * tmp[1]) / temp) / rg->forgetting_factor;

    // Calculate the prediction error for the new data point
    const double prediction_error = y - (x[0] * rg->coefficients[0] + x[1] * rg->coefficients[1]);

    // Update the regression coefficients using the RLS update equations
    rg->coefficients[0] += (rg->inverse_cov_matrix[0][0] * x[0] + rg->inverse_cov_matrix[0][1] * x[1]) * prediction_error;
    rg->coefficients[1] += (rg->inverse_cov_matrix[1][0] * x[0] + rg->inverse_cov_matrix[1][1]) * prediction_error;

    // Recalculate the residual sum of squares
    rg->residual_sum_squares = 0.0;
    for (size_t i = 0; i < rg->num_points; ++i) {
        const double error = rg->y[i] - (rg->coefficients[0] * rg->x[i] + rg->coefficients[1]);
        rg->residual_sum_squares += error * error;
    }
}

/**
 * @brief Calculates the current gradient (slope) of the fitted line.
 *
 * This function calculates the gradient (slope) of the linear regression model fitted to the data points
 * in the RunningGradient structure. The gradient represents the rate of change of the dependent variable
 * with respect to the independent variable.
 *
 * @param rg Pointer to the RunningGradient structure.
 * @return The gradient (slope) of the fitted line.
 */
inline double calculate_gradient(const RunningGradient * const rg) {
    assert(rg->num_points > 1 && "You must add more values into this object before calling this function.");
    return rg->coefficients[0];
}

/**
 * @brief Calculates the standard error of the gradient estimate.
 *
 * This function calculates the standard error of the gradient estimate in the linear regression model.
 * The standard error provides a measure of the accuracy of the gradient estimate.
 *
 * @param rg Pointer to the RunningGradient structure.
 * @return The standard error of the gradient estimate.
 */
double calculate_standard_error(const RunningGradient * const rg) {
    assert(rg->num_points > 2 && "You must add more values into this object before calling this function.");
    const double s = rg->residual_sum_squares / (rg->num_points - 2);
    const double adjust = 12.0 / (pow(rg->num_points, 3.0) - rg->num_points);
    return sqrt(s * adjust);
}

/**
 * @brief Calculates the normal cumulative distribution function.
 *
 * This function calculates the cumulative distribution function (CDF) of the normal distribution
 * for a given value, mean, and standard deviation. The CDF gives the probability that a normally
 * distributed random variable is less than or equal to the specified value.
 *
 * @param value The value to calculate the CDF for.
 * @param mean The mean of the normal distribution.
 * @param stddev The standard deviation of the normal distribution.
 * @return The cumulative probability for the given value.
 */
double normal_cdf(const double value, const double mean, const double stddev) {
    if (stddev == 0) {
        if (value < mean)
            return 0;
        else if (value > mean)
            return 1;
        else
            return 0.5;
    }
    const double normalized_value = (value - mean) / stddev;
    return 0.5 * erfc(-normalized_value / sqrt(2.0));
}

/**
 * @brief Calculates the probability that the gradient is less than a threshold.
 *
 * This function calculates the probability that the gradient (slope) of the linear regression model
 * in the RunningGradient structure is less than a specified threshold.
 *
 * @param rg Pointer to the RunningGradient structure.
 * @param thresh The threshold for the gradient.
 * @return The probability that the gradient is less than the threshold.
 */
double probability_gradient_less_than(const RunningGradient * const rg, const double thresh) {
    return normal_cdf(thresh, calculate_gradient(rg), calculate_standard_error(rg));
}

/**
 * @brief Calculates the probability that the gradient is greater than a threshold.
 *
 * This function calculates the probability that the gradient (slope) of the linear regression model
 * in the RunningGradient structure is greater than a specified threshold.
 *
 * @param rg Pointer to the RunningGradient structure.
 * @param thresh The threshold for the gradient.
 * @return The probability that the gradient is greater than the threshold.
 */
double probability_gradient_greater_than(const RunningGradient * const rg, const double thresh) {
    return 1.0 - probability_gradient_less_than(rg, thresh);
}

/**
 * @brief Finds the upper quantile in the data.
 *
 * This function calculates the quantile value below which a specified percentage of data points fall.
 * It uses the quickselect algorithm to efficiently find the quantile value without fully sorting the data.
 * The quantile value is determined by the formula: idx_upper = round((size - 1) * (1 - quantile)).
 * For example, a quantile of 0.9 will find the value below which 90% of the data points fall.
 *
 * @param data The data points.
 * @param size The number of data points.
 * @param quantile The quantile to find (must be in the range [0, 1]).
 * @return The value such that quantile percent of the values are greater than it.
 */
double find_upper_quantile(const double * const data, const size_t size, const double quantile) {
    assert(0 <= quantile && quantile <= 1.0 && "quantile must be in the range [0, 1]");
    assert(size > 0 && "The container must contain more than 0 elements.");

    double *sorted_data = (double *)malloc(size * sizeof(double));
    if (!sorted_data) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }
    for (size_t i = 0; i < size; ++i) {
        sorted_data[i] = data[i];
    }

    const size_t idx_upper = (size_t)round((size - 1) * (1 - quantile));
    const double upper_q = quickselect(sorted_data, 0, size - 1, idx_upper);

    free(sorted_data);
    return upper_q;
}

/**
 * @brief Calculates the probability of an increasing trend.
 *
 * This function calculates the probability that the gradient (slope) of the linear regression model
 * in the RunningGradient structure is greater than zero, indicating an increasing trend.
 *
 * @param rg Pointer to the RunningGradient structure.
 * @return The probability of an increasing trend.
 */
double calculate_probability(RunningGradient *rg) {
    if (rg->num_points > 2) {
        return probability_gradient_greater_than(rg, 0);
    }
    return 0.49;  // Default probability if not enough points
}

//#define DEBUG_LINEAR_GRADIENT

/**
 * @brief Compares gradients, incremental counts, and probabilities of increase for two parts of the data.
 *
 * This function calculates the gradients and probabilities of increase for two consecutive chunks of data (each of 15 points),
 * compares the total gradient increase, incremental counts, and probabilities, and determines which side has
 * a stronger increasing trend.
 *
 * @param data The array of data points.
 * @param start_index The starting index in the data array.
 * @param forgetting_factor The forgetting factor used in the RLS algorithm.
 * @return A structure containing the comparison result, including the dominant side and probabilities.
 */
GradientComparisonResult compare_gradient_parts(const double *data, size_t start_index, double forgetting_factor) {
    // Initialize the result structure with default values
    GradientComparisonResult result = {0.0, 0.0, 0, 0, 0.0, 0.0, UNDECIDED};
    
    // Initialize RunningGradient structures for the first (right) and second (left) parts
    RunningGradient rg_right_part;
    RunningGradient rg_left_part;

    init_running_gradient(&rg_right_part);
    init_running_gradient(&rg_left_part);

    size_t middle_index = start_index + 15;  // Middle of the window

    // Process the right side first (from middle_index to middle_index + 15)
    #ifdef DEBUG_LINEAR_GRADIENT
    printf("Processing the right side from %zu to %zu...\n", middle_index, middle_index + 15);
    #endif

    for (size_t i = 0; i < (WINDOW_SIZE / 2); ++i) {
        size_t index = middle_index + i;

        #ifdef DEBUG_LINEAR_GRADIENT
        printf("Adding data point %zu: %f to the right part\n", index, data[index]);
        #endif

        add_data_point(&rg_right_part, data[index]);

        if (rg_right_part.num_points > 1) {  // Ensure we have enough points to calculate a gradient
            const double gradient = calculate_gradient(&rg_right_part);

            #ifdef DEBUG_LINEAR_GRADIENT
            printf("Calculated gradient for right part at index %zu: %f\n", index, gradient);
            #endif

            result.total_gradient_first_part += gradient;

            if (gradient > 0.1f) { // TODO: HARDCODED. SHOULDNT BE HARDCODED. 
                result.increase_count_first_part++;

                #ifdef DEBUG_LINEAR_GRADIENT
                printf("Increase count for right part incremented. Current count: %zu\n", result.increase_count_first_part);
                #endif
            }
        }
    }

    // Calculate the probability of an increasing trend for the right part
    result.probability_increase_first_part = calculate_probability(&rg_right_part);

    #ifdef DEBUG_LINEAR_GRADIENT
    printf("Probability of increase for the right part: %f\n", result.probability_increase_first_part);
    #endif

    // Process the left side (from middle_index to middle_index - 15)
    #ifdef DEBUG_LINEAR_GRADIENT
    printf("Processing the left side from %zu to %zu...\n", middle_index, middle_index - 15);
    #endif

    for (size_t i = 0; i < (WINDOW_SIZE / 2); ++i) { // NOT HARDCODED ANYMORE
        size_t index = middle_index - i;

        #ifdef DEBUG_LINEAR_GRADIENT
        printf("Adding data point %zu: %f to the left part\n", index, data[index]);
        #endif

        add_data_point(&rg_left_part, data[index]);

        if (rg_left_part.num_points > 1) {  // Ensure we have enough points to calculate a gradient
            const double gradient = calculate_gradient(&rg_left_part);

            #ifdef DEBUG_LINEAR_GRADIENT
            printf("Calculated gradient for left part at index %zu: %f\n", index, gradient);
            #endif

            result.total_gradient_second_part += gradient;

            if (gradient > 0.1f) {
                result.increase_count_second_part++;

                #ifdef DEBUG_LINEAR_GRADIENT
                printf("Increase count for left part incremented. Current count: %zu\n", result.increase_count_second_part);
                #endif
            }
        }
    }

    // Calculate the probability of an increasing trend for the left part
    result.probability_increase_second_part = calculate_probability(&rg_left_part);

    #ifdef DEBUG_LINEAR_GRADIENT
    printf("Probability of increase for the left part: %f\n", result.probability_increase_second_part);
    #endif

    // Compare the total gradients and incremental counts to determine the dominant side
    #ifdef DEBUG_LINEAR_GRADIENT
    printf("Comparing the results...\n");
    #endif

    if (result.total_gradient_first_part > result.total_gradient_second_part && 
        result.increase_count_first_part > result.increase_count_second_part) {
        result.dominant_side = RIGHT_SIDE;

        #ifdef DEBUG_LINEAR_GRADIENT
        printf("The right side has a stronger increasing trend.\n");
        #endif
    } else if (result.total_gradient_second_part > result.total_gradient_first_part && 
               result.increase_count_second_part > result.increase_count_first_part) {
        result.dominant_side = LEFT_SIDE;

        #ifdef DEBUG_LINEAR_GRADIENT
        printf("The left side has a stronger increasing trend.\n");
        #endif
    } else {
        result.dominant_side = UNDECIDED;

        #ifdef DEBUG_LINEAR_GRADIENT
        printf("The trend is undecided or similar on both sides.\n");
        #endif
    }

    // Return the result structure with the calculated values
    return result;
}
