#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <assert.h>
#include "windowed_running_gradient.h"
#include "rls_analysis_parameters.h"


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

    for (uint16_t i = 0; i < WINDOW_SIZE; ++i) {
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
static inline uint16_t partition(double* arr, uint16_t left, uint16_t right, uint16_t pivotIndex) {
    double pivotValue = arr[pivotIndex];
    double temp = arr[pivotIndex];
    arr[pivotIndex] = arr[right];
    arr[right] = temp;  // Move pivot to end
    uint16_t storeIndex = left;
    for (uint16_t i = left; i < right; i++) {
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
static inline double quickselect(double* arr, uint16_t left, uint16_t right, uint16_t k) {
    if (left == right) {
        return arr[left];  // If the list contains only one element, return that element
    }
    uint16_t pivotIndex = left + (right - left) / 2;
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
 * prediction error. The inverse covariance matrix is used to _probabili the adjustment to the coefficients.
 *
 * RLS is particularly suitable for running gradient calculations because it efficiently updates the model parameters
 * without the need to reprocess all previous data points. This is essential for real-time applications where data is
 * continuously arriving, and the gradient needs to be updated dynamically.
 *
 * @param rg Pointer to the RunningGradient structure.
 * @param y The new data point value.
 */
void add_data_point(RunningGradient *const rg, const MqsRawDataPoint_t *data_point) {
    // Extract the phase angle from the data point
    double y = data_point->phaseAngle;

    // Check if the current window is not full
    if (rg->num_points < rg->max_points) {
        // Add new data point to the window
        rg->y[rg->num_points] = y;
        rg->x[rg->num_points] = (double)rg->num_points;
        rg->num_points++;
    } else {
        // Shift the window to the left and add new data point
        for (uint16_t i = 0; i < rg->max_points - 1; ++i) {
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
    for (uint16_t i = 0; i < rg->num_points; ++i) {
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
double find_upper_quantile(const double * const data, const uint16_t size, const double quantile) {
    assert(0 <= quantile && quantile <= 1.0 && "quantile must be in the range [0, 1]");
    assert(size > 0 && "The container must contain more than 0 elements.");

    double *sorted_data = (double *)malloc(size * sizeof(double));
    if (!sorted_data) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }
    for (uint16_t i = 0; i < size; ++i) {
        sorted_data[i] = data[i];
    }

    const uint16_t idx_upper = (uint16_t)round((size - 1) * (1 - quantile));
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
 * @brief Processes the gradient for a part of the window (left or right).
 *
 * This function adds data points to the running gradient, calculates the gradient for each point,
 * and updates the total gradient and the count of increasing gradients based on the globally configured threshold.
 *
 * @param rg Pointer to the RunningGradient structure.
 * @param data Pointer to the array of MqsRawDataPoint_t containing the data points.
 * @param start_index The starting index for processing data points.
 * @param end_index The ending index for processing data points.
 * @param total_gradient Pointer to the variable storing the total gradient for this part.
 * @param increase_count Pointer to the variable storing the count of increases in this part.
 */
void process_gradient(RunningGradient *rg, const MqsRawDataPoint_t *data, uint16_t start_index, uint16_t end_index, 
                      double *total_gradient, uint16_t *increase_count) {
    for (uint16_t i = start_index; i < end_index; ++i) {
        add_data_point(rg, &data[i]);

        if (rg->num_points > 1) {  // Ensure we have enough points to calculate a gradient
            double gradient = calculate_gradient(rg);
            *total_gradient += gradient;

            if (gradient > gradient_analysis_params.gradient_threshold) {  // Use globally defined threshold
                (*increase_count)++;
            }
        }
    }
}

/**
 * @brief Compares the total gradients and incremental counts to determine the dominant side.
 *
 * This function evaluates the total gradients and increase counts for both parts of the window (left and right)
 * to determine which side has the stronger trend, either increasing or decreasing. If both sides have gradients
 * that do not meet the minimum threshold, the result will be set as UNDECIDED.
 *
 * @param result Pointer to the GradientComparisonResult structure containing the gradient comparison data.
 * @return PeakPosition indicating whether the left side, right side, or neither side (undecided) has a dominant trend.
 */
PeakPosition determine_dominant_side(const GradientComparisonResult *result) {
    double min_total = gradient_analysis_params.minimum_gradient_total;

    // Check if both gradients are insignificant
    if ((result->total_gradient_first_part > -min_total && result->total_gradient_first_part < min_total) &&
        (result->total_gradient_second_part > -min_total && result->total_gradient_second_part < min_total)) {
        return UNDECIDED;  // Both sides have insignificant gradients
    }

    // Check for dominant increase
    if (result->increase_count_first_part > result->increase_count_second_part) {
        return LEFT_SIDE;  // Left side has a stronger increase
    }
    if (result->increase_count_second_part > result->increase_count_first_part) {
        return RIGHT_SIDE;  // Right side has a stronger increase
    }

    // Check for decreases if no dominant increase
    if (result->total_gradient_first_part < result->total_gradient_second_part) {
        return RIGHT_SIDE;  // Right side has a stronger decrease
    }
    if (result->total_gradient_second_part < result->total_gradient_first_part) {
        return LEFT_SIDE;   // Left side has a stronger decrease
    }

    return UNDECIDED;  // No clear dominant side
}

/**
 * @brief Compares the gradient trends between two parts of the window (left and right) to determine the dominant side.
 *
 * This function processes the gradients for the left and right sides of the window, calculates the total gradients and 
 * increase counts, and compares them to determine whether the left side, right side, or neither side has a dominant trend.
 *
 * @param data Pointer to the array of MqsRawDataPoint_t containing the data points.
 * @param start_index The starting index for processing the data window.
 * @param forgetting_factor The forgetting factor used in the gradient calculations.
 * @return GradientComparisonResult A structure containing the results of the comparison, including total gradients, 
 *                                  increase counts, probabilities, and the dominant side.
 */
GradientComparisonResult compare_gradient_parts(const MqsRawDataPoint_t *data, uint16_t start_index, double forgetting_factor) {
    GradientComparisonResult result = {0.0, 0.0, 0, 0, 0.0, 0.0, UNDECIDED};
    
    RunningGradient rg_left_part, rg_right_part;
    init_running_gradient(&rg_left_part);
    init_running_gradient(&rg_right_part);

    uint16_t middle_index = start_index + WINDOW_SIZE / 2;  // Middle of the window

    // Process the left side (first part)
    //printf("Processing the left side (first part) from %u to %u...\n", middle_index - (WINDOW_SIZE / 2), middle_index);
    process_gradient(&rg_left_part, data, middle_index - (WINDOW_SIZE / 2), middle_index, 
                     &result.total_gradient_first_part, &result.increase_count_first_part);
    result.probability_increase_first_part = calculate_probability(&rg_left_part);

    // Process the right side (second part)
    //printf("Processing the right side (second part) from %u to %u...\n", middle_index, middle_index + (WINDOW_SIZE / 2));
    process_gradient(&rg_right_part, data, middle_index, middle_index + (WINDOW_SIZE / 2), 
                     &result.total_gradient_second_part, &result.increase_count_second_part);
    result.probability_increase_second_part = calculate_probability(&rg_right_part);

    // Compare the results to determine the dominant side
    //printf("Comparing the results...\n");
    result.dominant_side = determine_dominant_side(&result);

    switch (result.dominant_side) {
        case LEFT_SIDE:
            printf("The left side has a stronger trend.\n");
            break;
        case RIGHT_SIDE:
            printf("The right side has a stronger trend.\n");
            break;
        case UNDECIDED:
        default:
            printf("The trend is undecided or similar on both sides.\n");
            break;
    }

    return result;
}

