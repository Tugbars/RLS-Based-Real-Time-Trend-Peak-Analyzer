#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include "windowed_running_gradient.h"
#include "running_cubic_gradient.h"
#include "rls_analysis_parameters.h"

#include "running_peak_analysis.h"


/**
 * @brief Initializes the RunningCubicGradient structure for cubic regression using Recursive Least Squares (RLS).
 *
 * This function initializes the `RunningCubicGradient` structure by setting all coefficients to zero, initializing the
 * inverse covariance matrix with high uncertainty, and preparing the structure to start receiving data points.
 *
 * @param rg Pointer to the RunningCubicGradient structure to be initialized.
 * @param forgetting_factor The forgetting factor used in the RLS algorithm. Should be a value typically close to 1,
 *                          which determines the weight given to newer data points.
 *
 * The inverse covariance matrix is initialized with large values on the diagonal to represent high uncertainty in the initial state.
 * The coefficients and residual sum of squares are initialized to zero.
 */
void init_running_cubic_gradient(RunningCubicGradient *rg, double forgetting_factor) {
    rg->num_points = 0;
    rg->max_points = CUBIC_RLS_WINDOW;
    
    // Initialize coefficients to zero
    for (int i = 0; i < 4; ++i) {
        rg->coefficients[i] = 0.0;
    }
    
    // Initialize residual sum of squares to zero
    rg->residual_sum_squares = 0.0;
    
    // Set the forgetting factor
    rg->forgetting_factor = forgetting_factor;
    
    // Initialize the inverse covariance matrix
    memset(rg->inverse_cov_matrix, 0, sizeof(rg->inverse_cov_matrix));
    rg->inverse_cov_matrix[0][0] = 1e9; // High uncertainty in x^3 coefficient
    rg->inverse_cov_matrix[1][1] = 1e9; // High uncertainty in x^2 coefficient
    rg->inverse_cov_matrix[2][2] = 1e9; // High uncertainty in x coefficient
    rg->inverse_cov_matrix[3][3] = 1e9; // High uncertainty in constant term
}

/**
 * @brief Adds a new data point to the RunningCubicGradient structure and updates the cubic regression model using Recursive Least Squares (RLS).
 *
 * This function adds a new data point to the RunningCubicGradient structure, which maintains a cubic model of the data points
 * using the Recursive Least Squares (RLS) algorithm. The inverse covariance matrix is updated efficiently using the Sherman-Morrison formula,
 * which is well-suited for rank-1 updates. The function also updates the regression coefficients and recalculates the residual sum of squares
 * based on the new data point.
 *
 * @param rg Pointer to the RunningCubicGradient structure.
 * @param y The new data point to be added.
 *
 * The Sherman-Morrison formula is used to update the inverse covariance matrix:
 * \f[
 * A^{-1} \leftarrow A^{-1} - \frac{A^{-1} x x^T A^{-1}}{1 + x^T A^{-1} x}
 * \f]
 * where \( x \) is the input vector derived from the new data point.
 *
 * The regression coefficients are updated as follows:
 * \f[
 * \theta \leftarrow \theta + \text{inverse\_cov\_matrix} \times \text{prediction\_error}
 * \f]
 *
 * Symmetry in the inverse covariance matrix is enforced during the update process to maintain numerical stability and correctness.
 */
void add_cubic_data_point(RunningCubicGradient *const rg, const MqsRawDataPoint_t *data_point) {
    double y = data_point->phaseAngle; // Extract the phase angle from the data point

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
    double x0 = pow((double)(rg->num_points - 1), 3);  // Cubic term
    double x1 = pow((double)(rg->num_points - 1), 2);  // Quadratic term
    double x2 = (double)(rg->num_points - 1);  // Linear term
    double x3 = 1.0;  // Constant term

    // Calculate A^{-1} * x
    double tmp[4];
    tmp[0] = rg->inverse_cov_matrix[0][0] * x0 + rg->inverse_cov_matrix[0][1] * x1 + rg->inverse_cov_matrix[0][2] * x2 + rg->inverse_cov_matrix[0][3] * x3;
    tmp[1] = rg->inverse_cov_matrix[1][0] * x0 + rg->inverse_cov_matrix[1][1] * x1 + rg->inverse_cov_matrix[1][2] * x2 + rg->inverse_cov_matrix[1][3] * x3;
    tmp[2] = rg->inverse_cov_matrix[2][0] * x0 + rg->inverse_cov_matrix[2][1] * x1 + rg->inverse_cov_matrix[2][2] * x2 + rg->inverse_cov_matrix[2][3] * x3;
    tmp[3] = rg->inverse_cov_matrix[3][0] * x0 + rg->inverse_cov_matrix[3][1] * x1 + rg->inverse_cov_matrix[3][2] * x2 + rg->inverse_cov_matrix[3][3] * x3;

    // Calculate the denominator (1 + V * A^{-1} * U)
    double denom = rg->forgetting_factor + (x0 * tmp[0] + x1 * tmp[1] + x2 * tmp[2] + x3 * tmp[3]);

    // Check for numerical stability
    if (fabs(denom) < 1e-10) {
        printf("Numerical stability issue: denom is too close to zero.\n");
    }

    // Update the inverse covariance matrix using the Woodbury identity

    // First row and first column
    rg->inverse_cov_matrix[0][0] = (rg->inverse_cov_matrix[0][0] - (tmp[0] * tmp[0]) / denom) / rg->forgetting_factor;
    rg->inverse_cov_matrix[0][1] = (rg->inverse_cov_matrix[0][1] - (tmp[0] * tmp[1]) / denom) / rg->forgetting_factor;
    rg->inverse_cov_matrix[0][2] = (rg->inverse_cov_matrix[0][2] - (tmp[0] * tmp[2]) / denom) / rg->forgetting_factor;
    rg->inverse_cov_matrix[0][3] = (rg->inverse_cov_matrix[0][3] - (tmp[0] * tmp[3]) / denom) / rg->forgetting_factor;

    // Second row and second column
    rg->inverse_cov_matrix[1][1] = (rg->inverse_cov_matrix[1][1] - (tmp[1] * tmp[1]) / denom) / rg->forgetting_factor;
    rg->inverse_cov_matrix[1][2] = (rg->inverse_cov_matrix[1][2] - (tmp[1] * tmp[2]) / denom) / rg->forgetting_factor;
    rg->inverse_cov_matrix[1][3] = (rg->inverse_cov_matrix[1][3] - (tmp[1] * tmp[3]) / denom) / rg->forgetting_factor;

    // Third row and third column
    rg->inverse_cov_matrix[2][2] = (rg->inverse_cov_matrix[2][2] - (tmp[2] * tmp[2]) / denom) / rg->forgetting_factor;
    rg->inverse_cov_matrix[2][3] = (rg->inverse_cov_matrix[2][3] - (tmp[2] * tmp[3]) / denom) / rg->forgetting_factor;

    // Fourth row and fourth column
    rg->inverse_cov_matrix[3][3] = (rg->inverse_cov_matrix[3][3] - (tmp[3] * tmp[3]) / denom) / rg->forgetting_factor;

    // Enforce symmetry in the covariance matrix
    rg->inverse_cov_matrix[1][0] = rg->inverse_cov_matrix[0][1];
    rg->inverse_cov_matrix[2][0] = rg->inverse_cov_matrix[0][2];
    rg->inverse_cov_matrix[2][1] = rg->inverse_cov_matrix[1][2];
    rg->inverse_cov_matrix[3][0] = rg->inverse_cov_matrix[0][3];
    rg->inverse_cov_matrix[3][1] = rg->inverse_cov_matrix[1][3];
    rg->inverse_cov_matrix[3][2] = rg->inverse_cov_matrix[2][3];

    // Calculate the prediction error
    double prediction_error = y - (x0 * rg->coefficients[0] + x1 * rg->coefficients[1] + x2 * rg->coefficients[2] + x3 * rg->coefficients[3]);

    // Update the regression coefficients
    rg->coefficients[0] += tmp[0] * prediction_error / denom;
    rg->coefficients[1] += tmp[1] * prediction_error / denom;
    rg->coefficients[2] += tmp[2] * prediction_error / denom;
    rg->coefficients[3] += tmp[3] * prediction_error / denom;

    // Recalculate the residual sum of squares
    rg->residual_sum_squares = 0.0;
    for (uint16_t i = 0; i < rg->num_points; ++i) {
        double error = rg->y[i] - (rg->coefficients[0] * pow(rg->x[i], 3) + rg->coefficients[1] * pow(rg->x[i], 2) + rg->coefficients[2] * rg->x[i] + rg->coefficients[3]);
        rg->residual_sum_squares += error * error;
    }
}

/**
 * @brief Calculates the first derivative (slope) of the cubic function at a specific point.
 *
 * This function computes the first derivative of the cubic equation represented by the
 * coefficients stored in the RunningCubicGradient structure at a given value of x.
 *
 * @param rg Pointer to the RunningCubicGradient structure containing the coefficients.
 * @param x The point at which to calculate the first derivative (slope).
 * @return The slope of the cubic function at the given point x.
 */
static double calculate_first_derivative(const RunningCubicGradient *rg, double x) {
    double a3 = rg->coefficients[0];  // Coefficient for x^3
    double a2 = rg->coefficients[1];  // Coefficient for x^2
    double a1 = rg->coefficients[2];  // Coefficient for x

    // Calculate the first derivative at the given point x
    double first_derivative = 3.0 * a3 * pow(x, 2) + 2.0 * a2 * x + a1;
    return first_derivative;
}

/**
 * @brief Calculates the second-order gradient (curvature) of the cubic function.
 *
 * This function computes the second-order gradient of the cubic equation represented by the
 * coefficients stored in the RunningCubicGradient structure. The second-order gradient is
 * the second derivative of the cubic function:
 *
 *     second_order_gradient = 6 * a3 * x + 2 * a2
 *
 * where:
 *     a3 is the coefficient of the x^3 term.
 *     a2 is the coefficient of the x^2 term.
 *
 * If you only want to calculate the curvature, independent of x, you can
 * simplify this calculation by evaluating it at any given point of interest.
 *
 * @param rg Pointer to the RunningCubicGradient structure containing the coefficients.
 * @return The second-order gradient (curvature) of the cubic function.
 */
double calculate_cubic_second_order_gradient(const RunningCubicGradient *rg) {
    double a3 = rg->coefficients[0];  // Coefficient for x^3
    double a2 = rg->coefficients[1];  // Coefficient for x^2

    // If you're evaluating at a particular x, you can include it, but for now:
    // We return a general expression for second-order without any x-dependence.
    double second_order_gradient = 6.0 * a3 + 2.0 * a2;
    return second_order_gradient;
}

/**
 * @brief Computes the second-order gradients (curvature) for the next 30 values in the given array using a cubic regression model.
 *
 * This function takes an array of doubles and a starting index, and then adds each of the next 30 values
 * to the RunningCubicGradient structure. After each addition, it computes the second-order gradient
 * (curvature) of the cubic model. The function returns an array containing these second-order gradients.
 *
 * The function can also print debug information if the `DEBUG_CUBIC_GRADIENT_CALC` macro is defined.
 *
 * @param second_order_gradients Array to store the computed second-order gradients.
 * @param values Array of double values from which the gradients are computed.
 * @param length Length of the array.
 * @param start_index Starting index from which to begin the gradient calculation.
 * @param forgetting_factor The forgetting factor for the RLS algorithm.
 */
void compute_cubic_second_order_gradients(double *second_order_gradients, const MqsRawDataPoint_t *values, uint16_t length, uint16_t start_index, double forgetting_factor) {
    // Ensure that the start index and the length allow for calculations
    if (start_index + CUBIC_RLS_WINDOW > length) {
        printf("Insufficient data to compute second-order gradients from the starting index.\n");
        return;
    }

    // Initialize the running cubic gradient structure
    RunningCubicGradient rg;
    init_running_cubic_gradient(&rg, forgetting_factor);

    // Loop through the next values starting from start_index
    for (uint16_t i = 0; i < CUBIC_RLS_WINDOW; ++i) {
        uint16_t current_index = start_index + i;
        const MqsRawDataPoint_t *current_value = &values[current_index];

        // Add the current value's phaseAngle to the running cubic gradient
        add_cubic_data_point(&rg, current_value);

        // Calculate the second-order gradient (curvature) only if we have at least 4 data points
        if (rg.num_points >= 4) {
            double second_order_gradient = calculate_cubic_second_order_gradient(&rg);
            second_order_gradients[i] = second_order_gradient;

            #ifdef DEBUG_CUBIC_GRADIENT_CALC
            // Print the current second-order gradient and its index
            printf("Second-order gradient at rg index %u after adding phase angle %.6f: %.6f\n", i, current_value->phaseAngle, second_order_gradient);
            #endif
        } else {
            second_order_gradients[i] = NAN; // Not enough points yet to calculate the gradient

            #ifdef DEBUG_CUBIC_GRADIENT_CALC
            // Print a message indicating insufficient data and its index
            printf("Not enough data points to calculate the gradient at rg index %u after adding phase angle %.6f.\n", i, current_value->phaseAngle);
            #endif
        }
    }
}

/**
 * @brief Computes the first-order gradients (slopes) for a given window size starting from a specified index.
 *
 * This function calculates the first-order gradients (slopes) of the cubic model at each point within a given window.
 * It uses the cubic regression model to derive the gradients by computing the first derivative at each point.
 *
 * @param first_order_gradients Array to store the computed first-order gradients.
 * @param values Array of double values representing the data.
 * @param length Length of the data array.
 * @param start_index Starting index from which to begin the gradient calculation.
 * @param window_size The number of points to include in the gradient calculation.
 * @param forgetting_factor The forgetting factor used in the RLS algorithm.
 */
void compute_cubic_first_order_gradients(double *first_order_gradients, const MqsRawDataPoint_t *values, uint16_t length, uint16_t start_index, uint16_t window_size, double forgetting_factor) {
    // Ensure that the start index and the window size allow for calculations
    if (start_index + window_size > length) {
        printf("Insufficient data to compute first-order gradients for the specified window size.\n");
        return;
    }

    // Initialize the running cubic gradient structure
    RunningCubicGradient rg;
    init_running_cubic_gradient(&rg, forgetting_factor);

    double cumulative_sum = 0.0;  // Variable to store the cumulative sum of first-order gradients

    // Loop through the specified window size
    for (uint16_t i = 0; i < window_size; ++i) {
        uint16_t current_index = start_index + i;
        const MqsRawDataPoint_t *current_value = &values[current_index];

        // Add the current value's phaseAngle to the running cubic gradient
        add_cubic_data_point(&rg, current_value);

        // Calculate the first-order gradient (slope) only if we have at least 4 data points
        if (rg.num_points >= 4) {
            double x_current = (double)(rg.num_points - 1);  // Current x value
            double first_order_gradient = calculate_first_derivative(&rg, x_current);
            first_order_gradients[i] = first_order_gradient;

            // Update the cumulative sum of first-order gradients
            cumulative_sum += first_order_gradient;

            // Debugging: Print the calculated first-order gradient
            printf("First-order gradient at rg index %u after adding phase angle %.6f: %.6f\n", start_index + i, current_value->phaseAngle, first_order_gradient);
            printf("Cumulative sum of first-order gradients up to index %u: %.6f\n", start_index + i, cumulative_sum);
        } else {
            first_order_gradients[i] = NAN; // Not enough points yet to calculate the gradient

            #ifdef DEBUG_CUBIC_GRADIENT_CALC
            // Print a message indicating insufficient data and its index
            printf("Not enough data points to calculate the gradient at rg index %u after adding phase angle %.6f.\n", i, current_value->phaseAngle);
            #endif
        }
    }

    #ifdef DEBUG_CUBIC_GRADIENT_CALC
    // Print the final cumulative sum after the loop is done
    printf("Final cumulative sum of first-order gradients: %.6f\n", cumulative_sum);
    #endif
}

/**
 * @brief Finds a consistent increasing trend in the given gradient array.
 *
 * This function scans through an array of gradient values to identify a consistent increasing trend.
 * It starts tracking when positive gradients are encountered and stops if more than one negative gradient is found consecutively.
 * The function returns a `GradientTrendIndices` struct containing the start and end indices of the consistent increase, 
 * a validity flag, and the cumulative sum of the gradients within this range.
 *
 * @param gradients The array of gradient values to analyze.
 * @param start_index The starting index in the original data array from which the gradient array begins.
 * @param window_size The number of elements in the gradient array to analyze.
 * @return GradientTrendIndices A structure containing the start index, end index, validity, and maximum sum of the consistent increasing trend.
 */
GradientTrendIndices find_consistent_increase(double *gradients, size_t start_index, size_t window_size) {
    // Initialize the structure to store information about the consistent increase.
    GradientTrendIndices increase_info = {0, 0, false, 0.0};
    bool tracking_increase = false;  // Flag to track if we're currently observing a consistent increase.
    double cumulative_sum = 0.0;  // Variable to accumulate the sum of positive gradients.
    int decrease_count = 0;  // Counter to track consecutive decreases.

    #ifdef CONSISTENT_TREND_DEBUG
    printf("Starting find_consistent_increase\n");
    #endif

    // Iterate through the gradient array within the specified window.
    for (size_t i = 0; i < window_size; ++i) {
        double gradient = gradients[i];

        #ifdef CONSISTENT_TREND_DEBUG
        printf("Index: %zu, Gradient: %.6f, Cumulative Sum: %.6f, Decrease Count: %d\n", start_index + i, gradient, cumulative_sum, decrease_count);
        #endif

        // Check if the current gradient is positive.
        if (gradient > 0) {
            // If we are not already tracking an increase, start now.
            if (!tracking_increase) {
                increase_info.start_index = start_index + i;  // Set the start index for the increase.
                tracking_increase = true;  // Begin tracking the increase.
                cumulative_sum = 0.0;  // Reset the cumulative sum.
                
                #ifdef CONSISTENT_TREND_DEBUG
                printf("Started tracking increase at index %zu\n", start_index + i);
                #endif
            }
            // Update the end index for the increase and add the gradient to the cumulative sum.
            increase_info.end_index = start_index + i;
            cumulative_sum += gradient;
            decrease_count = 0;  // Reset the decrease counter.
            
            #ifdef CONSISTENT_TREND_DEBUG
            printf("Tracking increase - Updated end_index to %zu, Cumulative Sum: %.6f\n", start_index + i, cumulative_sum);
            #endif

        } else if (gradient < 0) {  // If the current gradient is negative.
            decrease_count++;  // Increment the decrease counter.
            
            #ifdef CONSISTENT_TREND_DEBUG
            printf("Negative gradient found - Decrease Count: %d\n", decrease_count);
            #endif
            
            // If the number of consecutive decreases exceeds the allowed limit, stop tracking the increase.
            if (decrease_count > cubic_analysis_params.max_third_order_trend_decrease_count) {
                tracking_increase = false;  // Stop tracking the increase.
                
                #ifdef CONSISTENT_TREND_DEBUG
                printf("Stopped tracking increase due to consecutive decreases at index %zu\n", start_index + i);
                #endif
                break;  // Exit the loop early as the increase is no longer valid.
            }
        }
    }

    // If we have accumulated a positive sum of gradients, mark the increase as valid.
    if (cumulative_sum > 0) {
        increase_info.valid = true;
        increase_info.max_sum = cumulative_sum;
        
        #ifdef CONSISTENT_TREND_DEBUG
        printf("Valid increase found from index %zu to %zu with max sum %.6f\n", increase_info.start_index, increase_info.end_index, increase_info.max_sum);
        #endif
    } else {
        #ifdef CONSISTENT_TREND_DEBUG
        printf("No valid increase found.\n");
        #endif
    }

    return increase_info;  // Return the structure containing the details of the increase.
}


/**
 * @brief Finds a consistent decreasing trend in the given gradient array.
 *
 * This function scans through an array of gradient values to identify a consistent decreasing trend.
 * It starts tracking when negative gradients are encountered and stops if more than one positive gradient is found consecutively.
 * The function returns a `GradientTrendIndices` struct containing the start and end indices of the consistent decrease, 
 * a validity flag, and the cumulative sum of the gradients within this range.
 *
 * @param gradients The array of gradient values to analyze.
 * @param start_index The starting index in the original data array from which the gradient array begins.
 * @param window_size The number of elements in the gradient array to analyze.
 * @return GradientTrendIndices A structure containing the start index, end index, validity, and maximum sum of the consistent decreasing trend.
 */
GradientTrendIndices find_consistent_decrease(double *gradients, size_t start_index, size_t window_size) {
    GradientTrendIndices decrease_info = {0, 0, false, 0.0};
    bool tracking_decrease = false;
    double cumulative_sum = 0.0;
    int increase_count = 0;

    #ifdef CONSISTENT_TREND_DEBUG
    printf("Starting find_consistent_decrease\n");
    #endif

    for (size_t i = 0; i < window_size; ++i) {
        double gradient = gradients[i];
        
        #ifdef CONSISTENT_TREND_DEBUG
        printf("Index: %zu, Gradient: %.6f, Cumulative Sum: %.6f, Increase Count: %d\n", start_index + i, gradient, cumulative_sum, increase_count);
        #endif

        if (gradient < 0) {  // If the current gradient is negative
            if (!tracking_decrease) {
                // Start tracking the decrease
                decrease_info.start_index = start_index + i;
                tracking_decrease = true;
                cumulative_sum = 0.0;  // Reset cumulative sum when starting to track
                #ifdef CONSISTENT_TREND_DEBUG
                printf("Started tracking decrease at index %zu\n", start_index + i);
                #endif
            }
            decrease_info.end_index = start_index + i;
            cumulative_sum += gradient;
            increase_count = 0;  // Reset increase counter on negative gradient
            
            #ifdef CONSISTENT_TREND_DEBUG
            printf("Tracking decrease - Updated end_index to %zu, Cumulative Sum: %.6f\n", start_index + i, cumulative_sum);
            #endif

        } else if (gradient > 0) {  // If the current gradient is positive
            increase_count++;
            #ifdef CONSISTENT_TREND_DEBUG
            printf("Positive gradient found - Increase Count: %d\n", increase_count);
            #endif
            if (increase_count > cubic_analysis_params.max_third_order_trend_increase_count) {
                tracking_decrease = false;  // Stop tracking
                #ifdef CONSISTENT_TREND_DEBUG
                printf("Stopped tracking decrease due to consecutive increases at index %zu\n", start_index + i);
                #endif
            }
        }

        if (!tracking_decrease && gradient < 0 && increase_count > cubic_analysis_params.max_third_order_trend_increase_count) {
            // Restart tracking from this new decrease
            decrease_info.start_index = start_index + i;
            tracking_decrease = true;
            cumulative_sum = gradient;  // Reset cumulative sum with the new decrease
            decrease_info.end_index = start_index + i;
            increase_count = 0;  // Reset the increase counter
            #ifdef CONSISTENT_TREND_DEBUG
            printf("Restarted tracking decrease at index %zu, Cumulative Sum: %.6f\n", start_index + i, cumulative_sum);
            #endif
        }
    }

    if (cumulative_sum < 0) {  // Ensure there's a valid decrease before marking it valid
        decrease_info.valid = true;
        decrease_info.max_sum = cumulative_sum;
        #ifdef CONSISTENT_TREND_DEBUG
        printf("Valid decrease found from index %zu to %zu with max sum %.6f\n", decrease_info.start_index, decrease_info.end_index, decrease_info.max_sum);
        #endif
    } else {
        #ifdef CONSISTENT_TREND_DEBUG
        printf("No valid decrease found.\n");
        #endif
    }

    return decrease_info;
}

/**
 * @brief Tracks the gradient trends within a specified window using cubic regression.
 *
 * This function tracks the gradient trends of a given dataset within a specified window using cubic regression.
 * It calculates first-order gradients (slopes) at each point in the window, identifies regions of consistent increase,
 * and stores the maximum cumulative sum of gradients. It returns a structure containing the start and end indices
 * of the consistent increase, along with a flag indicating if the trend was valid and the maximum sum of gradients.
 *
 * @param first_order_gradients Array to store the computed first-order gradients for each point in the window.
 * @param values Array of double values representing the data points.
 * @param length The length of the values array.
 * @param start_index The starting index in the values array from which to begin the gradient calculation.
 * @param window_size The number of points to include in the gradient calculation.
 * @param forgetting_factor The forgetting factor used in the Recursive Least Squares (RLS) algorithm, typically close to 1,
 *        which determines the weight given to newer data points.
 * 
 * @return GradientTrendIndices A struct containing:
 * - start_index: The index where the consistent increase in gradients begins.
 * - end_index: The index where the consistent increase in gradients ends.
 * - valid: A boolean flag indicating whether a valid trend was detected.
 * - max_sum: The maximum cumulative sum of gradients during the consistent increase.
 *
 * ### Mathematical Explanation:
 * - **First-order Gradient Calculation:** 
 *   - The first-order gradient (slope) is calculated as the derivative of the cubic regression model:
 *   \[
 *   \text{first\_order\_gradient} = 3a_3x^2 + 2a_2x + a_1
 *   \]
 *   where \(a_3\), \(a_2\), and \(a_1\) are the coefficients of the cubic, quadratic, and linear terms, respectively.
 *
 * - **Cumulative Sum:** 
 *   - The function accumulates the sum of first-order gradients to track the overall trend within the window. A consistent positive
 *     cumulative sum indicates a strong upward trend.
 *
 * - **Trend Tracking Logic:**
 *   - The function identifies a trend as consistent if the first-order gradients are positive over a series of points.
 *   - It ignores one decrease in the trend to account for minor fluctuations.
 *   - If the cumulative sum is greater than the previous maximum sum, it updates the maximum cumulative sum and marks the trend as valid.
 */
GradientTrendResult track_gradient_trends_with_cubic_regression(const MqsRawDataPoint_t *values, uint16_t length, uint16_t start_index, uint16_t window_size, double forgetting_factor) {
    GradientTrendResult trend_result = {0}; // Initialize the struct with default values

    // Ensure that the start index and the window size allow for calculations
    if (start_index + window_size > length) {
        printf("Insufficient data to compute first-order gradients for the specified window size.\n");
        return trend_result;
    }

    // Array to store the first-order gradients
    double first_order_gradients[window_size];

    // Initialize the running cubic gradient structure
    RunningCubicGradient rg;
    init_running_cubic_gradient(&rg, forgetting_factor);

    // Loop through the specified window size and calculate first-order gradients
    for (uint16_t i = 0; i < window_size; ++i) {
        uint16_t current_index = start_index + i;
        const MqsRawDataPoint_t *current_value = &values[current_index];

        // Add the current value's phaseAngle to the running cubic gradient
        add_cubic_data_point(&rg, current_value);

        // Check if we have enough points to calculate the gradient
        if (rg.num_points >= 4) {
            double x_current = (double)(rg.num_points - 1);  // Current x value
            double first_order_gradient = calculate_first_derivative(&rg, x_current);
            first_order_gradients[i] = first_order_gradient;

            // Debugging: Print the calculated first-order gradient
            //printf("Calculated first-order gradient at index %u: %.6f\n", current_index, first_order_gradient);
        } else {
            first_order_gradients[i] = NAN; // Not enough points yet to calculate the gradient
            // Debugging: Print a message indicating insufficient data points
           // printf("Not enough data points to calculate the gradient at index %u.\n", current_index);
        }
    }

    // Debugging: Print a message before finding consistent increase
    //printf("Finding consistent increase...\n");

    // Find the consistent increase trend
    trend_result.increase_info = find_consistent_increase(first_order_gradients, start_index, window_size);

    // Debugging: Print a message before finding consistent decrease
    //printf("Finding consistent decrease...\n");

    // Find the consistent decrease trend
    trend_result.decrease_info = find_consistent_decrease(first_order_gradients, start_index, window_size);

    return trend_result;
}

/**
 * @brief Evaluates if the detected gradient trend is significant.
 *
 * This function analyzes the gradient trends detected by the `track_gradient_trends_with_cubic_regression` function
 * and determines whether the trend is significant. A trend is considered significant if the cumulative sum of gradients
 * exceeds a threshold (11.0) and the difference between the start and end indices is greater than 5.
 *
 * @param values Array of double values representing the data points.
 * @param length The length of the values array.
 * @param start_index The starting index in the values array from which to begin the gradient calculation.
 * @param window_size The number of points to include in the gradient calculation.
 * @param forgetting_factor The forgetting factor used in the Recursive Least Squares (RLS) algorithm.
 * 
 * @return bool Returns true if a significant gradient trend is detected, otherwise false.
 *
 * ### Logic Behind the Evaluation:
 * - **Threshold Evaluation:** 
 *   - The function checks if the maximum cumulative sum of gradients (calculated by `track_gradient_trends_with_cubic_regression`)
 *     exceeds a predefined threshold (11.0). This threshold is used to filter out minor trends that are not significant.
 * 
 * - **Index Difference Criterion:**
 *   - The function also verifies that the distance between the start and end indices of the trend is greater than 5 points. This ensures
 *     that the trend is sustained over a meaningful period.
 * 
 * - **Decision Logic:**
 *   - If both the threshold and index difference criteria are met, the function concludes that the trend is significant.
 */
PeakTrendAnalysisResult detect_significant_gradient_trends(const MqsRawDataPoint_t *values, uint16_t length, uint16_t start_index, uint16_t window_size, double forgetting_factor) {
    // Call the function to track gradient trends
    GradientTrendResult trend_result = track_gradient_trends_with_cubic_regression(values, length, start_index, window_size, forgetting_factor);

    // Initialize the result struct
    PeakTrendAnalysisResult result = {0};
    result.increase_info = trend_result.increase_info;
    result.decrease_info = trend_result.decrease_info;

    double significance_threshold = cubic_analysis_params.significance_threshold;
    uint16_t duration_threshold = cubic_analysis_params.duration_threshold;

    double average_increase = 0.0;
    double average_decrease = 0.0;

    // Calculate average increase if a valid increase trend is detected
    if (trend_result.increase_info.valid) {
        double sum_of_gradients = trend_result.increase_info.max_sum;
        uint16_t index_difference = trend_result.increase_info.end_index - trend_result.increase_info.start_index;
        average_increase = sum_of_gradients / (double)index_difference;

        // Evaluate significance based on global threshold and new average conditions
        if ((sum_of_gradients > significance_threshold && index_difference > duration_threshold) ||
            (index_difference >= peak_analysis_params.min_consistent_trend_count && average_increase > peak_analysis_params.min_average_increase)) {
            result.significant_increase = true;
        }
    }

    // Calculate average decrease if a valid decrease trend is detected
    if (trend_result.decrease_info.valid) {
        double sum_of_gradients = trend_result.decrease_info.max_sum;
        uint16_t index_difference = trend_result.decrease_info.end_index - trend_result.decrease_info.start_index;
        average_decrease = sum_of_gradients / (double)index_difference;

        // Evaluate significance based on global threshold and new average conditions
        if ((sum_of_gradients < -significance_threshold && index_difference > duration_threshold) ||
            (index_difference >= peak_analysis_params.min_consistent_trend_count && average_decrease < peak_analysis_params.min_average_decrease)) {
            result.significant_decrease = true;
        }
    }

    // Debugging: Print the final result summary
    printf("Trend detection complete. Significant Increase: %s, Significant Decrease: %s\n",
           result.significant_increase ? "Yes" : "No", result.significant_decrease ? "Yes" : "No");

    if (trend_result.increase_info.valid) {
        printf("Average Increase: %.6f over interval [%u - %u]\n", 
               average_increase, trend_result.increase_info.start_index, trend_result.increase_info.end_index);
    }

    if (trend_result.decrease_info.valid) {
        printf("Average Decrease: %.6f over interval [%u - %u]\n", 
               average_decrease, trend_result.decrease_info.start_index, trend_result.decrease_info.end_index);
    }

    return result;
}

//ilk average increase or decrease meselesini halletmemiz lazım. 
