#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include "trend_analysis.h"
#include "windowed_running_gradient.h"
#include "running_gradient_parameters.h"
#include "running_cubic_gradient.h"



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

void add_cubic_data_point(RunningCubicGradient *const rg, const double y) {
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
        fprintf(stderr, "Numerical stability issue: denom is too close to zero.\n");
        exit(EXIT_FAILURE);
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
    for (size_t i = 0; i < rg->num_points; ++i) {
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
 * @param values Array of double values.
 * @param length Length of the array.
 * @param start_index Starting index from which to begin the gradient calculation.
 * @param forgetting_factor The forgetting factor for the RLS algorithm.
 * @return An array of double containing the second-order gradients.
 */
void compute_cubic_second_order_gradients(double *second_order_gradients, const double *values, size_t length, size_t start_index, double forgetting_factor) {
    // Ensure that the start index and the length allow for calculations
    if (start_index + CUBIC_RLS_WINDOW > length) {
        printf("Insufficient data to compute second-order gradients from the starting index.\n");
        return;
    }

    // Initialize the running cubic gradient structure
    RunningCubicGradient rg;
    init_running_cubic_gradient(&rg, forgetting_factor);

    // Loop through the next values starting from start_index
    for (size_t i = 0; i < CUBIC_RLS_WINDOW; ++i) {
        size_t current_index = start_index + i;
        double current_value = values[current_index];

        // Add the current value to the running cubic gradient
        add_cubic_data_point(&rg, current_value);

        // Calculate the second-order gradient (curvature) only if we have at least 4 data points
        if (rg.num_points >= 4) {
            double second_order_gradient = calculate_cubic_second_order_gradient(&rg);
            second_order_gradients[i] = second_order_gradient;

            // Print the current second-order gradient and its index
            printf("Second-order gradient at rg index %zu after adding value %.6f: %.6f\n", i, current_value, second_order_gradient);
        } else {
            second_order_gradients[i] = NAN; // Not enough points yet to calculate the gradient

            // Print a message indicating insufficient data and its index
            printf("Not enough data points to calculate the gradient at rg index %zu after adding value %.6f.\n", i, current_value);
        }
    }
}

#define MINIMUM_REQUIRED_TREND_COUNT 2

/**
 * @brief Verifies the detected peak by checking for a consistent trend of increases on the left side
 * and decreases on the right side, with additional checks and adjustments for truncated data using cubic RLS.
 *
 * This function verifies whether a detected peak is a true peak by analyzing the second-order gradients
 * calculated using cubic regression. Specifically, it checks for a consistent increasing trend on the left
 * side and a decreasing trend on the right side of the peak. If either side of the peak has insufficient
 * data points (due to reaching the boundaries of the analysis window), the function attempts to shift the
 * analysis window to obtain more data and re-verifies the trend consistency.
 *
 * @param values Array of double values representing the data.
 * @param length Length of the data array.
 * @param second_order_gradients Array containing the precomputed second-order gradients.
 * @param peak_index The index of the detected peak within the second_order_gradients array.
 * @param start_index The starting index in the original data array corresponding to the first element of second_order_gradients.
 * @param forgetting_factor The forgetting factor used in the RLS algorithm.
 * @return bool True if the peak is verified based on the trend analysis, false otherwise.
 */
bool verify_cubic_peak(const double *values, size_t length, const double *second_order_gradients, size_t peak_index, size_t start_index, double forgetting_factor) {
    size_t left_trend_count = 0;
    size_t right_trend_count = 0;
    bool left_truncated = false;
    bool right_truncated = false;

    // Count increasing trends on the left side of the peak
    for (size_t i = peak_index; i > 0; --i) {
        if (second_order_gradients[i - 1] > 0) {
            left_trend_count++;
            if (left_trend_count >= MINIMUM_REQUIRED_TREND_COUNT) break;
        } else {
            if (i == 1) left_truncated = true;  // If we hit the start of the window
            break;
        }
    }

    // Count decreasing trends on the right side of the peak
    for (size_t i = peak_index; i < CUBIC_RLS_WINDOW - 1; ++i) {
        if (second_order_gradients[i + 1] < 0) {
            right_trend_count++;
            if (right_trend_count >= MINIMUM_REQUIRED_TREND_COUNT) break;
        } else {
            if (i == CUBIC_RLS_WINDOW - 2) right_truncated = true;  // If we hit the end of the window
            break;
        }
    }

    printf("Left trends count: %zu, Right trends count: %zu\n", left_trend_count, right_trend_count);

    // If verification was truncated on the left side, attempt to re-verify by shifting the window left
    if (left_truncated && left_trend_count < MINIMUM_REQUIRED_TREND_COUNT && start_index >= 10) {
        printf("Left truncation detected, shifting window 10 points to the left for re-verification...\n");
        size_t new_start_index = start_index - 10;
        double shifted_gradients[CUBIC_RLS_WINDOW];
        compute_cubic_second_order_gradients(shifted_gradients, values, length, new_start_index, forgetting_factor);

        // Recount increasing trends on the left side of the peak
        for (size_t i = 10; i > 0; --i) {
            if (shifted_gradients[i - 1] > 0) {
                left_trend_count++;
                if (left_trend_count >= MINIMUM_REQUIRED_TREND_COUNT) break;
            } else {
                break;
            }
        }
    }

    // If verification was truncated on the right side, attempt to re-verify by shifting the window right
    if (right_truncated && right_trend_count < MINIMUM_REQUIRED_TREND_COUNT && start_index + CUBIC_RLS_WINDOW <= length - 10) {
        printf("Right truncation detected, shifting window 10 points to the right for re-verification...\n");
        size_t new_start_index = start_index + 10;
        double shifted_gradients[CUBIC_RLS_WINDOW];
        compute_cubic_second_order_gradients(shifted_gradients, values, length, new_start_index, forgetting_factor);

        // Recount decreasing trends on the right side of the peak
        for (size_t i = 10; i < 20; ++i) {
            if (shifted_gradients[i] < 0) {
                right_trend_count++;
                if (right_trend_count >= MINIMUM_REQUIRED_TREND_COUNT) break;
            } else {
                break;
            }
        }
    }

    // Verify if the peak meets the criteria
    return left_trend_count >= MINIMUM_REQUIRED_TREND_COUNT && right_trend_count >= MINIMUM_REQUIRED_TREND_COUNT;
}

/**
 * @brief Finds and verifies a peak in the data using cubic RLS.
 *
 * This function first detects a peak using the `find_cubic_peak` function, and then verifies the peak
 * using the `verify_cubic_peak` function. The peak is considered valid only if it meets the trend
 * criteria on both sides of the peak.
 *
 * @param values Array of double values representing the data.
 * @param length Length of the data array.
 * @param start_index The start index in the data array from which to begin the gradient calculation.
 * @param forgetting_factor The forgetting factor used in the RLS algorithm.
 * @return CubicPeakAnalysisResult Structure containing the peak detection status and the peak index if found and verified.
 */
CubicPeakAnalysisResult find_and_verify_cubic_peak(const double *values, size_t length, size_t start_index, double forgetting_factor) {
    CubicPeakAnalysisResult result = { .peak_found = false, .peak_index = 0 };
    
    // Declare the second_order_gradients array
    double second_order_gradients[CUBIC_RLS_WINDOW];
    
    // Compute the second-order gradients
    compute_cubic_second_order_gradients(second_order_gradients, values, length, start_index, forgetting_factor);

    // Find and verify the peak based on the computed gradients
    for (size_t i = 1; i < CUBIC_RLS_WINDOW; ++i) {
        if (second_order_gradients[i - 1] > 0 && second_order_gradients[i] < 0) {
            // Temporarily set the peak index
            result.peak_index = start_index + i;

            // Verify the detected peak
            if (verify_cubic_peak(values, length, second_order_gradients, i, start_index, forgetting_factor)) {
                result.peak_found = true;
                printf("Verified peak found at index %zu\n", result.peak_index);
                break; // Exit the loop once a verified peak is found
            } else {
                printf("Peak at index %zu did not pass verification. Continuing search...\n", result.peak_index);
                result.peak_found = false; // Reset peak_found as the peak failed verification
            }
        }
    }

    if (!result.peak_found) {
        printf("No verified peak found in the specified window.\n");
    }

    return result;
}

