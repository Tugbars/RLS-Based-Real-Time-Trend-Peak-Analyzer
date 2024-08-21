#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include "trend_analysis.h"
#include "windowed_running_gradient.h"
#include "running_gradient_parameters.h"
#include "running_quadratic_gradient.h"

#define EXTEND_VALUE 15  // Should depend on the width of the peak where concavity occurs.

/**
 * @brief Structure to store running gradient data for a quadratic model.
 *
 * This structure is used to maintain the data for a quadratic model, which includes the coefficients for the quadratic term,
 * linear term, and intercept.
 */
typedef struct {
	size_t num_points; /**< Number of data points currently stored */
	size_t max_points; /**< Maximum number of data points that can be stored */
	double x[RLS_WINDOW]; /**< Array of x values */
	double y[RLS_WINDOW]; /**< Array of y values */
	double coefficients[3]; /**< Coefficients of the quadratic regression model (a2, a1, a0) */
	double residual_sum_squares; /**< Residual sum of squares for the regression */
	double inverse_cov_matrix[3][3]; /**< Inverse covariance matrix for RLS update */
	double forgetting_factor; /**< Forgetting factor for the RLS algorithm */
} RunningQuadraticGradient;

/**
 * @brief Initializes the RunningQuadraticGradient structure.
 *
 * This function sets the initial values for the coefficients, residual sum of squares,
 * inverse covariance matrix, and allocates memory for x and y arrays for the quadratic model.
 *
 * @param rg Pointer to the RunningQuadraticGradient structure to initialize.
 * @param forgetting_factor Forgetting factor for the RLS algorithm.
 */
void init_running_quadratic_gradient(RunningQuadraticGradient *rg, double forgetting_factor) {
	rg->num_points = 0;
	rg->max_points = RLS_WINDOW;
	
	// Initialize coefficients to zero
	rg->coefficients[0] = 0.0;
	rg->coefficients[1] = 0.0;
	rg->coefficients[2] = 0.0;
	
	// Initialize residual sum of squares to zero
	rg->residual_sum_squares = 0.0;
	
	// Set the forgetting factor
	rg->forgetting_factor = forgetting_factor;
	
	// Initialize the inverse covariance matrix
	memset(rg->inverse_cov_matrix, 0, sizeof(rg->inverse_cov_matrix));
	rg->inverse_cov_matrix[0][0] = 1e9; // High uncertainty in x^2 coefficient
	rg->inverse_cov_matrix[1][1] = 1e9; // High uncertainty in x coefficient
	rg->inverse_cov_matrix[2][2] = 1e9; // High uncertainty in constant term
}

/**
 * @brief Adds a new data point to the RunningQuadraticGradient structure and updates the quadratic regression model using Recursive Least Squares (RLS).
 *
 * This function adds a new data point to the RunningQuadraticGradient structure, which maintains a quadratic model of the data points
 * using the Recursive Least Squares (RLS) algorithm to update the regression coefficients.
 *
 * @param rg Pointer to the RunningQuadraticGradient structure.
 * @param y The new data point to be added.
 */
void add_quadratic_data_point(RunningQuadraticGradient *const rg, const double y) {
	// Check if the current window is not full
	if (rg->num_points < rg->max_points) {
		// Add new data point to the window
		rg->y[rg->num_points] = y;
		rg->x[rg->num_points] = (double)rg->num_points;
		rg->num_points++;
	} else {
		// Shift the window to the left and add new data point
		// This shift operation removes the oldest data point and adds the new one at the end.
		for (size_t i = 0; i < rg->max_points - 1; ++i) {
			rg->y[i] = rg->y[i + 1];
			rg->x[i] = rg->x[i + 1];
		}
		rg->y[rg->max_points - 1] = y;  // Add new data point
		rg->x[rg->max_points - 1] = (double)(rg->num_points);
		rg->num_points++;
	}

	// Define the input vector for the new data point
	// The input vector for quadratic fitting is x^2, x, and 1 (for the intercept).
	double x0 = pow((double)(rg->num_points - 1), 2); // Quadratic term
	double x1 = (double)(rg->num_points - 1); // Linear term
	double x2 = 1.0; // Constant term

	// Calculate the denominator term in the RLS update equations
	// This term helps in updating the inverse covariance matrix efficiently using the Sherman-Morrison formula.
	double temp = rg->forgetting_factor +
	              (x0 * rg->inverse_cov_matrix[0][0] + x1 * rg->inverse_cov_matrix[1][0] + x2 * rg->inverse_cov_matrix[2][0]) * x0 +
	              (x0 * rg->inverse_cov_matrix[0][1] + x1 * rg->inverse_cov_matrix[1][1] + x2 * rg->inverse_cov_matrix[2][1]) * x1 +
	              (x0 * rg->inverse_cov_matrix[0][2] + x1 * rg->inverse_cov_matrix[1][2] + x2 * rg->inverse_cov_matrix[2][2]) * x2;

	// Check for numerical stability
	if (fabs(temp) < 1e-10) {
		fprintf(stderr, "Numerical stability issue: temp is too close to zero.\n");
		exit(EXIT_FAILURE);
	}

	// Calculate the temporary vector used for updating the inverse covariance matrix
	double tmp0 = rg->inverse_cov_matrix[0][0] * x0 + rg->inverse_cov_matrix[0][1] * x1 + rg->inverse_cov_matrix[0][2] * x2;
	double tmp1 = rg->inverse_cov_matrix[1][0] * x0 + rg->inverse_cov_matrix[1][1] * x1 + rg->inverse_cov_matrix[1][2] * x2;
	double tmp2 = rg->inverse_cov_matrix[2][0] * x0 + rg->inverse_cov_matrix[2][1] * x1 + rg->inverse_cov_matrix[2][2] * x2;

	// Update the inverse covariance matrix using the Sherman-Morrison formula
	// The following lines update the inverse covariance matrix to incorporate the new data point.
	rg->inverse_cov_matrix[0][0] = (rg->inverse_cov_matrix[0][0] - (tmp0 * tmp0) / temp) / rg->forgetting_factor;
	rg->inverse_cov_matrix[0][1] = (rg->inverse_cov_matrix[0][1] - (tmp0 * tmp1) / temp) / rg->forgetting_factor;
	rg->inverse_cov_matrix[0][2] = (rg->inverse_cov_matrix[0][2] - (tmp0 * tmp2) / temp) / rg->forgetting_factor;
	rg->inverse_cov_matrix[1][0] = rg->inverse_cov_matrix[0][1]; // Enforcing symmetry
	rg->inverse_cov_matrix[1][1] = (rg->inverse_cov_matrix[1][1] - (tmp1 * tmp1) / temp) / rg->forgetting_factor;
	rg->inverse_cov_matrix[1][2] = (rg->inverse_cov_matrix[1][2] - (tmp1 * tmp2) / temp) / rg->forgetting_factor;
	rg->inverse_cov_matrix[2][0] = rg->inverse_cov_matrix[0][2]; // Enforcing symmetry
	rg->inverse_cov_matrix[2][1] = rg->inverse_cov_matrix[1][2]; // Enforcing symmetry
	rg->inverse_cov_matrix[2][2] = (rg->inverse_cov_matrix[2][2] - (tmp2 * tmp2) / temp) / rg->forgetting_factor;

	// Calculate the prediction error for the new data point
	// Prediction error is the difference between the actual data point and the value predicted by the quadratic model.
	double prediction_error = y - (x0 * rg->coefficients[0] + x1 * rg->coefficients[1] + x2 * rg->coefficients[2]);

	// Update the regression coefficients using the RLS update equations
	// The coefficients for the quadratic model are updated using the inverse covariance matrix and the prediction error.
	rg->coefficients[0] += (rg->inverse_cov_matrix[0][0] * x0 + rg->inverse_cov_matrix[0][1] * x1 + rg->inverse_cov_matrix[0][2] * x2) * prediction_error;
	rg->coefficients[1] += (rg->inverse_cov_matrix[1][0] * x0 + rg->inverse_cov_matrix[1][1] * x1 + rg->inverse_cov_matrix[1][2] * x2) * prediction_error;
	rg->coefficients[2] += (rg->inverse_cov_matrix[2][0] * x0 + rg->inverse_cov_matrix[2][1] * x1 + rg->inverse_cov_matrix[2][2] * x2) * prediction_error;

	// Recalculate the residual sum of squares
	// The residual sum of squares (RSS) is a measure of the discrepancy between the data and the estimation model.
	rg->residual_sum_squares = 0.0;
	for (size_t i = 0; i < rg->num_points; ++i) {
		double error = rg->y[i] - (rg->coefficients[0] * pow(rg->x[i], 2) + rg->coefficients[1] * rg->x[i] + rg->coefficients[2]);
		rg->residual_sum_squares += error * error;
	}
}

/**
 * @brief Calculates the slope (first derivative) of the quadratic function at a specific point.
 *
 * This function computes the slope of the quadratic equation represented by the
 * coefficients stored in the RunningQuadraticGradient structure at a given value of x.
 * The slope is derived from the first derivative of the quadratic function:
 *
 *     slope = 2 * a2 * x + a1
 *
 * where:
 *     a2 is the coefficient of the x^2 term,
 *     a1 is the coefficient of the x term,
 *     x is the point at which the slope is calculated.
 *
 * @param rg Pointer to the RunningQuadraticGradient structure containing the coefficients.
 * @param x The point at which to calculate the slope.
 * @return The slope (first derivative) of the quadratic function at the given point x.
 */
double calculate_slope_at_point(const RunningQuadraticGradient *rg, double x) {
	double a2 = rg->coefficients[0];  // Coefficient for x^2
	double a1 = rg->coefficients[1];  // Coefficient for x
	// Calculate the slope at the given point x
	double slope = 2.0 * a2 * x + a1;
	return slope;
}

/**
 * @brief Calculates the second-order gradient (curvature) of the quadratic function.
 *
 * This function computes the second-order gradient of the quadratic equation represented by the
 * coefficients stored in the RunningQuadraticGradient structure. The second-order gradient is
 * the second derivative of the quadratic function, which for a quadratic function is constant:
 *
 *     second_order_gradient = 2 * a2
 *
 * where:
 *     a2 is the coefficient of the x^2 term.
 *
 * @param rg Pointer to the RunningQuadraticGradient structure containing the coefficients.
 * @return The second-order gradient (curvature) of the quadratic function.
 */
double calculate_second_order_gradient(const RunningQuadraticGradient *rg) {
	double a2 = rg->coefficients[0];  // Coefficient for x^2
	double second_order_gradient = 2.0 * a2;
	return second_order_gradient;
}


/**
 * @brief Computes the second-order gradients for the next 30 values in the given array.
 *
 * This function takes an array of doubles and a starting index, and then adds each of the next 30 values
 * to the RunningQuadraticGradient structure. After each addition, it computes the second-order gradient
 * (curvature) of the quadratic model. The function returns an array containing these second-order gradients.
 *
 * @param values Array of double values.
 * @param length Length of the array.
 * @param start_index Starting index from which to begin the gradient calculation.
 * @param forgetting_factor The forgetting factor for the RLS algorithm.
 * @return An array of double containing the second-order gradients.
 */
double* compute_second_order_gradients(const double *values, size_t length, size_t start_index, double forgetting_factor) {
	// Ensure that the start index and the length allow for 30 calculations
	if (start_index + RLS_WINDOW > length) {
		printf( "Insufficient data to compute 30 second-order gradients from the starting index.\n");
		return NULL;
	}

	// Allocate memory for the result array
	double *second_order_gradients = (double *)malloc(RLS_WINDOW * sizeof(double));
	if (second_order_gradients == NULL) {
		fprintf(stderr, "Memory allocation failed.\n");
		return NULL;
	}

	// Initialize the running quadratic gradient structure
	RunningQuadraticGradient rg;
	init_running_quadratic_gradient(&rg, forgetting_factor);

	// Loop through the next 30 values starting from start_index
	for (size_t i = 0; i < RLS_WINDOW; ++i) {
		double current_value = values[start_index + i];

		// Add the current value to the running quadratic gradient
		add_quadratic_data_point(&rg, current_value);

		// Print the current value added
		printf("Added value: %.6f\n", current_value);

		// Calculate the second-order gradient only if we have at least 3 data points
		if (rg.num_points >= 3) {
			double second_order_gradient = calculate_second_order_gradient(&rg);
			second_order_gradients[i] = second_order_gradient;

			// Print the current second-order gradient
			printf("Second-order gradient after adding value %.6f: %.6f\n", current_value, second_order_gradient);
		} else {
			second_order_gradients[i] = NAN; // Not enough points yet to calculate the gradient

			// Print a message indicating insufficient data
			printf("Not enough data points to calculate the gradient after adding value %.6f.\n", current_value);
		}
	}

	return second_order_gradients;
}

double* compute_first_order_gradients(const double *values, size_t length, size_t start_index, double forgetting_factor) {
	// Ensure that the start index and the length allow for 30 calculations
	if (start_index + RLS_WINDOW > length) {
		fprintf(stderr, "Insufficient data to compute 30 first-order gradients from the starting index.\n");
		return NULL;
	}

	// Allocate memory for the result array
	double *first_order_gradients = (double *)malloc(RLS_WINDOW * sizeof(double));
	if (first_order_gradients == NULL) {
		fprintf(stderr, "Memory allocation failed.\n");
		return NULL;
	}

	// Initialize the running quadratic gradient structure
	RunningQuadraticGradient rg;
	init_running_quadratic_gradient(&rg, forgetting_factor);

	// Loop through the next 30 values starting from start_index
	for (size_t i = 0; i < RLS_WINDOW; ++i) {
		double current_value = values[start_index + i];

		// Add the current value to the running quadratic gradient
		add_quadratic_data_point(&rg, current_value);

		// Print the current value added
		printf("Added value: %.6f\n", current_value);

		// Calculate the first-order gradient only if we have at least 3 data points
		if (rg.num_points >= 2) {
			// Calculate the first-order gradient (slope) at the current x value
			double x_current = (double)(rg.num_points - 1); // Current x value
			double first_order_gradient = calculate_slope_at_point(&rg, x_current);
			first_order_gradients[i] = first_order_gradient;

			// Print the current first-order gradient
			printf("First-order gradient after adding value %.6f: %.6f\n", current_value, first_order_gradient);
		} else {
			first_order_gradients[i] = NAN; // Not enough points yet to calculate the gradient

			// Print a message indicating insufficient data
			printf("Not enough data points to calculate the gradient after adding value %.6f.\n", current_value);
		}
	}

	return first_order_gradients;
}

/**
 * @brief Structure to store the result of gradient analysis.
 *
 * This structure contains the position relative to the peak and the sum of first-order gradients.
 */
typedef struct {
    IndexRelativeToPeak position;
    double first_order_sum;
} GradientAnalysisResult;

/**
 * @brief Computes and verifies the gradients to determine the position relative to the peak.
 *
 * This function computes second-order and first-order gradients starting from the `start_index`. 
 * It then determines the position relative to the peak based on the second-order gradients and verifies 
 * this position by summing the first-order gradients.
 *
 * @param values Array of double values representing the data.
 * @param length Length of the data array.
 * @param start_index The starting index for gradient calculations.
 * @param target_index The index whose position relative to the peak is being determined.
 * @param forgetting_factor The forgetting factor used in the RLS algorithm.
 * @return A `GradientAnalysisResult` structure containing the position relative to the peak and the sum of first-order gradients.
 */
GradientAnalysisResult compute_and_verify_gradients(const double *values, size_t length, size_t start_index, size_t target_index, double forgetting_factor) {
    GradientAnalysisResult result;
    result.position = INDEX_POSITION_UNDETERMINED;
    result.first_order_sum = 0.0;

    // Compute second-order gradients starting from start_index
    double *second_order_gradients = compute_second_order_gradients(values, length, start_index, forgetting_factor);
    if (second_order_gradients == NULL) {
        printf("Failed to compute second-order gradients.\n");
        return result;
    }

    // Determine the specific index within the gradient array corresponding to target_index
    size_t gradient_index = target_index - start_index;
    double current_second_order_gradient = second_order_gradients[gradient_index];

    // Determine peak position based on second-order gradients
    if (fabs(current_second_order_gradient) < GRADIENT_THRESHOLD) {
        result.position = INDEX_AT_PEAK;
        printf("********Position determined: AT PEAK (Gradient: %.6f)\n", current_second_order_gradient);
    } else if (current_second_order_gradient > 0) {
        result.position = INDEX_LEFT_TO_PEAK;
        printf("*******Position determined: LEFT OF PEAK (Gradient: %.6f)\n", current_second_order_gradient);
    } else {
        result.position = INDEX_RIGHT_TO_PEAK;
        printf("*******Position determined: RIGHT OF PEAK (Gradient: %.6f)\n", current_second_order_gradient);
    }

    free(second_order_gradients);

    // Compute first-order gradients starting from start_index
    double *first_order_gradients = compute_first_order_gradients(values, length, start_index, forgetting_factor);
    if (first_order_gradients == NULL) {
        printf("Failed to compute first-order gradients.\n");
        return result;
    }

    if (result.position == INDEX_LEFT_TO_PEAK) {
    for (size_t i = 0; i <= gradient_index; ++i) {
        if (isnan(first_order_gradients[i])) {
            printf("NaN found in first-order gradient at index %zu, skipping...\n", i);
            continue;
        }
        result.first_order_sum += first_order_gradients[i];
        printf("Adding first-order gradient %.6f to sum. Current sum: %.6f\n", first_order_gradients[i], result.first_order_sum);
    }
    printf("Final sum of first-order gradients on the left side: %.6f\n", result.first_order_sum);
} else if (result.position == INDEX_RIGHT_TO_PEAK) {
    for (size_t i = gradient_index; i < RLS_WINDOW; ++i) {
        if (isnan(first_order_gradients[i])) {
            printf("NaN found in first-order gradient at index %zu, skipping...\n", i);
            continue;
        }
        result.first_order_sum += first_order_gradients[i];
        printf("Adding first-order gradient %.6f to sum. Current sum: %.6f\n", first_order_gradients[i], result.first_order_sum);
    }
    printf("Final sum of first-order gradients on the right side: %.6f\n", result.first_order_sum);
}

    free(first_order_gradients);

    return result;
}

/**
 * @brief Determines the position relative to the peak for a given index.
 *
 * This function determines the position relative to the peak for a specified `target_index` by 
 * computing and verifying the second-order and first-order gradients. The result is based on both the 
 * second-order gradients and the verification process using first-order gradients.
 *
 * @param values Array of double values representing the data.
 * @param length Length of the data array.
 * @param start_index The starting index for gradient calculations.
 * @param target_index The index whose position relative to the peak is being determined.
 * @param forgetting_factor The forgetting factor used in the RLS algorithm.
 * @return An `IndexRelativeToPeak` enum indicating the position relative to the peak (LEFT, RIGHT, AT, or UNDETERMINED).
 */
IndexRelativeToPeak determine_position_relative_to_peak(const double *values, size_t length, size_t start_index, size_t target_index, double forgetting_factor) {
    // Ensure the target index is within a valid range considering the RLS window size
    if (target_index >= length || target_index < start_index || start_index + RLS_WINDOW > length) {
        printf("Indices out of valid range for RLS_WINDOW size.\n");
        return INDEX_POSITION_UNDETERMINED;
    }

    printf("Analyzing gradients and verifying position...\n");

    // Compute and verify gradients
    GradientAnalysisResult result = compute_and_verify_gradients(values, length, start_index, target_index, forgetting_factor);

    // Final determination based on both first-order and second-order gradients
    if (result.position == INDEX_AT_PEAK) {
        printf("Final determination: AT PEAK\n");
        return INDEX_AT_PEAK;
    } else if (result.position == INDEX_LEFT_TO_PEAK && result.first_order_sum >= 0) {
        printf("Final determination: LEFT OF PEAK\n");
        return INDEX_LEFT_TO_PEAK;
    } else if (result.position == INDEX_RIGHT_TO_PEAK && result.first_order_sum <= 0) {
        printf("Final determination: RIGHT OF PEAK\n");
        return INDEX_RIGHT_TO_PEAK;
    } else {
        printf("Final determination: POSITION UNDETERMINED\n");
        return INDEX_POSITION_UNDETERMINED;
    }
}

#define MINIMUM_REQUIRED_TREND_COUNT 5

/**
 * @brief Verifies the detected peak by checking for a consistent trend of increases on the left side
 * and decreases on the right side, with additional checks and adjustments for truncated data.
 *
 * This function verifies whether a detected peak is a true peak by analyzing the second-order gradients.
 * Specifically, it checks for a consistent increasing trend on the left side and a decreasing trend on the right side of the peak.
 * If either side of the peak has insufficient data points (due to reaching the boundaries of the analysis window), 
 * the function attempts to shift the analysis window to obtain more data and re-verifies the trend consistency.
 *
 * @param values Array of double values representing the data.
 * @param length Length of the data array.
 * @param second_order_gradients Array of second-order gradients corresponding to the data.
 * @param peak_index The index of the detected peak within the second_order_gradients array.
 * @param start_index The starting index in the original data array corresponding to the first element of second_order_gradients.
 * @param forgetting_factor The forgetting factor used in the RLS algorithm.
 * @return bool True if the peak is verified based on the trend analysis, false otherwise.
 *
 * The function follows these steps:
 * 1. **Initial Trend Verification**:
 *    - Counts the number of data points on the left of the peak where the second-order gradient is positive, indicating an increasing trend.
 *    - Counts the number of data points on the right of the peak where the second-order gradient is negative, indicating a decreasing trend.
 *    - If the required number of consistent trends is found on both sides (as defined by `MINIMUM_REQUIRED_TREND_COUNT`), the peak is considered verified.
 * 
 * 2. **Handling Truncated Data**:
 *    - If the left side or right side analysis is truncated (i.e., reaches the boundary of the analysis window), the function attempts to re-verify by shifting the analysis window.
 * 
 * 3. **Re-verification by Expanding the Window**:
 *    - **Left Side Truncation**: If the left side was truncated and the number of increasing trends is insufficient, the function shifts the window 10 points to the left (if possible) and re-analyzes the left side.
 *    - **Right Side Truncation**: If the right side was truncated and the number of decreasing trends is insufficient, the function shifts the window 10 points to the right (if possible) and re-analyzes the right side.
 * 
 * 4. **Final Verification**:
 *    - After the re-verification, if the minimum required trends are found on both sides of the peak, the peak is confirmed as verified.
 *    - If either side fails the trend consistency check even after re-verification, the peak is not considered verified.
 */
bool verify_peak(const double *values, size_t length, const double *second_order_gradients, size_t peak_index, size_t start_index, double forgetting_factor) {
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
    for (size_t i = peak_index; i < RLS_WINDOW - 1; ++i) {
        if (second_order_gradients[i + 1] < 0) {
            right_trend_count++;
            if (right_trend_count >= MINIMUM_REQUIRED_TREND_COUNT) break;
        } else {
            if (i == RLS_WINDOW - 2) right_truncated = true;  // If we hit the end of the window
            break;
        }
    }

    printf("Left trends count: %zu, Right trends count: %zu\n", left_trend_count, right_trend_count);

    // If verification was truncated on the left side, attempt to re-verify by shifting the window left
    if (left_truncated && left_trend_count < MINIMUM_REQUIRED_TREND_COUNT && start_index >= 10) {
        printf("Left truncation detected, shifting window 10 points to the left for re-verification...\n");
        size_t new_start_index = start_index - 10;
        size_t new_end_index = new_start_index + 20;
        double *shifted_gradients = compute_second_order_gradients(values, length, new_start_index, forgetting_factor);

        if (shifted_gradients) {
            // Recount increasing trends on the left side of the peak
            for (size_t i = 10; i > 0; --i) {
                if (shifted_gradients[i - 1] > 0) {
                    left_trend_count++;
                    if (left_trend_count >= MINIMUM_REQUIRED_TREND_COUNT) break;
                } else {
                    break;
                }
            }
            free(shifted_gradients);
        }
    }

    // If verification was truncated on the right side, attempt to re-verify by shifting the window left
    if (right_truncated && right_trend_count < MINIMUM_REQUIRED_TREND_COUNT && start_index + RLS_WINDOW <= length - 10) {
        printf("Right truncation detected, shifting window 10 points to the left for re-verification...\n");
        size_t new_start_index = start_index + 10;
        size_t new_end_index = new_start_index + 20;
        double *shifted_gradients = compute_second_order_gradients(values, length, new_start_index, forgetting_factor);

        if (shifted_gradients) {
            // Recount decreasing trends on the right side of the peak
            for (size_t i = 10; i < 20; ++i) {
                if (shifted_gradients[i] < 0) {
                    right_trend_count++;
                    if (right_trend_count >= MINIMUM_REQUIRED_TREND_COUNT) break;
                } else {
                    break;
                }
            }
            free(shifted_gradients);
        }
    }

    // Verify if the peak meets the criteria
    return left_trend_count >= MINIMUM_REQUIRED_TREND_COUNT && right_trend_count >= MINIMUM_REQUIRED_TREND_COUNT;
}

/**
 * @brief Detects and verifies a peak in the data using the second-order gradient sign change.
 *
 * This function computes second-order gradients for a given interval and attempts to detect a peak. 
 * If a peak is found, it verifies whether the peak has the required increasing trend on the left 
 * and decreasing trend on the right. If the verification fails, it continues to search for the 
 * next potential peak.
 *
 * @param values Array of double values representing the data.
 * @param length Length of the data array.
 * @param start_index The start index for computing gradients.
 * @param forgetting_factor The forgetting factor for the RLS algorithm.
 * @return PeakAnalysisResult containing the peak detection status and index.
 */
PeakAnalysisResult find_and_verify_peak(const double *values, size_t length, size_t start_index, double forgetting_factor) {
    PeakAnalysisResult result = { .peak_found = false, .peak_index = 0 };

    // Compute second-order gradients starting from start_index
    double *second_order_gradients = compute_second_order_gradients(values, length, start_index, forgetting_factor);
    if (second_order_gradients == NULL) {
        printf("Failed to compute second-order gradients.\n");
        return result;
    }

    // Loop through the gradients to find the peak (sign change from positive to negative)
    for (size_t i = 1; i < RLS_WINDOW; ++i) {
        if (second_order_gradients[i - 1] > 0 && second_order_gradients[i] < 0) {
            result.peak_index = start_index + i;

            // Verify the detected peak
            if (verify_peak(values, length, second_order_gradients, i, start_index, forgetting_factor)) {
                result.peak_found = true;
                printf("Verified peak found at index %zu\n", result.peak_index);
                break;
            } else {
                printf("Peak at index %zu did not pass verification. Continuing search...\n", result.peak_index);
            }
        }
    }

    free(second_order_gradients);

    if (!result.peak_found) {
        printf("No verified peak found in the specified window.\n");
    }

    return result;
}

/**
 * @brief Performs initial concavity analysis by summing the second-order gradients for every 10 points in the RLS window.
 *
 * This function computes the sum of second-order gradients over every 10 points within the RLS window,
 * starting from the provided `start_index`. The results are returned in a `ConcavityAnalysisResult` structure.
 * The function also allows for re-initialization of the `RunningQuadraticGradient` structure after every 10 points,
 * depending on the value of the `reinitialize_after_each_segment` parameter.
 *
 * @param values Array of double values representing the data.
 * @param length Length of the data array.
 * @param start_index The starting index for gradient calculations.
 * @param forgetting_factor The forgetting factor used in the RLS algorithm.
 * @param reinitialize_after_each_segment Boolean flag to control whether to re-initialize after every 10 data points.
 * @return A `ConcavityAnalysisResult` structure containing the sums of second-order gradients for each 10-point segment.
 */
ConcavityAnalysisResult initial_concavity_analysis(const double *values, size_t length, size_t start_index, double forgetting_factor, bool reinitialize_after_each_segment) {
    ConcavityAnalysisResult result = { .sums = {0.0} };

    // Ensure that the start index and the length allow for the analysis
    if (start_index + RLS_WINDOW > length) {
        printf("Insufficient data for initial concavity analysis.\n");
        return result;
    }

    RunningQuadraticGradient rg;
    init_running_quadratic_gradient(&rg, forgetting_factor);

    // Loop through the RLS window in segments of 10 points
    for (size_t i = 0; i < RLS_WINDOW; i += 10) {
        if (reinitialize_after_each_segment) {
            // Re-initialize the running quadratic gradient structure
            init_running_quadratic_gradient(&rg, forgetting_factor);
        }

        double sum = 0.0;
        size_t valid_point_count = 0;

        // Process each segment of 10 points
        for (size_t j = i; j < i + 10 && j < RLS_WINDOW; ++j) {
            double current_value = values[start_index + j];
            add_quadratic_data_point(&rg, current_value);

            if (rg.num_points >= 3 && !isnan(rg.coefficients[0])) { // Ensure at least 3 points before calculating
                double second_order_gradient = calculate_second_order_gradient(&rg);
                sum += second_order_gradient;
                valid_point_count++;
            } else {
                printf("Skipping index %zu due to insufficient data points or NaN.\n", j);
            }
        }

        // Only record the sum if it was computed with at least one valid point
        if (valid_point_count > 0) {
            result.sums[i / 10] = sum;
            printf("Sum of second-order gradients from index %zu to %zu: %.6f\n", start_index + i, start_index + i + 9, sum);
        } else {
            result.sums[i / 10] = NAN;
            printf("Insufficient data to calculate sum for segment %zu.\n", i / 10 + 1);
        }
    }

    return result;
}


/**
 * @brief Analyzes the behavior of the three segments and returns the corresponding concavity pattern.
 *
 * This function checks the sum of second-order gradients for each of the three segments and determines 
 * the overall pattern of concavity across these segments. It returns an enum value representing the 
 * pattern identified. It also identifies potential and true peaks.
 *
 * @param concavity_result The `ConcavityAnalysisResult` containing the sums of the second-order gradients.
 * @param isPotentialPeak Pointer to a boolean that will be set to true if a potential peak is detected.
 * @param isTruePeak Pointer to a boolean that will be set to true if a true peak is detected.
 * @return A `ConcavityPattern` enum value representing the pattern of concavity across the segments.
 */
ConcavityPattern analyze_concavity_segments(const ConcavityAnalysisResult *concavity_result, bool *isPotentialPeak, bool *isTruePeak) {
    if (isnan(concavity_result->sums[0]) || isnan(concavity_result->sums[1]) || isnan(concavity_result->sums[2])) {
        printf("Insufficient data for determining concavity pattern.\n");
        *isPotentialPeak = false;
        *isTruePeak = false;
        return UNDETERMINED_PATTERN;
    }

    *isPotentialPeak = false;
    *isTruePeak = false;

    // Map increase (> 0.3) to 1 and decrease (<= 0.3) to 0
    int first_segment = concavity_result->sums[0] > 0.3 ? 1 : 0;
    int second_segment = concavity_result->sums[1] > 0.3 ? 1 : 0;
    int third_segment = concavity_result->sums[2] > 0.3 ? 1 : 0;

    // Create a binary pattern based on the segments
    int pattern = (first_segment << 2) | (second_segment << 1) | third_segment;

    // Switch on the binary pattern to determine the concavity pattern
    switch (pattern) {
        case 0b110: // 6: INCREASE, INCREASE, DECREASE
            *isPotentialPeak = true;
            return INCREASE_INCREASE_DECREASE;
        case 0b101: // 5: INCREASE, DECREASE, INCREASE
            *isPotentialPeak = true;
            return INCREASE_DECREASE_INCREASE;
        case 0b100: // 4: INCREASE, DECREASE, DECREASE
            *isPotentialPeak = true;
            *isTruePeak = true;
            return INCREASE_DECREASE_DECREASE;
        case 0b011: // 3: DECREASE, INCREASE, INCREASE
            return DECREASE_INCREASE_INCREASE;
        case 0b010: // 2: DECREASE, INCREASE, DECREASE
            *isPotentialPeak = true;
            return DECREASE_INCREASE_DECREASE;
        case 0b001: // 1: DECREASE, DECREASE, INCREASE
            return DECREASE_DECREASE_INCREASE;
        case 0b000: // 0: DECREASE, DECREASE, DECREASE
            *isTruePeak = true;
            return DECREASE_DECREASE_DECREASE;
        default:
            return UNDETERMINED_PATTERN;
    }
}

/**
 * @brief Iteratively moves the RLS window by 10 items to find a true peak.
 *
 * This function starts with an initial RLS window and moves it by 10 items until a true peak is found.
 * At each step, it performs an initial concavity analysis and checks for the presence of a true peak.
 * The function returns the starting index of the search and the end index where the true peak was found.
 *
 * @param values Array of double values representing the data.
 * @param length Length of the data array.
 * @param start_index The initial starting index for gradient calculations.
 * @param forgetting_factor The forgetting factor used in the RLS algorithm.
 * @param reinitialize_after_each_segment Boolean flag to control whether to re-initialize after every 10 data points.
 * @param start_idx Pointer to store the initial starting index of the search.
 * @param end_idx Pointer to store the ending index where the true peak was found.
 * @return bool True if a true peak was found, false otherwise.
 */
bool find_true_peak_with_concavity_analysis(const double *values, size_t length, size_t start_index, double forgetting_factor, bool reinitialize_after_each_segment, size_t *start_idx, size_t *end_idx) {
    *start_idx = start_index;
    bool isPotentialPeak = false;
    bool isTruePeak = false;

    while (start_index + RLS_WINDOW <= length) {
        // Perform initial concavity analysis
        ConcavityAnalysisResult concavity_result = initial_concavity_analysis(values, length, start_index, forgetting_factor, reinitialize_after_each_segment);

        // Analyze the concavity segments
        ConcavityPattern pattern = analyze_concavity_segments(&concavity_result, &isPotentialPeak, &isTruePeak);

        // Check if a true peak was found
        if (isTruePeak) {
            *end_idx = start_index + RLS_WINDOW - 1;
            printf("True peak found between indices %zu and %zu with pattern: %d\n", *start_idx, *end_idx, pattern);
            return true;
        }

        // Move the start index by 10 for the next window
        start_index += 10;
    }

    // No true peak was found
    *end_idx = start_index + RLS_WINDOW - 1;
    printf("No true peak found after examining all possible windows.\n");
    return false;
}

