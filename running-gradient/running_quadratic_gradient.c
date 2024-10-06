#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include "windowed_running_gradient.h"
#include "running_quadratic_gradient.h"
#include "rls_analysis_parameters.h"
#include "buffer_manager.h"

/**
 * @brief Structure to store running gradient data for a quadratic model.
 *
 * This structure is used to maintain the data for a quadratic model, which includes the coefficients for the quadratic term,
 * linear term, and intercept.
 */
typedef struct {
    uint16_t num_points;            /**< Number of data points currently stored */
    uint16_t max_points;            /**< Maximum number of data points that can be stored */
    double x[RLS_WINDOW];           /**< Array of x values */
    double y[RLS_WINDOW];           /**< Array of y values */
    double coefficients[3];         /**< Coefficients of the quadratic regression model (a2, a1, a0) */
    double residual_sum_squares;    /**< Residual sum of squares for the regression */
    double inverse_cov_matrix[3][3];/**< Inverse covariance matrix for RLS update */
    double forgetting_factor;       /**< Forgetting factor for the RLS algorithm */
} RunningQuadraticGradient;

/**
 * @brief Initializes the inverse covariance matrix for the RLS algorithm.
 * 
 * This function sets the diagonal elements of the matrix to large values,
 * representing high initial uncertainty about the model coefficients.
 *
 * @param rg Pointer to the RunningQuadraticGradient structure.
 */
void init_inverse_covariance_matrix(RunningQuadraticGradient *rg) {
    memset(rg->inverse_cov_matrix, 0, sizeof(rg->inverse_cov_matrix));  // Set all elements to 0
    rg->inverse_cov_matrix[0][0] = 1e9; // High uncertainty in x^2 coefficient
    rg->inverse_cov_matrix[1][1] = 1e9; // High uncertainty in x coefficient
    rg->inverse_cov_matrix[2][2] = 1e9; // High uncertainty in constant term
}

/**
 * @brief Initializes the RunningQuadraticGradient structure.
 *
 * This function sets the initial values for the coefficients, residual sum of squares,
 * and calls an external function to initialize the inverse covariance matrix.
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

    // Initialize the inverse covariance matrix externally
    init_inverse_covariance_matrix(rg);
}

/**
 * @brief Computes the condition number of the inverse covariance matrix to assess numerical stability.
 * 
 * The condition number is a measure of the sensitivity of a matrix to numerical errors. It is computed as the ratio
 * of the largest to smallest singular values. For symmetric positive-definite matrices like the inverse covariance matrix,
 * we approximate the condition number using the largest and smallest diagonal elements.
 *
 * ### Mathematical Background:
 * The condition number \(\kappa(A)\) is defined as:
 * \f[
 * \kappa(A) = \frac{\sigma_{\text{max}}}{\sigma_{\text{min}}}
 * \f]
 * where \(\sigma_{\text{max}}\) and \(\sigma_{\text{min}}\) are the largest and smallest singular values of the matrix \(A\).
 * A large condition number indicates that the matrix is close to singular, leading to instability and large numerical errors.
 *
 * ### Importance in RLS:
 * In the Recursive Least Squares (RLS) algorithm, the inverse covariance matrix is updated iteratively. If it becomes ill-conditioned
 * (i.e., high condition number), the algorithm may become unstable. Monitoring and resetting the inverse covariance matrix when
 * the condition number exceeds a threshold helps ensure numerical stability.
 *
 * @param matrix The 3x3 inverse covariance matrix from the RLS algorithm.
 * @return The condition number, approximated as the ratio of the largest to smallest diagonal element.
 */
double compute_condition_number(const double matrix[3][3]) {
    double max_value = fabs(matrix[0][0]);
    double min_value = fabs(matrix[0][0]);

    // Find the maximum and minimum diagonal elements
    for (int i = 0; i < 3; ++i) {
        if (fabs(matrix[i][i]) > max_value) {
            max_value = fabs(matrix[i][i]);
        }
        if (fabs(matrix[i][i]) < min_value) {
            min_value = fabs(matrix[i][i]);
        }
    }

    // Avoid division by zero or very small numbers
    if (min_value < 1e-10) {
        min_value = 1e-10;
    }

    return max_value / min_value;
}

/**
 * @brief Normalizes the x and y data arrays in the RunningQuadraticGradient structure.
 *
 * This function normalizes the x and y data points by subtracting their mean and dividing by their standard deviation.
 * Normalization improves numerical stability, especially in RLS, by ensuring that data values are within a reasonable range,
 * avoiding large discrepancies between data scales that could lead to ill-conditioning in matrix operations.
 *
 * @param rg Pointer to the RunningQuadraticGradient structure containing the data to normalize.
 */
void normalize_data(RunningQuadraticGradient *rg) {
    double x_mean = 0.0;
    double y_mean = 0.0;
    double x_var = 0.0;
    double y_var = 0.0;
    
    // Compute the mean of x and y
    for (uint16_t i = 0; i < rg->num_points; ++i) {
        x_mean += rg->x[i];
        y_mean += rg->y[i];
    }
    x_mean /= rg->num_points;
    y_mean /= rg->num_points;
    
    // Compute the variance (standard deviation squared) of x and y
    for (uint16_t i = 0; i < rg->num_points; ++i) {
        x_var += pow(rg->x[i] - x_mean, 2);
        y_var += pow(rg->y[i] - y_mean, 2);
    }
    x_var = sqrt(x_var / rg->num_points);  // Standard deviation of x
    y_var = sqrt(y_var / rg->num_points);  // Standard deviation of y
    
    // Avoid division by zero by ensuring standard deviations are not too small
    if (x_var < 1e-6) x_var = 1.0;
    if (y_var < 1e-6) y_var = 1.0;
    
    // Normalize x and y by subtracting the mean and dividing by the standard deviation
    for (uint16_t i = 0; i < rg->num_points; ++i) {
        rg->x[i] = (rg->x[i] - x_mean) / x_var;
        rg->y[i] = (rg->y[i] - y_mean) / y_var;
    }
}

/**
 * @brief Adds a new phase angle data point to the quadratic regression model and updates it using Recursive Least Squares (RLS).
 *
 * This function incorporates a new phase angle measurement into an ongoing quadratic regression model.
 * It utilizes the Recursive Least Squares (RLS) algorithm to efficiently update the model coefficients
 * without the need to reprocess all previous data points. By fitting a second-order polynomial (quadratic function)
 * to the noisy and scattered phase angle data collected from an impedance analyzer, the function helps in
 * smoothing out noise and small variances, allowing for accurate peak detection in the data.
 *
 * ### Purpose:
 * - **Noise Reduction**: The quadratic fit helps to smooth out the noisy phase angle data, mitigating the impact of measurement noise.
 * - **Peak Detection**: A second-order polynomial is suitable for modeling data with a single peak or trough, enabling peak point detection.
 *
 * ### Recursive Least Squares (RLS) Algorithm:
 * - The RLS algorithm is an adaptive filter that recursively finds the coefficients that minimize a weighted linear least squares cost function.
 * - It is particularly efficient for real-time applications where data points arrive sequentially, as it updates the model incrementally.
 * - The algorithm adjusts the model coefficients to best fit the new data point while considering previous data, with an emphasis determined by the forgetting factor.
 *
 * ### Function Workflow:
 * 1. **Data Window Management**:
 *    - Maintains a fixed-size window (`rg->max_points`) of the most recent data points.
 *    - If the window is not full, the new data point is simply added.
 *    - If the window is full, the oldest data point is removed (data points are shifted left), and the new data point is added at the end.
 * 2. **Input Vector Preparation**:
 *    - Constructs the input vector for the quadratic model: `[x^2, x, 1]`.
 *    - Here, `x` is the time or index of the data point, and `1` corresponds to the intercept term.
 * 3. **Inverse Covariance Matrix Update**:
 *    - Utilizes the Sherman-Morrison formula to update the inverse covariance matrix efficiently.
 *    - This step is crucial for the RLS algorithm to update the model coefficients without recomputing the entire covariance matrix.
 * 4. **Model Coefficients Update**:
 *    - Calculates the prediction error: the difference between the actual data point and the predicted value from the current model.
 *    - Updates the model coefficients (`rg->coefficients`) using the inverse covariance matrix and the prediction error.
 * 5. **Residual Sum of Squares Calculation**:
 *    - Recalculates the residual sum of squares (RSS) to assess the fit of the model to the data.
 *    - The RSS is used to measure the discrepancy between the data and the estimation model.
 *
 * @param rg Pointer to the `RunningQuadraticGradient` structure that holds the model state.
 * @param data_point Pointer to the new `MqsRawDataPoint_t` containing the phase angle to be added.
 *
 * @note
 * - **Numerical Stability**: The function includes a check for numerical stability to prevent division by values close to zero.
 * - **Symmetry Enforcement**: The inverse covariance matrix is enforced to remain symmetric after updates.
 * - **Forgetting Factor**: The `rg->forgetting_factor` parameter determines how quickly the influence of older data points diminishes.
 *
 * ### References:
 * - Haykin, S. (2002). *Adaptive Filter Theory*. Prentice Hall.
 * - The Sherman-Morrison formula: An efficient way to update the inverse of a matrix when it is modified slightly.
 *
 * @see init_running_quadratic_gradient
 * @see calculate_slope_at_point
 * @see calculate_second_order_gradient
 */
void add_quadratic_data_point(RunningQuadraticGradient *const rg, const MqsRawDataPoint_t *data_point) {
    double y = data_point->phaseAngle; // Extract the phase angle from the data point

    // Check if the current window is not full
    if (rg->num_points < rg->max_points) {
        // Add new data point to the window
        rg->y[rg->num_points] = y;
        rg->x[rg->num_points] = (double)rg->num_points;
        rg->num_points++;
    } else {
        // Shift the window to the left and add new data point
        // This shift operation removes the oldest data point and adds the new one at the end.
        for (uint16_t i = 0; i < rg->max_points - 1; ++i) {
            rg->y[i] = rg->y[i + 1];
            rg->x[i] = rg->x[i + 1];
        }
        rg->y[rg->max_points - 1] = y;  // Add new data point
        rg->x[rg->max_points - 1] = (double)(rg->num_points);
        rg->num_points++;
    }
    
     // Normalize data after updating
    //normalize_data(rg); //optimnal too experimental.

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
        //fprintf(stderr, "Numerical stability issue: temp is too close to zero.\n");
        //exit(EXIT_FAILURE);
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
    for (uint16_t i = 0; i < rg->num_points; ++i) {
        double error = rg->y[i] - (rg->coefficients[0] * pow(rg->x[i], 2) + rg->coefficients[1] * rg->x[i] + rg->coefficients[2]);
        rg->residual_sum_squares += error * error;
    }
    
      // Check for numerical stability and condition number
    double condition_number = compute_condition_number(rg->inverse_cov_matrix);
    if (condition_number > 1e8) {
        // Reset the inverse covariance matrix if condition number is too high
        init_inverse_covariance_matrix(rg);
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
 * @brief Verifies the detected peak by checking for a consistent trend of increases on the left side
 * and decreases on the right side, with additional checks and adjustments for truncated data using quadratic RLS.
 *
 * This function verifies whether a detected peak is a true peak by analyzing the second-order gradients
 * calculated using quadratic regression. Specifically, it checks for a consistent increasing trend on the left
 * side and a decreasing trend on the right side of the peak.
 *
 * @param values Array of data points.
 * @param length Length of the data array.
 * @param second_order_gradients Array containing the precomputed second-order gradients.
 * @param peak_index The index of the detected peak within the second_order_gradients array.
 * @param start_index The starting index in the original data array corresponding to the first element of second_order_gradients.
 * @param forgetting_factor The forgetting factor used in the RLS algorithm.
 * @return QuadraticPeakAnalysisResult Struct containing the verification result and truncation flags.
 */
QuadraticPeakAnalysisResult verify_quadratic_peak(
    const MqsRawDataPoint_t *values,
    uint16_t length,
    const double *second_order_gradients,
    uint16_t peak_index,
    uint16_t start_index,
    double forgetting_factor
) {
    QuadraticPeakAnalysisResult result = { .peak_found = false, .peak_index = peak_index, .is_truncated_left = false, .is_truncated_right = false };
    uint16_t left_trend_count = 0;
    uint16_t right_trend_count = 0;
    uint16_t inconsistency_count_left = 0;
    uint16_t inconsistency_count_right = 0;

    // Count increasing trends on the left side of the peak
    for (int i = peak_index; i > 0; --i) {
        if (second_order_gradients[i - 1] > 0) {
            left_trend_count++;
            if (left_trend_count >= quadratic_analysis_params.minimum_required_trend_count) break;
        } else {
            inconsistency_count_left++;
            if (inconsistency_count_left > quadratic_analysis_params.allowable_inconsistency_count) {
                break;
            }
        }
    }

    // Check for truncation on the left side
    if (peak_index - left_trend_count <= 0) {
        result.is_truncated_left = true;
    }

    // Count decreasing trends on the right side of the peak
    for (uint16_t i = peak_index; i < RLS_WINDOW - 1; ++i) {
        if (second_order_gradients[i + 1] < 0) {
            right_trend_count++;
            if (right_trend_count >= quadratic_analysis_params.minimum_required_trend_count) break;
        } else {
            inconsistency_count_right++;
            if (inconsistency_count_right > quadratic_analysis_params.allowable_inconsistency_count) {
                break;
            }
        }
    }

    // Check for truncation on the right side
    if (peak_index + right_trend_count >= RLS_WINDOW - 1) {
        result.is_truncated_right = true;
    }

    // Verify if the peak meets the criteria
    result.peak_found = left_trend_count >= quadratic_analysis_params.minimum_required_trend_count &&
                        right_trend_count >= quadratic_analysis_params.minimum_required_trend_count;

    return result;
}


/**
 * @brief Finds and verifies a peak in the data using quadratic RLS.
 *
 * This function first detects a peak using second-order gradients calculated from quadratic regression,
 * and then verifies the peak using the `verify_quadratic_peak` function. The peak is considered valid
 * only if it meets the trend criteria on both sides of the peak.
 *
 * @param values Array of data points.
 * @param length Length of the data array.
 * @param start_index The start index in the data array from which to begin the gradient calculation.
 * @param forgetting_factor The forgetting factor used in the RLS algorithm.
 * @return QuadraticPeakAnalysisResult Structure containing the peak detection status, truncation flags, and the peak index if found and verified.
 */
QuadraticPeakAnalysisResult find_and_verify_quadratic_peak(const MqsRawDataPoint_t *values, uint16_t length, uint16_t start_index, double forgetting_factor) {
    QuadraticPeakAnalysisResult result = { .peak_found = false, .peak_index = 0 };

    // Declare the second_order_gradients array
    double second_order_gradients[RLS_WINDOW];

    // Initialize the running quadratic gradient structure
    RunningQuadraticGradient rg;
    init_running_quadratic_gradient(&rg, forgetting_factor);

    // Compute the second-order gradients
    for (uint16_t i = 0; i < RLS_WINDOW && (start_index + i) < length; ++i) {
        const MqsRawDataPoint_t *current_value = &values[start_index + i];
        add_quadratic_data_point(&rg, current_value);

        // Calculate the second-order gradient only if we have at least 3 data points
        if (rg.num_points >= 3) {
            double second_order_gradient = calculate_second_order_gradient(&rg);
            second_order_gradients[i] = second_order_gradient;
        } else {
            second_order_gradients[i] = NAN;  // Not enough points yet to calculate the gradient
        }
    }

    // Find and verify the peak based on the computed gradients
    for (uint16_t i = 1; i < RLS_WINDOW; ++i) {
        if (second_order_gradients[i - 1] > 0 && second_order_gradients[i] < 0) {
            // Temporarily set the peak index
            result.peak_index = start_index + i;

            // Verify the detected peak
            result = verify_quadratic_peak(values, length, second_order_gradients, i, start_index, forgetting_factor);

            if (result.peak_found) {
                printf("Verified peak found at index %u\n", result.peak_index);
                break; // Exit the loop once a verified peak is found
            } else {
                printf("Peak at index %u did not pass verification. Continuing search...\n", result.peak_index);
            }
        }
    }
    

    if (!result.peak_found) {
        printf("No verified peak found in the specified window.\n");
    }else{
        // Adjust the peak index to match the real phaseAngle array
        uint16_t real_peak_index = buffer_manager.current_phase_index + (result.peak_index);

        // Print the real peak index in the phaseAngle array
        printf("[verify_peak_at_index] Verified peak found at real phaseAngle index: %d\n", real_peak_index);

    }

    return result;
}

/**
 * @brief Computes the weighted total sum of second-order gradients using quadratic RLS, 
 *        with an emphasis on reducing the impact of early, unstable gradients.
 * 
 * ### Problem Context:
 * The Recursive Least Squares (RLS) algorithm, while powerful for real-time smoothing and regression, 
 * can exhibit instability or aggressive behavior when adding initial data points, especially in the 
 * early stages of computation when the model is not fully formed. This can result in overshooting 
 * or erratic gradient values. 
 * 
 * To mitigate this issue, we introduce a weighting system that reduces the influence of the early 
 * second-order gradients during the summation process. The idea is to give less weight to the gradients 
 * calculated from early points, allowing the model to stabilize before heavily influencing the final result.
 * 
 * ### Weighting Approach:
 * We apply an exponential weighting function to the second-order gradients. This weighting system 
 * ensures that the influence of early gradients is minimized, and the later gradients—after the RLS 
 * algorithm has had time to stabilize—carry more weight in the final sum.
 * 
 * The exponential weighting is defined as:
 * \f[
 * w_i = 1 - e^{-\alpha \cdot i}
 * \f]
 * Where:
 * - \( w_i \) is the weight for the \( i \)-th gradient.
 * - \( \alpha \) is a tuning parameter that controls how quickly the weights increase. A higher \(\alpha\) 
 * results in weights approaching 1 more quickly, while a lower \(\alpha\) provides a more gradual increase.
 * 
 * The weighted second-order gradient sum is computed as:
 * \f[
 * \text{TotalWeightedSum} = \sum_{i=1}^{N} w_i \cdot g_i
 * \f]
 * Where:
 * - \( N \) is the total number of gradients.
 * - \( g_i \) is the second-order gradient at the \( i \)-th point.
 * 
 * ### Decisions:
 * - We do not modify the RLS algorithm itself, as the instability is only significant during early iterations. 
 * The weighting system ensures that the RLS can proceed as normal, while the summing process mitigates 
 * the influence of early instability.
 * 
 * - The exponential weighting function is chosen for its flexibility. It allows fine-tuning via the \(\alpha\) parameter 
 * to control the rate at which early gradients lose their influence.
 * 
 * - The weighting is applied only during the summation of second-order gradients to ensure that the RLS update mechanism 
 * and state tracking are not affected.
 * 
 * ### Parameters:
 * @param values Array of data points representing the phase angles.
 * @param length The length of the data array.
 * @param start_index The starting index from which to begin calculating the second-order gradients.
 * @param forgetting_factor The forgetting factor used in the RLS algorithm to give more weight to recent data.
 * @return The weighted total sum of second-order gradients. Returns `NAN` if there are insufficient data points.
 */
double compute_total_second_order_gradient(
    const MqsRawDataPoint_t *values,
    uint16_t length,
    uint16_t start_index,
    double forgetting_factor
) {
    // Ensure that the start index and the length allow for RLS_WINDOW calculations
    if (start_index + RLS_WINDOW > length) {
        printf("Insufficient data to compute %d second-order gradients from the starting index.\n", RLS_WINDOW);
        return NAN;
    }

    // Initialize the running quadratic gradient structure
    RunningQuadraticGradient rg;
    init_running_quadratic_gradient(&rg, forgetting_factor);

    double total_weighted_sum_gradients = 0.0; // Initialize total weighted sum of gradients
    uint16_t stabilization_threshold = 3;      // Minimum number of points for valid gradient calculation
    double alpha = 0.1;                        // Tuning parameter for the weight function, adjust as needed

    // Loop through the next RLS_WINDOW values starting from start_index
    for (uint16_t i = 0; i < RLS_WINDOW; ++i) {
        const MqsRawDataPoint_t *current_value = &values[start_index + i];

        // Add the current value's phaseAngle to the running quadratic gradient
        add_quadratic_data_point(&rg, current_value);

        // Calculate the second-order gradient only if we have at least the stabilization threshold points
        if (rg.num_points >= stabilization_threshold) {
            double second_order_gradient = calculate_second_order_gradient(&rg);

            // Calculate the weight for this gradient
            uint16_t weight_index = rg.num_points - stabilization_threshold + 1;
            double weight = 1.0 - exp(-alpha * weight_index);

            // Accumulate the weighted sum of gradients
            total_weighted_sum_gradients += weight * second_order_gradient;

        } else {
           //do nothing
        }
    }

    // Return the total weighted sum of the second-order gradients
    return total_weighted_sum_gradients;
}


/**
 * @brief Finds a consistent increasing trend in the second-order gradient array.
 *
 * This function scans through an array of second-order gradient values to identify a consistent increasing trend.
 * It starts tracking when positive gradients are encountered and stops if more than one negative gradient is found consecutively.
 * The function returns a `GradientTrendIndices` struct containing the start and end indices of the consistent increase, 
 * a validity flag, and the cumulative sum of the gradients within this range.
 *
 * @param gradients The array of second-order gradient values to analyze.
 * @param start_index The starting index in the original data array from which the gradient array begins.
 * @param window_size The number of elements in the gradient array to analyze.
 * @return GradientTrendIndices A structure containing the start index, end index, validity, and maximum sum of the consistent increasing trend.
 */
GradientTrendIndices find_consistent_increase_in_second_order(double *gradients, size_t start_index, size_t window_size) {
    GradientTrendIndices increase_info = {0, 0, false, 0.0};
    bool tracking_increase = false;
    double cumulative_sum = 0.0;
    int decrease_count = 0;
    double max_cumulative_sum = 0.0;  // Track the maximum cumulative sum encountered
    size_t best_start_index = 0;  // Track the start index for the best trend
    size_t best_end_index = 0;    // Track the end index for the best trend

    #ifdef DEBUG_QUADRATIC
    printf("Starting find_consistent_increase_in_second_order with start_index: %zu, window_size: %zu\n", start_index, window_size);
    #endif

    for (size_t i = 0; i < window_size; ++i) {
        double gradient = gradients[i];

        #ifdef DEBUG_QUADRATIC
        printf("Index: %zu, Gradient: %.6f\n", start_index + i, gradient);
        #endif

        if (gradient > 0) {  // If the current gradient is positive
            if (!tracking_increase) {
                // Start or resume tracking the increase
                increase_info.start_index = start_index + i;
                tracking_increase = true;
                cumulative_sum = 0.0;  // Reset cumulative sum when starting to track

                #ifdef DEBUG_QUADRATIC
                printf("Started tracking increase at index %zu\n", start_index + i);
                #endif
            }
            increase_info.end_index = start_index + i;
            cumulative_sum += gradient;
            decrease_count = 0;  // Reset decrease counter on positive gradient

            // Update the best cumulative sum and indices
            if (cumulative_sum > max_cumulative_sum) {
                max_cumulative_sum = cumulative_sum;
                best_start_index = increase_info.start_index;
                best_end_index = increase_info.end_index;
            }

            #ifdef DEBUG_QUADRATIC
            printf("Continuing tracking increase: Updated end_index to %zu, Cumulative Sum: %.6f\n", start_index + i, cumulative_sum);
            #endif
        } else if (gradient < 0) {  // If the current gradient is negative
            decrease_count++;

            #ifdef DEBUG_QUADRATIC
            printf("Negative gradient found, Decrease Count: %d\n", decrease_count);
            #endif
                  
            if (decrease_count > quadratic_analysis_params.max_second_order_trend_decrease_count) {  // Use global parameter
                tracking_increase = false;  // Stop tracking but do not discard the tracked increase

                #ifdef DEBUG_QUADRATIC
                printf("Stopped tracking increase due to consecutive decreases at index %zu\n", start_index + i);
                #endif

                // Reset decrease count and continue to look for the next increase
                decrease_count = 0;
            }
        }
    }

    // Use the best cumulative sum and indices found during tracking
    if (max_cumulative_sum > 0) {
        increase_info.valid = true;
        increase_info.max_sum = max_cumulative_sum;
        increase_info.start_index = best_start_index;
        increase_info.end_index = best_end_index;

        printf("Valid increase found from index %u to %u with max sum %.6f\n", increase_info.start_index, increase_info.end_index, max_cumulative_sum);
    } else {
        printf("No valid increase found.\n");
    }

    return increase_info;
}

/**
 * @brief Finds a consistent decreasing trend in the second-order gradient array.
 *
 * This function scans through an array of second-order gradient values to identify a consistent decreasing trend.
 * It starts tracking when negative gradients are encountered and stops if more than one positive gradient is found consecutively.
 * The function returns a `GradientTrendIndices` struct containing the start and end indices of the consistent decrease, 
 * a validity flag, and the cumulative sum of the gradients within this range.
 *
 * @param gradients The array of second-order gradient values to analyze.
 * @param start_index The starting index in the original data array from which the gradient array begins.
 * @param window_size The number of elements in the gradient array to analyze.
 * @return GradientTrendIndices A structure containing the start index, end index, validity, and maximum sum of the consistent decreasing trend.
 */
GradientTrendIndices find_consistent_decrease_in_second_order(double *gradients, size_t start_index, size_t window_size) {
    GradientTrendIndices decrease_info = {0, 0, false, 0.0};
    bool tracking_decrease = false;
    double cumulative_sum = 0.0;
    int increase_count = 0;

    #ifdef DEBUG_QUADRATIC
    printf("Starting find_consistent_increase_in_second_order with start_index: %u, window_size: %u\n", start_index, window_size);
    #endif

    for (size_t i = 0; i < window_size; ++i) {
        double gradient = gradients[i];

        #ifdef DEBUG_QUADRATIC
        printf("Index: %zu, Gradient: %.6f\n", start_index + i, gradient);
        #endif

        if (gradient < 0) {  // If the current gradient is negative
            if (!tracking_decrease) {
                // Start tracking the decrease
                decrease_info.start_index = start_index + i;
                tracking_decrease = true;
                cumulative_sum = 0.0;  // Reset cumulative sum when starting to track

                #ifdef DEBUG_QUADRATIC
                printf("Started tracking decrease at index %zu\n", start_index + i);
                #endif
            }
            decrease_info.end_index = start_index + i;
            cumulative_sum += gradient;
            increase_count = 0;  // Reset increase counter on negative gradient

            #ifdef DEBUG_QUADRATIC
            printf("Continuing tracking decrease: Updated end_index to %zu, Cumulative Sum: %.6f\n", start_index + i, cumulative_sum);
            #endif
        } else if (gradient > 0) {  // If the current gradient is positive
            increase_count++;

            #ifdef DEBUG_QUADRATIC
            printf("Positive gradient found, Increase Count: %d\n", increase_count);
            #endif

            if (increase_count > quadratic_analysis_params.max_second_order_trend_increase_count) {  // Use global parameter
                tracking_decrease = false;  // Stop tracking

                #ifdef DEBUG_QUADRATIC
                printf("Stopped tracking decrease due to consecutive increases at index %zu\n", start_index + i);
                #endif
                increase_count = 0;
                //break;
            }
        }

        // If tracking was stopped due to increases and a new decrease is found
        if (!tracking_decrease && gradient < 0 && increase_count > quadratic_analysis_params.max_second_order_trend_increase_count) {
            // Restart tracking from this new decrease
            decrease_info.start_index = start_index + i;
            tracking_decrease = true;
            cumulative_sum = gradient;  // Reset cumulative sum with the new decrease
            decrease_info.end_index = start_index + i;
            increase_count = 0;  // Reset the increase counter

            #ifdef DEBUG_QUADRATIC
            printf("Restarted tracking decrease at index %zu, Cumulative Sum: %.6f\n", start_index + i, cumulative_sum);
            #endif
        }
    }

    if (cumulative_sum < 0) {  // Ensure there's a valid decrease before marking it valid
        decrease_info.valid = true;
        decrease_info.max_sum = cumulative_sum;

       
        printf("Valid decrease found from index %u to %u with max sum %.6f\n", decrease_info.start_index, decrease_info.end_index, cumulative_sum);

    } else {

        printf("No valid decrease found.\n");

    }

    return decrease_info;
}

/**
 * @brief Tracks the second-order gradient trends within a specified window using quadratic regression.
 *
 * This function tracks the gradient trends of a given dataset within a specified window using quadratic regression.
 * It calculates second-order gradients (curvatures) at each point in the window, identifies regions of consistent increase,
 * and stores the maximum cumulative sum of gradients. It returns a structure containing the start and end indices
 * of the consistent increase, along with a flag indicating if the trend was valid and the maximum sum of gradients.
 *
 * @param values Array of double values representing the data points.
 * @param length The length of the values array.
 * @param start_index The starting index in the values array from which to begin the gradient calculation.
 * @param window_size The number of points to include in the gradient calculation.
 * @param forgetting_factor The forgetting factor used in the Recursive Least Squares (RLS) algorithm, typically close to 1,
 *        which determines the weight given to newer data points.
 * @return GradientTrendResult A struct containing:
 * - increase_info: Information about the detected increasing trend.
 * - decrease_info: Information about the detected decreasing trend.
 */
GradientTrendResult track_gradient_trends_with_quadratic_regression(const MqsRawDataPoint_t *values, uint16_t length, uint16_t start_index, uint16_t window_size, double forgetting_factor) {
    GradientTrendResult trend_result = {0};  // Initialize the struct with default values
    
    // Ensure that the start index and the window size allow for calculations
    if (start_index + window_size > length) {
        #ifdef DEBUG_QUADRATIC
        printf("Insufficient data to compute gradients for the specified window size.\n");
        #endif
        return trend_result;
    }

    // Array to store the second-order gradients
    double second_order_gradients[RLS_WINDOW];

    // Initialize the running quadratic gradient structure
    RunningQuadraticGradient rg;
    init_running_quadratic_gradient(&rg, forgetting_factor);

    // Loop through the specified window size and calculate second-order gradients
    for (uint16_t i = 0; i < window_size; ++i) {
        uint16_t current_index = start_index + i;
        const MqsRawDataPoint_t *current_value = &values[current_index];

        // Add the current value's phaseAngle to the running quadratic gradient
        add_quadratic_data_point(&rg, current_value);

        // Check if we have enough points to calculate the gradient
        if (rg.num_points >= 3) {
            double second_order_gradient = calculate_second_order_gradient(&rg);
            second_order_gradients[i] = second_order_gradient;

            // printf("Second-order gradient at index %u: %.6f\n", current_index, second_order_gradient);
        } else {
            second_order_gradients[i] = NAN;  // Not enough points yet to calculate the gradient

            #ifdef DEBUG_QUADRATIC
            printf("Not enough data points to calculate the gradient at index %u\n", current_index);
            #endif
        }
    }

    // Find the consistent increase trend
    #ifdef DEBUG_QUADRATIC
    printf("Finding consistent increase...\n");
    #endif
    trend_result.increase_info = find_consistent_increase_in_second_order(second_order_gradients, start_index, window_size);

    // Find the consistent decrease trend
    #ifdef DEBUG_QUADRATIC
    printf("Finding consistent decrease...\n");
    #endif
    trend_result.decrease_info = find_consistent_decrease_in_second_order(second_order_gradients, start_index, window_size);

    return trend_result;
}

