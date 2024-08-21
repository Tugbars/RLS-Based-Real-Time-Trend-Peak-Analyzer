#ifndef RUNNING_CUBIC_GRADIENT_H
#define RUNNING_CUBIC_GRADIENT_H

#include <stddef.h>  // For size_t

// Define the window size for the cubic RLS
#define CUBIC_RLS_WINDOW 30


typedef struct {
    bool peak_found;
    size_t peak_index;
} CubicPeakAnalysisResult;

/**
 * @struct RunningCubicGradient
 * @brief A structure to represent the state of a cubic regression model using Recursive Least Squares (RLS).
 *
 * This structure maintains the necessary data for performing cubic regression with RLS. It includes the coefficients
 * of the cubic model, the inverse covariance matrix for efficient coefficient updates, and other relevant data.
 */
typedef struct {
    double coefficients[4];  // Coefficients for x^3, x^2, x, and constant term
    double inverse_cov_matrix[4][4];  // 4x4 inverse covariance matrix
    double forgetting_factor;
    double y[CUBIC_RLS_WINDOW];  // Store y values
    double x[CUBIC_RLS_WINDOW];  // Store x values (could also calculate on the fly)
    double residual_sum_squares;
    size_t num_points;
    size_t max_points;
} RunningCubicGradient;


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
void compute_cubic_second_order_gradients(double *second_order_gradients, const double *values, size_t length, size_t start_index, double forgetting_factor);


double* compute_cubic_first_order_gradients(const double *values, size_t length, size_t start_index, double forgetting_factor);

CubicPeakAnalysisResult find_and_verify_cubic_peak(const double *values, size_t length, size_t start_index, double forgetting_factor);

#endif // RUNNING_CUBIC_GRADIENT_H
