#ifndef WINDOWED_RUNNING_GRADIENT_H
#define WINDOWED_RUNNING_GRADIENT_H

#include <stddef.h>
#include <stdint.h>
#include <stdbool.h>

#define WINDOW_SIZE 30
#define FORGETTING_FACTOR 0.5

/**
 * @struct RunningGradient
 * @brief A structure to represent the state of a running linear regression model using Recursive Least Squares (RLS).
 *
 * This structure maintains the necessary data for performing linear regression with RLS. It includes the coefficients
 * of the linear model, the residual sum of squares, the inverse covariance matrix for efficient coefficient updates,
 * and the arrays to store recent x and y values used in the regression.
 */
typedef struct {
    size_t num_points;                    /**< The number of data points currently in the window. */
    size_t max_points;                    /**< The maximum number of points allowed in the window (typically WINDOW_SIZE). */
    double coefficients[2];               /**< Coefficients for the linear model: slope and intercept. */
    double residual_sum_squares;          /**< The sum of squared residuals for the regression model. */
    double inverse_cov_matrix[2][2];      /**< 2x2 inverse covariance matrix for RLS updates. */
    double x[WINDOW_SIZE];                /**< Array to store x values used in the regression. */
    double y[WINDOW_SIZE];                /**< Array to store y values used in the regression. */
    double forgetting_factor;             /**< The forgetting factor for the RLS algorithm, typically close to 1. */
} RunningGradient;

/**
 * @enum PeakPosition
 * @brief An enumeration to represent the position relative to a peak.
 *
 * This enumeration is used to indicate whether a detected trend is on the left side, right side, or undecided relative to a peak.
 */
typedef enum {
    LEFT_SIDE = -1,                       /**< Indicates that the trend is on the left side of a peak. */
    RIGHT_SIDE = 1,                       /**< Indicates that the trend is on the right side of a peak. */
    UNDECIDED = 0,                         /**< Indicates that the trend's position relative to a peak is undecided. */
    ON_PEAK = 2                          // New enumeration value for when the trend is on the peak
} PeakPosition;

/**
 * @struct GradientComparisonResult
 * @brief A structure to store the results of comparing gradients between two parts of a dataset.
 *
 * This structure contains the total gradient, the count of increases, and the probability of an increasing trend 
 * for two parts of a dataset. It also includes information about which side (left or right) shows a stronger trend.
 */
typedef struct {
    double total_gradient_first_part;     /**< Total gradient for the first part (usually the right side) of the dataset. */
    double total_gradient_second_part;    /**< Total gradient for the second part (usually the left side) of the dataset. */
    size_t increase_count_first_part;     /**< Count of increases detected in the first part of the dataset. */
    size_t increase_count_second_part;    /**< Count of increases detected in the second part of the dataset. */
    double probability_increase_first_part; /**< Probability of an increasing trend in the first part of the dataset. */
    double probability_increase_second_part; /**< Probability of an increasing trend in the second part of the dataset. */
    PeakPosition dominant_side;           /**< Indicates which side (left or right) has a stronger increasing trend. */
} GradientComparisonResult;

GradientComparisonResult compare_gradient_parts(const double *data, size_t start_index, double forgetting_factor);
double calculate_gradient(const RunningGradient * const rg);

#endif // WINDOWED_RUNNING_GRADIENT_H
