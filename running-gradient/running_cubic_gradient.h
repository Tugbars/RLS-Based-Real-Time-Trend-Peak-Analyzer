#ifndef RUNNING_CUBIC_GRADIENT_H
#define RUNNING_CUBIC_GRADIENT_H

#include <stddef.h>  // For size_t

// Define the window size for the cubic RLS
#define CUBIC_RLS_WINDOW 30

/**
 * @struct CubicPeakAnalysisResult
 * @brief A structure to store the results of cubic peak analysis.
 *
 * This structure holds information about whether a peak was found during the cubic analysis and 
 * the index of that peak within the analyzed data.
 */
typedef struct {
    bool peak_found;     /**< Indicates if a peak was found during analysis. */
    size_t peak_index;   /**< The index of the detected peak within the data array. */
} CubicPeakAnalysisResult;

/**
 * @struct RunningCubicGradient
 * @brief A structure to represent the state of a cubic regression model using Recursive Least Squares (RLS).
 *
 * This structure maintains the necessary data for performing cubic regression with RLS. It includes the coefficients
 * of the cubic model, the inverse covariance matrix for efficient coefficient updates, and other relevant data.
 */
typedef struct {
    double coefficients[4];               /**< Coefficients for the cubic model: x^3, x^2, x, and constant term. */
    double inverse_cov_matrix[4][4];      /**< 4x4 inverse covariance matrix for RLS updates. */
    double forgetting_factor;             /**< The forgetting factor for the RLS algorithm, typically close to 1. */
    double y[CUBIC_RLS_WINDOW];           /**< Array to store y values used in the regression. */
    double x[CUBIC_RLS_WINDOW];           /**< Array to store x values used in the regression. */
    double residual_sum_squares;          /**< The sum of squared residuals for the regression model. */
    size_t num_points;                    /**< The number of data points currently in the window. */
    size_t max_points;                    /**< The maximum number of points allowed in the window (usually CUBIC_RLS_WINDOW). */
} RunningCubicGradient;

/**
 * @struct GradientTrendIndices
 * @brief A structure to store information about a detected trend (either increasing or decreasing).
 *
 * This structure contains the start and end indices of a detected trend, a validity flag, 
 * and the maximum cumulative sum of the gradients within this range.
 */
typedef struct {
    size_t start_index;    /**< The starting index of the detected trend in the data array. */
    size_t end_index;      /**< The ending index of the detected trend in the data array. */
    bool valid;            /**< A flag indicating whether the detected trend is valid. */
    double max_sum;        /**< The maximum cumulative sum of the gradients within the detected trend. */
} GradientTrendIndices;

/**
 * @struct GradientTrendResult
 * @brief A structure to store the results of tracking gradient trends.
 *
 * This structure contains information about both increasing and decreasing trends detected within a dataset.
 */
typedef struct {
    GradientTrendIndices increase_info;   /**< Information about the detected increasing trend. */
    GradientTrendIndices decrease_info;   /**< Information about the detected decreasing trend. */
} GradientTrendResult;

/**
 * @struct PeakTrendAnalysisResult
 * @brief A structure to store the results of analyzing trends for significant peaks.
 *
 * This structure contains detailed information about whether significant increasing or decreasing trends were found,
 * along with specific data about those trends.
 */
typedef struct {
    bool significant_increase;               /**< Flag indicating if a significant increasing trend was detected. */
    bool significant_decrease;               /**< Flag indicating if a significant decreasing trend was detected. */
    GradientTrendIndices increase_info;      /**< Detailed information about the significant increasing trend. */
    GradientTrendIndices decrease_info;      /**< Detailed information about the significant decreasing trend. */
} PeakTrendAnalysisResult;

void compute_cubic_second_order_gradients(double *second_order_gradients, const double *values, size_t length, size_t start_index, double forgetting_factor);

CubicPeakAnalysisResult find_and_verify_cubic_peak(const double *values, size_t length, size_t start_index, double forgetting_factor);

void compute_cubic_first_order_gradients(double *first_order_gradients, const double *values, size_t length, size_t start_index, size_t window_size, double forgetting_factor);

GradientTrendResult track_gradient_trends_with_cubic_regression(const double *values, size_t length, size_t start_index, size_t window_size, double forgetting_factor);

PeakTrendAnalysisResult detect_significant_gradient_trends(const double *values, size_t length, size_t start_index, size_t window_size, double forgetting_factor);

#endif // RUNNING_CUBIC_GRADIENT_H
