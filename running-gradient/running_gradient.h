#ifndef RUNNING_GRADIENT_H
#define RUNNING_GRADIENT_H

#include <stddef.h>

// Structure to hold the running gradient calculation data
typedef struct {
    size_t num_points;
    double inverse_cov_matrix[2][2];
    double coefficients[2];
    double residual_sum_squares;
} RunningGradient;

/**
 * @brief Initializes the RunningGradient structure.
 * 
 * @param rg Pointer to the RunningGradient structure.
 */
void init_running_gradient(RunningGradient *rg);

/**
 * @brief Partitions the array for quickselect.
 * 
 * @param arr The array.
 * @param left The left index.
 * @param right The right index.
 * @param pivotIndex The pivot index.
 * @return The index of the pivot after partition.
 */
size_t partition(double* arr, size_t left, size_t right, size_t pivotIndex);

/**
 * @brief Quickselect algorithm to find the k-th smallest element.
 * 
 * @param arr The array.
 * @param left The left index.
 * @param right The right index.
 * @param k The index to find.
 * @return The k-th smallest element.
 */
double quickselect(double* arr, size_t left, size_t right, size_t k);

/**
 * @brief Adds a new data point to the RunningGradient structure.
 * 
 * @param rg Pointer to the RunningGradient structure.
 * @param y The new data point value.
 */
void add_data_point(RunningGradient *rg, double y);

/**
 * @brief Calculates the current gradient (slope) of the fitted line.
 * 
 * @param rg Pointer to the RunningGradient structure.
 * @return The gradient (slope) of the fitted line.
 */
double calculate_gradient(const RunningGradient *rg);

/**
 * @brief Calculates the standard error of the gradient estimate.
 * 
 * @param rg Pointer to the RunningGradient structure.
 * @return The standard error of the gradient estimate.
 */
double calculate_standard_error(const RunningGradient *rg);

/**
 * @brief Calculates the normal cumulative distribution function.
 * 
 * @param value The value to calculate the CDF for.
 * @param mean The mean of the normal distribution.
 * @param stddev The standard deviation of the normal distribution.
 * @return The cumulative probability for the given value.
 */
double normal_cdf(double value, double mean, double stddev);

/**
 * @brief Calculates the probability that the gradient is less than a threshold.
 * 
 * @param rg Pointer to the RunningGradient structure.
 * @param thresh The threshold for the gradient.
 * @return The probability that the gradient is less than the threshold.
 */
double probability_gradient_less_than(const RunningGradient *rg, double thresh);

/**
 * @brief Calculates the probability that the gradient is greater than a threshold.
 * 
 * @param rg Pointer to the RunningGradient structure.
 * @param thresh The threshold for the gradient.
 * @return The probability that the gradient is greater than the threshold.
 */
double probability_gradient_greater_than(const RunningGradient *rg, double thresh);

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
double find_upper_quantile(const double* data, size_t size, double quantile);

/**
 * @brief Counts the number of steps without a noticeable decrease, robust to outliers.
 * 
 * This function counts the number of consecutive data points, starting from the end of the array,
 * that do not show a significant decreasing trend. It discards the upper quantile of data points 
 * to mitigate the impact of outliers. The function returns the count of such steps.
 * 
 * @param data The array of data points.
 * @param size The number of data points.
 * @param probability_of_decrease The probability threshold to consider a decrease.
 * @param quantile_discard The upper quantile of data points to discard (e.g., 0.10 for the top 10%).
 * @return The number of steps without a noticeable decrease.
 * 
 * @note This function is useful for analyzing trends in noisy data sets, where outliers can distort the overall trend.
 */
size_t count_steps_without_decrease_robust(const double* data, size_t size, double probability_of_decrease, double quantile_discard);

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
double probability_values_are_increasing_robust(const double* data, size_t size, double quantile_discard);

/**
 * @brief Processes the data in chunks of 10 and calculates various statistics for each chunk.
 * 
 * @param data The data points.
 * @param size The number of data points.
 */
void process_data_in_chunks(const double *data, size_t size);

#endif // RUNNING_GRADIENT_H
