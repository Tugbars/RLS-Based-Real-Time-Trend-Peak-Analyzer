#include "rls_analysis_parameters.h"
#include <stdint.h> // Include this to ensure uint16_t is defined
#include <stdio.h>

/**
 * @brief Global variable to store the parameters for cubic RLS analysis.
 *
 * This variable holds the configuration parameters that influence the behavior
 * of cubic Recursive Least Squares (RLS) analysis. It includes thresholds and
 * criteria used to determine significant trends in the data.
 */
cubic_rls_analysis_parameters cubic_analysis_params;

/**
 * @brief Global variable to store the parameters for quadratic RLS analysis.
 *
 * This variable holds the configuration parameters that influence the behavior
 * of quadratic Recursive Least Squares (RLS) analysis. It includes thresholds and
 * criteria used to identify and track trends in second-order gradients.
 */
quadratic_rls_analysis_parameters quadratic_analysis_params;


/**
 * @brief Global variable to store the parameters for peak analysis.
 */
on_peak_analysis_parameters peak_analysis_params;


/**
 * @brief Global variable to store the gradient analysis parameters.
 *
 * This variable holds the configuration parameters for gradient thresholding during analysis.
 * It determines the thresholds for increase detection and minimum total gradient required to flag significant trends.
 */
GradientAnalysisParams gradient_analysis_params;

/**
 * @brief Initializes the cubic RLS analysis parameters.
 *
 * This function sets up the configuration parameters for cubic RLS analysis.
 * These parameters determine how trends are detected and analyzed within the data.
 *
 * @param significance_thresh The threshold for determining if a gradient trend is significant.
 *                            A higher value means that only stronger trends will be considered significant.
 * @param duration_thresh The minimum number of consecutive points required for a trend to be considered.
 *                        This ensures that only sustained trends are analyzed.
 * @param min_trend_count The minimum number of consistent trends required for analysis.
 *                        This helps in filtering out noise and minor fluctuations.
 */
void init_cubic_rls_analysis_parameters(double significance_thresh, uint16_t duration_thresh, uint16_t min_trend_count, int max_decrease_count, int max_increase_count) {
    cubic_analysis_params.significance_threshold = significance_thresh;
    cubic_analysis_params.duration_threshold = duration_thresh;
    cubic_analysis_params.minimum_required_cubic_trend_count = min_trend_count;
    cubic_analysis_params.max_third_order_trend_decrease_count = max_decrease_count;
    cubic_analysis_params.max_third_order_trend_increase_count = max_increase_count;
}

/**
 * @brief Initializes the quadratic RLS analysis parameters.
 *
 * This function sets up the configuration parameters for quadratic RLS analysis.
 * These parameters influence how second-order gradients are analyzed and how trends are detected.
 *
 * @param min_gradient_sum The minimum cumulative sum of second-order gradients required for a trend to be considered valid.
 *                         This helps in identifying significant concave or convex regions.
 * @param max_decrease_count The maximum number of consecutive negative gradients allowed when tracking an increasing trend.
 *                           This allows for minor fluctuations without discarding a potentially valid trend.
 * @param max_increase_count The maximum number of consecutive positive gradients allowed when tracking a decreasing trend.
 *                           This helps to filter out noise when identifying a decreasing trend.
 */
void init_quadratic_rls_analysis_parameters(
    double centered_gradient_sum,
    int max_decrease_count,
    int max_increase_count,
    int min_trend_count,
    int allowable_inconsistency_count
) {
    quadratic_analysis_params.centered_gradient_sum = centered_gradient_sum;
    quadratic_analysis_params.max_second_order_trend_decrease_count = max_decrease_count;
    quadratic_analysis_params.max_second_order_trend_increase_count = max_increase_count;
    quadratic_analysis_params.minimum_required_trend_count = min_trend_count;
    quadratic_analysis_params.allowable_inconsistency_count = allowable_inconsistency_count;
}

/**
 * @brief Initializes the on-peak analysis parameters.
 */
void init_on_peak_analysis_parameters(double min_avg_increase, double min_avg_decrease, uint16_t min_trend_count) {
    peak_analysis_params.min_average_increase = min_avg_increase;
    peak_analysis_params.min_average_decrease = min_avg_decrease;
    peak_analysis_params.min_consistent_trend_count = min_trend_count;
}

/**
 * @brief Initializes the parameters for gradient analysis.
 *
 * @param gradient_thresh The threshold for the gradient to consider a significant trend.
 * @param min_gradient_total The minimum total gradient required to flag a significant trend.
 */
void init_gradient_analysis_params(double gradient_thresh, double min_gradient_total) {
    gradient_analysis_params.gradient_threshold = gradient_thresh;
    gradient_analysis_params.minimum_gradient_total = min_gradient_total;
}
