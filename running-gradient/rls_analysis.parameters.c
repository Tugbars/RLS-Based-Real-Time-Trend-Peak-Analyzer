#include "rls_analysis_parameters.h"

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
void init_cubic_rls_analysis_parameters(double significance_thresh, size_t duration_thresh, size_t min_trend_count, int max_decrease_count, int max_increase_count) {
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
void init_quadratic_rls_analysis_parameters(double min_gradient_sum, int max_decrease_count, int max_increase_count) {
    quadratic_analysis_params.minimum_second_order_gradient_sum = min_gradient_sum;
    quadratic_analysis_params.max_second_order_trend_decrease_count = max_decrease_count;
    quadratic_analysis_params.max_second_order_trend_increase_count = max_increase_count;
}
