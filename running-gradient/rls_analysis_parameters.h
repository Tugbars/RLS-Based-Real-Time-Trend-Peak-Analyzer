#ifndef RLS_ANALYSIS_PARAMETERS_H
#define RLS_ANALYSIS_PARAMETERS_H

#include <stdint.h>  // Include for uint16_t

// Struct to hold parameters for cubic RLS analysis
typedef struct {
    double significance_threshold;
    uint16_t duration_threshold;
    int minimum_required_cubic_trend_count;
    int max_third_order_trend_decrease_count;
    int max_third_order_trend_increase_count;
} cubic_rls_analysis_parameters;

// Struct to hold parameters for quadratic RLS analysis
typedef struct {
    double centered_gradient_sum;            // Threshold to consider peak centered (e.g., 1.0)
    int max_second_order_trend_decrease_count; // Max allowed decreases when tracking increase
    int max_second_order_trend_increase_count; // Max allowed increases when tracking decrease
    int minimum_required_trend_count;        // Minimum trends to consider peak valid (e.g., 5)
    int allowable_inconsistency_count;       // Allowed inconsistencies in trend detection (e.g., 2)
} quadratic_rls_analysis_parameters;

// Struct to hold parameters for peak analysis
typedef struct {
    double min_average_increase;         // Minimum average increase required to consider a significant trend
    double min_average_decrease;         // Minimum average decrease required to consider a significant trend
    uint16_t min_consistent_trend_count; // Minimum number of data points for a consistent trend
} on_peak_analysis_parameters;

// Struct to hold parameters for gradient analysis
typedef struct {
    double gradient_threshold;         // Threshold for gradient detection during analysis
    double minimum_gradient_total;     // Minimum total gradient required to flag significant trends
} GradientAnalysisParams;

// Declare the global variables for the parameters
extern cubic_rls_analysis_parameters cubic_analysis_params;
extern quadratic_rls_analysis_parameters quadratic_analysis_params;
extern on_peak_analysis_parameters peak_analysis_params;
extern GradientAnalysisParams gradient_analysis_params;

// Function to initialize the cubic RLS analysis parameters
void init_cubic_rls_analysis_parameters(double significance_thresh, uint16_t duration_thresh, uint16_t min_trend_count, int max_decrease_count, int max_increase_count);

// Function to initialize the quadratic RLS analysis parameters
void init_quadratic_rls_analysis_parameters(
    double centered_gradient_sum,
    int max_decrease_count,
    int max_increase_count,
    int min_trend_count,
    int allowable_inconsistency_count
);

// Function to initialize the on-peak analysis parameters
void init_on_peak_analysis_parameters(double min_avg_increase, double min_avg_decrease, uint16_t min_trend_count);

void init_gradient_analysis_params(double gradient_thresh, double min_gradient_total);

#endif // RLS_ANALYSIS_PARAMETERS_H

