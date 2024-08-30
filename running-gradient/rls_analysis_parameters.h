// rls_analysis_parameters.h

#ifndef RLS_ANALYSIS_PARAMETERS_H
#define RLS_ANALYSIS_PARAMETERS_H

#include <stddef.h>

// Struct to hold parameters for cubic RLS analysis
typedef struct {
    double significance_threshold;
    size_t duration_threshold;
    int minimum_required_cubic_trend_count;
    int max_third_order_trend_decrease_count;
    int max_third_order_trend_increase_count;
} cubic_rls_analysis_parameters;

// Struct to hold parameters for quadratic RLS analysis
typedef struct {
    double minimum_second_order_gradient_sum;
    int max_second_order_trend_decrease_count;
    int max_second_order_trend_increase_count;
} quadratic_rls_analysis_parameters;

// Declare the global variables for the parameters
extern cubic_rls_analysis_parameters cubic_analysis_params;
extern quadratic_rls_analysis_parameters quadratic_analysis_params;

// Function to initialize the cubic RLS analysis parameters
void init_cubic_rls_analysis_parameters(double significance_thresh, size_t duration_thresh, size_t min_trend_count, int max_decrease_count, int max_increase_count);

// Function to initialize the quadratic RLS analysis parameters
void init_quadratic_rls_analysis_parameters(double min_gradient_sum, int max_decrease_count, int max_increase_count);

#endif // RLS_ANALYSIS_PARAMETERS_H

