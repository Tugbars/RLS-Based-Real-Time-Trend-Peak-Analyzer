#include "running_gradient.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

void init_running_gradient(RunningGradient *rg) {
    rg->num_points = 0;
    rg->inverse_cov_matrix[0][0] = 1e6; rg->inverse_cov_matrix[0][1] = 0;
    rg->inverse_cov_matrix[1][0] = 0;   rg->inverse_cov_matrix[1][1] = 1e6;
    rg->coefficients[0] = 0;
    rg->coefficients[1] = 0;
    rg->residual_sum_squares = 0.0;
}

size_t partition(double* arr, size_t left, size_t right, size_t pivotIndex) {
    double pivotValue = arr[pivotIndex];
    double temp = arr[pivotIndex];
    arr[pivotIndex] = arr[right];
    arr[right] = temp;  // Move pivot to end
    size_t storeIndex = left;
    for (size_t i = left; i < right; i++) {
        if (arr[i] < pivotValue) {
            temp = arr[i];
            arr[i] = arr[storeIndex];
            arr[storeIndex] = temp;
            storeIndex++;
        }
    }
    temp = arr[storeIndex];
    arr[storeIndex] = arr[right];
    arr[right] = temp;  // Move pivot to its final place
    return storeIndex;
}

double quickselect(double* arr, size_t left, size_t right, size_t k) {
    if (left == right) {
        return arr[left];  // If the list contains only one element, return that element
    }
    size_t pivotIndex = left + (right - left) / 2;
    pivotIndex = partition(arr, left, right, pivotIndex);
    if (k == pivotIndex) {
        return arr[k];
    } else if (k < pivotIndex) {
        return quickselect(arr, left, pivotIndex - 1, k);
    } else {
        return quickselect(arr, pivotIndex + 1, right, k);
    }
}

void add_data_point(RunningGradient *rg, double y) {
    double x[2] = {rg->num_points, 1}; // Setup data point vector x

    // Do recursive least squares computations
    // Calculate a scalar value 'temp' to ensure numerical stability
    double temp = 1 + (x[0] * rg->inverse_cov_matrix[0][0] + x[1] * rg->inverse_cov_matrix[1][0]) * x[0] + 
                      (x[0] * rg->inverse_cov_matrix[0][1] + x[1] * rg->inverse_cov_matrix[1][1]) * x[1];

    // Compute temporary vector 'tmp' to adjust the inverse covariance matrix
    double tmp[2] = {rg->inverse_cov_matrix[0][0] * x[0] + rg->inverse_cov_matrix[0][1] * x[1], 
                     rg->inverse_cov_matrix[1][0] * x[0] + rg->inverse_cov_matrix[1][1] * x[1]};

    // Update the inverse covariance matrix by removing old data influence and adding new data point influence
    rg->inverse_cov_matrix[0][0] -= (tmp[0] * tmp[0]) / temp;
    rg->inverse_cov_matrix[0][1] -= (tmp[0] * tmp[1]) / temp;
    rg->inverse_cov_matrix[1][0] -= (tmp[1] * tmp[0]) / temp;
    rg->inverse_cov_matrix[1][1] -= (tmp[1] * tmp[1]) / temp;

    // Ensure the inverse covariance matrix remains symmetric
    rg->inverse_cov_matrix[0][0] = 0.5 * (rg->inverse_cov_matrix[0][0] + rg->inverse_cov_matrix[0][1]);
    rg->inverse_cov_matrix[1][1] = 0.5 * (rg->inverse_cov_matrix[1][0] + rg->inverse_cov_matrix[1][1]);

    // Update the coefficients (slope and intercept) of the linear model using the new data point
    rg->coefficients[0] += (rg->inverse_cov_matrix[0][0] * x[0] + rg->inverse_cov_matrix[0][1] * x[1]) * 
                           (y - (x[0] * rg->coefficients[0] + x[1] * rg->coefficients[1]));
    rg->coefficients[1] += (rg->inverse_cov_matrix[1][0] * x[0] + rg->inverse_cov_matrix[1][1] * x[1]) * 
                           (y - (x[0] * rg->coefficients[0] + x[1] * rg->coefficients[1]));

    // Update the residual sum of squares
    rg->residual_sum_squares += pow((y - (x[0] * rg->coefficients[0] + x[1] * rg->coefficients[1])), 2.0) * temp;

    // Increment the count of data points processed
    rg->num_points += 1;
}

double calculate_gradient(const RunningGradient *rg) {
    assert(rg->num_points > 1 && "You must add more values into this object before calling this function.");
    return rg->coefficients[0];
}

double calculate_standard_error(const RunningGradient *rg) {
    assert(rg->num_points > 2 && "You must add more values into this object before calling this function.");
    double s = rg->residual_sum_squares / (rg->num_points - 2);
    double adjust = 12.0 / (pow(rg->num_points, 3.0) - rg->num_points);
    return sqrt(s * adjust);
}

static double normal_cdf(double value, double mean, double stddev) {
    if (stddev == 0) {
        if (value < mean)
            return 0;
        else if (value > mean)
            return 1;
        else
            return 0.5;
    }
    value = (value - mean) / stddev;
    return 0.5 * erfc(-value / sqrt(2.0));
}

double probability_gradient_less_than(const RunningGradient *rg, double thresh) {
    return normal_cdf(thresh, calculate_gradient(rg), calculate_standard_error(rg));
}

double probability_gradient_greater_than(const RunningGradient *rg, double thresh) {
    return 1.0 - probability_gradient_less_than(rg, thresh);
}

double find_upper_quantile(const double* data, size_t size, double quantile) {
    assert(0 <= quantile && quantile <= 1.0 && "quantile must be in the range [0, 1]");
    assert(size > 0 && "The container must contain more than 0 elements.");

    double* sorted_data = (double*)malloc(size * sizeof(double));
    if (!sorted_data) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }
    for (size_t i = 0; i < size; ++i) {
        sorted_data[i] = data[i];
    }

    size_t idx_upper = (size_t)round((size - 1) * (1 - quantile));
    double upper_q = quickselect(sorted_data, 0, size - 1, idx_upper);

    free(sorted_data);
    return upper_q;
}

void process_data_in_chunks(const double *data, size_t size) {
    size_t chunk_size = 10;

    for (size_t i = 0; i < size; i += chunk_size) {
        size_t end = (i + chunk_size < size) ? i + chunk_size : size;

        RunningGradient rg;
        init_running_gradient(&rg);

        for (size_t j = i; j < end; ++j) {
            add_data_point(&rg, data[j]);
        }

        double gradient = calculate_gradient(&rg);
        double std_error = calculate_standard_error(&rg);
        double prob_less_than_0 = probability_gradient_less_than(&rg, 0);
        double prob_greater_than_0 = probability_gradient_greater_than(&rg, 0);

        printf("Chunk %zu-%zu: Gradient = %f, Std Error = %f, P(Grad < 0) = %f, P(Grad > 0) = %f\n", 
                i, end - 1, gradient, std_error, prob_less_than_0, prob_greater_than_0);
    }
}
