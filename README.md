# Running Gradient Calculation

## Overview

This implementation provides functions to perform running (online) gradient calculations using a recursive least squares (RLS) algorithm. This method allows for updating regression parameters as new data points are added without the need to reprocess the entire dataset. The primary focus is on maintaining the current state of the regression calculation, including the inverse covariance matrix and the regression coefficients.

## Why Use Recursive Least Squares (RLS)?

Recursive Least Squares (RLS) is a powerful method used for online parameter estimation in linear regression. Unlike the ordinary least squares (OLS) method, which requires the entire dataset to be available for processing, RLS updates the estimates of the regression parameters incrementally as new data points are received. This makes RLS particularly useful for real-time applications and systems where data is collected sequentially over time.

### Advantages of RLS:
- **Efficiency**: RLS updates the regression parameters in constant time, making it computationally efficient for large datasets.
- **Real-time Processing**: Ideal for scenarios where data arrives in a stream, such as sensor data or financial time series.
- **Memory Usage**: Requires less memory as it does not need to store the entire dataset.

## How It Uses RLS

The RLS algorithm maintains and updates the following key components:
- **Inverse Covariance Matrix**: Used to update the parameter estimates efficiently.
- **Regression Coefficients**: Represents the slope and intercept of the fitted line.
- **Residual Sum of Squares**: Keeps track of the sum of squared residuals for calculating standard error.

### Process for Each New Data Point:

1. **Data Point Representation**: The new data point is represented as a vector `x` where `x[0]` is the index of the data point (time step) and `x[1]` is always 1 (for the intercept term).

2. **Intermediate Calculations**:
   - **Temporary Scalar (`temp`)**: Ensures numerical stability and is used to adjust the inverse covariance matrix.
   - **Temporary Vector (`tmp`)**: Helps in updating the inverse covariance matrix.

3. **Update Inverse Covariance Matrix**:
   - Adjust the matrix to account for the new data point.
   - Ensure the matrix remains symmetric.

4. **Update Regression Coefficients**:
   - Update the slope and intercept based on the new data point.
   - Use the adjusted inverse covariance matrix to refine the estimates.

5. **Update Residual Sum of Squares**:
   - Calculate the contribution of the new data point to the residual sum of squares.

6. **Increment Data Point Count**:
   - Keep track of the number of data points processed.

## Components and Functionality

### RunningGradient Structure

The `RunningGradient` structure holds the state for the running gradient calculation. It includes:
- `num_points`: The number of data points processed.
- `inverse_cov_matrix`: The inverse of the covariance matrix (2x2 matrix).
- `coefficients`: The regression coefficients (slope and intercept).
- `residual_sum_squares`: The sum of the squared residuals.
