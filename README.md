# Running Gradient Calculation

## Overview

This implementation provides functions to perform running gradient calculations using a recursive least squares algorithm. This method allows for updating regression parameters as new data points are added without the need to reprocess the entire dataset. The primary focus is on maintaining the current state of the regression calculation, including the inverse covariance matrix and the regression coefficients.

## Why Use Recursive Least Squares (RLS)?

Recursive Least Squares (RLS) is a powerful method used for online parameter estimation in linear regression. Unlike the ordinary least squares (OLS) method, which requires the entire dataset to be available for processing, RLS updates the estimates of the regression parameters incrementally as new data points are received. This makes RLS particularly useful for real-time applications and systems where data is collected sequentially over time.

### Advantages of RLS

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

## Windowed Running Gradient

### Overview

To address the limitations of the standard RLS method, a windowed approach to performing running gradient calculations has been introduced. This approach limits the number of data points considered at any given time, preserving local features and improving peak detection in noisy data.

### Why Use a Windowed Running Gradient?

The windowed running gradient approach is designed to maintain the benefits of RLS while addressing its limitations in preserving local features. By focusing on a fixed number of recent data points, the windowed approach maintains local features, making it better suited for detecting sudden changes.

#### Advantages of the Windowed Approach

- **Local Feature Preservation**: By focusing on a fixed number of recent data points, the windowed approach maintains local features, making it better suited for detecting sudden changes.
- **Noise Reduction**: Helps in reducing the impact of older data, which may introduce noise if the underlying data distribution changes over time.
- **Adaptability**: More adaptable to changes in the data pattern, which is crucial for applications like financial data analysis and real-time monitoring systems.

### How It Uses a Windowed RLS

The windowed RLS algorithm maintains and updates the following key components within a fixed-size window:
- **Data Arrays (`x` and `y`)**: Store the most recent data points.
- **Regression Coefficients**: Represent the slope and intercept of the fitted line for the current window.
- **Residual Sum of Squares**: Tracks the sum of squared residuals for the data points within the window.

#### Process for Each New Data Point

1. **Data Point Addition**:
   - Add the new data point to the window. If the window is full, remove the oldest data point to make room for the new one.
   - Represent the data point as a vector `x`, where `x[0]` is the index within the window and `x[1]` is always 1.

2. **Intermediate Calculations**:
   - **Sum Variables**: Compute sums required for the least squares calculations within the window.
   - **Denominator Calculation**: Ensure numerical stability and adjust the regression coefficients.

3. **Update Regression Coefficients**:
   - Recompute the slope and intercept using only the data points within the current window.
   - Use the computed sums to efficiently update the coefficients.

4. **Update Residual Sum of Squares**:
   - Calculate the contribution of the new data point and remove the contribution of the oldest data point if the window is full.

5. **Increment Data Point Count**:
   - Track the number of data points processed within the current window.

## Summary

### Running Gradient

- **Method**: Recursive Least Squares (RLS)
- **Advantages**: Efficient, real-time processing, low memory usage
- **Limitations**: Over-smoothing with large datasets, potential masking of local features

### Windowed Running Gradient

- **Method**: Windowed Recursive Least Squares (Windowed RLS)
- **Advantages**: Preserves local features, reduces noise impact from older data, adaptable to data pattern changes
- **Limitations**: Fixed window size, potential for less stable long-term trends

This approach ensures you can choose the appropriate method based on your specific application requirements, balancing between long-term trend analysis and local feature preservation.
