#ifndef RUNNING_PEAK_ANALYSIS_H
#define RUNNING_PEAK_ANALYSIS_H

#include <stdint.h> // For uint16_t
#include "running_quadratic_gradient.h"
#include "mqs_def.h" // Assuming this is where MqsRawDataPoint_t is defined

// Struct to hold the result of segment analysis, including concavity analysis
typedef struct {
    bool isPotentialPeak;
    bool isTruePeak;
    PeakPosition nextDirection;
    ConcavityAnalysisOutput concavityOutput; // Updated to use ConcavityAnalysisOutput instead of ConcavityPattern
} SegmentAnalysisResult;

typedef struct {
    bool move_to_right;
    bool go_to_left;
    bool close_to_peak;  
    bool far_to_peak;    
    bool low_consistent_right; 
    bool low_consistent_left;  
    bool on_the_peak;   // Flag indicating the trend is exactly on the peak
} TrendDirectionFlags;

/**
 * @brief Analyzes the trend and concavity within a data segment to determine the direction of movement.
 * 
 * @param data Array of MqsRawDataPoint_t values representing the data points
 * @param size Total number of data points in the array
 * @param start_index Starting index within the data array where analysis should begin
 * @param forgetting_factor Forgetting factor that determines the weight of newer data points
 * @return SegmentAnalysisResult Structure containing information about the detected trends and concavity pattern
 */
SegmentAnalysisResult segment_trend_and_concavity_analysis(const MqsRawDataPoint_t *data, uint16_t window_size, double forgetting_factor);

void verify_peak_at_index(uint16_t buffer_start_index);

#endif // RUNNING_PEAK_ANALYSIS_H