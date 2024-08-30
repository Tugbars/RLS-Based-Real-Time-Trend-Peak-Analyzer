// running_peak_analysis.h

#ifndef RUNNING_PEAK_ANALYSIS_H
#define RUNNING_PEAK_ANALYSIS_H

#include <stddef.h> // For size_t
#include "running_quadratic_gradient.h"

typedef struct {
    bool isPotentialPeak;
    bool isTruePeak;
    PeakPosition nextDirection;
    ConcavityPattern pattern;
} SegmentAnalysisResult;

typedef struct {
    bool move_to_right;
    bool go_to_left;
    bool close_to_peak;  
    bool far_to_peak;    
    bool low_consistent_right; 
    bool low_consistent_left;  
    bool on_the_peak;   // New flag indicating the trend is exactly on the peak
} TrendDirectionFlags;

// Function to compare gradient parts, find direction, and then check initial concavity
SegmentAnalysisResult segment_trend_and_concavity_analysis(const double *data, size_t size, size_t start_index, double forgetting_factor);

#endif // RUNNING_PEAK_ANALYSIS_H
