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

// Buffer manager struct to handle sliding window analysis
typedef struct {
    int16_t current_phase_index;  // Tracks current index in the phaseAngles array
    int16_t current_buffer_index; // Tracks where we are in the MqsRawDataPoint_t buffer
    uint16_t buffer_size;         // Size of the buffer
    uint16_t window_size;         // Size of the sliding window (e.g., 30)
    MqsRawDataPoint_t* buffer;    // Pointer to the buffer
    int16_t mapping_start_index;  // Tracker variable for mapping start index
    int16_t mapping_end_index;    // Tracker variable for mapping end index
    double start_frequency;       // Starting frequency of the sweep (e.g., 11300 Hz)
    double frequency_increment;   // Frequency increment per step (e.g., 1 Hz)
} BufferManager;

// Extern declaration of the buffer manager for usage across different translation units
extern BufferManager buffer_manager;

/**
 * @brief Initializes the buffer manager with the provided parameters.
 * 
 * @param buffer Pointer to the buffer that holds the data
 * @param buffer_size Size of the buffer
 * @param window_size Size of the sliding window
 * @param start_index Starting index for the buffer manager
 * @param start_frequency Starting frequency in Hz
 * @param frequency_increment Frequency increment per step in Hz
 */
void init_buffer_manager(MqsRawDataPoint_t* buffer, uint16_t buffer_size, uint16_t window_size, int16_t start_index, double start_frequency, double frequency_increment);

/**
 * @brief Analyzes the trend and concavity within a data segment to determine the direction of movement.
 * 
 * @param data Array of MqsRawDataPoint_t values representing the data points
 * @param size Total number of data points in the array
 * @param start_index Starting index within the data array where analysis should begin
 * @param forgetting_factor Forgetting factor that determines the weight of newer data points
 * @return SegmentAnalysisResult Structure containing information about the detected trends and concavity pattern
 */
SegmentAnalysisResult segment_trend_and_concavity_analysis(const MqsRawDataPoint_t *data, uint16_t size, uint16_t start_index, double forgetting_factor);

/**
 * @brief Performs sliding window analysis with buffer management.
 * 
 * @param phaseAngles Pointer to the array of phase angles
 * @param phase_angle_size Size of the phase angle array
 */
void perform_sliding_window_analysis(const double* phaseAngles, uint16_t phase_angle_size);


#endif // RUNNING_PEAK_ANALYSIS_H