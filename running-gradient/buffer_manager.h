/**
 * @file buffer_manager.h
 * @brief Header file for managing the buffer used in sliding window analysis.
 *
 * This header defines the structures, global variables, and function prototypes
 * used for managing the buffer in sliding window analysis of phase angle data.
 */

#ifndef BUFFER_MANAGER_H
#define BUFFER_MANAGER_H

#include <stdint.h>
#include "mqs_def.h"  // Assuming this contains the definition for MqsRawDataPoint_t
#include "running_quadratic_gradient.h"  // For QuadraticPeakAnalysisResult

/** @brief Define the buffer size for the sliding window analysis. */
#define BUFFER_SIZE 200
//#define BUFFER_DEBUG

// Global indices for tracking analysis window
extern int analysis_start_index;  // Track the initial start index
extern int analysis_end_index;    // Track the final end index

/**
 * @struct BufferManager
 * @brief Structure to manage the buffer and indices for sliding window analysis.
 *
 * This structure holds the buffer data and manages indices for the current phase,
 * buffer positions, and frequency mapping in the analysis.
 */
typedef struct {
    int16_t current_phase_index;   /**< Current index in the phaseAngles array. */
    int16_t current_buffer_index;  /**< Current index in the buffer array. */
    uint16_t buffer_size;          /**< Size of the buffer. */
    uint16_t window_size;          /**< Size of the sliding window (e.g., 30). */
    MqsRawDataPoint_t* buffer;     /**< Pointer to the buffer array. */
    int16_t mapping_start_index;   /**< Tracker variable for mapping start index. */
    int16_t mapping_end_index;     /**< Tracker variable for mapping end index. */
    double start_frequency;        /**< Starting frequency of the sweep (e.g., 11300 Hz). */
    double frequency_increment;    /**< Frequency increment per step (e.g., 1 Hz). */
} BufferManager;

/** @brief Extern declaration of the buffer manager for usage across different files. */
extern BufferManager buffer_manager;

/**
 * @brief Initializes the BufferManager structure with starting values.
 *
 * @param buffer              Pointer to the buffer that holds the data.
 * @param buffer_size         Size of the buffer.
 * @param window_size         Size of the sliding window.
 * @param start_index         Starting index for phase data analysis.
 * @param start_frequency     The starting frequency for the phase index.
 * @param frequency_increment Increment value to calculate frequency for each phase index.
 */
void init_buffer_manager(
    MqsRawDataPoint_t* buffer,
    uint16_t buffer_size,
    uint16_t window_size,
    int16_t start_index,
    double start_frequency,
    double frequency_increment
);

/**
 * @brief Loads the initial window of phase angle data into the buffer.
 *
 * @param phaseAngles      Pointer to the array of phase angles.
 * @param phase_angle_size Size of the phase angle array.
 */
void load_initial_buffer(const double* phaseAngles, uint16_t phase_angle_size);

/**
 * @brief Updates the buffer with phaseAngle values and optionally sets impedance.
 *
 * @param phaseAngles        Pointer to the array of phase angles.
 * @param phase_index_start  The starting index of the phase angles to update from.
 * @param buffer_start_index The starting index of the buffer where the phase angles will be stored.
 */
void update_phaseAngle_to_buffer(
    const double* phaseAngles,
    int16_t phase_index_start,
    int16_t buffer_start_index
);

/**
 * @brief Updates the buffer to shift the analysis window in a specified direction.
 *
 * @param phaseAngles Pointer to the array of phase angles used in the analysis.
 * @param direction   Direction to move the buffer window. Accepts `LEFT_SIDE` or `RIGHT_SIDE`.
 * @param move_amount The number of positions to move the buffer window in the specified direction.
 */
void update_buffer_for_direction(const double* phaseAngles, int direction, int move_amount);

/**
 * @brief Handles the case when the segment analysis result is undecided.
 *
 * @param phaseAngles      Pointer to the array of phase angles.
 * @param phase_angle_size The total size of the phase angle array.
 */
void handle_undecided_case(const double* phaseAngles, uint16_t phase_angle_size);

void move_window_and_update_if_needed(const double* phaseAngles, int direction, int move_amount);

void print_analysis_interval(const double* phaseAngles, int total_size);


#endif // BUFFER_MANAGER_H