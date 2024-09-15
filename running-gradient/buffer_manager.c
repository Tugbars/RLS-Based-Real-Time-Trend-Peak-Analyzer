/**
 * @file buffer_manager.c
 * @brief Implementation of buffer management functions for sliding window analysis.
 *
 * This file contains the implementations of functions declared in `buffer_manager.h`,
 * including buffer initialization, updates, and management during sliding window analysis.
 */

#include <stdio.h>
#include <stdint.h>
#include "running_quadratic_gradient.h"
#include "rls_analysis_parameters.h"

#include "buffer_manager.h"
#include "windowed_running_gradient.h"

// Global variables
BufferManager buffer_manager;

/** @brief Define the buffer array used for sliding window analysis. */
MqsRawDataPoint_t buffer[BUFFER_SIZE];

/**
 * @brief Initializes the BufferManager struct with starting values.
 *
 * This function sets up the buffer manager by initializing its buffer, buffer size, window size, 
 * current phase index, and frequency parameters. The buffer index is set to the middle of the buffer.
 *
 * @param buffer A pointer to the phase data buffer to be managed.
 * @param buffer_size Size of the buffer.
 * @param window_size Size of the sliding window.
 * @param start_index Starting index for phase data analysis.
 * @param start_frequency The starting frequency for the phase index.
 * @param frequency_increment Increment value to calculate frequency for each phase index.
 *
 * @note This function prints out debugging information, including the initial state of the buffer manager.
 */
void init_buffer_manager(
    MqsRawDataPoint_t* buffer_ptr,
    uint16_t buffer_size,
    uint16_t window_size,
    int16_t start_index,
    double start_frequency,
    double frequency_increment
) {
    buffer_manager.buffer = buffer_ptr;
    buffer_manager.buffer_size = buffer_size;
    buffer_manager.window_size = window_size;
    buffer_manager.current_phase_index = start_index;
    buffer_manager.current_buffer_index = buffer_size / 2; // Start from the middle of the buffer

    buffer_manager.mapping_start_index = start_index;
    buffer_manager.mapping_end_index = start_index + window_size - 1;

    buffer_manager.start_frequency = start_frequency;
    buffer_manager.frequency_increment = frequency_increment;

    // Debugging: Print initialization state
    #ifdef BUFFER_DEBUG
    printf("Buffer Manager initialized:\n");
    printf("Buffer size: %d, Window size: %d, Start index: %d\n", buffer_size, window_size, start_index);
    printf("Start frequency: %.2f Hz, Frequency increment: %.2f Hz\n", start_frequency, frequency_increment);
    #endif
}

/**
 * @brief Loads the initial window of phase angle data into the buffer.
 *
 * @param phaseAngles      Pointer to the array of phase angles.
 * @param phase_angle_size Size of the phase angle array.
 */
void load_initial_buffer(const double* phaseAngles, uint16_t phase_angle_size) {
    #ifdef BUFFER_DEBUG
    printf("Loading initial window of data into the buffer.\n");
    #endif

    for (uint16_t i = 0; i < buffer_manager.window_size; ++i) {
        int16_t buffer_index = (buffer_manager.current_buffer_index + i) % buffer_manager.buffer_size;
        int16_t phase_index = buffer_manager.current_phase_index + i;

        // Ensure phase_index is within bounds
        if (phase_index >= phase_angle_size || phase_index < 0) {
            //printf("[load_initial_buffer] Phase index %d out of bounds. Setting phaseAngle to 0.0.\n", phase_index);
            buffer_manager.buffer[buffer_index].phaseAngle = 0.0;
        } else {
            buffer_manager.buffer[buffer_index].phaseAngle = phaseAngles[phase_index];
            buffer_manager.buffer[buffer_index].impedance = 0.0; // Initialize impedance as needed
            //printf("[load_initial_buffer] Buffer[%d] set to phaseAngle %.6f\n", buffer_index, phaseAngles[phase_index]);
        }
    }

    #ifdef BUFFER_DEBUG
    printf("Initial buffer loading complete.\n");
    #endif
}

/**
 * @brief Updates the buffer with phaseAngle values and optionally sets impedance.
 *
 * This function maps phase index values to buffer indices and updates the buffer with the corresponding
 * phaseAngle values. It also logs the updates made to the buffer.
 *
 * @param phaseAngles A pointer to the array of phase angles.
 * @param phase_index_start The starting index of the phase angles to update from.
 * @param buffer_start_index The starting index of the buffer where the phase angles will be stored.
 */
void update_phaseAngle_to_buffer(const double* phaseAngles, int16_t phase_index_start, int16_t buffer_start_index) {
    // Loop to update the buffer with new phaseAngle values
    for (uint16_t i = 0; i < buffer_manager.window_size; ++i) {
        int16_t phase_index = phase_index_start + i;
        int16_t buffer_index = (buffer_start_index + i) % buffer_manager.buffer_size;

        // Debugging: Log the details of buffer update before actual update
        //printf("[update_phaseAngle_to_buffer] [DEBUG] Before update -> Buffer index: %d, Phase index: %d, PhaseAngle: %.6f\n", buffer_index, phase_index, phaseAngles[phase_index]);

        // Copy phaseAngle and optionally set impedance
        buffer_manager.buffer[buffer_index].phaseAngle = phaseAngles[phase_index];
        buffer_manager.buffer[buffer_index].impedance = 0.0; // Mock impedance, update if needed

        // Debugging: Log the update in the buffer
        //printf("[update_phaseAngle_to_buffer] Buffer updated at phase index %d with phaseAngle %.6f (buffer index: %d)\n", phase_index, phaseAngles[phase_index], buffer_index);
    }
}

/**
 * @brief Updates the buffer to shift the analysis window in a specified direction.
 *
 * This function adjusts the buffer indices and updates the phase angle data within the buffer to move the analysis window
 * left or right by a specified amount. It handles circular buffer wrap-around and ensures that the phase indices remain
 * within the bounds of the phase angle array. The function is crucial for centering the peak within the analysis window
 * during sliding window analysis.
 *
 * @param phaseAngles Pointer to the array of phase angles used in the analysis.
 * @param direction   Direction to move the buffer window. Accepts `LEFT_SIDE` or `RIGHT_SIDE` as valid inputs.
 * @param move_amount The number of positions to move the buffer window in the specified direction.
 *
 * ### Function Workflow:
 * 1. **Initialize Indices**:
 *    - `phase_index_start`: Starting index of the current phase in the phase angle array.
 *    - `phase_index_end`: Ending index of the current phase window in the phase angle array.
 *    - `buffer_start_index`: Starting index in the buffer where the phase angles are stored.
 *    - Logs the initial indices before movement.
 * 2. **Calculate New Indices Based on Direction**:
 *    - If moving `LEFT_SIDE` (towards lower indices):
 *      - Subtracts `move_amount` from `phase_index_start` and `buffer_start_index`.
 *      - Updates `phase_index_end` accordingly.
 *      - Handles circular buffer wrap-around for `buffer_start_index`.
 *      - Logs the movement and new indices.
 *    - If moving `RIGHT_SIDE` (towards higher indices):
 *      - Adds `move_amount` to `phase_index_start` and `buffer_start_index`.
 *      - Updates `phase_index_end` accordingly.
 *      - Handles circular buffer wrap-around for `buffer_start_index`.
 *      - Logs the movement and new indices.
 *    - If direction is `UNDECIDED`, logs that no movement is performed and returns.
 * 3. **Bounds Checking**:
 *    - Ensures that `phase_index_start` is within the bounds of `phaseAngles`.
 *    - If out of bounds, logs an error message and aborts the operation.
 * 4. **Update Buffer with New Phase Angles**:
 *    - Calls `update_phaseAngle_to_buffer` to load the new phase angles into the buffer starting at `buffer_start_index`.
 * 5. **Update Buffer Manager Indices**:
 *    - Updates `buffer_manager.current_phase_index` with the new `phase_index_start`.
 *    - Updates `buffer_manager.current_buffer_index` with the new `buffer_start_index`.
 *    - Logs the updated indices after the operation.
 *
 * ### Important Notes:
 * - The function handles circular buffer wrap-around for `buffer_start_index` to ensure it stays within the buffer size.
 * - If `phase_index_start` goes out of bounds (negative or beyond the phase angle array size), the function aborts and logs an error.
 * - The buffer is updated in place, and the buffer manager's indices are adjusted to reflect the new positions.
 *
 *
 * @see update_phaseAngle_to_buffer
 * @see buffer_manager
 */
void update_buffer_for_direction(const double* phaseAngles, int direction, int move_amount) {
    int16_t phase_index_start = buffer_manager.current_phase_index;
    int16_t phase_index_end = phase_index_start + buffer_manager.window_size;
    int16_t buffer_start_index = buffer_manager.current_buffer_index;

    // Log initial indices before movement
    #ifdef BUFFER_DEBUG
    printf("[update_buffer] Initial phase index start: %d, end: %d, buffer start index: %d\n",
           phase_index_start, phase_index_end, buffer_start_index);
    #endif

    // Calculate the shift direction (LEFT or RIGHT)
    if (direction == LEFT_SIDE) {
        phase_index_start -= move_amount;
        phase_index_end = phase_index_start + buffer_manager.window_size;
        buffer_start_index -= move_amount;

        // Handle circular buffer wrap-around for buffer index
        if (buffer_start_index < 0) {
            buffer_start_index += buffer_manager.buffer_size;
        }
        #ifdef BUFFER_DEBUG
        printf("[update_buffer] Moving left by %d. New phase index range: [%d - %d], buffer start index: %d\n",
               move_amount, phase_index_start, phase_index_end, buffer_start_index);
        #endif
    } else if (direction == RIGHT_SIDE) {
        phase_index_start += move_amount;
        phase_index_end = phase_index_start + buffer_manager.window_size;
        buffer_start_index += move_amount;

        // Handle circular buffer wrap-around for buffer index
        if (buffer_start_index >= buffer_manager.buffer_size) {
            buffer_start_index -= buffer_manager.buffer_size;
        }
        #ifdef BUFFER_DEBUG
        printf("[update_buffer] Moving right by %d. New phase index range: [%d - %d], buffer start index: %d\n",
               move_amount, phase_index_start, phase_index_end, buffer_start_index);
        #endif
    } else {
        #ifdef BUFFER_DEBUG
        printf("[update_buffer] Direction undecided. No movement.\n");
        #endif
        return;
    }

    // Ensure phase index stays within bounds of phaseAngles
    if (phase_index_start < 0 || buffer_start_index > buffer_manager.buffer_size) {
        #ifdef BUFFER_DEBUG
        printf("[update_buffer] Out of bounds of phase angles (index: [%d - %d]). Operation aborted.\n",
               phase_index_start, phase_index_end);
        #endif
        return;
    }

    // Update the buffer with new phase angles
    update_phaseAngle_to_buffer(phaseAngles, phase_index_start, buffer_start_index);

    // Update the current phase and buffer indices
    buffer_manager.current_phase_index = phase_index_start;
    buffer_manager.current_buffer_index = buffer_start_index;

    // Log final indices after update
    //#ifdef BUFFER_DEBUG
    printf("[update_buffer] Updated buffer indices. Current phase index: %d, Buffer start index: %d\n",
           buffer_manager.current_phase_index, buffer_manager.current_buffer_index);
    //#endif
}


/**
 * @brief Handles the case when the segment analysis result is undecided.
 *
 * This function jumps forward in the phaseAngle array by the window size, discards the previous analysis,
 * and reloads the buffer from the new phase index. The new portion of the phaseAngle array is loaded into the buffer,
 * overwriting the previous values in the middle section of the buffer.
 *
 * @param phaseAngles A pointer to the array of phase angles.
 * @param phase_angle_size The total size of the phase angle array.
 */
void handle_undecided_case(const double* phaseAngles, uint16_t phase_angle_size) {
    // Flag to indicate entry into the function
    //printf("[handle_undecided_case] Entering undecided case handler.\n");

    // Jump forward by the window size in the phaseAngle array
    buffer_manager.current_phase_index += buffer_manager.window_size;

    // Ensure we don't exceed the phase angle size
    if (buffer_manager.current_phase_index + buffer_manager.window_size > phase_angle_size) {
        #ifdef BUFFER_DEBUG
        printf("[handle_undecided_case] Phase index exceeds phase angle size. Stopping analysis.\n");
        #endif
        return;
    }

    // Load new phase angle values from the phaseAngle array into the middle section of the buffer
    #ifdef BUFFER_DEBUG
    printf("[handle_undecided_case] Jumping to phase index: %d\n", buffer_manager.current_phase_index);
    #endif
    update_phaseAngle_to_buffer(phaseAngles, buffer_manager.current_phase_index, buffer_manager.current_buffer_index);

    // Log the buffer after loading new data
    #ifdef BUFFER_DEBUG
    printf("[handle_undecided_case] Loaded new phase angle values into the buffer. Starting analysis from buffer index %d.\n",
           buffer_manager.current_buffer_index);
    #endif
}

