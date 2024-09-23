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
int analysis_start_index = -1;  // Track the initial start index
int analysis_end_index = -1;    // Track the final end index
int buffer_shift_offset = 0;  // New global variable to track the shift relative to the buffer.


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
/**
 * @brief Initializes the BufferManager struct with starting values.
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

    // Initialize the analysis tracking variables
    if (analysis_start_index == -1) {  // Set the initial start index only once at the beginning
        analysis_start_index = start_index;
    }

    // Set the initial end index (will be updated dynamically)
    analysis_end_index = start_index + window_size - 1;
    buffer_shift_offset = 0;

    #ifdef BUFFER_DEBUG
    printf("Buffer Manager initialized:\n");
    printf("Buffer size: %d, Window size: %d, Start index: %d\n", buffer_size, window_size, start_index);
    printf("Start frequency: %.2f Hz, Frequency increment: %.2f Hz\n", start_frequency, frequency_increment);
    #endif
}

/**
 * @brief Loads the initial window of phase angle data into the buffer.
 *
 * This function updates both the start and end indices globally to reflect the phase angles loaded into the buffer.
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
            buffer_manager.buffer[buffer_index].phaseAngle = 0.0;
        } else {
            buffer_manager.buffer[buffer_index].phaseAngle = phaseAngles[phase_index];
            buffer_manager.buffer[buffer_index].impedance = 0.0; // Initialize impedance as needed
        }
    }

    // Update the global analysis indices
    if (analysis_start_index == -1) { // Set the start index only if it hasn't been set yet
        analysis_start_index = buffer_manager.current_phase_index;
    }

    // The end index is always the start index + window size - 1 after loading the initial window
    analysis_end_index = buffer_manager.current_phase_index + buffer_manager.window_size - 1;

    #ifdef BUFFER_DEBUG
    printf("Initial buffer loading complete.\n");
    printf("Analysis start index set to %d, end index set to %d\n", analysis_start_index, analysis_end_index);
    #endif
}

//move_amount mantığı yerleştirilmeli. 
//bir önceki phase Shit'e göre farkından move amount ölçülmeli. 
//bu sadece move amountu
// update_phaseAngle_toBuffer starts from the moved phaseIndex and updates to the buffer which has its index is also moved by the move_amount in the functipn. 


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
void update_phaseAngle_to_buffer(const double* phaseAngles, int16_t phase_index_start, int16_t buffer_start_index) { //  MesSweepAddDataPoint(currentRawSweep, real, imaginary); ile aynı olmalı.
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
        buffer_start_index -= move_amount;

        // Handle circular buffer wrap-around for buffer index
        if (buffer_start_index < 0) {
            buffer_start_index += buffer_manager.buffer_size;
        }

        // Only update the start index when sliding left
        if (phase_index_start < analysis_start_index) {
            analysis_start_index = phase_index_start;  // Move start index leftwards
        }

        #ifdef BUFFER_DEBUG
        printf("[update_buffer] Moving left by %d. New phase index range: [%d - %d], buffer start index: %d\n",
               move_amount, phase_index_start, phase_index_end, buffer_start_index);
        #endif
    } else if (direction == RIGHT_SIDE) {
        phase_index_start += move_amount;
        buffer_start_index += move_amount;

        // Handle circular buffer wrap-around for buffer index
        if (buffer_start_index >= buffer_manager.buffer_size) {
            buffer_start_index -= buffer_manager.buffer_size;
        }

        // Only update the end index when sliding right
        if ((phase_index_start + buffer_manager.window_size - 1) > analysis_end_index) {
            analysis_end_index = phase_index_start + buffer_manager.window_size - 1;  // Move end index rightwards
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
    
    
    // Update the current phase and buffer indices
    buffer_manager.current_phase_index = phase_index_start;
    buffer_manager.current_buffer_index = buffer_start_index;
    
    // Log final indices after update
    //#ifdef BUFFER_DEBUG
    printf("[update_buffer] Updated buffer indices. Current phase index: %d, Buffer start index: %d\n",
           buffer_manager.current_phase_index, buffer_manager.current_buffer_index);
    //#endif

    // Update the buffer with new phase angles
    update_phaseAngle_to_buffer(phaseAngles, phase_index_start, buffer_start_index);
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

    // Update the global analysis tracking indices
    // Since we're jumping forward, both start and end indices need to be updated
    analysis_start_index = buffer_manager.current_phase_index;
    analysis_end_index = buffer_manager.current_phase_index + buffer_manager.window_size - 1;

    // Load new phase angle values from the phaseAngle array into the middle section of the buffer
    #ifdef BUFFER_DEBUG
    printf("[handle_undecided_case] Jumping to phase index: %d\n", buffer_manager.current_phase_index);
    #endif
    
        // Log the buffer after loading new data
    #ifdef BUFFER_DEBUG
    printf("[handle_undecided_case] Loaded new phase angle values into the buffer. Starting analysis from buffer index %d.\n",
           buffer_manager.current_buffer_index);
    printf("Updated analysis interval: Start index = %d, End index = %d\n", analysis_start_index, analysis_end_index);
    #endif
    
    update_phaseAngle_to_buffer(phaseAngles, buffer_manager.current_phase_index, buffer_manager.current_buffer_index);
}

/**
 * @brief Moves the analysis window by a specified number of positions left or right and updates the buffer only if necessary.
 *
 * This function moves the current window either left or right by the specified `move_amount` positions and checks whether the values
 * for the new window are already available in the buffer by comparing the `analysis_start_index` and 
 * `analysis_end_index` against the new window range. If the values are already present, the window is simply 
 * moved without reloading data from `phaseAngles`. Otherwise, the missing data is loaded into the buffer.
 *
 * @param phaseAngles Pointer to the array of phase angles used in the analysis.
 * @param direction   Direction to move the window. Either `LEFT_SIDE` or `RIGHT_SIDE`.
 * @param move_amount Number of positions to move the window in the specified direction.
 */
void move_window_and_update_if_needed(const double* phaseAngles, int direction, int move_amount) {
    int16_t new_phase_index_start, new_phase_index_end;
    int16_t buffer_start_index = buffer_manager.current_buffer_index;

    // Calculate the new start and end indices based on the direction
    if (direction == LEFT_SIDE) {
        new_phase_index_start = buffer_manager.current_phase_index - move_amount;
        new_phase_index_end = new_phase_index_start + buffer_manager.window_size - 1;

        // Handle buffer wrap-around if moving left
        if (buffer_start_index - move_amount < 0) {
            buffer_start_index = buffer_manager.buffer_size + (buffer_start_index - move_amount);
        } else {
            buffer_start_index -= move_amount;
        }
    } else if (direction == RIGHT_SIDE) {
        new_phase_index_start = buffer_manager.current_phase_index + move_amount;
        new_phase_index_end = new_phase_index_start + buffer_manager.window_size - 1;

        // Handle buffer wrap-around if moving right
        buffer_start_index = (buffer_start_index + move_amount) % buffer_manager.buffer_size;

        // Only update the analysis_end_index if the new end index exceeds the current one
        if (new_phase_index_end > analysis_end_index) {
            analysis_end_index = new_phase_index_end;
        }
    } else {
        return;  // Invalid direction, no action
    }

    // Check if the new window is already covered by the global tracking indices
    if (new_phase_index_start >= analysis_start_index && new_phase_index_end <= analysis_end_index) {
        // Track how much we have shifted the buffer relative to the original buffer start
        buffer_shift_offset += move_amount;  // Track the shift without buffer update
        buffer_shift_offset %= buffer_manager.buffer_size;  // Wrap it within buffer bounds

        printf("[move_window_and_update_if_needed] Window already in buffer. Moving without update.\n");
    } else {
        // Reset the shift offset if we are reloading the buffer
        buffer_shift_offset = 0;

        // The new window needs to load values from `phaseAngles`
        printf("[move_window_and_update_if_needed] Loading new window values into the buffer.\n");
        update_buffer_for_direction(phaseAngles, direction, move_amount);

        // Update the global tracking indices
        analysis_start_index = new_phase_index_start;
        if (new_phase_index_end > analysis_end_index) {
            analysis_end_index = new_phase_index_end;
        }
    }

    // Update the buffer manager with the new indices
    buffer_manager.current_phase_index = new_phase_index_start;
    buffer_manager.current_buffer_index = buffer_start_index;

    // Debug: Print the updated indices
    printf("[move_window_and_update_if_needed] Updated window to phase index [%d - %d]\n", 
           new_phase_index_start, new_phase_index_end);
    printf("Updated buffer index to %d\n", buffer_manager.current_buffer_index);
}

/**
 * @brief Prints the values within the tracked interval of phase indices from the beginning to the end.
 *
 * This function will be called once the state machine has finished executing to print the values
 * within the sliding window interval of phase indices from `analysis_start_index` to `analysis_end_index`.
 */
void print_analysis_interval(const double* phaseAngles, int total_size) {
    if (analysis_start_index == -1 || analysis_end_index == -1) {
        printf("No analysis interval found.\n");
        return;
    }

    // Ensure the interval is within valid range
    if (analysis_start_index < 0 || analysis_end_index >= total_size) {
        printf("Invalid analysis interval. Start or end index out of bounds.\n");
        return;
    }

    printf("\n[RESULT] Sliding Window Analysis Interval:\n");
    printf("Analysis start index: %d, Analysis end index: %d\n", analysis_start_index, analysis_end_index);

    printf("\n[RESULT] Buffer Values in the Interval:\n");

    // Ensure the buffer index is correctly calculated for circular wrapping
    for (int phase_index = analysis_start_index; phase_index <= analysis_end_index && phase_index < total_size; ++phase_index) {
        // Calculate the corresponding buffer index relative to the window position
        int relative_position_in_window = phase_index - analysis_start_index;

        // Adjust buffer index by the shift offset to account for window movement without update
        int adjusted_buffer_index = (buffer_manager.current_buffer_index + relative_position_in_window - buffer_shift_offset) % buffer_manager.buffer_size;

        // Handle negative buffer index wrap-around
        if (adjusted_buffer_index < 0) {
            adjusted_buffer_index += buffer_manager.buffer_size;
        }

        // Print the buffer values
        printf("Phase Index %d -> Buffer Index %d: Phase Angle: %.6f, Impedance: %.6f\n",
               phase_index, adjusted_buffer_index,
               buffer_manager.buffer[adjusted_buffer_index].phaseAngle,
               buffer_manager.buffer[adjusted_buffer_index].impedance);
    }

    printf("--------------------------------------\n");
}
