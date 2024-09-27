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
#include "sliding_window_analysis.h"

// Global variables
BufferUpdateInfo buffer_update_info = {false, 0, 0, 0};

BufferManager buffer_manager;
int analysis_start_index = -1;  // Track the initial start index
int analysis_end_index = -1;    // Track the final end index
int buffer_shift_offset = 0;  // New global variable to track the shift relative to the buffer.

void AdptSweepAddDataPoint(const double* phaseAngles, int16_t phase_index_start, int16_t buffer_start_index, uint16_t move_amount);


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
    buffer_manager.current_phase_index = start_index;                                                 //STARTING FREQUENCY - MY OWN MODIFIER.
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
    //printf("Buffer Manager initialized:\n");
    //printf("Buffer size: %d, Window size: %d, Start index: %d\n", buffer_size, window_size, start_index);
    //printf("Start frequency: %.2f Hz, Frequency increment: %.2f Hz\n", start_frequency, frequency_increment);
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

    // Use AdptSweepAddDataPoint to load the initial window of data
    AdptSweepAddDataPoint(phaseAngles, buffer_manager.current_phase_index, buffer_manager.current_buffer_index, buffer_manager.window_size);  // Pass the full window size as move_amount

    // Ensure that phase indices are within bounds after updating
    for (uint16_t i = 0; i < buffer_manager.window_size; ++i) {
        int16_t phase_index = buffer_manager.current_phase_index + i;
        if (phase_index >= phase_angle_size || phase_index < 0) {
            int16_t buffer_index = (buffer_manager.current_buffer_index + i) % buffer_manager.buffer_size;
            buffer_manager.buffer[buffer_index].phaseAngle = 0.0; // Set out-of-bounds phase angle to 0.0
        }
    }

    // Update the global analysis indices
    if (analysis_start_index == -1) { // Set the start index only if it hasn't been set yet
        analysis_start_index = buffer_manager.current_phase_index;
    }

    // The end index is always the start index + window size - 1 after loading the initial window
    analysis_end_index = buffer_manager.current_phase_index + buffer_manager.window_size - 1;

   // #ifdef BUFFER_DEBUG
    //printf("Initial buffer loading complete.\n");
    printf("[BUFFER] Analysis start index set to %d, end index set to %d\n", analysis_start_index, analysis_end_index);
    //#endif
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
void AdptSweepAddDataPoint(const double* phaseAngles, int16_t phase_index_start, int16_t buffer_start_index, uint16_t move_amount) {
    // Loop to update the buffer with new phaseAngle values
    for (uint16_t i = 0; i < move_amount; ++i) {
        int16_t phase_index = phase_index_start + i;
        int16_t buffer_index = (buffer_start_index + i) % buffer_manager.buffer_size;

        // Debugging: Log the details of buffer update before actual update
        #ifdef BUFFER_DEBUG
        printf("[AdptSweepAddDataPoint] [DEBUG] Before update -> Buffer index: %d, Phase index: %d, PhaseAngle: %.6f\n", buffer_index, phase_index, phaseAngles[phase_index]);
        #endif

        // Copy phaseAngle and optionally set impedance
        buffer_manager.buffer[buffer_index].phaseAngle = phaseAngles[phase_index];
        //buffer_manager.buffer[buffer_index].impedance = 0.0; // Mock impedance, update if needed

        // Debugging: Log the update in the buffer
        #ifdef BUFFER_DEBUG
        printf("[AdptSweepAddDataPoint] Buffer updated at phase index %d with phaseAngle %.6f (buffer index: %d)\n", phase_index, phaseAngles[phase_index], buffer_index);
        #endif
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
 *    - Calls `AdptSweepAddDataPoint` to load the new phase angles into the buffer starting at `buffer_start_index`.
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
 * @see AdptSweepAddDataPoint
 * @see buffer_manager
 */
void update_buffer_for_direction(const double* phaseAngles, int direction, int move_amount) {
    int16_t buffer_start_index;
    int16_t phase_index_start;

    // Log initial indices before movement
    #ifdef BUFFER_DEBUG
    printf("[update_buffer] Initial phase index start: %d, buffer start index: %d, move amount: %d\n",
           buffer_manager.current_phase_index, buffer_manager.current_buffer_index, move_amount);
    #endif

    if (direction == LEFT_SIDE) {
        // Moving left
        phase_index_start = buffer_manager.current_phase_index - move_amount;
        buffer_start_index = buffer_manager.current_buffer_index - move_amount;
        
         // Check if the buffer indices are out of bounds
        if (buffer_start_index < 0 || buffer_start_index + buffer_manager.window_size > buffer_manager.buffer_size) {
            printf("[update_buffer] Buffer index out of bounds. Setting error flag.\n");
            set_boundary_error_flag(1);  // Set the boundary error flag
            return;  // Exit the function to stop further processing
        }

        if (buffer_start_index < 0) {
            buffer_start_index += buffer_manager.buffer_size;
        }

        // Debugging: Log the movement and new indices
        #ifdef BUFFER_DEBUG
        printf("[update_buffer] Moving left by %d. New phase index start: %d, new buffer start index: %d\n",
               move_amount, phase_index_start, buffer_start_index);
        #endif

        // Update indices
        buffer_manager.current_phase_index = phase_index_start;
        buffer_manager.current_buffer_index = buffer_start_index;

        // Update analysis start index if necessary
        if (phase_index_start < analysis_start_index) {
            analysis_start_index = phase_index_start;
        }

    } else if (direction == RIGHT_SIDE) {
        // Moving right
        phase_index_start = buffer_manager.current_phase_index + buffer_manager.window_size;
        buffer_start_index = (buffer_manager.current_buffer_index + buffer_manager.window_size) % buffer_manager.buffer_size;
        
         // Check if the buffer indices are out of bounds
        if (buffer_start_index < 0 || buffer_start_index + buffer_manager.window_size > buffer_manager.buffer_size) {
            printf("[update_buffer] Buffer index out of bounds. Setting error flag.\n");
            set_boundary_error_flag(1);  // Set the boundary error flag
            return;  // Exit the function to stop further processing
        }


        // Debugging: Log the movement and new indices
        #ifdef BUFFER_DEBUG
        printf("[update_buffer] Moving right by %d. New phase index start: %d, new buffer start index: %d\n",
               move_amount, phase_index_start, buffer_start_index);
        #endif

        // Update indices
        buffer_manager.current_phase_index += move_amount;
        buffer_manager.current_buffer_index = (buffer_manager.current_buffer_index + move_amount) % buffer_manager.buffer_size;

        // Update analysis end index if necessary
        int16_t new_end_index = buffer_manager.current_phase_index + buffer_manager.window_size - 1;
        if (new_end_index > analysis_end_index) {
            analysis_end_index = new_end_index;
        }

    } else {
        // Invalid direction
        #ifdef BUFFER_DEBUG
        printf("[update_buffer] Direction undecided. No movement.\n");
        #endif
        return;
    }

    // Instead of calling AdptSweepAddDataPoint, set up buffer update info
    buffer_update_info.needs_update = true;
    buffer_update_info.phase_index_start = phase_index_start;
    buffer_update_info.buffer_start_index = buffer_start_index;
    buffer_update_info.move_amount = move_amount;

    // Reset buffer_shift_offset since we've updated the buffer
    buffer_shift_offset = 0;

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
 * overwriting the previous values in the buffer.
 *
 * @param phaseAngles A pointer to the array of phase angles.
 * @param phase_angle_size The total size of the phase angle array.
 */
void handle_undecided_case(const double* phaseAngles, uint16_t phase_angle_size) {
    #ifdef BUFFER_DEBUG
    printf("[handle_undecided_case] Entering undecided case handler.\n");
    #endif

    // Jump forward by the window size in the phaseAngle array
    buffer_manager.current_phase_index += buffer_manager.window_size;

    // Ensure we don't exceed the phase angle size
    if (buffer_manager.current_phase_index + buffer_manager.window_size > phase_angle_size) {
        #ifdef BUFFER_DEBUG
        printf("[handle_undecided_case] Phase index exceeds phase angle size. Stopping analysis.\n");
        #endif
        return;
    }

    // Reset buffer_shift_offset since we're discarding previous analysis
    buffer_shift_offset = 0;

    // Update global analysis tracking indices
    analysis_start_index = buffer_manager.current_phase_index;
    analysis_end_index = buffer_manager.current_phase_index + buffer_manager.window_size - 1;

    //#ifdef BUFFER_DEBUG
    printf("[handle_undecided_case] Jumping to phase index: %d\n", buffer_manager.current_phase_index);
    //#endif

    // Set up buffer update info instead of calling AdptSweepAddDataPoint
    buffer_update_info.needs_update = true;
    buffer_update_info.phase_index_start = buffer_manager.current_phase_index;
    buffer_update_info.buffer_start_index = buffer_manager.current_buffer_index;
    buffer_update_info.move_amount = buffer_manager.window_size;
}

/**
 * @brief Handles the case when the segment analysis result is negative undecided.
 *
 * This function jumps backward in the phaseAngle array by the window size, discards the previous analysis,
 * and reloads the buffer from the new phase index. The new portion of the phaseAngle array is loaded into the buffer,
 * overwriting the previous values in the buffer.
 *
 * @param phaseAngles A pointer to the array of phase angles.
 * @param phase_angle_size The total size of the phase angle array.
 */
void handle_negative_undecided_case(const double* phaseAngles, uint16_t phase_angle_size) {
    #ifdef BUFFER_DEBUG
    printf("[handle_negative_undecided_case] Entering negative undecided case handler.\n");
    #endif

    // Jump backward by the window size in the phaseAngle array
    buffer_manager.current_phase_index -= buffer_manager.window_size;

    // Ensure we don't go below 0 in the phase angle size
    if (buffer_manager.current_phase_index < 0) {
        buffer_manager.current_phase_index = 0;
        #ifdef BUFFER_DEBUG
        printf("[handle_negative_undecided_case] Phase index went below 0. Stopping analysis.\n");
        #endif
        return;
    }

    // Reset buffer_shift_offset since we're discarding previous analysis
    buffer_shift_offset = 0;

    // Update global analysis tracking indices
    analysis_start_index = buffer_manager.current_phase_index;
    analysis_end_index = buffer_manager.current_phase_index + buffer_manager.window_size - 1;

    #ifdef BUFFER_DEBUG
    printf("[handle_negative_undecided_case] Jumping to phase index: %d\n", buffer_manager.current_phase_index);
    #endif

    // Set up buffer update info instead of calling AdptSweepAddDataPoint
    buffer_update_info.needs_update = true;
    buffer_update_info.phase_index_start = buffer_manager.current_phase_index;
    buffer_update_info.buffer_start_index = buffer_manager.current_buffer_index;
    buffer_update_info.move_amount = buffer_manager.window_size;
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
    int16_t new_phase_index_start = buffer_manager.current_phase_index + (direction == RIGHT_SIDE ? move_amount : -move_amount);
    int16_t new_phase_index_end = new_phase_index_start + buffer_manager.window_size - 1;

    bool need_buffer_update = false;

    //#ifdef BUFFER_DEBUG
    printf("[move_window_and_update_if_needed] Attempting to move window. Direction: %d, Move amount: %d\n", direction, move_amount);
    //#endif

    if (direction == LEFT_SIDE) {
        if (new_phase_index_start < analysis_start_index) {
            need_buffer_update = true;
        }
    } else if (direction == RIGHT_SIDE) {
        if (new_phase_index_end > analysis_end_index) {
            need_buffer_update = true;
        }
    } else {
        return;  // Invalid direction
    }

    if (need_buffer_update) {
        // Set up buffer update info
        update_buffer_for_direction(phaseAngles, direction, move_amount);
    } else {
        // Adjust indices without updating the buffer
        if (direction == LEFT_SIDE) {
            buffer_manager.current_phase_index = new_phase_index_start;
            buffer_manager.current_buffer_index = (buffer_manager.current_buffer_index - move_amount + buffer_manager.buffer_size) % buffer_manager.buffer_size;
        } else {
            buffer_manager.current_phase_index = new_phase_index_start;
            buffer_manager.current_buffer_index = (buffer_manager.current_buffer_index + move_amount) % buffer_manager.buffer_size;
        }

        // Update buffer_shift_offset
        //buffer_shift_offset = (buffer_shift_offset + move_amount) % buffer_manager.buffer_size;

        //#ifdef BUFFER_DEBUG
        printf("[move_window_and_update_if_needed] Window already in buffer. Shifted by %d. New buffer_shift_offset: %d\n", move_amount, buffer_shift_offset);
        //#endif
    }

    // Update analysis indices if necessary
    if (new_phase_index_start < analysis_start_index) {
        analysis_start_index = new_phase_index_start;
    }
    if (new_phase_index_end > analysis_end_index) {
        analysis_end_index = new_phase_index_end;
    }

    // Debug: Print the updated indices
    //#ifdef BUFFER_DEBUG
    printf("[move_window_and_update_if_needed] Updated window to phase index [%d - %d]\n", 
           new_phase_index_start, new_phase_index_end);
    printf("Updated buffer index to %d\n", buffer_manager.current_buffer_index);
    //#endif
}

/**
 * @brief Performs a buffer update if one is required.
 *
 * This function checks whether a buffer update is necessary based on the `buffer_update_info` structure,
 * which contains information about pending buffer updates. If a buffer update is required (indicated by the
 * `needs_update` flag), the function proceeds to update the phase angles in the buffer by calling 
 * `AdptSweepAddDataPoint`. Once the update is completed, the `needs_update` flag is reset to prevent
 * redundant updates.
 *
 * ### Why This Function is Needed:
 * In the sliding window analysis, the buffer is constantly shifting to analyze different sections of the
 * phase angle data. However, to avoid unnecessary buffer updates or synchronization issues between logic and
 * buffer management, the buffer updates are deferred and handled only when required. This function decouples
 * the buffer update logic from the main analysis flow, ensuring that the buffer is updated asynchronously,
 * only when changes are pending.
 *
 * ### Use Case:
 * This function is typically called after the state machine logic determines that a window shift has occurred
 * or a buffer needs to be updated based on changes in the sliding window. For instance:
 * - After determining a new direction of analysis in the sliding window, the state machine will check if
 *   a buffer update is needed and call this function to handle it.
 * - If phase angles have moved out of the current buffer range, this function ensures that the necessary
 *   new data is loaded into the buffer without interrupting the analysis logic.
 * 
 * By separating the update check and the actual buffer modification, this function helps maintain a clean
 * separation of concerns, ensuring that the analysis and buffer updates happen in an orderly, non-blocking fashion.
 * 
 * @param phaseAngles Pointer to the array of phase angles that are used in the analysis.
 */
void perform_buffer_update_if_needed(const double* phaseAngles) {
    if (buffer_update_info.needs_update) {
        AdptSweepAddDataPoint(
            phaseAngles,
            buffer_update_info.phase_index_start,
            buffer_update_info.buffer_start_index,
            buffer_update_info.move_amount
        );

        // Reset the update flag
        buffer_update_info.needs_update = false;
    }
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

    for (int phase_index = analysis_start_index; phase_index <= analysis_end_index && phase_index < total_size; ++phase_index) {
        int relative_position = phase_index - buffer_manager.current_phase_index;
        int adjusted_buffer_index = (buffer_manager.current_buffer_index + relative_position - buffer_shift_offset + buffer_manager.buffer_size) % buffer_manager.buffer_size;

        // Debugging: Log the index mappings
        #ifdef BUFFER_DEBUG
        printf("[print_analysis_interval] Phase Index: %d -> Buffer Index: %d\n", phase_index, adjusted_buffer_index);
        #endif

        printf("Phase Index %d -> Buffer Index %d: Phase Angle: %.6f, Impedance: %.6f\n",
               phase_index, adjusted_buffer_index,
               buffer_manager.buffer[adjusted_buffer_index].phaseAngle,
               buffer_manager.buffer[adjusted_buffer_index].impedance);
    }

    printf("--------------------------------------\n");
}
