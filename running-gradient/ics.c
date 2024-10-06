// ics.c - Mock implementation

#include "ics.h"
#include "sliding_window_analysis.h" // Include any necessary headers
#include <stdio.h>

#include "mes_buffers.h"

// Externally defined variables
extern const double phaseAngles[];
extern size_t phase_angle_size;
extern BufferUpdateInfo buffer_update_info;
extern MqsRawDataPoint_t* buffer_manager_buffer;

// Global variables to store the parameters
float_t startFrequency = 0.0f;
float_t frequencyIncrement = 0.0f;
uint16_t numberOfSamples = 0;

// Simulated hardware status
static bool isIdle = true;
static bool isSweepAborted = false;
static HalIcsSweepSampleCb_t sampleCb = NULL;

// Function to set start frequency
void HalIcsSetStartFreq(const float_t freq) {
    startFrequency = freq;
}

// Function to set frequency increment
void HalIcsSetFreqInc(const float_t increment) {
    frequencyIncrement = increment;
}

// Function to set number of samples
void HalIcsSetNoOfSamples(const uint16_t number) { // move_amount
    numberOfSamples = number;
}

// Helper function to simulate real and imaginary values
void simulate_real_imaginary(double phaseAngle, double magnitude, int16_t* real_int16, int16_t* imaginary_int16) {
    double radians = phaseAngle * M_PI / 180.0;
    double real = magnitude * cos(radians);
    double imaginary = magnitude * sin(radians);

    // Ensure the values fit in int16_t range
    if (real > 32767.0) real = 32767.0;
    if (real < -32768.0) real = -32768.0;
    if (imaginary > 32767.0) imaginary = 32767.0;
    if (imaginary < -32768.0) imaginary = -32768.0;

    *real_int16 = (int16_t)real;
    *imaginary_int16 = (int16_t)imaginary;
}

// Mock HalIcsStartFreqSweep function (phaseAngle as real part, 1.0 as imaginary)
bool HalIcsStartFreqSweep(HalIcsSweepSampleCb_t cb) {
    // Check if the system is currently idle. If it's not, we can't start a new sweep.
    sampleCb = cb;  // Store the callback function for later use

    // Retrieve the number of samples (buffer_update_info.move_amount) and the starting index
    uint16_t num_samples = numberOfSamples;  // Number of samples to process
    int16_t phase_index_start = (int16_t)startFrequency;  // Starting index for the sweep

    // Debug: Print the start of the sweep
    printf("Starting frequency sweep with %d samples from phase index %d\n", num_samples, phase_index_start);

    // Ensure we do not exceed the bounds of the phaseAngles array
    if (phase_index_start + num_samples > phase_angle_size) {
        // Adjust the number of samples to fit within the bounds of the phase_angle_size array
        num_samples = phase_angle_size - phase_index_start;
        //printf("Adjusted number of samples to prevent array overrun: %d samples\n", num_samples);
    }

    // Loop through the number of samples and simulate data collection
    for (uint16_t i = 0; i < num_samples; ++i) {
        // Calculate the current phase index
        int16_t phase_index = phase_index_start + i;

        // Fetch the phase angle for the current sample
        double phaseAngle = phaseAngles[phase_index];

        // Calculate the current buffer index (buffer starts at buffer_update_info.buffer_start_index)
        int16_t buffer_index = (buffer_update_info.buffer_start_index + i);
        
        // Debug: Print phase angle and buffer index
        //printf("Sample %d: Phase index %d, PhaseAngle = %.6f, Buffer index = %d\n", 
        //       i + 1, phase_index, phaseAngle, buffer_index);

        // Call the callback function, passing the phaseAngle as the real value and 1.0 as the imaginary value
        // The last argument 'isLast' is true only for the last sample in the loop
        cb(phaseAngle, 1.0, (i == num_samples - 1));

        // Optionally, store phaseAngle directly in the buffer (uncomment if needed)
        // buffer_manager_buffer[buffer_index].phaseAngle = phaseAngle;
    }

    // Debug: Sweep process completed
    printf("Frequency sweep completed with %d samples.\n", num_samples);

    // Set the system back to idle state after the sweep completes
    isIdle = true;

    return true;
}


