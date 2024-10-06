// ics.h

#ifndef ICS_H
#define ICS_H

#include <stdint.h>
#include <stdbool.h>
#include <math.h>

// Typedef for the callback function
//typedef void (*HalIcsSweepSampleCb_t)(const int16_t real, const int16_t imaginary, const bool isLast);

// Typedef for the callback function (now using double for real and imaginary values)
typedef void (*HalIcsSweepSampleCb_t)(const double real, const double imaginary, const bool isLast);


// Function prototypes
void HalIcsSetStartFreq(const float_t freq);
void HalIcsSetFreqInc(const float_t increment);
void HalIcsSetNoOfSamples(const uint16_t number);
bool HalIcsStartFreqSweep(HalIcsSweepSampleCb_t cb);

// Global variables (extern if needed)
extern float_t startFrequency;
extern float_t frequencyIncrement;
extern uint16_t numberOfSamples;

#endif // ICS_H
