#ifndef MES_BUFFERS_H
#define MES_BUFFERS_H

#include "mqs_def.h"
#include <stdint.h>

typedef struct {
	const MqsRawDataSweep_t* setup;
	MqsRawDataPoint_t* data;
	uint_fast16_t count;
	uint_fast16_t peak;
} MesSweep_t;

// Declaration of buffers and related variables
extern uint8_t rawBuffer[];
extern MqsRawDataSet_t* rawData;
extern MesSweep_t rawBaseSweep;
extern MesSweep_t rawAfterExposureSweep;
extern MesSweep_t* currentRawSweep;

extern uint8_t filterBuffer[];
extern MqsRawDataSet_t* filterData;
extern MesSweep_t filterBaseSweep;
extern MesSweep_t filterAfterExposureSweep;
extern MesSweep_t* currentFilterSweep;

void MesSweepApplySetupDefaults(MqsRawDataSweep_t* const setup);
void MesSweepSetupAndClear(MesSweep_t* const sweep, const MqsRawDataSweep_t* const setup, MqsRawDataPoint_t* const data);

#endif // MES_BUFFERS_H
