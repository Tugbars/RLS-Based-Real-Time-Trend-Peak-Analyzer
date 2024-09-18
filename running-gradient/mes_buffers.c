#include "mes_buffers.h"
#include "mqs_def.h"
#include <string.h>

// Definition of buffers and related variables
uint8_t rawBuffer[MQS_RAW_DATA_SET_SIZE_MAX];
MqsRawDataSet_t* rawData = (MqsRawDataSet_t*)rawBuffer;
MesSweep_t rawBaseSweep;
MesSweep_t rawAfterExposureSweep;
MesSweep_t* currentRawSweep = &rawBaseSweep;

uint8_t filterBuffer[MQS_RAW_DATA_SET_SIZE_MAX];
MqsRawDataSet_t* filterData = (MqsRawDataSet_t*)filterBuffer;
MesSweep_t filterBaseSweep;
MesSweep_t filterAfterExposureSweep;
MesSweep_t* currentFilterSweep = &filterBaseSweep;

//buradan alıyor. 
void MesSweepApplySetupDefaults(MqsRawDataSweep_t* const setup)
{
	setup->startFrequency = 10400;
	setup->frequencyIncrement = 1;
	//printf("%f, ", setup->frequencyIncrement);
	setup->dataCount = 360;
}


void MesSweepSetupAndClear(MesSweep_t* const sweep, const MqsRawDataSweep_t* const setup, MqsRawDataPoint_t* const data)
{
    memset(sweep, 0x00, sizeof(MesSweep_t));
	sweep->setup = setup;
	sweep->data = data;
}
