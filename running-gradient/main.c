#include <stdio.h>
#include <time.h>
#include <stdlib.h>

#include "buffer_manager.h"
#include "mes_buffers.h"

#include "windowed_running_gradient.h"
#include "running_cubic_gradient.h"
#include "running_peak_analysis.h"
#include "running_quadratic_gradient.h"
#include "rls_analysis_parameters.h"
#include "sliding_window_analysis.h"

void myCallbackFunction(void) {
    if (boundaryErrorOccurred) {
        printf("Callback executed: State machine interrupted due to boundary error.\n");
        boundaryErrorOccurred = false;  // Reset the persistent flag after handling the error
    } else {
        printf("Callback executed: State machine returned to SWP_WAITING.\n");
    }
}

void PrepareBaseSweep(MesSweep_t *const sweep, MqsRawDataSet_t *const data)
{
  MesSweepApplySetupDefaults(&data->base);   // sets the parameters up such as increment, starting frequency etc.
  data->base.startFrequency = 210;
  MesSweepSetupAndClear(sweep, &data->base, data->data); // memset
}

int main() {

    PrepareBaseSweep(&rawBaseSweep, rawData); 
    currentRawSweep = &rawBaseSweep;
    
    int start_index = 210;
    startSlidingWindowAnalysis(currentRawSweep, start_index, myCallbackFunction); 
    
	return 0;
}
