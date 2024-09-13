#include "sliding_window_analysis.h"
#include <stdio.h>

/******************************************************************************/
/* Global Variables */
/******************************************************************************/

typedef union {
    uint8_t value;
    struct {
        uint8_t isPeakFound : 1;
        uint8_t isUndecided : 1;
        uint8_t isSweepRequested : 1;
        uint8_t isSweepDone : 1;
        uint8_t reserved : 4;  // Unused bits for future use
    };
} SwpStatus_t;

SlidingWindowAnalysisContext ctx; // Sliding window analysis context
SwpStatus_t currentStatus = { .value = 0x0 };  // Initialize all flags to 0
SwpState_t currentState = SWP_WAITING;  // Initialize to SWP_WAITING state

/******************************************************************************/
/* Function Prototypes (Internal Functions) */
/******************************************************************************/

// Function Declarations for OnEntry and OnExit functions
static void OnEntryInitialAnalysis(void);
static void OnEntrySegmentAnalysis(void);
static void OnEntryUpdateBufferDirection(void);
static void OnEntryUndecidedTrendCase(void);
static void OnEntryPeakFindingAnalysis(void);
static void OnEntryWaiting(void);

static void OnExitInitialAnalysis(void);
static void OnExitSegmentAnalysis(void);
static void OnExitUpdateBufferDirection(void);
static void OnExitUndecidedTrendCase(void);
static void OnExitPeakFindingAnalysis(void);
static void OnExitWaiting(void);

static void initBufferManager(void);
static SwpState_t NextState(SwpState_t state);

/******************************************************************************/
/* State Function Definitions */
/******************************************************************************/
typedef struct {
    void (*onEntry)(void);
    void (*onExit)(void);
    bool isComplete;
} StateFuncs_t;

static StateFuncs_t STATE_FUNCS[SWP_STATE_LAST] = {
    {OnEntryInitialAnalysis, OnExitInitialAnalysis, false},  // SWP_INITIAL_ANALYSIS
    {OnEntrySegmentAnalysis, OnExitSegmentAnalysis, false},  // SWP_SEGMENT_ANALYSIS
    {OnEntryUpdateBufferDirection, OnExitUpdateBufferDirection, false},  // SWP_UPDATE_BUFFER_DIRECTION
    {OnEntryUndecidedTrendCase, OnExitUndecidedTrendCase, false},  // SWP_UNDECIDED_TREND_CASE
    {OnEntryPeakFindingAnalysis, OnExitPeakFindingAnalysis, false},  // SWP_PEAK_FINDING_ANALYSIS
    {OnEntryWaiting, OnExitWaiting, false}  // SWP_WAITING
};

/******************************************************************************/
/* Sliding Window Analysis Functions */
/******************************************************************************/
static void initBufferManager(void) {
    // Initialize the buffer manager with appropriate parameters
    // Arguments:
    // buffer          -> The buffer that will store the phase and impedance values (array of MqsRawDataPoint_t)
    // BUFFER_SIZE     -> Total size of the buffer (maximum number of data points the buffer can hold)
    // 30              -> Sliding window size (number of data points to analyze in each window)
    // 0               -> Starting index in the phaseAngles array (this could be any index where you want to start the analysis)
    // 11300.0         -> Starting frequency for the frequency sweep (in Hz, e.g., starting at 11300 Hz)
    // 1.0             -> Frequency increment per step (in Hz, e.g., increment by 1 Hz for each data point)
    
    init_buffer_manager(buffer, BUFFER_SIZE, WINDOW_SIZE, 220, 11300.0, 1.0);

    // Debugging: Print initialization state
    printf("[DEBUG] Buffer Manager initialized:\n");
    printf("Buffer size: %d, Window size: %d, Start index: %d\n", BUFFER_SIZE, WINDOW_SIZE, 140);
}

void startSlidingWindowAnalysis(const double* phaseAngles, uint16_t phase_angle_size, Callback_t callback) {
    ctx.phaseAngles = phaseAngles;
    ctx.phase_angle_size = phase_angle_size;
    ctx.callback = callback;

    // Initialize the buffer manager
    initBufferManager();
    
        // Enable sweep request
    currentStatus.isSweepRequested = 1;  // Start sweep request

    // Set the initial state to waiting
    currentState = SWP_WAITING;

    // Set the initial state and start the state machine
    SwpProcessStateChange();  // Start the state machine
}

/******************************************************************************/
/* State Machine Entry Functions */
/******************************************************************************/
static void OnEntryInitialAnalysis(void) {
    printf("→→→→→Entering SWP_INITIAL_ANALYSIS state.\n");

    load_initial_buffer(ctx.phaseAngles, ctx.phase_angle_size);
    
      // Debugging buffer state after load
    printf("[DEBUG] Buffer Manager State after initial loading:\n");
    printf("Current buffer index: %d\n", buffer_manager.current_buffer_index);
    printf("Current phase index: %d\n", buffer_manager.current_phase_index);
    printf("Window size: %d\n", buffer_manager.window_size);
    printf("Buffer size: %d\n", buffer_manager.buffer_size);
    currentStatus.isSweepRequested = true;
    STATE_FUNCS[SWP_INITIAL_ANALYSIS].isComplete = true;
    SwpProcessStateChange();  // Move to next state
}

static void OnEntrySegmentAnalysis(void) {
    printf("→→→→→Entering SWP_SEGMENT_ANALYSIS state.\n");

    SegmentAnalysisResult result = segment_trend_and_concavity_analysis(
        &buffer_manager.buffer[buffer_manager.current_buffer_index],
        buffer_manager.window_size,
        0.5
    );

    ctx.direction = result.nextDirection;

    if (ctx.direction == ON_PEAK) {
        currentStatus.isPeakFound = 1;
    } else if (ctx.direction == UNDECIDED) {
        currentStatus.isUndecided = 1;
    }

    STATE_FUNCS[SWP_SEGMENT_ANALYSIS].isComplete = true;
    SwpProcessStateChange();
}

static void OnEntryUpdateBufferDirection(void) {
    printf("→→→→→Entering SWP_UPDATE_BUFFER_DIRECTION state.\n");

    update_buffer_for_direction(ctx.phaseAngles, ctx.direction);

    STATE_FUNCS[SWP_UPDATE_BUFFER_DIRECTION].isComplete = true;
    SwpProcessStateChange();  // Move to next state
}

static void OnEntryUndecidedTrendCase(void) {
    printf("→→→→→→Entering SWP_UNDECIDED_TREND_CASE state.\n");

    handle_undecided_case(ctx.phaseAngles, ctx.phase_angle_size);

    STATE_FUNCS[SWP_UNDECIDED_TREND_CASE].isComplete = true;
    SwpProcessStateChange();
}

static void OnEntryPeakFindingAnalysis(void) {
    printf("----→Entering SWP_PEAK_FINDING_ANALYSIS state.\n");

    verify_peak_at_index(buffer_manager.current_buffer_index);

    STATE_FUNCS[SWP_PEAK_FINDING_ANALYSIS].isComplete = true;
    currentStatus.isSweepDone = 1;
    SwpProcessStateChange();
}

static void OnEntryWaiting(void) {
    printf("→→→→→Entering SWP_WAITING state.\n");

    if (currentStatus.isSweepRequested) {
        STATE_FUNCS[SWP_WAITING].isComplete = true;
    }
}

/******************************************************************************/
/* Exit Functions */
/******************************************************************************/
static void OnExitInitialAnalysis(void) {
    printf("Exiting SWP_INITIAL_ANALYSIS state.\n");
}

static void OnExitSegmentAnalysis(void) {
    printf("Exiting SWP_SEGMENT_ANALYSIS state.\n");
}

static void OnExitUpdateBufferDirection(void) {
    printf("Exiting SWP_UPDATE_BUFFER_DIRECTION state.\n");
}

static void OnExitUndecidedTrendCase(void) {
    printf("Exiting SWP_UNDECIDED_TREND_CASE state.\n");
}

static void OnExitPeakFindingAnalysis(void) {
    printf("Exiting SWP_PEAK_FINDING_ANALYSIS state.\n");
}

static void OnExitWaiting(void) {
    printf("Exiting SWP_WAITING state.\n");
}

/******************************************************************************/
/* State Machine Process and Transitions */
/******************************************************************************/
void SwpProcessStateChange(void) {
    bool stateChanged;
    
    // The state machine should process until no further state changes occur
    do {
        // Determine the next state from the current state
        SwpState_t nextState = NextState(currentState);  
        stateChanged = (nextState != currentState);

        // Reset the current status flags before moving to the next state
        currentStatus.value = 0x0;

        if (stateChanged) {
            printf("[DEBUG] Transitioning from state %d to state %d\n", currentState, nextState);

            // Exit the current state (if applicable)
            if (STATE_FUNCS[currentState].onExit != NULL) {
                STATE_FUNCS[currentState].onExit();
            }

            // Move to the new state
            currentState = nextState;

            // Reset the completion flag for the new state
            STATE_FUNCS[currentState].isComplete = false;

            // Enter the new state
            if (STATE_FUNCS[currentState].onEntry != NULL) {
                STATE_FUNCS[currentState].onEntry();
            }
        }

        // If the new state is complete, proceed to the next state
        if (STATE_FUNCS[currentState].isComplete) {
            printf("[DEBUG] State %d is complete. Proceeding to next state...\n", currentState);
        }

    } while (stateChanged && STATE_FUNCS[currentState].isComplete);

    // Handle the callback execution when we reach the waiting state
    if (currentState == SWP_WAITING && ctx.callback) {
        printf("[DEBUG] Executing callback as we have reached SWP_WAITING.\n");
        ctx.callback();  // Execute the callback function
        ctx.callback = NULL;  // Clear the callback to avoid re-execution
    }
}

/******************************************************************************/
/* Helper Functions for State Transitioning */
/******************************************************************************/

static SwpState_t NextState(SwpState_t state) {
    switch (state) {
        case SWP_INITIAL_ANALYSIS:
            if (currentStatus.isPeakFound) {
                return SWP_PEAK_FINDING_ANALYSIS;
            }
            if (!currentStatus.isPeakFound) {
                return SWP_SEGMENT_ANALYSIS;
            }
            if (!currentStatus.isUndecided) {
                return SWP_UNDECIDED_TREND_CASE;
            }
            return SWP_INITIAL_ANALYSIS;

        case SWP_SEGMENT_ANALYSIS:
            if (currentStatus.isUndecided) {
                return SWP_UNDECIDED_TREND_CASE;
            }
            if (currentStatus.isPeakFound) {
                return SWP_PEAK_FINDING_ANALYSIS;
            }
            return SWP_UPDATE_BUFFER_DIRECTION;

        case SWP_UPDATE_BUFFER_DIRECTION:
            return SWP_SEGMENT_ANALYSIS;

        case SWP_UNDECIDED_TREND_CASE:
            return SWP_SEGMENT_ANALYSIS;

        case SWP_PEAK_FINDING_ANALYSIS:
            if (currentStatus.isSweepDone) {
                return SWP_WAITING;
            }
            return SWP_PEAK_FINDING_ANALYSIS;

        case SWP_WAITING:
            if (currentStatus.isSweepRequested) {
                return SWP_INITIAL_ANALYSIS;
            }
            return SWP_WAITING;

        default:
            return SWP_WAITING;
    }
}
//220'de o dediğim problem verify oluyor. 
