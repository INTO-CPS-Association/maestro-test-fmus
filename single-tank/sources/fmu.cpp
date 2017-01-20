/*
 * fmi.cpp
 *
 *  Created on: Aug 21, 2015
 *      Author: parallels
 */

#include<stdio.h>

/*
 * fmi2Functions.c
 *
 *  Created on: May 22, 2015
 *      Author: kel
 */

#include "fmi2Functions.h"

#include <string>

#define VAR_COUNT 4
#define REAL_LEVEL_ID 0
#define REAL_LEVEL_MAX_ID 1
#define REAL_LEVEL_MIN_ID 2
#define BOOL_VALVE_ID 3

struct fmuState
{
	std::string name;
	double time;
	double reals[VAR_COUNT];
	int ints[VAR_COUNT];
	bool bools[VAR_COUNT];
	std::string strings[VAR_COUNT];
};

// ---------------------------------------------------------------------------
// FMI functions
// ---------------------------------------------------------------------------
extern "C" fmi2Component fmi2Instantiate(fmi2String instanceName, fmi2Type fmuType, fmi2String fmuGUID,
		fmi2String fmuResourceLocation, const fmi2CallbackFunctions *functions, fmi2Boolean visible,
		fmi2Boolean loggingOn)
{

	fmuState* state = new fmuState();
	state->name = std::string(instanceName);
	state->time = 0.0;

	for (int i = 0; i < VAR_COUNT; i++)
	{
		state->reals[i] = 0;
		state->ints[i] = 0;
		state->bools[i] = 0;
		state->strings[i] = std::string("");
	}

	return (fmi2Component) state;
}

extern "C" fmi2Status fmi2SetupExperiment(fmi2Component c, fmi2Boolean toleranceDefined, fmi2Real tolerance,
		fmi2Real startTime, fmi2Boolean stopTimeDefined, fmi2Real stopTime)
{
	return fmi2OK;
}

extern "C" fmi2Status fmi2EnterInitializationMode(fmi2Component c)
{
	return fmi2OK;
}

extern "C" fmi2Status fmi2ExitInitializationMode(fmi2Component c)
{
	return fmi2OK;
}

extern "C" fmi2Status fmi2Terminate(fmi2Component c)
{
	return fmi2OK;
}

extern "C" fmi2Status fmi2Reset(fmi2Component c)
{
	return fmi2OK;
}

extern "C" void fmi2FreeInstance(fmi2Component c)
{
}

// ---------------------------------------------------------------------------
// FMI functions: class methods not depending of a specific model instance
// ---------------------------------------------------------------------------

extern "C" const char* fmi2GetVersion()
{
	return fmi2Version;
}

extern "C" const char* fmi2GetTypesPlatform()
{
	return fmi2TypesPlatform;
}

// ---------------------------------------------------------------------------
// FMI functions: logging control, setters and getters for Real, Integer,
// Boolean, String
// ---------------------------------------------------------------------------

extern "C" fmi2Status fmi2SetDebugLogging(fmi2Component c, fmi2Boolean loggingOn, size_t nCategories,
		const fmi2String categories[])
{
	return fmi2OK;
}

extern "C" fmi2Status fmi2GetReal(fmi2Component c, const fmi2ValueReference vr[], size_t nvr, fmi2Real value[])
{
	for (int i = 0; i < nvr; i++)
	{
		value[i] = ((fmuState*) c)->reals[vr[i]];
	}
	return fmi2OK;
}

extern "C" fmi2Status fmi2GetInteger(fmi2Component c, const fmi2ValueReference vr[], size_t nvr, fmi2Integer value[])
{
	for (int i = 0; i < nvr; i++)
	{
		value[i] = ((fmuState*) c)->ints[vr[i]];
	}
	return fmi2OK;
}

extern "C" fmi2Status fmi2GetBoolean(fmi2Component c, const fmi2ValueReference vr[], size_t nvr, fmi2Boolean value[])
{
	for (int i = 0; i < nvr; i++)
	{
		value[i] = ((fmuState*) c)->bools[vr[i]];
	}
	return fmi2OK;
}

extern "C" fmi2Status fmi2GetString(fmi2Component c, const fmi2ValueReference vr[], size_t nvr, fmi2String value[])
{
	for (int i = 0; i < nvr; i++)
	{

		std::string s = ((fmuState*) c)->strings[vr[i]];

		char *b = new char[s.length() + 1]
		{ };
		std::copy(s.c_str(), s.c_str() + s.length(), b);

		value[i] = b;
	}
	return fmi2OK;

	return fmi2OK;
}

extern "C" fmi2Status fmi2SetReal(fmi2Component c, const fmi2ValueReference vr[], size_t nvr, const fmi2Real value[])
{
	for (int i = 0; i < nvr; i++)
	{
		((fmuState*) c)->reals[vr[i]] = value[i];
	}
	return fmi2OK;
}

extern "C" fmi2Status fmi2SetInteger(fmi2Component c, const fmi2ValueReference vr[], size_t nvr,
		const fmi2Integer value[])
{
	for (int i = 0; i < nvr; i++)
	{
		((fmuState*) c)->ints[vr[i]] = value[i];
	}
	return fmi2OK;
}

extern "C" fmi2Status fmi2SetBoolean(fmi2Component c, const fmi2ValueReference vr[], size_t nvr,
		const fmi2Boolean value[])
{
	for (int i = 0; i < nvr; i++)
	{
		((fmuState*) c)->bools[vr[i]] = value[i];
	}
	return fmi2OK;
}

extern "C" fmi2Status fmi2SetString(fmi2Component c, const fmi2ValueReference vr[], size_t nvr,
		const fmi2String value[])
{
	for (int i = 0; i < nvr; i++)
	{
		((fmuState*) c)->strings[vr[i]] = std::string(value[i]);
	}
	return fmi2OK;
}

extern "C" fmi2Status fmi2GetFMUstate(fmi2Component c, fmi2FMUstate* FMUstate)
{
	return fmi2OK;
}
extern "C" fmi2Status fmi2SetFMUstate(fmi2Component c, fmi2FMUstate FMUstate)
{
	return fmi2OK;
}
extern "C" fmi2Status fmi2FreeFMUstate(fmi2Component c, fmi2FMUstate* FMUstate)
{
	return fmi2OK;
}

extern "C" fmi2Status fmi2SerializedFMUstateSize(fmi2Component c, fmi2FMUstate FMUstate, size_t *size)
{
	return fmi2Error;
}
extern "C" fmi2Status fmi2SerializeFMUstate(fmi2Component c, fmi2FMUstate FMUstate, fmi2Byte serializedState[],
		size_t size)
{
	return fmi2Error;
}
extern "C" fmi2Status fmi2DeSerializeFMUstate(fmi2Component c, const fmi2Byte serializedState[], size_t size,
		fmi2FMUstate* FMUstate)
{
	return fmi2Error;
}

extern "C" fmi2Status fmi2GetDirectionalDerivative(fmi2Component c, const fmi2ValueReference vUnknown_ref[],
		size_t nUnknown, const fmi2ValueReference vKnown_ref[], size_t nKnown, const fmi2Real dvKnown[],
		fmi2Real dvUnknown[])
{
	return fmi2Error;
}

// ---------------------------------------------------------------------------
// Functions for FMI for Co-Simulation
// ---------------------------------------------------------------------------
#ifdef FMI_COSIMULATION
/* Simulating the slave */
extern "C" fmi2Status fmi2SetRealInputDerivatives(fmi2Component c, const fmi2ValueReference vr[], size_t nvr,
		const fmi2Integer order[], const fmi2Real value[])
{
	return fmi2Error;
}

extern "C" fmi2Status fmi2GetRealOutputDerivatives(fmi2Component c, const fmi2ValueReference vr[], size_t nvr,
		const fmi2Integer order[], fmi2Real value[])
{
	return fmi2Error;
}

extern "C" fmi2Status fmi2CancelStep(fmi2Component c)
{
	return fmi2Error;
}

//#define CONTROLLER
extern "C" fmi2Status fmi2DoStep(fmi2Component c, fmi2Real currentCommunicationPoint, fmi2Real communicationStepSize,
		fmi2Boolean noSetFMUStatePriorToCurrentPoint)
{
	fmuState* state = ((fmuState*) c);

	state->time += communicationStepSize;

	double* level = &state->reals[REAL_LEVEL_ID];

#ifdef CONTROLLER

	bool* valveState = &state->bools[BOOL_VALVE_ID];

	if (*level >= state->reals[REAL_LEVEL_MAX_ID])
	{
		*valveState = fmi2True;
	}

	if (*level <= state->reals[REAL_LEVEL_MIN_ID])
	{
		*valveState = fmi2False;
	}

#else
	if (state->bools[BOOL_VALVE_ID] == fmi2True)
	{
		*level = (*level) - 1;
	} else
	{
		*level = (*level) + 1;
	}

#endif

	return fmi2OK;
}

extern "C" fmi2Status fmi2GetStatus(fmi2Component c, const fmi2StatusKind s, fmi2Status *value)
{
	return fmi2OK;
}

extern "C" fmi2Status fmi2GetRealStatus(fmi2Component c, const fmi2StatusKind s, fmi2Real *value)
{
	return fmi2OK;
}

extern "C" fmi2Status fmi2GetIntegerStatus(fmi2Component c, const fmi2StatusKind s, fmi2Integer *value)
{
	return fmi2OK;
}

extern "C" fmi2Status fmi2GetBooleanStatus(fmi2Component c, const fmi2StatusKind s, fmi2Boolean *value)
{
	return fmi2OK;
}

extern "C" fmi2Status fmi2GetStringStatus(fmi2Component c, const fmi2StatusKind s, fmi2String *value)
{
	return fmi2OK;
}

/* INTO cps specific*/
extern "C" fmi2Status fmi2GetMaxStepsize(fmi2Component c, fmi2Real* size)
{
	return fmi2OK;
}

// ---------------------------------------------------------------------------
// Functions for FMI2 for Model Exchange
// ---------------------------------------------------------------------------
#else
/* Enter and exit the different modes */
fmi2Status fmi2EnterEventMode(fmi2Component c)
{
	return fmi2Error;
}

fmi2Status fmi2NewDiscreteStates(fmi2Component c, fmi2EventInfo *eventInfo)
{
	return fmi2Error;
}

fmi2Status fmi2EnterContinuousTimeMode(fmi2Component c)
{
	return fmi2Error;
}

fmi2Status fmi2CompletedIntegratorStep(fmi2Component c, fmi2Boolean noSetFMUStatePriorToCurrentPoint,
		fmi2Boolean *enterEventMode, fmi2Boolean *terminateSimulation)
{
	return fmi2Error;
}

/* Providing independent variables and re-initialization of caching */
fmi2Status fmi2SetTime(fmi2Component c, fmi2Real time)
{
	return fmi2Error;
}

fmi2Status fmi2SetContinuousStates(fmi2Component c, const fmi2Real x[], size_t nx)
{
	return fmi2Error;
}

/* Evaluation of the model equations */
fmi2Status fmi2GetDerivatives(fmi2Component c, fmi2Real derivatives[], size_t nx)
{
	return fmi2Error;
}

fmi2Status fmi2GetEventIndicators(fmi2Component c, fmi2Real eventIndicators[], size_t ni)
{
	return fmi2Error;
}

fmi2Status fmi2GetContinuousStates(fmi2Component c, fmi2Real states[], size_t nx)
{
	return fmi2Error;
}

fmi2Status fmi2GetNominalsOfContinuousStates(fmi2Component c, fmi2Real x_nominal[], size_t nx)
{
	return fmi2Error;
}
#endif // Model Exchange

