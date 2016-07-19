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

#include <string.h>

//#include <map>

//static std::map <int, ExternalClient> clients;
#include <vector>
#include <sstream>
#include <cmath>

const fmi2CallbackFunctions *g_functions;
static uintptr_t state = 0;
static uintptr_t expectedState = 0;
std::string* name;

#define SSTR( x ) dynamic_cast< std::ostringstream & >( \
        ( std::ostringstream() << std::dec << x ) ).str()

template<class T>
static void log(const fmi2CallbackFunctions *functions, fmi2ComponentEnvironment componentEnvironment,
		fmi2String instanceName, fmi2Status status, fmi2String category, fmi2String message, T arg)
{
	if (functions != NULL && functions->logger != NULL)
	{
		functions->logger(componentEnvironment, instanceName, status, category, message, arg);
	}
}

static void notimplemented(fmi2Component c, fmi2String message)
{
	std::string base("Not implemented: %s");
	std::string m(message);
	if (g_functions != NULL)
	{
		log(g_functions, (void*) 2, "", fmi2Error, "error", (base + m).c_str(), "");
	}
}

//static void message(fmi2String message)
//{
//	std::string m(message);
//	if (g_functions != NULL)
//	{
//		log(g_functions, (void*) 1, "Info", fmi2OK, "Message: %s", m.c_str());
//	}
//}
template<class T>
static void fmiprintf(fmi2String message, T arg)
{
	if (g_functions != NULL)
	{
		log(g_functions, (void*) 2, name->c_str(), fmi2OK, "logAll", message, arg);
	}
}

// ---------------------------------------------------------------------------
// FMI functions
// ---------------------------------------------------------------------------
extern "C" fmi2Component fmi2Instantiate(fmi2String instanceName, fmi2Type fmuType, fmi2String fmuGUID,
		fmi2String fmuResourceLocation, const fmi2CallbackFunctions *functions, fmi2Boolean visible,
		fmi2Boolean loggingOn)
{
	name = new std::string(instanceName);
	g_functions = functions;
	fmiprintf("instantiating rollback-test %s\n", "");
	return (void*) 2;
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
	for(int i =0; i < nvr ; i++)
	{
		value[i] = 0;
	}
	return fmi2OK;
}

extern "C" fmi2Status fmi2GetInteger(fmi2Component c, const fmi2ValueReference vr[], size_t nvr, fmi2Integer value[])
{
	for(int i =0; i < nvr ; i++)
	{
		value[i] = 0;
	}
	return fmi2OK;
}

extern "C" fmi2Status fmi2GetBoolean(fmi2Component c, const fmi2ValueReference vr[], size_t nvr, fmi2Boolean value[])
{
	for(int i =0; i < nvr ; i++)
	{
		value[i] = false;
	}
	return fmi2OK;
}

extern "C" fmi2Status fmi2GetString(fmi2Component c, const fmi2ValueReference vr[], size_t nvr, fmi2String value[])
{
	return fmi2Warning;
}

extern "C" fmi2Status fmi2SetReal(fmi2Component c, const fmi2ValueReference vr[], size_t nvr, const fmi2Real value[])
{
	return fmi2OK;
}

extern "C" fmi2Status fmi2SetInteger(fmi2Component c, const fmi2ValueReference vr[], size_t nvr,
		const fmi2Integer value[])
{
	return fmi2OK;
}

extern "C" fmi2Status fmi2SetBoolean(fmi2Component c, const fmi2ValueReference vr[], size_t nvr,
		const fmi2Boolean value[])
{
	return fmi2OK;
}

extern "C" fmi2Status fmi2SetString(fmi2Component c, const fmi2ValueReference vr[], size_t nvr,
		const fmi2String value[])
{
	return fmi2OK;
}

extern "C" fmi2Status fmi2GetFMUstate(fmi2Component c, fmi2FMUstate* FMUstate)
{
	*FMUstate = (void*) state;
	fmiprintf("returning state: %lu", state);
	return fmi2OK;
}
extern "C" fmi2Status fmi2SetFMUstate(fmi2Component c, fmi2FMUstate FMUstate)
{
	state = (uintptr_t) FMUstate;
	fmiprintf("Changed state to: %lu", state);
	return fmi2OK;
}
extern "C" fmi2Status fmi2FreeFMUstate(fmi2Component c, fmi2FMUstate* FMUstate)
{
	fmiprintf("freeing state: %lu", state);
	return fmi2OK;
}

extern "C" fmi2Status fmi2SerializedFMUstateSize(fmi2Component c, fmi2FMUstate FMUstate, size_t *size)
{
	notimplemented(c, "fmi2SerializedFMUstateSize");
	return fmi2Error;
}
extern "C" fmi2Status fmi2SerializeFMUstate(fmi2Component c, fmi2FMUstate FMUstate, fmi2Byte serializedState[],
		size_t size)
{
	notimplemented(c, "fmi2SerializeFMUstate");
	return fmi2Error;
}
extern "C" fmi2Status fmi2DeSerializeFMUstate(fmi2Component c, const fmi2Byte serializedState[], size_t size,
		fmi2FMUstate* FMUstate)
{
	notimplemented(c, "fmi2DeSerializeFMUstate");
	return fmi2Error;
}

extern "C" fmi2Status fmi2GetDirectionalDerivative(fmi2Component c, const fmi2ValueReference vUnknown_ref[],
		size_t nUnknown, const fmi2ValueReference vKnown_ref[], size_t nKnown, const fmi2Real dvKnown[],
		fmi2Real dvUnknown[])
{
	notimplemented(c, "fmi2GetDirectionalDerivative");
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
	notimplemented(c, "fmi2SetRealInputDerivatives");
	return fmi2Error;
}

extern "C" fmi2Status fmi2GetRealOutputDerivatives(fmi2Component c, const fmi2ValueReference vr[], size_t nvr,
		const fmi2Integer order[], fmi2Real value[])
{
	notimplemented(c, "fmi2GetRealOutputDerivatives");
	return fmi2Error;
}

extern "C" fmi2Status fmi2CancelStep(fmi2Component c)
{
	notimplemented(c, "fmi2CancelStep");
	return fmi2Error;
}

float epsilon = 0.000000001;

static bool set = false;

static bool first = true;
static double currentTime = 0;

static double nextSuccessTime = 0;

static fmi2Real communicationStepSizeHistory[2] =
{ 0 };

extern "C" fmi2Status fmi2DoStep(fmi2Component c, fmi2Real currentCommunicationPoint, fmi2Real communicationStepSize,
		fmi2Boolean noSetFMUStatePriorToCurrentPoint)
{
//	fmiprintf("doStep current state is: %lu\n", state);
	if (g_functions != NULL && g_functions->logger != NULL)
	{
		g_functions->logger((void*) 2, name->c_str(), fmi2OK, "logAll",
				"doStep curPoin %f, size %f, current state is: %lu", currentCommunicationPoint, communicationStepSize,
				state);
	}

	fmi2Status ret = fmi2Error;
	communicationStepSizeHistory[1] = communicationStepSizeHistory[0];
	communicationStepSizeHistory[0] = communicationStepSize;

	double timeOnSuccess = currentCommunicationPoint + communicationStepSize;

	if (first)
	{
		nextSuccessTime = timeOnSuccess;
		communicationStepSizeHistory[0] = communicationStepSize;
		communicationStepSizeHistory[1] = communicationStepSize;
		first = false;
	}

	g_functions->logger((void*) 2, name->c_str(), fmi2OK, "logAll", "doStep checking time accepts %f, proposed %f",
			nextSuccessTime, timeOnSuccess, state);

	if (std::abs(nextSuccessTime - timeOnSuccess) < epsilon)
	{
		currentTime = timeOnSuccess;
		if (expectedState != state)
		{
			fmiprintf("Stepping from the wrong state%s", "");
			ret = fmi2Warning;
		} else
		{
//			nextValidTime = 0; //disable
			ret = fmi2OK;
		}

	} else
	{
		fmiprintf("Stepping discard%s", "");
		currentTime = nextSuccessTime;
		ret = fmi2Discard;
	}

//pre-calc next time

	if (ret == fmi2OK)
	{
		if (!set)
		{
			fmiprintf("Setting next /2 step using %f", communicationStepSizeHistory[1]);
			nextSuccessTime = timeOnSuccess + (communicationStepSizeHistory[0] / 2);
			set = true;
		} else
		{
			fmiprintf("Setting next /2 step using %f --", communicationStepSizeHistory[0]);
			nextSuccessTime = timeOnSuccess + communicationStepSizeHistory[1]; //normal step
			set = false;
		}
		expectedState++;
	}

	state++;
	g_functions->logger((void*) 2, name->c_str(), ret, "logAll", "Current time is now %f and state %lu", currentTime,
			state);

	g_functions->logger((void*) 2, name->c_str(), ret, "logAll", "Next accept time is %f and state %lu",
			nextSuccessTime, expectedState);

	return ret;
}

extern "C" fmi2Status fmi2GetStatus(fmi2Component c, const fmi2StatusKind s, fmi2Status *value)
{
	return fmi2OK;
}

extern "C" fmi2Status fmi2GetRealStatus(fmi2Component c, const fmi2StatusKind s, fmi2Real *value)
{
	if (s == fmi2StatusKind::fmi2LastSuccessfulTime)
	{
		fmiprintf("fmi2GetRealStatus time is now %f", currentTime);
		*value = currentTime;
	}
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

