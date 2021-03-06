/* Implementation of util.h */

#include "util.h"
#include "result.h"
#include "dsblock.h"
#ifndef ONLY_INCLUDE_INLINE_INTEGRATION
#include "integration.h"
#endif
#include "fmiFunctions_fwd.h"


#include "adymosim.h"
#ifndef ONLY_INCLUDE_INLINE_INTEGRATION
#include <cvode/cvode.h>
/* Sundials */
#endif

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>


/* ----------------- local variables ----------------- */

/* when compiling as a single complilation unit, this is already defined */
#ifndef FMU_SOURCE_SINGLE_UNIT
extern Component* globalComponent2;
#endif

/* ----------------- function definitions ----------------- */

extern void ModelicaMessage(const char *string);
extern void ModelicaError(const char *string);

/* -------------------------------------------------------- */
DYMOLA_STATIC void util_logger(Component* comp, FMIString instanceName, FMIStatus status,
							   FMIString category, FMIString message, ...)
{
	char buf[4096];
	va_list ap;
	int capacity;
	int n = 0;

	if (comp->loggingOn == FMIFalse && (status < FMIWarning)) {
		return;
	}

	va_start(ap, message);
	capacity = sizeof(buf) - 1;
	if (comp->logbufp > comp->logbuf) {
		/* add buffered messages */
#if defined(_MSC_VER) && _MSC_VER>=1400
		n = vsnprintf_s(buf, capacity, _TRUNCATE, comp->logbuf, ap);
#else
		buf[capacity]=0;
		n = vsnprintf(buf, capacity, comp->logbuf, ap);
#endif
		if (n >= capacity || n <= 0) {
			goto done;
		}
	}

#if defined(_MSC_VER) && _MSC_VER>=1400
	vsnprintf_s(buf + n, capacity - n, _TRUNCATE, message, ap);
#else
	buf[capacity]=0;
	vsnprintf(buf + n, capacity - n, message, ap);
#endif
	va_end(ap);
#ifdef FMI_2
	comp->functions->logger(comp->functions->componentEnvironment, instanceName, status, category, buf);
#else
	comp->functions.logger(comp, instanceName, status, category, buf);
#endif

done:
	comp->logbufp = comp->logbuf;
}

/* -------------------------------------------------------- */
DYMOLA_STATIC void util_buflogger(Component* comp, FMIString instanceName, FMIStatus status,
								  FMIString category, FMIString message, ...)
{
	va_list ap;
	int capacity;
	int n;

	if (comp->loggingOn == FMIFalse && (status < FMIWarning)) {
		return;
	}

	capacity = (int) (sizeof(comp->logbuf) - (comp->logbufp - comp->logbuf) - 1);
	if (capacity <= 0) {
		goto fail;
	}

	va_start(ap, message);
	n = vsnprintf(comp->logbufp, capacity, message, ap);
	if (n >= capacity || n <= 0) {
		goto fail;
	}
	va_end(ap);
	comp->logbufp += n;
	return;

fail:
	comp->logbufp = comp->logbuf;
}

/* -------------------------------------------------------- */
DYMOLA_STATIC FMIString util_strdup(const FMICallbackFunctions *functions, FMIString s)
{
	static const FMIString nullString = "<NULL>";
	static const int maxLen = 1024;
	FMIString pos = s;
	int len;
	char* newS;

	if (s == NULL) {
		s = nullString;
	}
	len = (int) strlen(s);
	if (len > maxLen) {
		len = maxLen;
	}

	
	newS = (FMIChar*) functions->allocateMemory(len + 1, sizeof(FMIChar));
	if (newS == NULL) {
		return NULL;
	}
	strncpy(newS, s, len);
	return newS;
}
#ifndef ONLY_INCLUDE_INLINE_INTEGRATION
/* ------------------------------------------------------ */
DYMOLA_STATIC int util_check_flag(void *flagvalue, char *funcname, int opt, Component* comp)
{
	int *errflag;

	/* Check if SUNDIALS function returned NULL pointer - no memory allocated */
	if (opt == 0 && flagvalue == NULL) {
		util_logger(comp, comp->instanceName, FMIError, "",
			"SUNDIALS_ERROR: %s() failed - returned NULL pointer", funcname);
		return 1;
	}

	/* Check if flag < 0 */
	else if (opt == 1) {
		errflag = (int *) flagvalue;
		if (*errflag < 0) {
			util_logger(comp, comp->instanceName, FMIError, "",
				"SUNDIALS_ERROR: %s() failed with flag = %s",
				funcname, CVodeGetReturnFlagName(*errflag));
			return 1;
		}
	}

	return 0;
}
#endif /* ONLY_INCLUDE_INLINE_INTEGRATION */

/* ------------------------------------------------------ */
DYMOLA_STATIC int GetDymolaOneIteration(struct DYNInstanceData*);
DYMOLA_STATIC void SetDymolaOneIteration(struct DYNInstanceData*, int val);
DYMOLA_STATIC void SetDymolaJacobianPointers(struct DYNInstanceData*did_, double * QJacobian_,double * QBJacobian_,double * QCJacobian_,double * QDJacobian_,int QJacobianN_,
	int QJacobianNU_,int QJacobianNY_,double * QJacobianSparse_,int * QJacobianSparseR_,int * QJacobianSparseC_,int QJacobianNZ_);
DYMOLA_STATIC int util_refresh_cache(Component* comp, int idemand, const char* label, FMIBoolean* iterationConverged)
{
#ifdef FMI_2
	switch (comp->mStatus) {
	    case modelInstantiated:
			if (GetDymolaOneIteration(comp->did) == 0) {
				SetDymolaOneIteration(comp->did, 5);
			}
			break;
		case modelInitializationMode:
			if (GetDymolaOneIteration(comp->did) == 0) {
				SetDymolaOneIteration(comp->did, 5);
			}
			idemand=iDemandStart;
			break;
		case modelEventMode:
		case modelEventMode2:
			if(comp->firstEventCall){
				SetDymolaOneIteration(comp->did, 2);
				comp->firstEventCall=FMIFalse;
			}else if (GetDymolaOneIteration(comp->did) == 0) {
				SetDymolaOneIteration(comp->did, 5);
			}
			idemand=iDemandEventHandling;
			break;
		default:
			//assert(GetDymolaOneIteration(comp->did) == 0);
			break;
	}
#endif
	comp->QiErr=0;
	if (comp->did) {
		globalComponent2=comp;
		dsblock_tid(&idemand, &comp->icall, &comp->time, comp->states, 0,             
			comp->inputs, comp->parameters, 0, 0, comp->derivatives,       
			comp->outputs, comp->auxiliary,                                
			comp->crossingFunctions, comp->duser, comp->iuser, (void**) comp->sParameters, comp->did, &comp->QiErr, 0);
		globalComponent2=0;
	} else {
		dsblock_(&idemand, &comp->icall, &comp->time, comp->states, 0,             
			comp->inputs, comp->parameters, 0, 0, comp->derivatives,       
			comp->outputs, comp->auxiliary,                                
			comp->crossingFunctions, comp->duser, comp->iuser, (void**) comp->sParameters, &comp->QiErr);
	}
	/*Only check convergence criteria for functions that need it i.e. event handling*/
	if(iterationConverged){
		*iterationConverged = (GetDymolaOneIteration(comp->did)  == 0) ? FMITrue : FMIFalse;
	}
	/*Always reset*/
	SetDymolaOneIteration(comp->did, 0);
	comp->valWasSet = 0;
	if (comp->QiErr == 0) {
		return 0;
	} 
	
	if (label != NULL) {
		util_logger(comp, comp->instanceName, FMIError, "",                     
			"%s: dsblock_ failed, QiErr = %d", label, comp->QiErr);
	}
	return comp->QiErr;
}

/* ------------------------------------------------------ */
DYMOLA_STATIC FMIStatus util_error(Component* comp)
{	
	/*Termination should be done exteranly when error is thrown*/
	switch(comp->mStatus) {
	    case modelInstantiated:
			/* internally state "error" and "terminated" are equivalent */
			comp->mStatus = modelTerminated;
			break;
		case modelTerminated:
			break;
		default:
			if (comp->isCoSim == FMITrue) {
#ifndef FMU_SKIP_CO_SIM
				fmiTerminateSlave_(comp);
#endif
			} else {
#ifdef FMI_2
				fmiTerminateModel_(comp);
#elif !defined(FMU_SKIP_MODEL_EXCHANGE)
				fmiTerminate_(comp);
#endif
			}
			break;
	}
	return FMIError;
}

/* ------------------------------------------------------ */
extern void InitializeDymosimStruct(struct BasicDDymosimStruct* basicD, struct BasicIDymosimStruct* basicI);
extern void SetDymolaEventOptional(struct DYNInstanceData*did_, int val);
extern int GetDymolaHaveEventIterated(struct DYNInstanceData*did_);
DYMOLA_STATIC FMIStatus util_initialize_model(FMIComponent c, FMIBoolean toleranceControlled, FMIReal relativeTolerance, FMIBoolean complete)
{
	Component* comp = (Component*) c;
	FMIStatus status = FMIOK;
	int QiErr = 0;
	/* Start of integration */
	int idemand = iDemandStart;
	comp->terminationByModel = FMIFalse;
/* #ifdef _DEBUG fails on dSPACE platforms where _DEBUG is always defined, 0 or 1. */
#if _DEBUG
	/* catch numerical exceptions */
	/*_EM_INVALID_ is not compatible with the underlying handling of NaN */
	_controlfp(0, _EM_ZERODIVIDE | _EM_OVERFLOW);
#endif
	if (comp->mStatus != modelInstantiated) {
        util_logger(comp, comp->instanceName, FMIWarning, "", "model cannot be initialized in current state(%d)", comp->mStatus);
		return FMIWarning;
	}

	/* Ignore toleranceControlled for external integration according to recommendation from Hans Olsson:
	The tolerance for the equations to be solved at events are based on the crossing function
	tolerance (to avoid chattering), which is normally stricter than relativeTolerance.
	Hence, relativeTolerance is normally not significant. */
	if (toleranceControlled == FMITrue) {
		util_logger(comp, comp->instanceName, FMIWarning, "", "fmiInitialize: tolerance control not supported");
	}
	InitializeDymosimStruct(comp->dstruct, comp->istruct);
#ifdef OVERIDE_DYMOLA_EVENT_OPTIONAL
	SetDymolaEventOptional(comp->did,  1);
#else
  #ifdef ONLY_INCLUDE_INLINE_INTEGRATION
	if(comp->isCoSim){
		SetDymolaEventOptional(comp->did,  1);
	}else{
		SetDymolaEventOptional(comp->did,  0);
	}
  #else
	if(comp->isCoSim && comp->istruct->mInlineIntegration){
			SetDymolaEventOptional(comp->did,  1);
	}else{
			SetDymolaEventOptional(comp->did,  0);
	}
  #endif
#endif
	comp->icall = iDemandStart - 1;
	if (complete ==FMIFalse) {
		SetDymolaOneIteration(comp->did, 2);
	} else {
		SetDymolaOneIteration(comp->did, 0);
	}
	QiErr = util_refresh_cache(comp, idemand, NULL, NULL);
	if (QiErr != 0) {
		status = FMIError;
		util_logger(comp, comp->instanceName, status, "",
			"fmiInitialize: dsblock_ failed, QiErr = %d", QiErr);
		util_logger(comp, comp->instanceName, status, "",
			"Unless otherwise indicated by error messages, possible errors are (non-exhaustive):\n"
			"The model references external data that is not present on the target machine, at least not with the same location.\n",
			QiErr);
		return util_error(comp);
	} 

	if (comp->storeResult == FMITrue) {
		result_setup(comp);
	}

	return status;
}

/* ------------------------------------------------------ */
DYMOLA_STATIC FMIStatus util_event_update(FMIComponent c, FMIBoolean intermediateResults,
#ifdef FMI_2
/* needs another argument since not in eventInfo for FMI 2*/
FMIBoolean* terminateSolution
#else
fmiEventInfo* eventInfo
#endif
)
{
	Component* comp = (Component*) c;
	FMIStatus status = FMIOK;
	int QiErr = 0;
	FMIBoolean converged = 0;

	/* make sure idemand up to auxilliary variables are computed prior to starting event iteration */
	/* (necessary for e.g. co-simulating Modelica.Electrical.Machines.Examples.SynchronousInductionMachines.SMEE_Generator) */
	{
		static IDemandLevel idemand[] = {iDemandVariable, iDemandEventHandling};
		int i;
		for (i = 0; i < 2; i++) {
			/* TODO: figure out why intermediate results are not retrieved if dsblock_ is first called with idemand == iDemandVariable */
			if (i == 0) {
				if (intermediateResults == FMITrue) {
					continue;
				}
			} else {
				/* configure actual event iteration */
				if (intermediateResults == FMITrue)
				{
					if (comp->eventIterationOnGoing) {
						SetDymolaOneIteration(comp->did, 3);
					} else {
						SetDymolaOneIteration(comp->did, 2);
						comp->eventIterationOnGoing = 1;
					}
				} else {
					if (comp->eventIterationOnGoing) {
						SetDymolaOneIteration(comp->did, 4);
					} else {
						SetDymolaOneIteration(comp->did, 0);
					}
				}
				comp->icall = 0;
			}
			QiErr = util_refresh_cache(comp, idemand[i], NULL, &converged);
#ifdef FMI_2
			*terminateSolution = FMIFalse;
#else
			eventInfo->terminateSimulation = FMIFalse;
#endif
			if (QiErr>=-995 && QiErr<=-990) QiErr=0; /* Ignore special cases for now */
			if (QiErr != 0) {
				if (QiErr == -999) {
					util_logger(comp, comp->instanceName, FMIOK, "",
						"event updating: simulation terminated by model");
#ifdef FMI_2
					*terminateSolution = FMITrue;
#else
					eventInfo->terminateSimulation = FMITrue;
#endif
					return FMIOK;
				} else {
					util_logger(comp, comp->instanceName, FMIError, "",
						"event updating: dsblock_ failed, QiErr = %d", QiErr);
					return util_error(comp);
				}
			}
		}

#ifdef FMI_2
		/* for FMI 2 we only use this function internally for co-simulation and then always expect converged */
		assert(intermediateResults == FMIFalse && converged == FMITrue);
#else
		if (converged == FMIFalse)
		{
			/* more iterations needed */
			assert(intermediateResults == FMITrue);
			eventInfo->iterationConverged = FMIFalse;
		} else {
			eventInfo->iterationConverged = FMITrue;
			comp->eventIterationOnGoing = 0;
		}
#endif
	}

	/* for FMI 2 we only use this internally for co-simulation and do not use these values */
#ifndef FMI_2
	/* always fmiFalse since we hide states that might be exchanged */
	eventInfo->stateValueReferencesChanged = FMIFalse;
	/* TODO: introduce a flag in dymosim that can be checked for this purpose (for now it is faster to
	to fetch states each time than to conclude the truth value from copy + compare) */
	eventInfo->stateValuesChanged = (comp->nStates > 0) ? FMITrue : FMIFalse;
	eventInfo->upcomingTimeEvent = (comp->dstruct->mNextTimeEvent < (1.0E37 - 1)) ? FMITrue : FMIFalse;
	if (eventInfo->upcomingTimeEvent == FMITrue) {
		eventInfo->nextEventTime = comp->dstruct->mNextTimeEvent;
	}
#endif
	return status;
}
DYMOLA_STATIC FMIStatus util_initialize_slave(FMIComponent c,
	                                          FMIReal      relativeTolerance,
											  FMIReal      tStart,
											  FMIBoolean   StopTimeDefined,
											  FMIReal      tStop)
{
	FMIStatus status;
	FMIReal relTol = relativeTolerance;
#ifndef FMI_2
	FMIEventInfo eventInfo;
#endif
	Component* comp = (Component*) c;
	fmiSetTime_(c, tStart);

#ifdef FMI_2
	status = util_initialize_model(c, FMIFalse, 0, FMIFalse);
#else
	status = fmiInitialize_(c, FMIFalse, 0, &eventInfo);
#endif /* FMI_2 */
	if (status != FMIOK) {
		return status;
	}
	comp->StopTimeDefined = StopTimeDefined;
	comp->tStop = tStop;
	comp->valWasSet = 0;
#ifndef ONLY_INCLUDE_INLINE_INTEGRATION
	if(!comp->istruct->mInlineIntegration){
		if (integration_setup(c, relTol) != integrationOK) {
			return FMIFatal;
		}
	}
#endif
	return FMIOK;
}

#ifdef FMI_2
/* ------------------------------------------------------ */
DYMOLA_STATIC FMIStatus util_exit_model_initialization_mode(FMIComponent c, const char* label, ModelStatus nextStatus)
{
	Component* comp = (Component*) c;
	int QiErr;
	
	if (comp->mStatus != modelInitializationMode) {
		util_logger(comp, comp->instanceName, FMIWarning, "", "%s: may only called in initialization mode", label);
		return FMIWarning;
	}

	SetDymolaOneIteration(comp->did, 4);
	QiErr = util_refresh_cache(comp, iDemandStart, NULL, NULL);
	if (QiErr != 0) {
		return util_error(comp);
	}
	/* reset */
	SetDymolaOneIteration(comp->did, 0);

	comp->mStatus = nextStatus;
	return FMIOK;
}
#endif /* FMI_2 */

