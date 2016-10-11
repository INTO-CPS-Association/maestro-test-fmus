/* dsmodel.h */

/*
 * Copyright (C) 2005-2008 Dynasim AB.
 * All rights reserved.
 *
 */
/* For using model as dll */
/*
   Basically defines an ode-interface:
     der(x)=f(t,x,u,p)  <- dymCalculateDerivatives
	 y=g(t,x,u,p)       <- dymCalculateOutputs
	 
	 Trigger event when sign(h(t,x,u,p)) changes <- dymCalculateCrossingFunctions
*/
/* Procedure for using:*/
/* 1. Compile dsmodel.c to a dll with DYMOSIM_DLL_EXPORTS defined */
/* 2. Open the dll or link with this header and DYMOSIM_DLL_IMPORTS defined*/
/* 3. The work-procedure:

      m=dymAllocateModel();
	  dymChangeTime(m, ...);
	  dymChangeStates(m, ...);
	  dymChangeInputs(m, ...);
	  dymChangeParameters(m, ...);
	  [Optionally for HILS: dymSetFixedStepSize(m, ...)]

	  dymInitialize(m);
	  dymGetStates(m, ...);  Changes due to initialization 
	  
	  During each step:
	     [Optionally at start of step: dymChangeParameters(m, ...);]
		 dymChangeTime(m,...);
		 dymChangeStates(m,...);
		 if (directfeedthrough) dymChangeInputs(m, ...);
		 
		 dymCalculateOutputs(m, ...);
		 if (!directfeedthrough) dymChangeInputs(m, ...);
		 dymCalculateDerivatives(m, ...);

		 When checking for events at end of a succesful step:
		   dymCalculateCrossingFunctions(m, ...);
		   For accurate crossing detection, interpolate time/states/inputs and use:
		      dymChangeTime(m, ...);
			  dymChangeStates(m, ...);
			  dymChangeInputs(m, ...);
			  dymCalculateCrossingFunctions(m, ...);

		   When an event is detected [use dymEventHandling(m,1)] or at least an event is possible [use dymEventHandling(m,0)].
		      switch(dymEventHandling(m, ...)) {
			  case 0:break; Nothing happened.
			  case 1:goto failedSimulation;
			  case 2:goto terminateSimulation;
			  default:
			    dymGetStates(m,...); And reset solver for next step even if states were not changed.
		   }

  Simulation is completed - either succesfully or not:

   terminateSimulation:
      dymTerminateAndFreeModel(m);
      m=0;
      return;
   failedSimulation:
      dymFreeModel(m);
	  m=0;
	  return;

*/


#ifndef DSMODEL_H
#define DSMODEL_H

#ifdef DYMOSIM_DLL_EXPORTS
#define DYMOSIM_EXPORT __declspec(dllexport) 
#elif defined(DYMOSIM_DLL_IMPORTS)
#define DYMOSIM_EXPORT __declspec(dllimport) 
#else
#define DYMOSIM_EXPORT 
#endif

struct DsModel;
DYMOSIM_EXPORT struct DsModel* dymAllocateModel();
DYMOSIM_EXPORT void dymFreeModel(struct DsModel*);
DYMOSIM_EXPORT long dymNrParameters(struct DsModel*);
DYMOSIM_EXPORT long dymNrStates(struct DsModel*);
DYMOSIM_EXPORT long dymNrInputs(struct DsModel*);
DYMOSIM_EXPORT long dymNrOutputs(struct DsModel*);
DYMOSIM_EXPORT long dymNrCrossingFunctions(struct DsModel*); /* Includes one for time event */
DYMOSIM_EXPORT int dymHasDirectFeedThrough(struct DsModel*); /* True if model has direct feed-through */

/* Variables can also be changed using these helpers. Nothing is changed until the next computation */
DYMOSIM_EXPORT void dymChangeTime(struct DsModel*,double time);
DYMOSIM_EXPORT void dymChangeStates(struct DsModel*,const double*newStates);
DYMOSIM_EXPORT void dymChangeInputs(struct DsModel*,const double*newInputs);
DYMOSIM_EXPORT void dymChangeParameters(struct DsModel*,const double*newParameters);

/* Routines below return 0 for ok and non-zero for error */

/* Computes initial values for states at the given time. This guarantees that initialization has been performed. */
DYMOSIM_EXPORT int dymInitialize(struct DsModel*); 

/* Compute values after changing time/states/inputs/(parameters) */
DYMOSIM_EXPORT int dymCalculateDerivatives(struct DsModel*,double*derivatives);
DYMOSIM_EXPORT int dymCalculateOutputs(struct DsModel*,double*outputs);

DYMOSIM_EXPORT void dymGetStates(struct DsModel*,double*states); 
/* Returns current values of states*/
/* States may change during initialization and at events due to reinit */

DYMOSIM_EXPORT void dymGetParameters(struct DsModel*,double*parameters); 
/* For reading the values of parameters. It will later support reading defaults from model */

/* For event handling */
DYMOSIM_EXPORT int dymCalculateCrossingFunctions(struct DsModel*,double*indicators,double*nextTime); /* Next-time is scalar and may be nil */
DYMOSIM_EXPORT int dymEventHandling(struct DsModel*,int knownEvent); 
/* if you have ensured that this routine is only called when you have detected events for the model set knownEvent=1, 
   otherwise set knownEvent=0 */
/* Incorrectly setting knownEvent=1 can cause slowdowns */
/* 0 - no event, 1 - error, 2 - terminate simulation, 3 - normal event. Please use dymGetStates afterwards in case return value>=2 */

DYMOSIM_EXPORT int dymTerminateAndFreeModel(struct DsModel*); /* Use after successful simulation */

/* For setting up HILS-integration. Call before Initialize, and only call if you will be using a fixed step size algorithm. */
DYMOSIM_EXPORT int dymSetFixedStepSize(struct DsModel*,double fixedStep);

/* For better handling of direct feed through. Use these instead of NrInputs,NrOutputs,ChangeInputs,CalculateOutputs, and HasDirectFeedThrough */
/* Currently, Dymola does not provide so fine-grained resolution and thus all inputs and outputs are either seen as direct or not */
DYMOSIM_EXPORT long dymNrDirectInputs(struct DsModel*);
DYMOSIM_EXPORT long dymNrDirectOutputs(struct DsModel*);
DYMOSIM_EXPORT long dymNrNonDirectInputs(struct DsModel*);
DYMOSIM_EXPORT long dymNrNonDirectOutputs(struct DsModel*);
DYMOSIM_EXPORT void dymChangeDirectInputs(struct DsModel*,const double*newInputs);
DYMOSIM_EXPORT void dymChangeNonDirectInputs(struct DsModel*,const double*newInputs);

DYMOSIM_EXPORT int dymCalculateDirectOutputs(struct DsModel*,double*outputs);
DYMOSIM_EXPORT int dymCalculateNonDirectOutputs(struct DsModel*,double*outputs);

/* For generating Jacobians */
DYMOSIM_EXPORT int dymHasAnalyticJacobian(struct DsModel*); /* 1 if true */

DYMOSIM_EXPORT int dymComputeJacobianABCD(struct DsModel*,double*A,double*B,double*C,double*D);
/* Linearize the model (analytic or good numeric). 0 result if ok. */
DYMOSIM_EXPORT int dymComputeJacobianA(struct DsModel*,double*A);
/* Linearize the model for solver (analytic or good numeric). 0 result if ok. */

#if !defined(DYMOSIM_DLL_IMPORTS)
/*   Include interface code */
#include "dsmodelDll.c"
#endif
#endif
