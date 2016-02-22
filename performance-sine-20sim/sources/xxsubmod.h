/**********************************************************
 * This file is generated by 20-sim ANSI-C Code Generator
 *
 *  file:  src\xxsubmod.h
 *  subm:  Sine
 *  model: Sinus Integrator Model
 *  expmt: Sinus Integrator Model
 *  date:  September 22, 2015
 *  time:  9:16:27 am
 *  user:  INTO-CPS
 *  from:  20-sim 4.5 Professional Single
 *  build: 4.5.1.5561
 **********************************************************/

/* This file describes the model functions
   that are supplied for computation.

   The model itself is the xxmodel.c file
*/

#ifndef XX_SUBMOD_H
#define XX_SUBMOD_H

/* Our own include files */
#include "xxmodel.h"



/* The submodel functions */
void XXInitializeSubmodel (XXDouble t);
void XXCalculateSubmodel (XXDouble t);
void XXTerminateSubmodel (XXDouble t);

#endif
