/**********************************************************
 * This file is generated by 20-sim ANSI-C Code Generator
 *
 *  file:  src\xxmodel.c
 *  model: MassSpringDamper
 *  expmt: MassSpringDamper
 *  date:  April 3, 2017
 *  time:  11:26:28 AM
 *  user:  INTO-CPS
 *  from:  20-sim 4.6 Professional Single
 *  build: 4.6.3.7711
 **********************************************************/

/* This file contains the actual model variables and equations */

/* Note: Alias variables are the result of full optimization
   of the model in 20-sim. As a result, only the real variables
   are used in the model for speed. The user may also include
   the alias variables by adding them to the end of the array:

   XXDouble xx_variables[NUMBER_VARIABLES + NUMBER_ALIAS_VARIABLES + 1];
   XXCharacter *xx_variable_names[] = {
     VARIABLE_NAMES, ALIAS_VARIABLE_NAMES, NULL
   };

   and calculate them directly after the output equations:

   void XXCalculateOutput (void)
   {
     OUTPUT_EQUATIONS
     ALIAS_EQUATIONS
   }
*/

/* system include files */
#include <stdlib.h>
#include <math.h>

/* 20-sim include files */
#include "xxmodel.h"
#include "xxfuncs.h"


#if (3 > 8192) && defined _MSC_VER
#pragma optimize("", off)
#endif
/* this method is called before calculation is possible */
void XXModelInitialize (xx_ModelInstance* model_instance)
{
	/* set the parameters */
	xx_P[0] = 1.0;		/* c1 */
	xx_P[1] = 1.0;		/* d1 */
	xx_P[2] = 1.0;		/* m1 */


	/* set the initial values */
	xx_I[1] = 1.0;		/* x1_initial */
	xx_I[0] = 0.0;		/* v1_initial */


	/* set the states */
	xx_s[1] = xx_I[1];		/* x1 */
	xx_s[0] = xx_I[0];		/* v1 */


}
#if (3 > 8192) && defined _MSC_VER
#pragma optimize("", on)
#endif

/* This function calculates the initial equations of the model.
 * These equations are calculated before anything else
 */
void XXCalculateInitial (xx_ModelInstance* model_instance)
{

}

/* This function calculates the static equations of the model.
 * These equations are only dependent from parameters and constants
 */
void XXCalculateStatic (xx_ModelInstance* model_instance)
{

}

/* This function calculates the input equations of the model.
 * These equations are dynamic equations that must not change
 * in calls from the integration method (like random and delay).
 */
void XXCalculateInput (xx_ModelInstance* model_instance)
{

}

/* This function calculates the dynamic equations of the model.
 * These equations are called from the integration method
 * to calculate the new model rates (that are then integrated).
 */
void XXCalculateDynamic (xx_ModelInstance* model_instance)
{
	/* additional code for v1; */
	xx_R[1] = xx_s[0];

	/* v1_dot = (1 / m1) * (((-c1 * x1) - d1 * v1) + fk); */
	xx_R[0] = (1.0 / xx_P[2]) * (((-xx_P[0] * xx_s[1]) - xx_P[1] * xx_s[0]) + xx_V[0]);


	/* increment the step counter */
	model_instance->steps++;
}

/* This function calculates the output equations of the model.
 * These equations are not needed for calculation of the rates
 * and are kept separate to make the dynamic set of equations smaller.
 * These dynamic equations are called often more than one time for each
 * integration step that is taken. This makes model computation much faster.
 */
void XXCalculateOutput (xx_ModelInstance* model_instance)
{

}

/* This function calculates the final equations of the model.
 * These equations are calculated after all the calculations
 * are performed
 */
void XXCalculateFinal (xx_ModelInstance* model_instance)
{

}

/* this method is called after all calculations are performed */
void XXModelTerminate(xx_ModelInstance* model_instance)
{
}


