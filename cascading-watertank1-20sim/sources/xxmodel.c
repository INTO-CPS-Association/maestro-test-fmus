/**********************************************************
 * This file is generated by 20-sim ANSI-C Code Generator
 *
 *  file:  src\xxmodel.c
 *  model: WaterTankCascade
 *  expmt: WaterTankCascade
 *  date:  October 30, 2015
 *  time:  12:29:45 pm
 *  user:  INTO-CPS
 *  from:  20-sim 4.5 Professional Single
 *  build: 4.5.1.5561
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

/* the global variables */
XXDouble xx_start_time = 0.0;
XXDouble xx_finish_time = 12.0;
XXDouble xx_step_size = 0.01;
XXDouble xx_time = 0.0;
XXInteger xx_steps = 0;
XXBoolean xx_initialize = XXTRUE;
XXBoolean xx_major = XXTRUE;
XXBoolean xx_stop_simulation = XXFALSE;

/* the variable arrays */
XXDouble xx_MEMORY[0 + 6 + 1 + 16 + 1 + 1 + 1];
XXDouble* xx_C = xx_MEMORY;		/* constants */
XXDouble* xx_P = xx_MEMORY + 0;			/* parameters */
XXDouble* xx_I = xx_MEMORY + 0 + 6;		/* initial values */
XXDouble* xx_V = xx_MEMORY + 0 + 6 + 1;		/* variables */
XXDouble* xx_s = xx_MEMORY + 0 + 6 + 1 + 16;		/* states */
XXDouble *xx_R = xx_MEMORY + 0 + 6 + 1 + 16 + 1;		/* rates (or new states) */

/* the names of the variables as used in the arrays above
   uncomment this part if these names are needed
XXCharacter *xx_parameter_names[] = {
	"FlowSource\\phi",
	"tank1\\Tank\\area",
	"tank1\\Tank\\gravity",
	"tank1\\Tank\\liquid_density",
	"valve1\\R\\r",
	"valve1\\Valve\\leakage_value"
,	NULL
};
XXCharacter *xx_initial_value_names[] = {
	"tank1\\Tank\\volume_initial"
,	NULL
};
XXCharacter *xx_variable_names[] = {
	"FlowSource\\p.phi",
	"FlowSource\\inflow",
	"tank1\\Tank\\p.e",
	"tank1\\Tank\\level",
	"valve1\\R\\p.f",
	"valve1\\R\\waterFlow",
	"valve1\\Valve\\powerIn.phi",
	"valve1\\Valve\\powerOut.p",
	"valve1\\waterOut.e",
	"valve1\\valveControl",
	"valve1\\leakage",
	"waterOut.e",
	"waterOut.f",
	"valveControl",
	"leakage",
	"level"
,	NULL
};
XXCharacter *xx_state_names[] = {
	"tank1\\Tank\\volume"
,	NULL
};
XXCharacter *xx_rate_names[] = {
	"tank1\\Tank\\p.f"
,	NULL
};
*/

/* this method is called before calculation is possible */
void XXModelInitialize (void)
{
	/* set the parameters */
	xx_P[0] = 1.0;		/* FlowSource\phi {m3/s} */
	xx_P[1] = 1.0;		/* tank1\Tank\area */
	xx_P[2] = 9.81;		/* tank1\Tank\gravity */
	xx_P[3] = 1.0;		/* tank1\Tank\liquid_density */
	xx_P[4] = 5.0;		/* valve1\R\r */
	xx_P[5] = 0.1;		/* valve1\Valve\leakage_value */


	/* set the initial values */
	xx_I[0] = 0.0;		/* tank1\Tank\volume_initial */


	/* set the states */
	xx_s[0] = xx_I[0];		/* tank1\Tank\volume */


}

/* This function calculates the initial equations of the model.
 * These equations are calculated before anything else
 */
void XXCalculateInitial (void)
{


}

/* This function calculates the static equations of the model.
 * These equations are only dependent from parameters and constants
 */
void XXCalculateStatic (void)
{
	/* FlowSource\p.phi = FlowSource\phi; */
	xx_V[0] = xx_P[0];

	/* FlowSource\inflow = if FlowSource\p.phi < 0... ; */
	xx_V[1] = (xx_V[0] < 0.0) ? 
		/* 0 */
		0.0
	:
		/* 0.1 * abs (FlowSource\p.phi) */
		(0.1 * XXAbsolute (xx_V[0]))
	;

}

/* This function calculates the input equations of the model.
 * These equations are dynamic equations that must not change
 * in calls from the integration method (like random and delay).
 */
void XXCalculateInput (void)
{

}

/* This function calculates the dynamic equations of the model.
 * These equations are called from the integration method
 * to calculate the new model rates (that are then integrated).
 */
void XXCalculateDynamic (void)
{
	/* tank1\Tank\p.e = (tank1\Tank\gravity * tank1\Tank\volume) * tank1\Tank\liquid_density; */
	xx_V[2] = (xx_P[2] * xx_s[0]) * xx_P[3];

	/* tank1\Tank\level = tank1\Tank\volume / tank1\Tank\area; */
	xx_V[3] = xx_s[0] / xx_P[1];

	/* valve1\valveControl = valveControl; */
	xx_V[9] = xx_V[13];

	/* valve1\leakage = leakage; */
	xx_V[10] = xx_V[14];

	/* level = tank1\Tank\level; */
	xx_V[15] = xx_V[3];

	/* valve1\Valve\powerOut.p = if (valve1\valveControl <> 0.0)... ; */
	xx_V[7] = (xx_V[9] != 0.0) ? 
		/* tank1\Tank\p.e */
		xx_V[2]
	:
		/* if valve1\leakage...  */
		(xx_V[10] ? 
			/* (tank1\Tank\p.e * valve1\Valve\leakage_value) */
			(xx_V[2] * xx_P[5])
		:
			/* 0.0 */
			0.0
		)
	;

	/* valve1\R\p.f = valve1\Valve\powerOut.p / valve1\R\r; */
	xx_V[4] = xx_V[7] / xx_P[4];

	/* valve1\Valve\powerIn.phi = if (valve1\valveControl <> 0.0)... ; */
	xx_V[6] = (xx_V[9] != 0.0) ? 
		/* valve1\R\p.f */
		xx_V[4]
	:
		/* if valve1\leakage...  */
		(xx_V[10] ? 
			/* (valve1\R\p.f * valve1\Valve\leakage_value) */
			(xx_V[4] * xx_P[5])
		:
			/* 0.0 */
			0.0
		)
	;

	/* waterOut.f = valve1\R\p.f; */
	xx_V[12] = xx_V[4];

	/* tank1\Tank\p.f = FlowSource\p.phi - valve1\Valve\powerIn.phi; */
	xx_R[0] = xx_V[0] - xx_V[6];


	/* increment the step counter */
	xx_steps++;
}

/* This function calculates the output equations of the model.
 * These equations are not needed for calculation of the rates
 * and are kept separate to make the dynamic set of equations smaller.
 * These dynamic equations are called often more than one time for each
 * integration step that is taken. This makes model computation much faster.
 */
void XXCalculateOutput (void)
{
	/* valve1\waterOut.e = waterOut.e; */
	xx_V[8] = xx_V[11];

	/* valve1\R\waterFlow = 0.1 * valve1\R\p.f; */
	xx_V[5] = 0.1 * xx_V[4];

}

/* This function calculates the final equations of the model.
 * These equations are calculated after all the calculations
 * are performed
 */
void XXCalculateFinal (void)
{

}

/* this method is called after all calculations are performed */
void XXModelTerminate(void)
{
}


