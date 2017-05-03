/**********************************************************
 * This file is generated by 20-sim ANSI-C Code Generator
 *
 *  file:  src\xxfuncs.c
 *  model: MassSpringDamper
 *  expmt: MassSpringDamper
 *  date:  April 3, 2017
 *  time:  11:26:43 AM
 *  user:  INTO-CPS
 *  from:  20-sim 4.6 Professional Single
 *  build: 4.6.3.7711
 **********************************************************/

/* This file contains support functions for several SIDOPS functions

   For flexibility, ANSI-C is created, and typedefs are used
   for integers and doubles, see the xxfuncs.h file for more
   information on these types.

   This means that all used functions follow the ANSI definition.

   Please check the math.h file of your particular compiler
   to see if this is indeed the case. Otherwise, you might have
   to adapt the used functions below to obtain the same behavior.

*/

/* The system include files */
#include <stdlib.h>
#include <math.h>

/* Our own include files */
#include "xxfuncs.h"

/* Constants that are used in our functions below */

typedef union
{
	double m_double;
	const char* m_char;
}str2dbl;

XXDouble XXString2Double(const char* argument)
{
	str2dbl myConversion;
	myConversion.m_char = argument;
	return myConversion.m_double;

}

const char* XXDouble2String(XXDouble argument)
{
	str2dbl myConversion;
	myConversion.m_double = argument;
	return myConversion.m_char;
}

/* The 20-sim SIDOPS functions */
/* 20-sim stubs. Implement them yourself if needed */

