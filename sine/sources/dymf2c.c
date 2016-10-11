/* dymf2c.c

   Files needed for c-files transformed from Fortran to C
   with the f2c compiler

   Author : Martin Otter, DLR.
   Release: 1997-11-26: implemented and tested
*/
/*
 * Copyright (C) 1997-2001 Dynasim AB.
 * All rights reserved.
 *
 */

/* DLL-Status: No static data used/defined in file */
#if !defined(NO_FILE)
#include <stdio.h>
#endif
#include <math.h>
#include "dymutil.h"
#include "f2c.h"
#include "localeless.h"
#include "libdssetup.h"

#define LEN_TEXTLINE 200

#if !defined(RT)
LIBDS_API void aprintLine(int type, const char *str);
#endif
static int CopyString(const char message[], int lstr,
                       char *textline, int lstrMax) {

   /* copy Fortran string "message" into textline

      -> message   Message as string (not NULL terminated)
      -> lstr      Length of message
      <- textline  String in which message has to be copied
      -> lstrMax   Length of textline
      <- RETURN    textline is filled upto textline[ic-1],
                   where ic is returned.
   */
      int  i, ic;

   /* remove trailing blanks in message */
      ic = Dymola_min(lstr, lstrMax-2);

   /* copy message into textline */
      for (i=0; i<ic && message[i]!='\0'; i++) {
         textline[i] = message[i];
      }
	  ic=i;
      while ( ic > 1 && message[ic-1] == ' ' ) --ic;
      textline[ic] = '\0';
      return ic;
}
      
DYMOLA_STATIC void DymosimMessageSeverity(int severity, const char*message) {
#if defined(RT) || defined(DynSimStruct)
	DymosimMessage(message);
#else
	if (severity==0) {
		DymosimMessage(message);
	} else {
		/* Skipping pre-amble */
		aprintLine(severity, message);
	}
#endif
}
DYMOLA_STATIC void dymosimmessageSev_(int severity, const char message[], int lstr) {
   /* Print Fortran string */
      char textline[LEN_TEXTLINE];

      CopyString(message, lstr, textline, LEN_TEXTLINE);

      DymosimMessageSeverity(severity, textline);
}
DYMOLA_STATIC void dymosimmessage_(const char message[], int lstr) {
	dymosimmessageSev_(0, message, lstr);
}


DYMOLA_STATIC void dymosimmessagedoubleSev_(int severity, const char message[], const double *d, int lstr) {

   /* Print Fortran string and double number */
      char textline[LEN_TEXTLINE];
      int ic;

      ic = CopyString(message, lstr, textline, LEN_TEXTLINE-20);
      sprintf(&textline[ic], " %G", *d);

      DymosimMessageSeverity(severity, textline);
}
DYMOLA_STATIC void dymosimmessagedouble_(const char message[], const double *d, int lstr) {
	dymosimmessagedoubleSev_(0, message, d, lstr);
}


DYMOLA_STATIC void dymosimmessagedouble2Sev_(int severity, const char message[], const double *d, const char message2[], 
                            int lstr, int lstr2) {

   /* Print Fortran string, double number, Fortran string */
      char textline[LEN_TEXTLINE];
      int ic;

      ic = CopyString(message, lstr, textline, LEN_TEXTLINE-20);
      ic = ic + sprintf(&textline[ic], " %G", *d);
      ic = CopyString(message2, lstr2, &textline[ic], LEN_TEXTLINE-ic);
      DymosimMessageSeverity(severity, textline);
}
DYMOLA_STATIC void dymosimmessagedouble2_(const char message[], const double *d, const char message2[], 
                            int lstr, int lstr2) {
	dymosimmessagedouble2Sev_(0, message, d, message2, lstr, lstr2);
}



DYMOLA_STATIC void dymosimmessagedoublesSev_(int severity, const char message[], const double d[], const int *nd, int lstr) {

   /* Print Fortran string and double vector */
      char textline[LEN_TEXTLINE];
      int ic, i;
	  int onLine=0;

      ic = CopyString(message, lstr, textline, LEN_TEXTLINE);
	  ic = ic+CopyString("{", 1L, textline+ic, LEN_TEXTLINE-ic);
      for (i=0; i<*nd; i++) {
		  if (i>0) 
			  ic = ic+CopyString(",", 1L, textline+ic, LEN_TEXTLINE-ic);
		  if (onLine>=5) {
			  DymosimMessageSeverity(severity, textline);
			  ic = 0;
			  ic = ic+CopyString(" ", 1L, textline, LEN_TEXTLINE);
			  onLine = 0;
		  }
          if ( ic+20 >= LEN_TEXTLINE ) {
			  break; /* Error stop*/
		  }

          ic = ic + sprintf(&textline[ic], " %G", d[i]);
		  ++onLine;
      }
	  ic = ic+CopyString(" }", 2L, textline+ic, LEN_TEXTLINE-ic);

      DymosimMessageSeverity(severity, textline);
}
DYMOLA_STATIC void dymosimmessagedoubles_(const char message[], const double d[], const int *nd, int lstr) {
	dymosimmessagedoublesSev_(0, message, d, nd, lstr);
}


DYMOLA_STATIC void DymosimMessageIntSev_(int severity, const char* message, const int *i, int lstr) {
   /* Print Fortran string and int number */
      char textline[LEN_TEXTLINE];
      int ic;

      ic = CopyString(message, lstr, textline, LEN_TEXTLINE-15);
      sprintf(&textline[ic], " %ld", *i);

      DymosimMessageSeverity(severity, textline);
}
DYMOLA_STATIC void DymosimMessageInt_(const char* message, const int *i, int lstr) {
	DymosimMessageIntSev_(0, message, i, lstr);
}


/* ------------------------------------------------------------------------- */
/* functions needed from f2c-compiled code                                   */
/* ------------------------------------------------------------------------- */

#ifndef GODESS
DYMOLA_STATIC double d_sign(const doublereal *a, const doublereal *b)
{
double x;
x = (*a >= 0 ? *a : - *a);
return( *b >= 0 ? x : -x);
}

#endif


DYMOLA_STATIC integer i_sign(const integer *a, const integer *b)
{
integer x;
x = (*a >= 0 ? *a : - *a);
return( *b >= 0 ? x : -x);
}


#ifndef GODESS
DYMOLA_STATIC double pow_di(const doublereal *ap, const integer *bp)
{
double pow, x;
integer n;

pow = 1;
x = *ap;
n = *bp;

if(n != 0)
	{
	if(n < 0)
		{
		n = -n;
		x = 1/x;
		}
	for( ; ; )
		{
		if(n & 01)
			pow *= x;
		if(n >>= 1)
			x *= x;
		else
			break;
		}
	}
return(pow);
}
#endif


DYMOLA_STATIC integer pow_ii(const integer *ap, const integer *bp)
{
	integer pow, x, n;

	x = *ap;
	n = *bp;

	if (n <= 0) {
		if (n == 0 || x == 1)
			return 1;
		if (x != -1)
			return x == 0 ? 1/x : 0;
		n = -n;
		}
	for(pow = 1; ; )
		{
		if(n & 01)
			pow *= x;
		if(n >>= 1)
			x *= x;
		else
			break;
		}
	return(pow);
	}



DYMOLA_STATIC integer i_indx(const char *a, const char *b, ftnlen la, ftnlen lb)
{
ftnlen i, n;
const char *s, *t, *bend;

n = la - lb + 1;
bend = b + lb;

for(i = 0 ; i < n ; ++i)
	{
	s = a + i;
	t = b;
	while(t < bend)
		if(*s++ != *t++)
			goto no;
	return(i+1);
	no: ;
	}
return(0);
}


DYMOLA_STATIC integer i_len(const char *s, ftnlen n)
{
return(n);
}

#ifndef GODESS
DYMOLA_STATIC int s_cat(char *lp, const char *const rpp[], const ftnlen rnp[], const ftnlen *np, ftnlen ll)
{
ftnlen i, n, nc;
const char *f__rp;

n = (int)*np;
for(i = 0 ; i < n ; ++i)
	{
	nc = ll;
	if(rnp[i] < nc)
		nc = rnp[i];
	ll -= nc;
	f__rp = rpp[i];
	while(--nc >= 0)
		*lp++ = *f__rp++;
	}
while(--ll >= 0)
	*lp++ = ' ';
	return 0;
}



/* compare two strings */
DYMOLA_STATIC integer s_cmp(register const char *a, register const char *b, ftnlen la, ftnlen lb)
{
register const char *aend, *bend;
aend = a + la;
bend = b + lb;

if(la <= lb)
	{
	while(a < aend)
		if(*a != *b)
			return( *a - *b );
		else
			{ ++a; ++b; }

	while(b < bend)
		if(*b != ' ')
			return( ' ' - *b );
		else	++b;
	}

else
	{
	while(b < bend)
		if(*a == *b)
			{ ++a; ++b; }
		else
			return( *a - *b );
	while(a < aend)
		if(*a != ' ')
			return(*a - ' ');
		else	++a;
	}
return(0);
}



/* assign strings:  a = b */
DYMOLA_STATIC int s_copy(register char *a, register const char *b, ftnlen la, ftnlen lb)
{
register char *aend;
register const char *bend;

aend = a + la;

if(la <= lb)
	while(a < aend)
		*a++ = *b++;

else
	{
	bend = b + lb;
	while(b < bend)
		*a++ = *b++;
	while(a < aend)
		*a++ = ' ';
	}
 return 0;
}
#endif


DYMOLA_STATIC double pow_dd(const doublereal *ap, const doublereal *bp)
{
return(pow(*ap, *bp) );
}



#ifndef GODESS
#define log10e 0.43429448190325182765

DYMOLA_STATIC double d_lg10(const doublereal *x)
{
return( log10e * log(*x) );
}


DYMOLA_STATIC integer i_dnnt(const doublereal *x)
{
return (integer) ( (*x)>=0 ?
	floor(*x + .5) : -floor(.5 - *x) );
}
#endif

#if defined(RT)

/* Special wrapper functions of certain Lapack functions for code export */
DYMOLA_STATIC int dgetrf_(const integer*m,const integer*n,doublereal *a,const integer*lda,integer *ipiv,integer *info);
DYMOLA_STATIC int dgetrs_(const char*trans,const integer*n,const integer*nrhs,doublereal *a,const integer*lda,integer *ipiv,doublereal *b,const integer*ldb,integer *info, 
	ftnlen trans_len);
DYMOLA_STATIC int dgeqpf_(const integer*m,const integer*n,doublereal *a,const integer*lda,integer *jpvt,doublereal *tau,doublereal *work,integer *info);
DYMOLA_STATIC doublereal dlange_(char *norm,const integer*m,const integer*n,doublereal *a,const integer*lda,doublereal *work,ftnlen norm_len);
DYMOLA_STATIC int dgecon_(char *norm,const integer*n,doublereal *a,const integer*lda,doublereal *anorm,doublereal *rcond,doublereal *work,integer* iwork,integer* info,
	 ftnlen norm_len);
DYMOLA_STATIC int dgemm_(const char*transa,const char*transb,const integer*m,const integer*n,const integer * k,const doublereal*alpha,doublereal * a,const integer*lda,doublereal * b,const integer*ldb, 
	const doublereal*beta,doublereal * c__, const integer*ldc, ftnlen transa_len, ftnlen transb_len);
DYMOLA_STATIC int dtrsm_(const char*side,const char*uplo,const char*transa,const char*diag,const integer*m,const integer*n,const doublereal*alpha,doublereal*a,const integer*lda, 
		   doublereal*b,const integer*ldb,ftnlen side_len,ftnlen uplo_len,ftnlen transa_len,ftnlen diag_len);
DYMOLA_STATIC int dtrsv_(const char*uplo,const char*trans,const char * diag,const integer*n,doublereal * a,const integer*lda,doublereal *x,const integer*incx, ftnlen uplo_len, 
	ftnlen trans_len, ftnlen diag_len);
DYMOLA_STATIC logical (lsame_)(const char*ca,const char*cb,ftnlen ca_len,ftnlen cb_len);

DYMOLA_STATIC int dgetrf__RT(const integer*m,const integer*n,doublereal *a,const integer*lda,integer *ipiv,integer *info)
{	
	return dgetrf_(m,n,a,lda,ipiv,info);
}
DYMOLA_STATIC int dgetrs__RT(const char*trans,const integer*n,const integer*nrhs,doublereal *a,const integer*lda,integer *ipiv,doublereal *b,const integer*ldb,integer *info)
{
	return dgetrs_(trans,n,nrhs,a,lda,ipiv,b,ldb,info,0);
}
DYMOLA_STATIC int dgeqpf__RT(const integer*m,const integer*n,doublereal *a,const integer*lda,integer *jpvt,doublereal *tau,doublereal *work,integer *info)
{
	return dgeqpf_(m,n,a,lda,jpvt,tau,work,info);
}
DYMOLA_STATIC doublereal dlange__RT(char *norm,const integer*m,const integer*n,doublereal *a,const integer*lda,doublereal *work)
{
	return dlange_(norm,m,n,a,lda,work,0);
}
DYMOLA_STATIC int dgecon__RT(char *norm,const integer*n,doublereal *a,const integer*lda,doublereal *anorm,doublereal *rcond,doublereal *work,integer* iwork,integer* info)
{
	return dgecon_(norm,n,a,lda,anorm,rcond,work,iwork,info,0);
}
DYMOLA_STATIC int dgemm__RT(const char*transa,const char*transb,const integer*m,const integer*n,const integer * k,const doublereal*alpha,doublereal * a,const integer*lda,doublereal * b,const integer*ldb, 
	const doublereal*beta,doublereal * c__, const integer*ldc)
{
	return dgemm_(transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c__,ldc,0,0);
}
DYMOLA_STATIC int dtrsm__RT(const char*side,const char*uplo,const char*transa,const char*diag,const integer*m,const integer*n,const doublereal*alpha,doublereal*a,const integer*lda, 
		   doublereal*b,const integer*ldb)
{
	return dtrsm_(side,uplo,transa,diag,m,n,alpha,a,lda,b,ldb,0,0,0,0);
}
DYMOLA_STATIC int dtrsv__RT(const char*uplo,const char*trans,const char * diag,const integer*n,doublereal * a,const integer*lda,doublereal *x,const integer*incx)
{
	return dtrsv_(uplo,trans,diag,n,a,lda,x,incx,0,0,0);
}

DYMOLA_STATIC logical lsame__RT(const char*ca,const char*cb)
{
	return lsame_(ca, cb, 0, 0);
}

#endif /* defined(RT) */
