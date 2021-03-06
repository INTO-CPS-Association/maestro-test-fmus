#ifndef MATRIXOP_H
#define MATRIXOP_H
/*  * Declarations of matrix functions in Dymola
 *
 * Copyright (C) 1997-2001 Dynasim AB.
 * All rights reserved.
 *
 * Author: Hans Olsson Dynasim AB, 1999
 * Version: 1.4, 1999-09-24*/
/* */
#include "matrixop1.h"
#include <math.h>

#if !defined(DSE_STRUCT)
 #if defined(DS_EMBEDDED)
  #define DSE_STRUCT sts->
 #else
  #define DSE_STRUCT
 #endif
#endif
#ifndef DYNX
#if _MSC_VER>=1700
#define DYNX(s,i) (*(s+i))
#else
#define DYNX(s,i) s[i]
#endif
#endif

#if defined(DS_EMBEDDED)
#include <dsembedded.h>
#include <dse_dymosim.c>
#endif

#if defined(Matlab5) || defined(Matlab51) || defined(SimStruct) || defined(DYM2DS) || defined(FMU_SOURCE_CODE_EXPORT)
#if !defined(DynSimStruct)
#define DynSimStruct 1
/* Needed for dsutil.h in Modelica mode. */
#endif
#endif

#ifndef DYNassert_static
#define DYNassert_static(e) do { enum { DYN_assert_static_=1/(e) }; } while(0)
#endif

#define Modelica

#include <userdefs.h>
#include <stddef.h>
#include <stdarg.h>
#include <math.h>

#include "dsblock.h"
#include "dsutil.h"
#include "ModelicaUtilities.h"

/* Also present in moutil.c*/
#define IF   (
#define THEN ) ? (
#define ELSE ) :
#define AND &&
#define OR ||
#define NOT !
#define true 1
#define false 0

#if !defined(NDEBUG) && (defined(DYMOSIM) || defined(BUILDFMU))
#undef Assert
#undef AssertModelica
#define Assert(b,x) ((!(b))?AssertModelicaF(b,#b,x):(void)(0))
#define AssertModelica(b,bs,x) ((!(b))?AssertModelicaF(b,bs,x):(void)(0))
#define AssertModelica3(b,bs,x,l) ((!(b))?AssertModelicaF3(b,bs,x,l):(void)(0))
#else
#undef Assert
#undef AssertModelica
#define Assert(b,x)
#define AssertModelica(b,bs,x)
#define AssertModelica3(b,bs,x,l)
#endif

#if !defined(RT) && !defined(DYN_MULTINSTANCE)
Real Realbuffer[
#ifndef DYNREALBUFFER
	7000000
#else
	DYNREALBUFFER
#endif
];
Integer Integerbuffer[50000];
SizeType Sizebuffer[50000];
String Stringbuffer[10000];
char simplestring[100000];
#endif

#if !defined(RT)
Real RealbufferNon[
#ifndef DYNREALBUFFER
	70000
#else
	DYNREALBUFFER/10
#endif
];
Integer IntegerbufferNon[50000];
SizeType SizebufferNon[50000];
String StringbufferNon[10000];
char simplestringNon[1000];
DYMOLA_STATIC Real* EndRealbufferNon=RealbufferNon+sizeof(RealbufferNon)/sizeof(*RealbufferNon);
DYMOLA_STATIC Integer* EndIntegerbufferNon=IntegerbufferNon+sizeof(IntegerbufferNon)/sizeof(*IntegerbufferNon);
DYMOLA_STATIC SizeType* EndSizebufferNon=SizebufferNon+sizeof(SizebufferNon)/sizeof(*SizebufferNon);
DYMOLA_STATIC String* EndStringbufferNon=StringbufferNon+sizeof(StringbufferNon)/sizeof(*StringbufferNon);
DYMOLA_STATIC char* EndsimplestringNon=simplestringNon+sizeof(simplestringNon)/sizeof(*simplestringNon);
#endif

#if !defined(RT) && !defined(DYN_MULTINSTANCE)
#if (defined(_OPENMP) && !defined(DISABLE_DYMOLA_OPENMP))
DYMOLA_STATIC Real* EndRealbuffer=0;
#pragma omp threadprivate(EndRealbuffer)
DYMOLA_STATIC Integer* EndIntegerbuffer=0;
#pragma omp threadprivate(EndIntegerbuffer)
DYMOLA_STATIC SizeType* EndSizebuffer=0;
#pragma omp threadprivate(EndSizebuffer)
DYMOLA_STATIC String* EndStringbuffer=0;
#pragma omp threadprivate(EndStringbuffer)
DYMOLA_STATIC char* Endsimplestring=0;
#pragma omp threadprivate(Endsimplestring)
DYMOLA_STATIC Real* EndRealbuffer2=Realbuffer+sizeof(Realbuffer)/sizeof(*Realbuffer);
DYMOLA_STATIC Integer* EndIntegerbuffer2=Integerbuffer+sizeof(Integerbuffer)/sizeof(*Integerbuffer);
DYMOLA_STATIC SizeType* EndSizebuffer2=Sizebuffer+sizeof(Sizebuffer)/sizeof(*Sizebuffer);
DYMOLA_STATIC String* EndStringbuffer2=Stringbuffer+sizeof(Stringbuffer)/sizeof(*Stringbuffer);
DYMOLA_STATIC char* Endsimplestring2=simplestring+sizeof(simplestring)/sizeof(*simplestring);
#else
DYMOLA_STATIC Real* EndRealbuffer=Realbuffer+sizeof(Realbuffer)/sizeof(*Realbuffer);
DYMOLA_STATIC Integer* EndIntegerbuffer=Integerbuffer+sizeof(Integerbuffer)/sizeof(*Integerbuffer);
DYMOLA_STATIC SizeType* EndSizebuffer=Sizebuffer+sizeof(Sizebuffer)/sizeof(*Sizebuffer);
DYMOLA_STATIC String* EndStringbuffer=Stringbuffer+sizeof(Stringbuffer)/sizeof(*Stringbuffer);
DYMOLA_STATIC char* Endsimplestring=simplestring+sizeof(simplestring)/sizeof(*simplestring);
#endif
#endif

DYMOLA_STATIC void* RecordElement(const RealArray a,...) {
  SizeType index;
  va_list ap;
  va_start(ap,a);
  index=FindIndex(a.ndims,a.dims,ap);
  va_end(ap);
  return (void*)(a.data+index);
}
DYMOLA_STATIC Integer RecordSize(const RealArray a,SizeType i) {
	return RealSize(a,i);
}
DYMOLA_STATIC RealArray RecordLeading(int i,int j,RealArray a) {
	return RealLeading(i,j,a);
}
DYMOLA_STATIC IntegerArray RecordSizes(const RealArray a) {
	IntegerArray b;
	b=RealSizes(a);
	Assert(b.ndims>0 && b.dims[0]>0,"Error in record logic");
	b.dims[0]-=1;
	return b;
}
DYMOLA_STATIC Real *RealTemp(SizeType r);
DYMOLA_STATIC Integer *IntegerTemp(SizeType r);
DYMOLA_STATIC SizeType *SizeTemp(SizeType r);
DYMOLA_STATIC String * StringTemp(SizeType r);
DYMOLA_STATIC Integer RealMatchingSizes(const RealArray a,const RealArray b);
DYMOLA_STATIC Integer IntegerMatchingSizes(const IntegerArray a,const IntegerArray b);
DYMOLA_STATIC Integer StringMatchingSizes(const StringArray a,const StringArray b);
DYMOLA_STATIC RealArray RecordArrayRealSlice(const RealArray a,size_t off) {
  RealArray temp;
  SizeType i,len;
  temp.ndims=a.ndims-1;
  temp.dims=SizeTemp(temp.ndims);
  for(i=0;i<temp.ndims;++i) temp.dims[i]=a.dims[i];
  len=RealNrElements(temp);
  temp.data=RealTemp(len);
  for(i=0;i<len;++i) temp.data[i]=*((Real*)(((char*)(a.data+i*a.dims[a.ndims-1]))+off));
  return temp;
}
DYMOLA_STATIC String Enum2String2(Integer x,Integer leftJustified,Integer minwidth,StringArray names) {
	char buf[100];
	char*ret;
	Assert(minwidth<40,"");
	Assert(names.ndims==1 && x>=1 && x<=names.dims[0],"");
	sprintf(buf,leftJustified?"%-*s":"%*s",(int)minwidth,names.data[x-1]);
	ret=StringAllocate(strlen(buf));
	strcpy(ret,buf);
	return ret;
}
DYMOLA_STATIC IntegerArray RecordArrayIntegerSlice(const RealArray a,size_t off) {
  IntegerArray temp;
  SizeType i,len;
  temp.ndims=a.ndims-1;
  temp.dims=SizeTemp(temp.ndims);
  for(i=0;i<temp.ndims;++i) temp.dims[i]=a.dims[i];
  len=IntegerNrElements(temp);
  temp.data=IntegerTemp(len);
  for(i=0;i<len;++i) temp.data[i]=*((Integer*)(((char*)(a.data+i*a.dims[a.ndims-1]))+off));
  return temp;
}
DYMOLA_STATIC StringArray RecordArrayStringSlice(const RealArray a,size_t off) {
  StringArray temp;
  SizeType i,len;
  temp.ndims=a.ndims-1;
  temp.dims=SizeTemp(temp.ndims);
  for(i=0;i<temp.ndims;++i) temp.dims[i]=a.dims[i];
  len=StringNrElements(temp);
  temp.data=StringTemp(len);
  for(i=0;i<len;++i) temp.data[i]=*((String*)(((char*)(a.data+i*a.dims[a.ndims-1]))+off));
  return temp;
}

DYMOLA_STATIC RealArray RecordArrayRealArraySlice(const RealArray a,size_t off) {
  RealArray temp;
  RealArray elem,elemi;
  SizeType i,j,lenA,lenTemp,lenElem;
  Assert(RealNrElements(a)>0,"Can only make record slices of non-empty arrays");
  elem=*(RealArray*)((char*)(a.data+0)+off);
  temp.ndims=a.ndims-1+elem.ndims;
  temp.dims=SizeTemp(temp.ndims);
  lenElem=RealNrElements(elem);
  for(i=0;i<temp.ndims;++i) if (i<a.ndims-1) temp.dims[i]=a.dims[i]; else temp.dims[i]=elem.dims[i-(a.ndims-1)];
  lenTemp=RealNrElements(temp);
  lenA=RealNrElements(a)/a.dims[a.ndims-1];
  Assert(lenTemp==lenA*lenElem,"Consistency");
  temp.data=RealTemp(lenTemp);
  for(i=0;i<lenA;++i) {
	  elemi=*(RealArray*)((char*)(a.data+i*a.dims[a.ndims-1])+off);
	  Assert(RealMatchingSizes(elemi,elem),"Can only slice homogenous record arrays");
	  for(j=0;j<lenElem;++j) temp.data[i*lenElem+j]=elemi.data[j];
  }
  return temp;
}
DYMOLA_STATIC IntegerArray RecordArrayIntegerArraySlice(const RealArray a,size_t off) {
  IntegerArray temp;
  IntegerArray elem,elemi;
  SizeType i,j,lenA,lenTemp,lenElem;
  Assert(RealNrElements(a)>0,"Can only make record slices of non-empty arrays");
  elem=*(IntegerArray*)((char*)(a.data+0)+off);
  temp.ndims=a.ndims-1+elem.ndims;
  temp.dims=SizeTemp(temp.ndims);
  lenElem=IntegerNrElements(elem);
  for(i=0;i<temp.ndims;++i) if (i<a.ndims-1) temp.dims[i]=a.dims[i]; else temp.dims[i]=elem.dims[i-(a.ndims-1)];
  lenTemp=IntegerNrElements(temp);
  lenA=RealNrElements(a)/a.dims[a.ndims-1];
  Assert(lenTemp==lenA*lenElem,"Consistency");
  temp.data=IntegerTemp(lenTemp);
  for(i=0;i<lenA;++i) {
	  elemi=*(IntegerArray*)((char*)(a.data+i*a.dims[a.ndims-1])+off);
	  Assert(IntegerMatchingSizes(elemi,elem),"Can only slice homogenous record arrays");
	  for(j=0;j<lenElem;++j) temp.data[i*lenElem+j]=elemi.data[j];
  }
  return temp;
}
DYMOLA_STATIC StringArray RecordArrayStringArraySlice(const RealArray a,size_t off) {
  StringArray temp;
  StringArray elem,elemi;
  SizeType i,j,lenA,lenTemp,lenElem;
  Assert(RealNrElements(a)>0,"Can only make record slices of non-empty arrays");
  elem=*(StringArray*)((char*)(a.data+0)+off);
  temp.ndims=a.ndims-1+elem.ndims;
  temp.dims=SizeTemp(temp.ndims);
  lenElem=StringNrElements(elem);
  for(i=0;i<temp.ndims;++i) if (i<a.ndims-1) temp.dims[i]=a.dims[i]; else temp.dims[i]=elem.dims[i-(a.ndims-1)];
  lenTemp=StringNrElements(temp);
  lenA=RealNrElements(a)/a.dims[a.ndims-1];
  Assert(lenTemp==lenA*lenElem,"Consistency");
  temp.data=StringTemp(lenTemp);
  for(i=0;i<lenA;++i) {
	  elemi=*(StringArray*)((char*)(a.data+i*a.dims[a.ndims-1])+off);
	  Assert(StringMatchingSizes(elemi,elem),"Can only slice homogenous record arrays");
	  for(j=0;j<lenElem;++j) temp.data[i*lenElem+j]=elemi.data[j];
  }
  return temp;
}

static double invmodDymola(double x,double y) {double quot=x/y;return (floor(quot));}
/* static double invmodDymola(double x,double y) {double quot=x/y;return (quot>=0)?floor(quot):ceil(quot);} */

static double abs_inverse(double x,double y, const char*ys) {	
	Assert(x>=0, StringAdd("The inverse of the function abs is not defined for negative values.\nFor: ",ys));
	return (y>=0 ? x :-x);
}

static double powIntExp_inverse(double x, int exponent, double y, const char*ys) {
	double xabs;
	double yabs;
	double invExp;
	Assert(exponent!=0, StringAdd("The inverse of the function y^e is not unique for e = 0.\nFor: ", ys));
	invExp = 1.0/exponent;
	xabs = (x < 0.0 ? -x : x);
	if (exponent % 2 == 0) {
		/* even exponent */
		Assert(x >=0, StringAdd("The inverse of the function y^e is undefined for negative values if e is even.\nFor: ",ys));
		yabs = pow(xabs, invExp);
		return abs_inverse(yabs, y, ys);
	}
	yabs = pow(xabs, invExp);
	return (x < 0.0 ? -yabs : yabs);
}

static double sin_inverse(double x,double y, const char*ys) {
	double pi;
	double n;
	Assert(x >=-1 && x <=1, StringAdd("The inverse of the function sin is not defined for values greater than one.\nFor: ", ys));
	pi = 4*atan(1);
	n = invmodDymola((y+pi/2), 2*pi);
	return ((y - n*2*pi) > pi/2 ? pi - asin(x) : asin(x)) + n*2*pi;
}

static double cos_inverse(double x,double y, const char*ys) {
	double pi;
	double n;
	Assert(x >=-1 && x <=1, StringAdd("The inverse of the function cos is not defined for values greater than one.\nFor: ", ys));
	pi = 4*atan(1);
	n = invmodDymola(y, 2*pi);
	return (((y - n*2*pi) > pi) ? 2*pi - acos(x) : acos(x)) + n*2*pi;
}

static double tan_inverse(double x,double y, const char*ys) {
	double pi;
	double n;
	pi = 4*atan(1);
	n = invmodDymola((y+pi/2), pi);
	return atan(x) + n*pi;
}

static double semiLinear(double md,double p,double n) {
	if (md>0) return md*p;
	else return md*n;
}
DYMOLA_STATIC void terminate(const char*message) {
	if (*(GlobalErrorPointer())!=-999) {
		DymosimMessage("Simulation successfully terminated");
		DymosimMessage((char*)message);
		*(GlobalErrorPointer()) = -999;
	}
}
DYMOLA_STATIC void resetState(void) {
	*(GlobalErrorPointer()) = -998;
}




#define SetRealVectorElement(v,a,i) (DYNX((a).data,(i)-1)=(v))
#define SetIntegerVectorElement(v,a,i) (DYNX((a).data,(i)-1)=(v))
#define SetStringVectorElement(v,a,i) (DYNX((a).data,(i)-1)=(v))
#define RealVectorElement(a,i) (DYNX((a).data,(i)-1))
#define IntegerVectorElement(a,i) (DYNX((a).data,(i)-1))
#define StringVectorElement(a,i) (DYNX((a).data,(i)-1))


/* Routines for moutil.h */
#if !defined(LINALG_PACK)
#define LINALG_PACK
#define DeclareRealMatrix(Mat,r,c) RealArray Mat=RealTemporaryMatrix(r,c);int Mat##__m=r
#define DeclareRealVector(Mat,r) RealArray Mat=RealTemporaryVector(r)

/* Transposed to be compatible with Fortran */
#define DeclareStaticRealMatrix(Mat,r,c) \
static Real Mat##__internal[(r==0 || c==0)? 1 : r*c] ZERO_INITIALIZED;\
static SizeType Mat##__dims[2]={r,c};\
static RealArray Mat={2,Mat##__dims,Mat##__internal};\
static int Mat##__m=r

#define DeclareDidRealMatrix(Mat,r,c,hv) \
Real *Mat##__internal=&(hv);\
static SizeType Mat##__dims[2]={r,c};\
RealArray Mat={2,Mat##__dims,Mat##__internal};\
static int Mat##__m=r

#define MatrixNrRows(Mat) (Mat##__m,Mat.dims[0])
#define MatrixNrColumns(Mat) (Mat##__m,Mat.dims[1])

#define DeclareStaticRealVector(Mat,r) \
static Real Mat##__internal[r] ZERO_INITIALIZED;\
static SizeType Mat##__dims[1]={r};\
static RealArray Mat={1,Mat##__dims,Mat##__internal}

#define DeclareDidRealVector(Mat,r,hv) \
Real *Mat##__internal=&(hv);\
static SizeType Mat##__dims[1]={r};\
RealArray Mat={1,Mat##__dims,Mat##__internal}

#define DeclareDidRealVector1(Mat,r,hv) \
Real *Mat##__internal=&(hv);\
static SizeType Mat##__dims[1]={r};\
RealArray Mat={1,Mat##__dims,0}

#define DeclareDidRealVector2(Mat,r,hv) \
	Mat.data=Mat##__internal

#define DeclareDidRealMatrix1(Mat,r,c,hv) \
Real *Mat##__internal=&(hv);\
static SizeType Mat##__dims[2]={r,c};\
RealArray Mat={2,Mat##__dims,0};\
static int Mat##__m=r

#define DeclareDidRealMatrix2(Mat,r,c,hv) \
	Mat.data=Mat##__internal

#define SetMatrix(Mat,i,j,val) DYNX(Mat##__internal,(i-1)+(j-1)*Mat##__m)=val

#define SetMatrixLeading(Mat,i,j,ld,val) DYNX(Mat##__internal,(i-1)+(j-1)*ld)=val

#define GetMatrix(Mat,i,j) DYNX(Mat##__internal, (i-1)+(j-1)*Mat##__m)

#define GetVector(Mat,i) DYNX(Mat##__internal,(i-1))
#define SetVector(Mat,i,val) DYNX(Mat##__internal,(i-1))=val
#define SetVectorSlice(Mat,i,val) CopySlice(Mat##__internal+i-1,val)

/* Special case for Factored matrices. We only consider square factored matrices */
#define DeclareStaticFactoredRealMatrix(Mat,n) \
 DeclareStaticRealMatrix(Mat,n,n);\
 static Real Mat##__DWork[n*n+6*n] ZERO_INITIALIZED; \
 static int Mat##__IWork[2*n+1];\
 static int Mat##__FactoredV=0;int* Mat##__FactoredP=&(Mat##__FactoredV)

#define DeclareDidFactoredRealMatrix1(Mat,n,hv,nhv,iv,niv) \
 DeclareDidRealMatrix1(Mat,n,n,hv);\
 Real* Mat##__DWork=(&hv)+(n*n); \
 int* Mat##__IWork=&(iv);\
 int* Mat##__FactoredP=&(iv)+(2*n+1)

#define DeclareDidFactoredRealMatrix2(Mat,n,hv,nhv,iv,niv) \
 DeclareDidRealMatrix2(Mat,n,n,hv);DYNassert_static((n*n+(n*n+6*n))==(nhv) && (2*n+2)==(niv))

#define DeclareFactoredRealMatrix(Mat,n) \
  DeclareRealMatrix(Mat,n,n);\
  Real *Mat##__DWork=RealTemp(n*n+6*n);\
  int *Mat##__IWork=IntegerTemp(2*n+1);\
  int Mat##__FactoredV=0;int* Mat##__FactoredP=&(Mat##__FactoredV)

#define NeedFactor(Mat) (*(Mat##__FactoredP)==0)
#define SetNeedFactor(Mat) *(Mat##__FactoredP)=0

#define LinearSystemOfEquations(A, b, x, n) \
    DeclareFactoredRealMatrix(A, n); \
    DeclareRealVector(b, n); \
    DeclareRealVector(x, n); \
	int factoredForParametersV_=0;int*factoredForParametersP_=&factoredForParametersV_;\
	int factoredForEventsV_=0;int*factoredForEventsP_=&factoredForEventsV_;\
	const int NewParametersJac=1;\
    A##__Mark = PushMark();\
    RealFillAssign(b,0.0); \
    SetNeedFactor(A); 

#define StaticLinearSystemOfEquations(A,b,x,n) \
	DeclareStaticFactoredRealMatrix(A, n); \
	DeclareStaticRealVector(b, n); \
	DeclareStaticRealVector(x, n); \
	static int factoredForParametersV_=0;int*factoredForParametersP_=&factoredForParametersV_;\
	static int factoredForEventsV_=0;int*factoredForEventsP_=&factoredForEventsV_;\
	const int NewParametersJac=*(factoredForParametersP_) != dymolaParametersNr_;\
	RealFillAssign(b,0.0); \
        if (Init_) SetNeedFactor(A); 

#define DidLinearSystemOfEquations(A,b,x,n,hv,nhv,iv,niv) \
	DeclareDidFactoredRealMatrix1(A, n, hv, (n*n*2+6*n), iv, (2*n+2)); \
	DeclareDidRealVector1(b, n, (&hv)[(n*n*2+6*n)]); \
	DeclareDidRealVector1(x, n, (&hv)[(n*n*2+6*n+n)]); \
	int* factoredForParametersP_=&(iv)+(2*n+2);\
	int* factoredForEventsP_=&(iv)+(2*n+3);\
	int NewParametersJac=*(factoredForParametersP_) != dymolaParametersNr_;\
	DYNassert_static( (n*n+(n*n+8*n))==nhv && (2*n+4)==niv);\
	DeclareDidFactoredRealMatrix2(A, n, hv, (n*n*2+6*n), iv, (2*n+2)); \
	DeclareDidRealVector2(b, n, (&hv)[(n*n*2+6*n)]); \
	DeclareDidRealVector2(x, n, (&hv)[(n*n*2+6*n+n)]); \
	RealFillAssign(b,0.0); \
        if (Init_) SetNeedFactor(A); 

#if 0
#define OverdeterminedLinearSystemOfEquations(A,b,x,m,n) \
	static Real A##__internal[m*n] ZERO_INITIALIZED;\
	static SizeType A##__dims[2]={m,n};\
	static RealArray A={2,A##__dims,A##__internal};\
	static int A##__m=m;\
    static Real A##__DWork[n*n+6*n] ZERO_INITIALIZED; \
    static int A##__IWork[2*n+1];\
    static int A##__FactoredV=0;int* Mat##__FactoredP=Mat##__FactoredV;\
	DeclareStaticRealVector(b, m); \
	DeclareStaticRealVector(x, n); \
	static int factoredForParametersV_=0;int*factoredForParametersP_=&factoredForParametersV_;\
	static int factoredForEventsV_=0;int*factoredForEventsP_=&factoredForEventsV_;\
	const int NewParametersJac=*(factoredForParametersP_) != dymolaParametersNr_;\
	RealFillAssign(b,0.0); \
        if (Init_) SetNeedFactor(A); 

#define OverdeterminedDidLinearSystemOfEquations(A,b,x,m,n,hv,nhv,iv,niv) \
	Real A##__internal=&(hv);\
	static SizeType A##__dims[2]={m,n};\
	RealArray A={2,A##__dims,0};\
	static int A##__m=m;\
    Real* A##__DWork=&(hv)+m*n; \
    int* A##__IWork=&(iv);\
    int* Mat##__FactoredP=&(iv)+(2*n+1);\
	DeclareDidRealVector1(b, m, (&hv)[m*n+n*n+6*n]); \
	DeclareDidRealVector1(x, n, (&hv)[m*n+n*n+6*n+m]); \
	int*factoredForParametersP_=&(iv)+(2*n+2);\
	int*factoredForEventsP_=&(iv)+(2*n+3);\
	const int NewParametersJac=*(factoredForParametersP_) != dymolaParametersNr_;
	DYNassert_static(m*n+n*n+6*n+m+n==nhv && (2*n+4)==niv);\
	RealFillAssign(b,0.0); A.data=A##__internal;DeclareDidRealVector2(b, m, (&hv)[m*n+n*n+6*n]);DeclareDidRealVector2(x, n, (&hv)[m*n+n*n+6*n+m]);\
        if (Init_) SetNeedFactor(A); 
#endif

static void ReportDummySelection(const char*const varnames_[],const double*Pvar,double t,int m,int n) {
	int i,j,first=1;
	m=n-m;
	for(i=0;i<n;++i) {
		int selected=0;
		for(j=0;j<m;++j) selected=selected || Pvar[j*n+i];
		if (selected) {
			char str[1000];
			if (first)
				sprintf(str,"Selected at %g:\n %s.stateSelect=StateSelect.always",t,varnames_[i]);
			else
				sprintf(str,", %s.stateSelect=StateSelect.always",varnames_[i]);
			DymosimMessage(str);
			first=0;
		}
	}
}

#define DummyLinear(Jac,b,m,n,Pvar,var,stateP) \
	{\
		LIBDS_API int CheckP(const double*,double*,int,int,int*,const double*);\
		LIBDS_API void ComputeP(const double*,double*,int,int,int,int*,double*);\
		DeclareStaticRealMatrix(JacCopy_,m,n);\
		int i_,j_;\
		static int Pivots_[3*n];\
		if (Init_ || DYNEvent && CheckP(&GetMatrix(Jac,1,1),&(GetMatrix(JacCopy_,1,1)),m,n,Pivots_,Pvar)) {\
			RealArray varP_=var;\
			ComputeP(&GetMatrix(Jac,1,1),&(GetMatrix(JacCopy_,1,1)),m,m,n,Pivots_,Pvar);\
			{int m1_=n-m,n_=n,uno_=1;double d1_=1.0,d2_=0.0;\
				dgemv_("T",&n_,&m1_,&d1_,Pvar,&n_,varP_.data,&uno_,\
				&d2_,stateP,&uno_,1);}\
			Release();\
			for(i_=m+1;i_<=n;i_++) SetVector(b,i_,(stateP)[i_-(m+1)]);\
			if (PrintEvent&(1<<13)) ReportDummySelection(varnames_,Pvar,DYNTime,m,n);\
			for(i_=m+1;i_<=n;i_++) for(j_=1;j_<=n;j_++) SetMatrix(Jac,i_,j_,(Pvar)[(i_-(m+1))*n+(j_-1)]);\
			AnyREvent_=true;\
		} else if (CheckP(&GetMatrix(Jac,1,1),&(GetMatrix(JacCopy_,1,1)),m,n,Pivots_,Pvar)) triggerStepEvent_=1;\
		if ((PrintEvent&(1<<14))&&terminal()) ReportDummySelection(varnames_,Pvar,DYNTime,m,n);\
	}

#define DummyLinearDid(Jac,b,m,n,Pvar,var,stateP, hv, nhv, iv, niv) \
	{\
		LIBDS_API int CheckP(const double*,double*,int,int,int*,const double*);\
		LIBDS_API void ComputeP(const double*,double*,int,int,int,int*,double*);\
		DeclareDidRealMatrix1(JacCopy_,m,n,hv);\
		int i_,j_;\
		int*Pivots_=&(iv);\
		DYNassert_static(m*n==nhv && (3*n)==niv);\
		if (Init_ || DYNEvent && CheckP(&GetMatrix(Jac,1,1),&(GetMatrix(JacCopy_,1,1)),m,n,Pivots_,Pvar)) {\
			RealArray varP_=var;\
			ComputeP(&GetMatrix(Jac,1,1),&(GetMatrix(JacCopy_,1,1)),m,m,n,Pivots_,Pvar);\
			{int m1_=n-m,n_=n,uno_=1;double d1_=1.0,d2_=0.0;\
				dgemv_("T",&n_,&m1_,&d1_,Pvar,&n_,varP_.data,&uno_,\
				&d2_,stateP,&uno_,1);}\
			Release();\
			for(i_=m+1;i_<=n;i_++) SetVector(b,i_,(stateP)[i_-(m+1)]);\
			if (PrintEvent&(1<<13)) ReportDummySelection(varnames_,Pvar,DYNTime,m,n);\
			for(i_=m+1;i_<=n;i_++) for(j_=1;j_<=n;j_++) SetMatrix(Jac,i_,j_,(Pvar)[(i_-(m+1))*n+(j_-1)]);\
			AnyREvent_=true;\
		} else if (CheckP(&GetMatrix(Jac,1,1),&(GetMatrix(JacCopy_,1,1)),m,n,Pivots_,Pvar)) triggerStepEvent_=1;\
		if ((PrintEvent&(1<<14))&&terminal()) ReportDummySelection(varnames_,Pvar,DYNTime,m,n);\
	}

#define DummyLinearVar(Jac,b,m1,m2,n,Pvar,var,stateP) \
	{\
		LIBDS_API void DummyLinearVar1(int,double*,double*,double*,int,int,int,int*,double*);\
        LIBDS_API void DummyLinearVar2a(double*,double*,double*,int,int,int,int*,double*,const double*,double*);\
		LIBDS_API void DummyLinearVar2b(double*,double*,double*,int,int,int,int*,double*,double*);\
		DeclareStaticRealMatrix(Jac2_,n,m2);\
		static int Pivots_[1+3*n+m2];\
		DummyLinearVar1(DYNEvent,&GetMatrix(Jac,1,1),&GetMatrix(Jac2_,1,1),&GetVector(b,1),m1,m2,n,Pivots_,Pvar);\
		if (DYNEvent) {\
			RealArray varP_=var;\
			DummyLinearVar2a(&GetMatrix(Jac,1,1),&GetMatrix(Jac2_,1,1),&GetVector(b,1),m1,m2,n,Pivots_,Pvar,varP_.data,stateP);\
		    Release();\
			if (PrintEvent&(1<<13)) ReportDummySelection(varnames_,Pvar,DYNTime,m1,n);\
			AnyREvent_=true;\
		} else DummyLinearVar2b(&GetMatrix(Jac,1,1),&GetMatrix(Jac2_,1,1),&GetVector(b,1),m1,m2,n,Pivots_,Pvar,stateP);\
	}

#define DummyLinearVarDid(Jac,b,m1,m2,n,Pvar,var,stateP, hv, nhv, iv, niv) \
	{\
		LIBDS_API void DummyLinearVar1(int,double*,double*,double*,int,int,int,int*,double*);\
        LIBDS_API void DummyLinearVar2a(double*,double*,double*,int,int,int,int*,double*,const double*,double*);\
		LIBDS_API void DummyLinearVar2b(double*,double*,double*,int,int,int,int*,double*,double*);\
		DeclareDidRealMatrix(Jac2_,n,m2,hv);\
		int*Pivots_=&(iv);\
		DYNassert_static( (n*m2)==nhv && (1+3*n+m2)==niv );\
		DummyLinearVar1(DYNEvent,&GetMatrix(Jac,1,1),&GetMatrix(Jac2_,1,1),&GetVector(b,1),m1,m2,n,Pivots_,Pvar);\
		if (DYNEvent) {\
			RealArray varP_=var;\
			DummyLinearVar2a(&GetMatrix(Jac,1,1),&GetMatrix(Jac2_,1,1),&GetVector(b,1),m1,m2,n,Pivots_,Pvar,varP_.data,stateP);\
		    Release();\
			if (PrintEvent&(1<<13)) ReportDummySelection(varnames_,Pvar,DYNTime,m1,n);\
			AnyREvent_=true;\
		} else DummyLinearVar2b(&GetMatrix(Jac,1,1),&GetMatrix(Jac2_,1,1),&GetVector(b,1),m1,m2,n,Pivots_,Pvar,stateP);\
	}

#define OverdeterminedDummyLinear(Jac,b,l,m,m1,m2,n,Pvar,var,stateP) \
	{\
		LIBDS_API void CreateJacFullLin(double*,double*,const double*,int,int,int,const int*);\
		LIBDS_API void ComputePfull(const double*,double*,int,int,int,int,int,int*,int*,double*,const double*,const double*);\
		LIBDS_API int CheckPfull(const double*,double*,int,int,int,int,int,int*,const int*);\
		DeclareStaticRealMatrix(JacCopy_,m,n);\
		static int IWork_[4*m2+2*m+2*n];\
		if (DYNEvent) {\
			RealArray varP_=var;\
			ComputePfull(&(GetMatrix(Jac,1,1)),&(GetMatrix(JacCopy_,1,1)),l,m,n,m1,m2,IWork_+(m+n),IWork_,Pvar,varP_.data,stateP);\
			Release();\
			if (PrintEvent&(1<<13)) ReportDummySelection(varnames_,Pvar,DYNTime,m,n);\
			AnyREvent_=true;\
		} else if (CheckPfull(&(GetMatrix(Jac,1,1)),&(GetMatrix(JacCopy_,1,1)),l,m,n,m1,m2,IWork_+(m+n),IWork_)) triggerStepEvent_=1;\
		CreateJacFullLin(&(GetMatrix(Jac,1,1)),&(GetVector(b,1)),stateP,l,m,n,IWork_);\
	}

#define OverdeterminedDummyLinearDid(Jac,b,l,m,m1,m2,n,Pvar,var,stateP, hv, nhv, iv, niv) \
	{\
		LIBDS_API void CreateJacFullLin(double*,double*,const double*,int,int,int,const int*);\
		LIBDS_API void ComputePfull(const double*,double*,int,int,int,int,int,int*,int*,double*,const double*,const double*);\
		LIBDS_API int CheckPfull(const double*,double*,int,int,int,int,int,int*,const int*);\
		DeclareDidRealMatrix(JacCopy_,m,n,hv);\
		int*IWork_=&iv;\
		DYNassert_static( (m*n)==nhv && (4*m2+2*m+2*n)==niv);\
		if (DYNEvent) {\
			RealArray varP_=var;\
			ComputePfull(&(GetMatrix(Jac,1,1)),&(GetMatrix(JacCopy_,1,1)),l,m,n,m1,m2,IWork_+(m+n),IWork_,Pvar,varP_.data,stateP);\
			Release();\
			if (PrintEvent&(1<<13)) ReportDummySelection(varnames_,Pvar,DYNTime,m,n);\
			AnyREvent_=true;\
		} else if (CheckPfull(&(GetMatrix(Jac,1,1)),&(GetMatrix(JacCopy_,1,1)),l,m,n,m1,m2,IWork_+(m+n),IWork_)) triggerStepEvent_=1;\
		CreateJacFullLin(&(GetMatrix(Jac,1,1)),&(GetVector(b,1)),stateP,l,m,n,IWork_);\
	}

#define OverdeterminedLinear(Jac,b,l,m,m1,m2,n) \
	{\
		LIBDS_API void CreateJacFullLin(double*,double*,const double*,int,int,int,const int*);\
		LIBDS_API void ComputePfull(const double*,double*,int,int,int,int,int,int*,int*,double*,const double*,const double*);\
		LIBDS_API int CheckPfull(const double*,double*,int,int,int,int,int,int*,const int*);\
		DeclareStaticRealMatrix(JacCopy_,m,n);\
		static int IWork_[4*m2+2*m+2*n];\
		if (DYNEvent) {\
			ComputePfull(&(GetMatrix(Jac,1,1)),&(GetMatrix(JacCopy_,1,1)),l,m,n,m1,m2,IWork_+(m+n),IWork_,0,0,0);\
			if (PrintEvent&(1<<13)) ReportDummySelection(varnames_,Pvar,DYNTime,m,n);\
			AnyREvent_=true;\
		} else if (CheckPfull(&(GetMatrix(Jac,1,1)),&(GetMatrix(JacCopy_,1,1)),l,m,n,m1,m2,IWork_+(m+n),IWork_)) triggerStepEvent_=1;\
		CreateJacFullLin(&(GetMatrix(Jac,1,1)),&(GetVector(b,1)),0,l,m,n,IWork_);\
	}

	
#define OverdeterminedLinearDid(Jac,b,l,m,m1,m2,n, hv, nhv, iv, niv) \
	{\
		LIBDS_API void CreateJacFullLin(double*,double*,const double*,int,int,int,const int*);\
		LIBDS_API void ComputePfull(const double*,double*,int,int,int,int,int,int*,int*,double*,const double*,const double*);\
		LIBDS_API int CheckPfull(const double*,double*,int,int,int,int,int,int*,const int*);\
		DeclareDidRealMatrix(JacCopy_,m,n,hv);\
		int*IWork_=&iv;\
		DYNassert_static(m*n==nhv && 4*m2+2*m+2*n==niv);\
		if (DYNEvent) {\
			ComputePfull(&(GetMatrix(Jac,1,1)),&(GetMatrix(JacCopy_,1,1)),l,m,n,m1,m2,IWork_+(m+n),IWork_,0,0,0);\
			if (PrintEvent&(1<<13)) ReportDummySelection(varnames_,Pvar,DYNTime,m,n);\
			AnyREvent_=true;\
		} else if (CheckPfull(&(GetMatrix(Jac,1,1)),&(GetMatrix(JacCopy_,1,1)),l,m,n,m1,m2,IWork_+(m+n),IWork_)) triggerStepEvent_=1;\
		CreateJacFullLin(&(GetMatrix(Jac,1,1)),&(GetVector(b,1)),0,l,m,n,IWork_);\
	}

LIBDS_API int daxpy_(const int*,double*,double*,const int*,double*,const int*);
LIBDS_API int dgemv_(const char*,const int*,const int*,const double*,double*,const int*,double*,const int*,
			  const double*,double*,const int*,int);

#define MatrixZeros(J) RealFillAssign(J,0.0)
#define SolveLinearSystemOfEquations(A, b, x, num) \
    {double d1_,d2_;int i1_,i2_,i3_,i4_,Ierror_;\
	i1_ = A.dims[0]; \
    i2_ = A.dims[1]; \
	*(factoredForParametersP_) = dymolaParametersNr_; \
    i4_ = ! TryLUFactorization_ | (EquilibMatrix_ ? 2 : 0); \
	d1_ = -1;\
    d2_ = 1;\
	i3_ = 1;\
    if(1)dgemv_("N",&i2_,&i2_,&d1_,&(GetMatrix(A,1,1)),&i1_,&(GetVector(x,1)),&i3_,\
		&d2_,&(GetVector(b,1)),&i3_,1);\
	i3_ = num;\
    dymli3_(&i3_, &i4_, &(GetMatrix(A,1,1)), &i1_, &i2_, &(GetVector(b,1)), \
      &DYNTime, &DYNEvent, &PrintEvent, A##__DWork, A##__IWork,(A##__FactoredP),varnames_,&Ierror_);\
    i3_ = 1;\
	if(1)daxpy_(&i2_, &d2_, &(GetVector(b,1)), &i3_,&(GetVector(x,1)), &i3_);else RealAssign(x,b);\
	if (Ierror_!=0) {*(DSE_STRUCT QiErr)=Ierror_;DYN_GOTOLEAVE;}\
}

static void dymli3MI_(int* sysnr, int *fact, double *a,
					 int lda, int n,
		   double* b, double *Time, int *Event, int* PrintEvent,
		   double* dwork, int* iwork, int *factor, const char*const*varnames,int* ierr,
		   double*A0,double*b0,double*AR,double*ARwork,int*ARiwork,
		   int rrow,int rcol) {
	int i,j;
	int arFactor=0;
	if ((*factor&3)!=0) {
		*ierr=0;
		for(j=0;j<rcol;++j) {
			for(i=0;i<rrow;++i) b0[i]=a[i+j*lda]-A0[i+j*rrow];
			for(i=rrow;i<n;++i) b0[i]=0;
			dymli3_(sysnr,fact,a,&lda,&n,b0,Time,Event,PrintEvent,dwork,iwork,factor,varnames,ierr);
			if (*ierr!=0) {
				break;
			}
			for(i=0;i<rrow && i<rcol;++i) AR[i+j*rcol]=b0[i];
			for(i=rrow;i<rcol;++i) AR[i+j*rcol]=0;
			AR[j+j*rcol]+=1.0;
		}
		if (*ierr==0) {
			for(i=0;i<n;++i) b0[i]=b[i];
			dymli3_(sysnr,fact,a,&lda,&n,b0,Time,Event,PrintEvent,dwork,iwork,factor,varnames,ierr);
		}
		if (*ierr==0) dymli3_(sysnr,fact,AR,&rcol,&rcol,b0,Time,Event,PrintEvent,ARwork,ARiwork,&arFactor,0,ierr);
		if (*ierr==0) for(i=0;i<n;++i) b0[i+n]=b[i];
		if (*ierr==0) for(i=0;i<rrow;++i) {double d=0;for(j=0;j<rcol;++j)d+=(a[i+j*lda]-A0[i+j*rrow])*b0[j];b0[i+n]-=d;}
		if (*ierr==0) dymli3_(sysnr,fact,a,&lda,&n,b0+n,Time,Event,PrintEvent,dwork,iwork,factor,varnames,ierr);
		if (*ierr==0) for(i=0;i<n;++i) b[i]=b0[i+n];
		if (*ierr) {
			*factor=0;
		}
	}
	if ((*factor&3)==0) {
		dymli3_(sysnr,fact,a,&lda,&n,b,Time,Event,PrintEvent,dwork,iwork,factor,varnames,ierr);
		for(j=0;j<rcol;++j) for(i=0;i<rrow;++i) A0[i+j*rrow]=a[i+j*lda];
	}
}

#define SolveLinearSystemOfEquationsMI(A, b, x, num, rows, rrow, rcol) \
    {double d1_,d2_;int i1_,i2_,i3_,i4_,Ierror_;\
	static double b0[2*rows];\
	static double A0[rrow*rcol];\
	DeclareStaticFactoredRealMatrix(Ar,rcol);\
	i1_ = DYNX(A.dims,0); \
    i2_ = DYNX(A.dims,1); \
	*(factoredForParametersP_) = dymolaParametersNr_; \
    i4_ = ! TryLUFactorization_ | (EquilibMatrix_ ? 2 : 0); \
	d1_ = -1;\
    d2_ = 1;\
	i3_ = 1;\
    if(1)dgemv_("N",&i2_,&i2_,&d1_,&(GetMatrix(A,1,1)),&i1_,&(GetVector(x,1)),&i3_,\
		&d2_,&(GetVector(b,1)),&i3_,1);\
	i3_ = num;\
    dymli3MI_(&i3_, &i4_, &(GetMatrix(A,1,1)), i1_, i2_, &(GetVector(b,1)), \
	&DYNTime, &DYNEvent, &PrintEvent, A##__DWork, A##__IWork,(A##___FactoredP),varnames_,&Ierror_,\
	A0,b0,&(GetMatrix(Ar,1,1)),Ar##__DWork,Ar##__IWork, rrow,rcol\
	 );\
    i3_ = 1;\
	if(1)daxpy_(&i2_, &d2_, &(GetVector(b,1)), &i3_,&(GetVector(x,1)), &i3_);else RealAssign(x,b);\
	if (Ierror_!=0) {*QiErr=Ierror_;DYN_GOTOLEAVE;}\
}

#define DidSolveLinearSystemOfEquationsMI(A, b, x, num, rows, rrow, rcol, hv, nhv, iv, niv) \
    {double d1_,d2_;int i1_,i2_,i3_,i4_,Ierror_;\
	double*b0=hv;\
	double*A0=&(hv)+2*rows;\
	DeclareDidFactoredRealMatrix1(Ar,rcol, (&hv)[2*rows+rrow*rcol], nhv-(2*rows+rrow*rcol), iv, 2*rcol+2);\
	DeclareDidFactoredRealMatrix2(Ar,rcol, (&hv)[2*rows+rrow*rcol], nhv-(2*rows+rrow*rcol), iv, niv);\
	i1_ = DYNX(A.dims,0); \
    i2_ = DYNX(A.dims,1); \
	*(factoredForParametersP_) = dymolaParametersNr_; \
    i4_ = ! TryLUFactorization_ | (EquilibMatrix_ ? 2 : 0); \
	d1_ = -1;\
    d2_ = 1;\
	i3_ = 1;\
    if(1)dgemv_("N",&i2_,&i2_,&d1_,&(GetMatrix(A,1,1)),&i1_,&(GetVector(x,1)),&i3_,\
		&d2_,&(GetVector(b,1)),&i3_,1);\
	i3_ = num;\
    dymli3MI_(&i3_, &i4_, &(GetMatrix(A,1,1)), i1_, i2_, &(GetVector(b,1)), \
	&DYNTime, &DYNEvent, &PrintEvent, A##__DWork, A##__IWork,(A##___FactoredP),varnames_,&Ierror_,\
	A0,b0,&(GetMatrix(Ar,1,1)),Ar##__DWork,Ar##__IWork, rrow,rcol\
	 );\
    i3_ = 1;\
	if(1)daxpy_(&i2_, &d2_, &(GetVector(b,1)), &i3_,&(GetVector(x,1)), &i3_);else RealAssign(x,b);\
	if (Ierror_!=0) {*QiErr=Ierror_;DYN_GOTOLEAVE;}\
}

static void SolveScalarError(const char*xs,const char*bs,const char*As,double bval_,double Aval_,int alwaysError) {
	char s[1000];
	sprintf(s,"Error: %s for %.400s = (%.200s)/(%.200s) = %g/%g",
		(alwaysError&1)?"Scalar system is always singular":"Singular inconsistent scalar system",xs,bs,As,bval_,Aval_);
	if (alwaysError&2)
		DymosimError(s);
	else
		DymosimMessage(s);
}
static void SolveScalarWarning(const char*xs,double x) {
	char s[1000];
	sprintf(s,"Singular scalar system. Using minimum norm solution for %.400s = %g\n",xs,x);DymosimMessage(s);
}

#define SolveScalarLinear(A, As, b, bs, x, xs)\
{double Aval_=A;double bval_=b;if (Aval_!=0) x=bval_/Aval_; else if (bval_!=0) {SolveScalarError(xs,bs,As,bval_,Aval_,2);} else if(PrintEvent&(1<<12)) \
{ SolveScalarWarning(xs,x);}}

#define SolveScalarLinearMixed(A, As, b, bs, x, xs)\
{double Aval_=A;double bval_=b;if (Aval_!=0) x=bval_/Aval_; else if (bval_!=0) {SolveScalarError(xs,bs,As,bval_,Aval_,0);MixedFailFlag_=1;} else if(PrintEvent&(1<<12)) \
{ SolveScalarWarning(xs,x);}}

#define SolveScalarLinearParametric(A, As, b, bs, x, xs)\
{double Aval_=A;double bval_=b;if (Aval_!=0) x=bval_/Aval_; else {SolveScalarError(xs,bs,As,bval_,Aval_,3);}}

#define SolveScalarLinearMixedParametric(A, As, b, bs, x, xs)\
{double Aval_=A;double bval_=b;if (Aval_!=0) x=bval_/Aval_; else {SolveScalarError(xs,bs,As,bval_,Aval_,1);MixedFailFlag_=1;}}


#define SolveLinearSystemOfEquationsMixed(A, b, x, num) \
    {double d1_,d2_;int i1_,i2_,i3_,i4_,Ierror_;\
	i1_ = DYNX(A.dims,0); \
    i2_ = DYNX(A.dims,1); \
	*(factoredForParametersP_) = dymolaParametersNr_; \
    i4_ = ! TryLUFactorization_ | (EquilibMatrix_ ? 2 : 0); \
	d1_ = -1;\
    d2_ = 1;\
	i3_ = 1;\
    if(1)dgemv_("N",&i2_,&i2_,&d1_,&(GetMatrix(A,1,1)),&i1_,&(GetVector(x,1)),&i3_,\
		&d2_,&(GetVector(b,1)),&i3_,1);\
	i3_ = num;\
    dymli3_(&i3_, &i4_, &(GetMatrix(A,1,1)), &i1_, &i2_, &(GetVector(b,1)), \
      &DYNTime, &DYNEvent, &PrintEvent, A##__DWork, A##__IWork,(A##__FactoredP),varnames_,&Ierror_);\
    i3_ = 1;\
	if(1)daxpy_(&i1_, &d2_, &(GetVector(b,1)), &i3_,&(GetVector(x,1)), &i3_);else RealAssign(x,b);\
	if (Ierror_!=0) {MixedFailFlag_=Ierror_;Ierror_=0;}\
}

#define SolveLinearSystemOfEquationsMIMixed(A, b, x, num, rows, rrow, rcol) \
    {double d1_,d2_;int i1_,i2_,i3_,i4_,Ierror_;\
	static double b0[2*rows];\
	static double A0[rrow*rcol];\
	DeclareStaticFactoredRealMatrix(Ar,rcol);\
	i1_ = DYNX(A.dims,0); \
    i2_ = DYNX(A.dims,1); \
	*(factoredForParametersP_) = dymolaParametersNr_; \
    i4_ = ! TryLUFactorization_ | (EquilibMatrix_ ? 2 : 0); \
	d1_ = -1;\
    d2_ = 1;\
	i3_ = 1;\
    if(1)dgemv_("N",&i2_,&i2_,&d1_,&(GetMatrix(A,1,1)),&i1_,&(GetVector(x,1)),&i3_,\
		&d2_,&(GetVector(b,1)),&i3_,1);\
	i3_ = num;\
	dymli3MI_(&i3_, &i4_, &(GetMatrix(A,1,1)), i1_, i2_, &(GetVector(b,1)), \
	&DYNTime, &DYNEvent, &PrintEvent, A##__DWork, A##__IWork,(A##__FactoredP),varnames_,&Ierror_,\
	A0,b0,&(GetMatrix(Ar,1,1)),Ar##__DWork,Ar##__IWork, rrow,rcol\
	 );\
    i3_ = 1;\
	if(1)daxpy_(&i1_, &d2_, &(GetVector(b,1)), &i3_,&(GetVector(x,1)), &i3_);else RealAssign(x,b);\
	if (Ierror_!=0) {MixedFailFlag_=Ierror_;Ierror_=0;}\
}


#define DidSolveLinearSystemOfEquationsMIMixed(A, b, x, num, rows, rrow, rcol, hv, nhv, iv, niv) \
    {double d1_,d2_;int i1_,i2_,i3_,i4_,Ierror_;\
	double*b0=&hv;\
	double*A0=&hv+2*rows\
	DeclareDidFactoredRealMatrix1(Ar,rcol, (&hv)[2*rows+rrow*rcol], nhv -(2*rows+rrow*rcol), iv, niv);\
	DeclareDidFactoredRealMatrix2(Ar,rcol, (&hv)[2*rows+rrow*rcol], nhv -(2*rows+rrow*rcol), iv, niv);\
	i1_ = DYNX(A.dims,0); \
    i2_ = DYNX(A.dims,1); \
	*(factoredForParametersP_) = dymolaParametersNr_; \
    i4_ = ! TryLUFactorization_ | (EquilibMatrix_ ? 2 : 0); \
	d1_ = -1;\
    d2_ = 1;\
	i3_ = 1;\
    if(1)dgemv_("N",&i2_,&i2_,&d1_,&(GetMatrix(A,1,1)),&i1_,&(GetVector(x,1)),&i3_,\
		&d2_,&(GetVector(b,1)),&i3_,1);\
	i3_ = num;\
	dymli3MI_(&i3_, &i4_, &(GetMatrix(A,1,1)), i1_, i2_, &(GetVector(b,1)), \
	&DYNTime, &DYNEvent, &PrintEvent, A##__DWork, A##__IWork,(A##__FactoredP),varnames_,&Ierror_,\
	A0,b0,&(GetMatrix(Ar,1,1)),Ar##__DWork,Ar##__IWork, rrow,rcol\
	 );\
    i3_ = 1;\
	if(1)daxpy_(&i1_, &d2_, &(GetVector(b,1)), &i3_,&(GetVector(x,1)), &i3_);else RealAssign(x,b);\
	if (Ierror_!=0) {MixedFailFlag_=Ierror_;Ierror_=0;}\
}

#define EndLinearSystemOfEquations(A) PopMark(A##__Mark);Release()
#define EndStaticLinearSystemOfEquations(A) 
#endif

#ifdef NDEBUG
#define VerifyIndex(i,m,s) (i)
#else
static unsigned int VerifyIndex(unsigned int i,unsigned int maxVal,const char*s) {
	Assert(i<maxVal,s);
	return i;
}
#endif

/* Major pop, e.g. at end of integration */

DYMOLA_STATIC void PopAllMarks() {
  PopMarkMarks();
  PopModelContext();
}
static double spliceFunction2(double,double);
static double spliceFunction2der(double x,double deltax,double dx,double ddeltax);
static double divideDymola(double x,double y);
static double remainderDymola(double x,double y);
static double modulusDymola(double x,double y);
static 
#if defined(_MSC_VER) && _MSC_VER>=1200
__inline
#elif __GNUC__
__inline
#endif
double spliceFunctionF(double pos,double neg,double x,double deltax) {
	const double d=spliceFunction2(x,deltax);
	return d*(d>0 ? pos : 0.0)+(1-d)*(d<1 ? neg : 0.0);
}
DYMOLA_STATIC double spliceFunctionF_der(double pos,double neg,double x,double deltax,double pos_der,double neg_der,double x_der,double deltax_der) {
	const double d=spliceFunction2(x,deltax);
	return (d*(d>0 ? (pos_der) : 0.0)+(1-d)*(d<1 ? (neg_der) : 0.0)\
    + (d>0 && d<1 ? spliceFunction2der((x),(deltax),(x_der),(deltax_der))*((pos)-(neg)): 0.0));

}
  static void ModelicaInternal_mkdir(const char* directoryName);
  static void ModelicaInternal_rmdir(const char* directoryName);
  static int  ModelicaInternal_stat(const char* name);
  static void ModelicaInternal_rename(const char* oldName, const char* newName) ;
  static void ModelicaInternal_removeFile(const char* file);
  static void ModelicaInternal_copyFile(const char* oldFile, const char* newFile);
  static void ModelicaInternal_readDirectory(const char* directory, int nFiles, const char* files[]);
  static int  ModelicaInternal_getNumberOfFiles(const char* directory);
  static const char* ModelicaInternal_fullPathName(const char* name);
  static const char* ModelicaInternal_temporaryFileName();
  static void ModelicaInternal_print(const char* string, const char* fileName);
  static int  ModelicaInternal_countLines(const char* fileName);
  static void ModelicaInternal_readFile(const char* fileName, const char* string[], size_t nLines);
  static const char* ModelicaInternal_readLine(const char* fileName, int lineNumber, int* endOfFile);
  static void ModelicaInternal_chdir(const char* directoryName);
  static const char* ModelicaInternal_getcwd(int dummy);
  static void ModelicaInternal_getenv(const char* name, int convertToSlash, const char**,int* exist);
  static void ModelicaInternal_setenv(const char* name, const char* value, int convertFromSlash);

  static const char* ModelicaStrings_substring(const char* string, int startIndex, int endIndex);
  static int ModelicaStrings_length(const char* string);
  static int ModelicaStrings_compare(const char* string1, const char* string2, int caseSensitive);
  static int ModelicaStrings_skipWhiteSpace(const char* string, int i);
  static void ModelicaStrings_scanIdentifier(const char* string, int startIndex, int* nextIndex, const char** identifier);
  static void ModelicaStrings_scanInteger(const char* string, int startIndex, int unsignedNumber, 
                                        int* nextIndex, int* integerNumber);
  static void ModelicaStrings_scanReal(const char* string, int startIndex, int unsignedNumber,
                                     int* nextIndex, double* number);
  static void ModelicaStrings_scanString(const char* string, int startIndex,
                                       int* nextIndex, const char** result);
  static void ModelicaStreams_closeFile(const char* fileName);

  struct ModelicaInternal_readLine_M_struct {
  const char*   string0_0_0member;
  int   endOfFile0_0_0member;
};
  DYMOLA_STATIC struct ModelicaInternal_readLine_M_struct ModelicaInternal_readLine_M(const char*  fileName0_0, int  lineNumber0_0) {
  {
    const char*   string0_0;
    int   endOfFile0_0;
    string0_0="";
    endOfFile0_0=0;
     {
      string0_0 = (ModelicaInternal_readLine)(fileName0_0, lineNumber0_0, &endOfFile0_0);
      }
    {
      struct ModelicaInternal_readLine_M_struct out_;
      out_.string0_0_0member = string0_0;
      out_.endOfFile0_0_0member = endOfFile0_0;
      return out_;
    }
  }}
  struct ModelicaInternal_getenv_M_struct {
  const char*   content0_0_0member;
  int   exist0_0_0member;
};
  DYMOLA_STATIC struct ModelicaInternal_getenv_M_struct ModelicaInternal_getenv_M(const char*  
  name0_0, int  convertToSlash0_0) {
  PushContext("ModelicaInternal_getenv")
  {
    const char*   content0_0;
    int   exist0_0;
    content0_0="";
    exist0_0=0;
    {
       (ModelicaInternal_getenv)(name0_0, convertToSlash0_0, &content0_0, & 
        exist0_0);
      }
	PopContext()
    {
      struct ModelicaInternal_getenv_M_struct out_;
      out_.content0_0_0member = content0_0;
      out_.exist0_0_0member = exist0_0;
      return out_;
    }
  }}

  DYMOLA_STATIC StringArray    ModelicaInternal_readDirectory_M(const char*  directory0_0, int  
  nNames0_0) {
  {
    StringArray    names0_0;
    MarkObject retmark_ = PushMark();
    names0_0=StringTemporary( 1, nNames0_0);
    retmark_ = PushMark();
    StringFillAssign( names0_0, "");
    /* Start of real code */
    {
      (ModelicaInternal_readDirectory)(directory0_0, nNames0_0, names0_0.data);
      }
    /* Output section */
    PopMark(retmark_);
    return names0_0;
  }}
  struct ModelicaStrings_scanInteger_M_struct {
	 int nextIndex0_0_0member;
	 int number0_0_0member;
  };
  DYMOLA_STATIC struct ModelicaStrings_scanInteger_M_struct ModelicaStrings_scanInteger_M
		(const char*string0_0, int startTokenIndex0_0, int unsigned0_0) {
	struct ModelicaStrings_scanInteger_M_struct out;
	ModelicaStrings_scanInteger(string0_0,startTokenIndex0_0, unsigned0_0, &out.nextIndex0_0_0member, &out.number0_0_0member);
    return out;
}
struct ModelicaStrings_scanReal_M_struct {
	 int nextIndex0_0_0member;
	 double number0_0_0member;
  };
  DYMOLA_STATIC struct ModelicaStrings_scanReal_M_struct ModelicaStrings_scanReal_M
		(const char*string0_0, int startTokenIndex0_0, int unsigned0_0) {
	struct ModelicaStrings_scanReal_M_struct out;
	ModelicaStrings_scanReal(string0_0,startTokenIndex0_0, unsigned0_0, &out.nextIndex0_0_0member, &out.number0_0_0member);
    return out;
}
struct ModelicaStrings_scanString_M_struct {
	 int nextIndex0_0_0member;
	 const char* string20_0_0member;
  };
  DYMOLA_STATIC struct ModelicaStrings_scanString_M_struct ModelicaStrings_scanString_M
		(const char*string0_0, int startTokenIndex0_0) {
	struct ModelicaStrings_scanString_M_struct out;
	ModelicaStrings_scanString(string0_0,startTokenIndex0_0, &out.nextIndex0_0_0member, &out.string20_0_0member);
    return out;
}
struct ModelicaStrings_scanIdentifier_M_struct {
	 int nextIndex0_0_0member;
	 const char* identifier0_0_0member;
  };
  DYMOLA_STATIC struct ModelicaStrings_scanIdentifier_M_struct ModelicaStrings_scanIdentifier_M
		(const char*string0_0, int startTokenIndex0_0) {
	struct ModelicaStrings_scanIdentifier_M_struct out;
	ModelicaStrings_scanIdentifier(string0_0,startTokenIndex0_0, &out.nextIndex0_0_0member, &out.identifier0_0_0member);
    return out;
}

  static void OutputCoverageResult(const char*const classNamesForCoverage_[],
								 const char*const componentNameStrings_[],
								 const int componentNameStringMap_[],
								 const char*const expressionsStringsForCoverage_[],
								 const int firstExprIndexForCoverage_[],
                                 int resultFalse_[], int resultTrue_[], 
								 int nclasses)
{
	int i = 0;
	int j = 0;
#if 0
	DymosimMessage("Result of coverage analysis.");
	DymosimMessage("The following conditions are constant as indicated.");
	for (i = 0; i < nclasses; ++i) {
		Dymola_bool anyfound = false;
		int startindex = firstExprIndexForCoverage_[i];
		int stop = firstExprIndexForCoverage_[i+1];
		for (j = startindex; j < stop; ++j) {
			Dymola_bool isConstant = false;
			Dymola_bool value = false;
			if (resultFalse_[j] != resultTrue_[j]) {
				isConstant = true;
				if (resultFalse_[j] == 1) 
					value = false;
				else 
					value = true;
			}
			if (isConstant) {
				if (!anyfound) {
					anyfound = true;
					DymosimMessage(classNamesForCoverage_[i]);
				}
				if (!value) {
					DymosimMessage("false  ");
					DymosimMessage(expressionsStringsForCoverage_[j]);
				}
				else {
					DymosimMessage("true  ");
					DymosimMessage(expressionsStringsForCoverage_[j]);
				}
			}
		}
	}
#endif
#if !defined(DYMOLA_DSPACE) && !defined(NO_FILE)
	{	
		FILE *f=fopen("Coverage.mo","w");
		if (f) {
			int noutputs = 0;
			fprintf(f,"/* Result of coverage analysis. */\n\n");
			for (i = 0; i < nclasses; ++i) {
				int startindex = firstExprIndexForCoverage_[i];
				int stop = firstExprIndexForCoverage_[i+1];
				for (j = startindex; j < stop; ++j) {
					Dymola_bool isConstant = false;
					Dymola_bool value = false;
					if (resultFalse_[j] != resultTrue_[j]) {
						isConstant = true;
						if (resultFalse_[j] == 1) 
							value = false;
						else 
							value = true;
					}
					if (isConstant) {
						if (noutputs == 0) {
							fprintf(f,"{");
						} else {
							fprintf(f,",\n");
						}
						++noutputs;

						fprintf(f,"ModelManagement.Check.Coverage.CoverElement(");
						fprintf(f,"\"%s\"", classNamesForCoverage_[i]);
						fprintf(f,", \"%s\"", componentNameStrings_[componentNameStringMap_[j]]);
						fprintf(f,", \n     \"%s\"", expressionsStringsForCoverage_[j]);

						if (!value) fprintf(f,", \"false\"");
						else fprintf(f,", \"true\"");
						fprintf(f,")");
					} else {
						if (noutputs == 0) {
							fprintf(f,"{");
						} else {
							fprintf(f,",\n");
						}
						++noutputs;

						fprintf(f,"ModelManagement.Check.Coverage.CoverElement(");
						fprintf(f,"\"%s\"", classNamesForCoverage_[i]);
						fprintf(f,", \"%s\"", componentNameStrings_[componentNameStringMap_[j]]);
						fprintf(f,", \n     \"%s\"", expressionsStringsForCoverage_[j]);

						fprintf(f,", \"covered\")");
						

					}
				}
			}
			if (noutputs==0) {
				fprintf(f,"fill(ModelManagement.Check.Coverage.CoverElement(\"\",\"\",\"\",\"\"),0)\n");
			} else {
				fprintf(f,"}\n");
			}
		}
		fclose(f);
	}
#endif
}

struct spatialDistribution_M_struct {
  double out00_0_0member;
  double out10_0_0member;
};
DYMOLA_STATIC int spatialDistribution_M1(int sz,double x, int positiveVelocity,
	RealArray initialPoints, RealArray initialValues) 
{
	int res,i,n,firstTime=1;
	double u0=0,d1,d2;
	Assert(initialPoints.ndims==1,"Initial Points should be a vector");
	Assert(initialValues.ndims==1,"Initial Values should be a vector");
	Assert(initialPoints.dims[0]==initialValues.dims[0],"Initial Points and Values should have size");
	n=initialValues.dims[0];
	if (n>0) {
		
		for(i=n-1;i>=0;--i) {
			double xVal,uVal;
			xVal=x+initialPoints.data[i];
			uVal=initialValues.data[i];
			if (firstTime) {
				res=transportIni(sz, xVal, uVal);
			} else {
				if (i==0 && initialPoints.data[i]==0 && n>1 && initialPoints.data[i+1]>1e-13)
					transportFunction(uVal, uVal, xVal+1e-13, 1, res, &d1, &d2);
				transportFunction(uVal, uVal, xVal, 1, res, &d1, &d2);
			}
			firstTime=0;
		}
	} else {
		res=transportIni(sz, x, u0);
	}
	return res;
}
DYMOLA_STATIC struct spatialDistribution_M_struct spatialDistribution_M2(int id,double in0, double in1, double x, int positiveVelocity) 
{
	struct spatialDistribution_M_struct  res;
	transportFunction(in0, in1, -x, positiveVelocity, id, &(res.out00_0_0member), &(res.out10_0_0member));
	return res;
}
#define spatialDistribution_M(in0, in1,x,pv,ip,iv,nr) \
  (delayID[nr]= (delayID[nr]==0) ? 1+spatialDistribution_M1(4000,x,pv,ip,iv): delayID[nr], \
    spatialDistribution_M2(delayID[nr]-1,in0,in1,x,pv))

#if defined(LIBDS_DLL) || defined(LIBDS_NODLL_API) || (defined(_OPENMP) && !defined(DISABLE_DYMOLA_OPENMP)) || defined(DYN_MULTINSTANCE)
#ifndef MATRIXOP_IMP
#define MATRIXOP_IMP
#include "matrixop.c"
#endif
#endif
#endif
