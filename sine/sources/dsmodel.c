#define DYN_MULTINSTANCE 1 
#include "dsmodel_fmu.h" 
#include "dsutil.h" 
#if defined(FMU_SOURCE_CODE_EXPORT)
#define DYMOSIM_RT_IMP
#endif
#ifdef DYMOSIM_RT_IMP
/*
 * Copyright (C) 1997-2001 Dynasim AB.
 * All rights reserved.
 *
 */

#if defined(DYMOLA_DSPACE)
/* Realtime compilation. DB 1998-10-05 */
#define NO_FILE
#include <dsdefs.h>
#else
#include <stdio.h>
#endif
#include "assumption.h"

#include <stdlib.h>
#include <errno.h> 
#include <string.h>
#include <ctype.h>
#include <math.h>
#include "f2c.h"
#include "amat.h"
#include "float.h"
#if !defined(MOTIF)
#include "sprwat.h"
#endif
#include "localeless.h"

#if defined(DYMOLAB) || defined(OPENGL)
extern int maxRows();
#else
static integer maxRows() {return 0;}
#endif
#if !defined(NO_FILE) && !defined(DYMOLAB) && !defined(DYMOSIM) && !defined(FMU_SOURCE_CODE_EXPORT) && !defined(Matlab5) && !defined(Matlab51) && !defined(Matlab4)
#include "supportunicodefiles.h"
#endif

#define ALLOW_REOPEN
#ifndef LIBDS_DLL
#define NAMATERROR 1500
DYMOLA_STATIC char amatError[NAMATERROR] = { "No message." };
DYMOLA_STATIC char* amatErrorFunction() {
	return amatError;
}
                                  /* utility storage on which a message     */
                                  /* may be placed                          */
#endif


#if !defined(NO_FILE) || defined(DYM2CCUR)

#define RET_MZERO 1
#define RET_REQ   2
#define RET_EOF   3
#define RET_ERR   4

static int amatGetCHeader (AmatGetFile *file, Amatrix *matrix, Amatrix matrixReq);

static int amatGetFHeader (AmatGetFile *file, Amatrix *matrix, char name[], integer lname);
static int amatGetFMatrix (AmatGetFile *file, Amatrix *matrix);
static int amatGetFNMatrix(AmatGetFile *file, Amatrix *matrix);
static int amatGetFText   (AmatGetFile *file, Amatrix *matrix);

static int amatGetBHeader (AmatGetFile *file, Amatrix *matrix, char name[], integer lname, int *prec);
static int amatGetBMatrix (AmatGetFile *file, Amatrix *matrix, AmatType typeOnFile, int prec,int removeSpace);
static int amatGetBIMatrix(AmatGetFile *file, Amatrix *matrix);
static int amatGetBRMatrix(AmatGetFile *file, Amatrix *matrix);
static int amatGetBDMatrix(AmatGetFile *file, Amatrix *matrix);
static int amatGetBText   (AmatGetFile *file, Amatrix *matrix, int prec,int removeSpace);

DYMOLA_STATIC int (*amatStopRequest)(void)=0;
static int amatGetSkip    (AmatGetFile *file);


/* Error messages */
   static const char *const amatg_table[] = {
      /* 0 */   "Error reading matrix header from file \"%.400s\":\n%.400s\n"
      /* 1 */  ,"Wrong row dimension (= %d) of matrix \"%.400s\" on file \"%.400s\".\n"
                "Dimension must be >= 0.\n"
      /* 2 */  ,"Wrong column dimension (= %d) of matrix \"%.400s\" on file \"%.400s\".\n" 
                "Dimension must be >= 0.\n"
      /* 3 */  ,"Wrong precision flag (= %d) of numeric matrix \"%.400s\" on file \"%.400s\"\n"
                "Flag must be 0 (= 64-bit) or 1 (32-bit).\n"
      /* 4 */  ,"Wrong type flag (= %d) of matrix \"%.400s\" on file \"%.400s\".\n"
                "The type flag must be 0 (= numeric matrix) or 1 (= character matrix).\n"
      /* 5 */  ,"Error when reading numeric matrix \"%.400s(%d,%d)\" from file\n"
                "\"%.400s\", since it is not possible to allocate storage for\n"
                "%d numeric elements.\n"
      /* 6 */  ,"Error when reading numeric data of matrix \"%.400s(%d,%d)\" from file\n"
                "\"%.400s\"\n"
      /* 7 */  ,"Not possible to open file \"%.400s\": %.400s\n"
      /* 8 */  ,"Wrong precision flag (= %d) of character matrix \"%.400s(%d,%d)\"\n"
                "on file \"%.400s\".\n"
                "Flag must be 0 (= 64-bit) or 5 (= 8-bit unsigned integer).\n"
     /*  9 */  ,"Error reading first character from file \"%.400s\":\n%.400s\n"
     /* 10 */  ,"Error should never occur when reading file \"%.400s\"\n"
                "(FileFormat = %d).\n"
     /* 11 */  ,"Error should never occur when reading file \"%.400s\"\n"
                "(matrix.type = %d).\n"
     /* 12 */  ,"Wrong error number\n"
     /* 13 */  ,"Error when reading text matrix \"%.400s(%d,%d)\" from file \"%.400s\",\n"
                "since it is not possible to allocate storage for\n"
                "%d characters.\n"
     /* 14 */  ,"Error reading data of text matrix \"%.40s(%d,%d)\" from file \"%.400s\":\n"
                "%.400s\n"
     /* 15 */  ,"Not possible to allocate enough storage in order to open\n"
                "file \"%.400s\" for reading of matrices:\n%.400s\n"
     /* 16 */  ,"Not possible to read matrices from file \"%.400s\",\n"
                "because the version number  (= \"%.400s\") does not correspond to\n"
                "the required version number (= \"%.400s\").\n"
                "Note, that the version number is stored in the second row\n"
                "of the first matrix \"Aclass\" on file.\n"
     /* 17 */  ,"Not possible to read matrices from file \"%.400s\",\n"
                "because matrix \"Aclass\" has less than 3 rows.\n"
                "(in the first row the class-name, in the second row\n"
                "the version-number and in the third row the description text\n"
                "must be stored)\n"
     /* 18 */  ,"Error when reading matrix from file \"%.400s\",\n"
                "because it is not possible to allocate storage for\n"
                "the matrix name (length = %d)\n:%.400s\n"
     /* 19 */  ,"Not possible to read matrices from file \"%.400s\",\n"
                "because the class name (= \"%.400s\") does not correspond to\n"
                "the required class name (= \"%.400s\").\n"
                "Note, that the class name is stored in the first row\n"
                "of the first matrix \"Aclass\" on file.\n"
     /* 20 */  ,"Cannot read matrix from file \"%.400s\".\n"
	            "The file has incorrect format (note that binary files must be saved with '-v4' in MATLAB).\n"
     /* 21 */  ,"Error reading format and version information in first line of \n"
                "file \"%.400s\": \"#%.400s\" expected, but \"#%c\" found.\n"
     /* 22 */  ,"Error scanning the following declaration on file \"%.400s\":\n"
                "%.400s"
                "Expected is: <type> <name> (<row_dimension>, <column_dimension>)\n"
                "where <type> is \"integer\", \"real\", \"double\" or \"char\" and\n"
                "<name> is the name of the matrix.\n"
     /* 23 */  ,"Not possible to allocate enough storage in order to read\n"
                "matrix from file \"%.400s\":\n%.400s\n"
     /* 24 */  ,"Error scanning the following declaration on file \"%.400s\":\n"
                "%.400s"
                "Name of matrix is wrong, should be \"%.400s\".\n"
     /* 25 */  ,"Error when reading matrix \"%.400s(%d,%d)\" from file \"%.400s\",\n"
                "because matrix should be of type \"integer\", \"real\" or \"double\"\n"
                "and not of type \"%.400s\".\n"
     /* 26 */  ,"Error when reading matrix \"%.400s(%d,%d)\" from file \"%.400s\",\n"
                "because matrix should be of type \"char\" and not of type \"%.400s\".\n"
     /* 27 */  ,"Error when reading matrix \"%.400s(%d,%d)\" from file \"%.400s\",\n"
                "because row dimension = %d required.\n"
     /* 28 */  ,"Error when reading matrix \"%.400s(%d,%d)\" from file \"%.400s\",\n"
                "because column dimension = %d required.\n"
     /* 29 */  ,"Error when reading matrix \"%.400s(%d,%d)\" from file \"%.400s\",\n"
                "because number of matrix elements is zero.\n"
     /* 30 */  ,"Error should never occur when reading matrix \"%.400s(%d,%d)\"\n"
                "from file \"%.400s\" (matrix->type = %d)\n"
     /* 31 */  ,"Error when reading matrix \"%.400s\" from file \"%.400s\",\n"
                "because matrix-name = \"%.400s\" required.\n"
     /* 32 */  ,"Programming error when reading matrix \"%.400s(%d,%d)\" from file \"%.400s\".\n"
                "Function amatGetMatrix is called with storage supplied for the\n"
                "matrix data (matrix.data != NULL), but row dimension, column dimension\n"
                "and matrix type are not all predefined\n"
                "(matrix.nrow = %d, matrix.ncol = %d, matrix.type = %.400s)\n"
     /* 33 */  ,"Programming error detected when reading file \"%.400s\":\n"
                "Function amatReadAll is called, without providing storage for\n"
                "array \"%.400s\", (although %d number of elements required).\n"
     /* 34 */  ,"Matrix \"%.400s\" is not expected on file \"%.400s\".\n"
     /* 35 */  ,"Programming error detected when reading file \"%.400s\":\n"
                "Function amatReadAll is called with m[%d].matrix.nrow = %d (!= 0),\n"
                "although m[%d].rowDim > 0. This is not allowed.\n"
     /* 36 */  ,"Programming error detected when reading file \"%.400s\":\n"
                "Function amatReadAll is called with m[%d].matrix.ncol = %d (!= 0),\n"
                "although m[%d].colDim > 0. This is not allowed.\n"
     /* 37 */  ,"Error when reading matrix \"%.400s(%d,%d)\" from file \"%.400s\",\n"
                "because matrix \"%.400s\" is not stored in this file (which is required)\n"
     /* 38 */  ,"Error when reading matrix \"%.100s(%d,%d)\" from file \"%.400s\",\n"
                "because %.100s dimension does not agree with %.100s dimension of\n"
                "matrix \"%.100s(%d,%d)\".\n"
     /* 39 */  ,"Error when reading matrix \"%.400s(%d,%d)\" from file \"%.400s\",\n"
                "because row dimension <= %d required.\n"
     /* 40 */  ,"Error when reading matrix \"%.400s(%d,%d)\" from file \"%.400s\",\n"
                "because column dimension <= %d required.\n"
     /* 41 */  ,"Matrix \"%.400s\" has to be present on file \"%.400s\"\n"
                "which is not the case.\n"
     /* 42 */  ,"Programming error detected when reading file \"%.400s\":\n"
                "Index rowIndex/colIndex (= %d) stored in %d-th matrix \"%.400s\"\n"
                "is not in the required range 0 <= Index < %d.\n"
     /* 43 */  ,"Not possible to read matrices from file \"%.400s\",\n"
                "because the first matrix on file (=\"%.400s\")\n"
                "has not the name \"class\" or \"Aclass\".\n"
	 /* 44 */  ,"Reading of file \"%.400s\" stopped by user.\n"
   };



int amatRead (const char *fileName, Amatrix *matrix){ 

   /* Read single matrix from a file. */
      AmatGetFile  file;

      if ( amatGetOpen(fileName, "NoClass", "0.0", &file ) != 0 ) return 1;
      if ( amatGetMatrix(&file, matrix) != 0 ) {
         amatGetClose(&file);
         return 2;
      }
      amatGetClose(&file);
      return 0;
}



int amatReadAll (const char *fileName, const char *classReq, const char *versionReq, 
                 AclassRead c[], integer dim_c) {

   /* Read all matrices from file. */

   /* Declarations */
      AmatGetFile  file;
      int          err;

   /* Open file and check class name and version number */
      if ( (err=amatGetOpen(fileName, classReq, versionReq, &file)) != 0 ) return err;

   /* Read all matrices from file */
      err=amatReadAll2(&file, c, dim_c);
      return err;
}


static void amatUnexpectedNull(void) {
	sprintfC(amatError,"Unexpected null-pointer");
}
int amatReadAll2 (AmatGetFile *file, AclassRead c[], integer dim_c) {

   /* Read all matrices from file. */

   /* Declarations */
      char        *name=NULL, *str1, *str2;
      integer      i, j, ndim, found;
      int          err, eof;
      const char        *fileName;
	  if (!file) {amatUnexpectedNull();return RET_ERR;}
	  fileName = file->name;

   /* Check that storage is provided for all fixed dimensioned matrices */
      for (i=0; i<dim_c; i++) {
         ndim = Dymola_abs(c[i].matrix->nrow) * Dymola_abs(c[i].matrix->ncol);
         if ( ndim > 0 && c[i].matrix->data.v == NULL ) {
            if ( c[i].matrix->name != NULL ) {
               name = c[i].matrix->name;
            } else {
               name = "????";
            }
            goto error33;
         }

         if ( c[i].rowDim > 0 && c[i].matrix->nrow != 0 ) {j=i; goto error35;}
         if ( c[i].colDim > 0 && c[i].matrix->ncol != 0 ) {j=i; goto error36;} 

         c[i].onFile = FALSE_;
      }

   /* Read all matrices from file */
      while ( (name=amatGetName(file,&eof)) != NULL ) {
         /* Find matrix name */
            found = FALSE_;
            for (i=0; i<dim_c; i++) {
               if ( strcmp(c[i].matrix->name,name) == 0) {
                  found = TRUE_;
                  j     = i;
                  break;
               }
            }
            if ( !found ) goto error34;

         /* Read matrix j from file */
            if( (err = amatGetMatrix(file, c[j].matrix)) > 1 ) {
                if ( err==3 ) err=4; else err=5;
                goto errorExit;
            }
            c[j].onFile = TRUE_;

         /* Free allocated storage */
            free(name);
            name = NULL;  

      } /* End of  while */

      if ( !eof ) {err=4; goto errorExit;}


   /* Check additional semantic requirements */
      for (i=0; i<dim_c; i++) {
         if ( c[i].matrix->nrow == 0  ||  c[i].matrix->ncol == 0 ) c[i].onFile = FALSE_;
         if ( c[i].onFile ) {
            /* matrix i was read from file */
               if ( c[i].rowDim > 0 ) {
                  /* row dimension of matrix i must be equal to other dimension */
                     j = c[i].rowIndex;
                     if ( j < 0 || j >= dim_c ) goto error42;
                     if ( !c[j].onFile ) goto error37;
                     if ( c[i].rowDim==1 && c[i].matrix->nrow!=c[j].matrix->nrow ) {
                        str1 = "row"; str2 = "row"; goto error38;
                     }
                     if ( c[i].rowDim==2 && c[i].matrix->nrow!=c[j].matrix->ncol ) {
                        str1 = "row"; str2 = "column"; goto error38;
                     }
               }
               if ( c[i].colDim > 0 ) {
                  /* column dimension of matrix i must be equal to other dimension */
                     j = c[i].colIndex;
                     if ( j < 0 || j >= dim_c ) goto error42;
                     if ( !c[j].onFile ) goto error37;
                     if ( c[i].colDim==1 && c[i].matrix->ncol!=c[j].matrix->nrow ) {
                        str1 = "column"; str2 = "row"; goto error38;
                     }
                     if ( c[i].colDim==2 && c[i].matrix->ncol!=c[j].matrix->ncol ) {
                        str1 = "column"; str2 = "column"; goto error38;
                     }
               }
         } else if ( c[i].req ) {
            /* matrix i not on file as required */
               goto error41;
         }
      }


   /* Close reading of file */
      amatGetClose(file);
      return 0;


   /* Error handling */
      error33   : sprintfC(amatError, amatg_table[33], fileName, name, (int) ndim);
                  goto errorClose;

      error34   : sprintfC(amatError, amatg_table[34], name, fileName);
                  goto errorClose;

      error35   : sprintfC(amatError, amatg_table[35], fileName, (int) j, 
                          (int) c[j].matrix->nrow, (int) j);
                  goto errorClose;

      error36   : sprintfC(amatError, amatg_table[36], fileName, (int) j, 
                          (int) c[j].matrix->ncol, (int) j);
                  goto errorClose;

      error37   : sprintfC(amatError, amatg_table[37], c[i].matrix->name,
                          c[i].matrix->nrow, c[i].matrix->ncol, fileName, c[j].matrix->name);
                  goto errorClose;

      error38   : sprintfC(amatError, amatg_table[38], c[i].matrix->name,
                          c[i].matrix->nrow, c[i].matrix->ncol, fileName, str1, str2,
                          c[j].matrix->name, c[j].matrix->nrow, c[j].matrix->ncol);
                  goto errorClose;

      error41   : sprintfC(amatError, amatg_table[41], c[i].matrix->name, fileName);
                  goto errorClose;

      error42   : sprintfC(amatError, amatg_table[42], fileName, (int) j, (int) i, 
                          c[i].matrix->name, (int) dim_c);
                  goto errorClose;
      
      errorClose: err = 5;
      errorExit : amatGetClose(file);
                  free(name);
                  return err;
}



int amatGetOpen (const char *fileName, const char *classReq, const char *versionReq, AmatGetFile *file) {

   /* Open file for reading of matrices */

   /* Declarations */
      FILE        *fp;
      int          c_first;
      size_t       len;
      const char        *vers = "1";
      int          ic;
      Amatrix     *classInfo;
      int          err, i1, i2;
      const char        *str;

   /* Initialize file */
	  if (!file || !fileName)  {amatUnexpectedNull(); return RET_ERR;}
      file->name    = NULL;
      file->fp      = NULL;
      file->format  = amatASCII;
      file->bstruct = binNormal;
      file->header  = FALSE_;
      file->prec    = 0;
      classInfo = &(file->classInfo);
	  if (!classInfo) {amatUnexpectedNull(); return 13;}
      amatInit(classInfo);

   /* Open file, from which input data must be read */ 
	  if ( (fp = fopen(fileName, "r")) == NULL ) goto error7;
      file->fp = fp;

   /* Read first character */
      if ( (c_first = getc(fp)) == EOF ) goto error9;

   /* Determine file format (formatted or binary) */
      if ( c_first == '#' ) {
         file->format = amatASCII;
         /* read next character and check version number */
            ic = getc(fp);
            if ( ic == EOF ) goto error0;
            if ( ic != vers[0] ) goto error21;
      } else { 
         file->format = amatBinary;
#ifndef UNIX
         /* binary file must be opened in "binary modus" */
            fclose(fp);
            if ( (fp = fopen(fileName, "rb")) == NULL ) goto error7;
            file->fp = fp;
#else
         ungetc(c_first, fp);
#endif
      }

   /* Allocate storage for file-name and store it in "file" */
      len = strlen(fileName) + 1;
      file->name = (char *) malloc( len*sizeof(char) );
      if ( file->name == NULL ) goto error15;
      strcpy( file->name, fileName );

   /* If no class description should be read, terminate function */
      if ( classReq != NULL && strcmp(classReq,"NoClass") == 0 ) return 0;

   /* Read class information (= matrix "classInfo") */

   /* Note: the classInfo->name is 0 and will therefore be allocated in amatGetMatrix */
      classInfo->name = NULL;
      classInfo->type = charMatrix;
      if ( (err=amatGetMatrix(file, classInfo)) != 0 ) {
		 if (classInfo->name) free(classInfo->name);
         amatGetClose(file);
         if (err == RET_MZERO)  return 1;
         if (err == RET_REQ  )  return 1;
         if (err == RET_EOF  )  return 4;
         return 5;
      }

   /* Check that name of matrix is "class" or "Aclass" */
      i1 = strcmp(classInfo->name, "class");
      i2 = strcmp(classInfo->name, "Aclass");
      if ( i1 != 0  &&  i2 != 0 ) goto error43;

	  /* De allocate classInfo and replace with static data; consistent with amatGetClose */
	  free(classInfo->name);
	  if (i1==0) classInfo->name="class";
	  else if (i2==0) classInfo->name="Aclass";
	  else classInfo->name="";

   /* Check class name */
      if ( classInfo->nrow < 3 ) goto error17;
      if ( classReq != NULL ) {
         if ( strcmp(classInfo->data.c[0],classReq) != 0 ) {str=classReq; goto error19;}
      } else {
         str = "";
         goto error19;
      }

   /* Check version number */
      if ( versionReq != NULL ) {
         if ( strcmp(classInfo->data.c[1],versionReq) != 0 ) goto error16;
      }

   /* Inquire bstruct flag if available */
      if ( classInfo->nrow >= 4 ) {
         if ( strcmp(classInfo->data.c[3],"binTrans") == 0 ) {
            file->bstruct = binTrans;
         }
      }
                        
   return 0;


   /* Error handling */
      error0 : if ( feof(fp) ) {
                  sprintfC(amatError, amatg_table[0], fileName, "End-Of-File reached."); 
                  amatGetClose(file);
                  return 4;
               } else { 
                  sprintfC(amatError, amatg_table[0], fileName, strerror(errno));
                  amatGetClose(file);
                  return 5;
               }

      error7 : sprintfC(amatError, amatg_table[7], fileName, strerror(errno));
               return 5;

      error9 : if ( feof(fp) ) {
                  sprintfC(amatError, amatg_table[9], fileName, "End-Of-File reached."); 
               } else { 
                  sprintfC(amatError, amatg_table[9], fileName, strerror(errno));
               }
               amatGetClose(file);
               return 5;

      error15: sprintfC(amatError, amatg_table[15], fileName, strerror(errno));
               fclose(fp);
               return 5;

      error16: if (!classInfo) {amatUnexpectedNull(); return 13;}
               sprintfC(amatError, amatg_table[16], fileName, classInfo->data.c[1], versionReq);
               if ( classInfo->data.c != NULL ) {
                   amatTextDel(classInfo->data.c, classInfo->nrow);
                   classInfo->data.c = NULL;
			   }
               amatGetClose(file);
               return 3;

      error17: if (!classInfo) {amatUnexpectedNull(); return 13;}
               sprintfC(amatError, amatg_table[17], fileName);
               if ( classInfo->data.c != NULL ) {
                  amatTextDel(classInfo->data.c, classInfo->nrow);
                  classInfo->data.c = NULL;
               }
               goto END;

      error19: if (!classInfo) {amatUnexpectedNull(); return 13;}
               sprintfC(amatError, amatg_table[19], fileName, classInfo->data.c[0], str);
               if ( classInfo->data.c != NULL ) {
                  amatTextDel(classInfo->data.c, classInfo->nrow);
                  classInfo->data.c = NULL;
               }
               amatGetClose(file);
               return 2;

      error21: sprintfC(amatError, amatg_table[21], fileName, vers, ic);
               goto END;

      error43: if (!classInfo) {amatUnexpectedNull(); return 13;}
               sprintfC(amatError, amatg_table[43], fileName, classInfo->name);
               amatGetClose(file);
               return 2;

      END    : amatGetClose(file);
               return 5;

} /* End of  amatGetOpen */





void amatGetClose (AmatGetFile *file) {

   /* Close file */
      Amatrix *mat;

      if ( file && file->fp != NULL ) {
		 fclose (file->fp);   
         free   (file->name);
         /* file->classInfo.name is a static string, must not be freed */
            mat = &(file->classInfo);
            amatTextDel(mat->data.c, mat->nrow);
            amatInit(mat);
         file->fp   = NULL;
         file->name = NULL;
      }

} /* End of  amatGetClose */



char *amatGetClass (AmatGetFile *file) {
   
   /* Get class name from file */
      char    *p;
	  if (!file)  {amatUnexpectedNull(); return 0;}
      p = file->classInfo.data.c[0];
      file->classInfo.data.c[0] = NULL;
      return p;
}



char *amatGetDescr (AmatGetFile *file) {
   
   /* Get description string from file */
      char    *p;
	  if (!file)  {amatUnexpectedNull(); return 0;}
      p = file->classInfo.data.c[2];
      file->classInfo.data.c[2] = NULL;
      return p;
}




char *amatGetVersion (AmatGetFile *file) {
   
   /* Get version string from file */
      char *p;
	  if (!file)  {amatUnexpectedNull(); return 0;}
      p = file->classInfo.data.c[1];
      file->classInfo.data.c[1] = NULL;
      return p;
}




#define LINE_LEN 200
#define NAME_LEN 100
#define TYPE_LEN  10 


int amatGetMatrix (AmatGetFile *file, Amatrix *matrixReq) {

   return amatGetMatrixP (file, matrixReq, voidMatrix);

} /* End of  amatGetMatrix */


int amatGetMatrixP3 (AmatGetFile *file, Amatrix *matrixReq, AmatType typeReq, const fpos_t*dataCont,int removeSpace) {

   /* Read numeric or text matrix from file.
      File may be in binary or formatted form */

   /* Declaration */
      Amatrix *matrix;                  /* matrix read from file */
      AmatType typeOnFile;
      char     name[NAME_LEN+1];
      int      err=0, prec;
      size_t   len;
	  integer  nRowPositive;
      integer  lname = NAME_LEN;
      char     *ch;

	  if (!file || !matrixReq)  {amatUnexpectedNull(); return RET_ERR;}
   /* Read matrix declaration */
      matrix = &file->matrix;
      if ( file->header ) {
         /* Matrix declaration already read from file */
            file->header = FALSE_;
            matrix->name = matrixReq->name;
            matrix->data = matrixReq->data;
            prec         = file->prec;

      } else { 
         /* Initialize utility matrix */
             matrix->name = matrixReq->name;
             matrix->nrow = 0;
             matrix->ncol = 0;
             matrix->type = voidMatrix;
             matrix->data = matrixReq->data;

         /* Read header from file */
             switch ( file->format ) {
                case amatASCII : if ( (err=amatGetFHeader(file, matrix, name, lname)) != 0 ) goto END;
                                 break;
                case amatBinary: if ( (err=amatGetBHeader(file, matrix, name, lname, &prec)) != 0 ) goto END;
                                 break;
                default        : goto error10;
             }
      }
     

   /* Check header */
      if (dataCont==0) {
		  if ( (err=amatGetCHeader(file, matrix, *matrixReq)) > 1 ) goto END;
	  }


   /* Set actual matrix type */
      typeOnFile = matrix->type;
      switch ( matrixReq->type ) {
         case voidMatrix   :  if ( typeOnFile == realMatrix || typeOnFile == doubleMatrix ) {
                                 if ( typeReq == realMatrix ) {
                                    matrix->type = realMatrix;
                                 } else if ( typeReq == doubleMatrix ) {
                                    matrix->type = doubleMatrix;
                                 }
                              }
                              break;
         case integerMatrix:  matrix->type = integerMatrix; break;
         case realMatrix   :  matrix->type = realMatrix   ; break;
		 case doubleMatrix :
			 if (typeOnFile == realMatrix) {
				 /* reduce according to type on file */ 
				 matrix->type = realMatrix;
			 } else {
				 matrix->type = doubleMatrix;
			 }
			 break;
         case charMatrix   :  matrix->type = charMatrix   ; break;
		 case realRowMatrix : matrix->type = realRowMatrix; break;
		 case doubleRowMatrix : 
			 if (typeOnFile == realMatrix) {
				 /* reduce according to type on file */ 
				 matrix->type = realRowMatrix;
			 } else {
				 matrix->type = doubleRowMatrix;
			 }
			 break;
         default           :  matrix->type = matrixReq->type; goto error30;
      }
      if ( err > 1 ) goto END; else err=0;  /* "err" value from call to amatGetCHeader */
	  
   /* Allocate storage for data */
      if ( matrix->data.v == NULL ) {
		 nRowPositive =matrix->nrow; 		 /* Guard against null-size allocation, also for charMatrix */
		 if (nRowPositive==0) nRowPositive=1;
         len = (size_t) (  nRowPositive* matrix->ncol );
         switch ( matrix->type ) {
            case integerMatrix: matrix->data.i = (integer    *) malloc( len*sizeof(integer   ) ); break;
            case realMatrix   : matrix->data.r = (real       *) malloc( len*sizeof(real      ) ); break;
            case doubleMatrix : matrix->data.d = (doublereal *) malloc( len*sizeof(doublereal) ); break;
			case charMatrix   : matrix->data.c = (char      **) malloc( nRowPositive*sizeof(char *) );	matrix->nrowread = 0; break;
			case realRowMatrix: /* Fall thru */
			case doubleRowMatrix:
				{
					int stopped=0;
					if (maxRows()>0) {
						/* This case is only on a special flag. */
						matrix->nrowallocated = -matrix->nrow-1; /* Actual number of rows */
						if (matrix->nrow>maxRows())
							matrix->nrow=maxRows(); /* Useful number of rows */
					} else {
						if (matrix->nrow<=10) 
							matrix->nrowallocated = 10; /* allocate some more */
						else 
							matrix->nrowallocated = matrix->nrow; /* Since we always allocate at least 10 no need for guard */
					}

					if (matrix->type == realRowMatrix) {
						if (matrix->nrowallocated<0)
							matrix->data.rrow = (real**) malloc(matrix->nrow*sizeof(real*));
						else
							matrix->data.rrow = (real **) malloc(matrix->nrowallocated*sizeof(real*));
						matrix->nrowread = 0;
						if (matrix->data.rrow!=0) {
							int i;
							for(i=0;i<matrix->nrow;i++) {
								if (amatStopRequest!=0 && (*amatStopRequest)()) {
									matrix->data.rrow[i] = 0;
									stopped=1;
								} else
									matrix->data.rrow[i] = (real*) malloc(matrix->ncol*sizeof(real));
								if (matrix->data.rrow[i]==0) {
									for(;i>=0;i--) free(matrix->data.rrow[i]);
									free(matrix->data.rrow);
									matrix->data.rrow=0;
									break;
								}
							}
						}
					} else {
						/* doubleRowMatrix */
						if (matrix->nrowallocated<0)
							matrix->data.drow = (doublereal**) malloc(matrix->nrow*sizeof(double*));
						else
							matrix->data.drow = (doublereal **) malloc(matrix->nrowallocated*sizeof(double*));
						matrix->nrowread = 0;
						if (matrix->data.drow!=0) {
							int i;
							for(i=0;i<matrix->nrow;i++) {
								if (amatStopRequest!=0 && (*amatStopRequest)()) {
									matrix->data.drow[i] = 0;
									stopped=1;
								} else
									matrix->data.drow[i] = (doublereal*) malloc(matrix->ncol*sizeof(double));
								if (matrix->data.drow[i]==0) {
									for(;i>=0;i--) free(matrix->data.drow[i]);
									free(matrix->data.drow);
									matrix->data.drow=0;
									break;
								}
							}
						}
					}
					if (stopped) goto error44;
					break;
				}
            default           : goto error30;
         }
         if ( matrix->data.v == NULL ) goto error5;

      } else if (dataCont) {
		  int realRow = (matrixReq->type == realRowMatrix && matrix->type == realRowMatrix);
		  if (realRow || (matrixReq->type == doubleRowMatrix && matrix->type == doubleRowMatrix)) {
 			  len = matrix->nrow;
			  if (matrixReq->nrowallocated<0) {
				  /* Do not update nrowallocated yet */
				  integer deltaShift=matrix->nrow-(matrixReq->nrowread-matrixReq->nrow);
				  matrix->nrowallocated = -matrix->nrow-1;
				  if (maxRows()>0 && matrix->nrow>maxRows())
					  matrix->nrow=maxRows();
				  deltaShift-=matrix->nrow;
				  if (realRow) {
					  for(;deltaShift>0;--deltaShift) {
						  integer i;
						  real*d2=matrixReq->data.rrow[0];
						  for(i=1;i<matrixReq->nrow;++i) matrixReq->data.rrow[i-1]=matrixReq->data.rrow[i];
						  matrixReq->data.rrow[matrixReq->nrow-1]=d2;
					  }
					  if (matrix->nrow>matrixReq->nrow) {
						  matrixReq->data.rrow = (real**)realloc(matrixReq->data.rrow,matrix->nrow*sizeof(real*));
						  /* Note: Should only be able to shrink since matrix->nrow will decrease if greater than maxRows but not increase */
					  }
				  } else {
					  /* doubleRowMatrix */
					  for(;deltaShift>0;--deltaShift) {
						  integer i;
						  double*d2=matrixReq->data.drow[0];
						  for(i=1;i<matrixReq->nrow;++i) matrixReq->data.drow[i-1]=matrixReq->data.drow[i];
						  matrixReq->data.drow[matrixReq->nrow-1]=d2;
					  }
					  if (matrix->nrow>matrixReq->nrow) {
						  matrixReq->data.drow = (doublereal**)realloc(matrixReq->data.drow,matrix->nrow*sizeof(double*));
						  /* Note: Should only be able to shrink since matrix->nrow will decrease if greater than maxRows but not increase */
					  }
				  }
			  } else {
				  if (matrix->nrow > matrixReq->nrowallocated) {
					  integer toAllocate=Dymola_max(matrix->nrow,matrixReq->nrowallocated*2);
					  if (realRow) {
						  real**rr=0;
						  rr = (real**)realloc(matrixReq->data.rrow, toAllocate*sizeof(real*));
						  if (rr) {
							  matrixReq->data.rrow = rr;
							  matrixReq->nrowallocated = matrix->nrowallocated = toAllocate;
						  } else {
							  /* Cannot grow. Do nothing more. */
							  goto error5;
						  }
					  } else {
						  /* doubleRowMatrix */
						  doublereal**dd=0;
						  dd = (doublereal**)realloc(matrixReq->data.drow, toAllocate*sizeof(double*));
						  if (dd) {
							  matrixReq->data.drow = dd;
							  matrixReq->nrowallocated = matrix->nrowallocated = toAllocate;
						  } else {
							  /* Cannot grow. Do nothing more. */
							  goto error5;
						  }
					  }
				  } else {
					  matrix->nrowallocated = matrixReq->nrowallocated;
				  }
			  }
			  if (realRow) {
				  matrix->data.rrow = matrixReq->data.rrow;
				  matrix->nrowread = matrixReq->nrowread;
				  if (matrix->nrow > matrixReq->nrow && matrix->data.rrow) {
					  int i;
					  for(i=matrixReq->nrow;i<matrix->nrow;++i) {
						  matrix->data.rrow[i] = (real*) malloc(matrix->ncol*sizeof(real));
						  if (matrix->data.rrow[i]==0) {
							  for(;i>=matrixReq->nrow;i--) 
								  free(matrix->data.rrow[i]);
							  matrix->nrow=matrixReq->nrow; /* Has not grown */
							  goto error5;
						  }
					  }
				  }
			  } else {
				  /* doubleRowMatrix */
				  matrix->data.drow = matrixReq->data.drow;
				  matrix->nrowread = matrixReq->nrowread;
				  if (matrix->nrow > matrixReq->nrow && matrix->data.drow) {
					  int i;
					  for(i=matrixReq->nrow;i<matrix->nrow;++i) {
						  matrix->data.drow[i] = (doublereal*) malloc(matrix->ncol*sizeof(double));
						  if (matrix->data.drow[i]==0) {
							  for(;i>=matrixReq->nrow;i--) 
								  free(matrix->data.drow[i]);
							  matrix->nrow=matrixReq->nrow; /* Has not grown */
							  goto error5;
						  }
					  }
				  }
			  }
			  if (matrix->data.v == NULL) goto error5;
			  if (fsetpos(file->fp, dataCont)) goto error5; /* Or somewhere else?*/
		  } else {
         /* Check, whether nrow, ncol and type are user-defined */
            if ( (matrixReq->nrow == 0) || (matrixReq->ncol == 0 && matrixReq->type != charMatrix) 
              || (matrixReq->type == voidMatrix) ) goto error32;
		  }
	  }


   /* switch according to file format */
      switch ( file->format ) {
         case amatASCII : err = amatGetFMatrix(file, matrix);
                          break;
         case amatBinary: err = amatGetBMatrix(file, matrix, typeOnFile, prec,removeSpace);
                          break;
         default        : goto error10;
      }


   /* Copy matrix from file into output argument */
   END:
      switch ( err ) {
         case 0 : /* No error occured */
                     if ( matrixReq->name == NULL ) matrixReq->name = matrix->name;
                     matrixReq->nrow = matrix->nrow;
                     matrixReq->ncol = matrix->ncol;
                     matrixReq->type = matrix->type;
					 matrixReq->nrowread = matrix->nrowread;
					 matrixReq->nrowallocated = matrix->nrowallocated;
                     if ( matrixReq->data.v == NULL ) matrixReq->data = matrix->data;
                     return 0;

         case 1 : /* Matrix has zero elements */
                     if ( matrixReq->name == NULL ) matrixReq->name = matrix->name;
                     matrixReq->nrow = 0;
                     matrixReq->ncol = 0;
					 matrixReq->nrowread = matrix->nrowread;
					 matrixReq->nrowallocated = matrix->nrowallocated;
                     matrixReq->type = matrix->type;
                     return RET_MZERO;

         default: /* Error, remove allocated storage */
                     if ( matrixReq->name   == NULL ) free(matrix->name);
                     if ( matrixReq->data.v == NULL ) free(matrix->data.v);
                     if ( err == RET_REQ ) return RET_REQ;
                     if ( err == RET_EOF ) return RET_EOF;
                     return RET_ERR;
      }


   /* Error handling */
      error5 : sprintfC(amatError, amatg_table[5], matrix->name, (int) matrix->nrow, 
                       (int) matrix->ncol, file->name, (int) len);
               err = RET_ERR;
               goto END;

      error10: sprintfC(amatError, amatg_table[10], file->name, (int) file->format);
               err = RET_ERR;
               goto END;

      error30: sprintfC(amatError, amatg_table[30], matrix->name, (int) matrix->nrow, 
                       (int) matrix->ncol, file->name, (int) matrix->type);
               err = RET_ERR;
               goto END;
 
      error32: switch (matrixReq->type) {
                  case voidMatrix   : ch = "voidMatrix"   ; break;
                  case integerMatrix: ch = "integerMatrix"; break;
                  case realMatrix   : ch = "realMatrix"   ; break;
                  case doubleMatrix : ch = "doubleMatrix" ; break;
                  case charMatrix   : ch = "charMatrix"   ; break;
				  case realRowMatrix : ch ="realRowMatrix"; break;
				  case doubleRowMatrix : ch ="doubleRowMatrix"; break;
                  default           : ch = "unknown"      ; break;
               }
               sprintfC(amatError, amatg_table[32], matrix->name, (int) matrix->nrow, 
                       (int) matrix->ncol, file->name, (int) matrixReq->nrow,
                       (int) matrixReq->ncol, ch);
               err = RET_ERR;
			   goto END;
      error44:
			   sprintfC(amatError, amatg_table[44], file->name);
			   err = RET_ERR;
               goto END;
 
} /* End of  amatGetMatrix */
int amatGetMatrixP  (AmatGetFile *file, Amatrix *matrixReq, AmatType typeReq) {
	return amatGetMatrixP3(file,matrixReq,typeReq,0,1);
}
int amatGetMatrixP2  (AmatGetFile *file, Amatrix *matrixReq, AmatType typeReq,const fpos_t*dataCont) {
	return amatGetMatrixP3(file,matrixReq,typeReq,dataCont,1);
}


char *amatGetName (AmatGetFile *file, int *eof) {

   /* Get name of next matrix from file.
      File may be in binary or formatted form */

   /* Declaration */
      Amatrix *matrix;                  /* matrix read from file */
      char     name[NAME_LEN+1];
      int      err, prec=0;
      integer  lname = NAME_LEN;

   /* Initialize utility matrix */
	  if (!file || !eof)  {amatUnexpectedNull(); return 0;}
      *eof = FALSE_;
      file->header   = FALSE_;
      matrix         = &file->matrix;
      matrix->name   = NULL;
      matrix->nrow   = 0;
      matrix->ncol   = 0;
      matrix->type   = voidMatrix;
      matrix->data.v = NULL;

   /* Read header from file */
      switch ( file->format ) {
         case amatASCII : if ( (err=amatGetFHeader(file, matrix, name, lname)) != 0 ) {
                             if ( err == RET_EOF ) *eof = TRUE_;
                             return NULL;
                          }
                          break;
         case amatBinary: if ( (err=amatGetBHeader(file, matrix, name, lname, &prec)) != 0 ) {
                             if ( err == RET_EOF ) *eof = TRUE_;
                             return NULL;
                          }
                          break;
         default        : sprintfC(amatError, amatg_table[10], file->name, (int) file->format);
                          return NULL;
      }

   /* Return matrix name */
      file->header = TRUE_;
      file->prec   = prec;
      return matrix->name;

} /* End of  amatGetName */



static int amatGetFMatrix (AmatGetFile *file, Amatrix *matrix) {

   /* Read numeric or text matrix in formatted format from file */

   /* Read data */
	if (!file||!matrix)  {amatUnexpectedNull(); return RET_ERR;}
      switch ( matrix->type ) {
	     case realRowMatrix:
	     case doubleRowMatrix:
         case integerMatrix:
         case realMatrix   : 
         case doubleMatrix : return amatGetFNMatrix (file, matrix);
         case charMatrix   : return amatGetFText    (file, matrix);
      }

   /* Error handling */
      sprintfC(amatError, amatg_table[11], file->name, (int) matrix->type);
      return RET_ERR;

} /* End of  amatGetFMatrix */


static int amatGetFHeader (AmatGetFile *file, Amatrix *matrix, char name[], integer lname){

   /* Read matrix header from file in formatted format and 
      check matrix attributes
      <- RETURN: = 0: No error occured.
                 = 1: The matrix on file has zero dimensions (nrow*ncol=0). 
                 = 2: The matrix on file does not have the required attributes.
                 = 3: End-of-File reached (file was already at EOF when calling this
                      function).
                 = 4: Another error occured.
   */
    
   /* Declarations */
      FILE    *fp;
      char     line[LINE_LEN], type[TYPE_LEN+1];
      char     fmt[50];
      size_t   len;
      int      err;
	  	
	  if (!file||!matrix)  {amatUnexpectedNull(); return RET_ERR;}

   /* Build sscanf format string */
      sprintfC(fmt, "%%%ds %%%d[^( 	](%%d,%%d)",(int) TYPE_LEN, (int) lname);
	  /* Gives: % <TYPE_LEN> s % <lname> [^( ....](%d, %d) */

   /* Read next input line and scan declaration */
      fp = file->fp;
      if ( (err=amatGetSkip(file)) != 0 ) return err;
      if ( fgets(line,LINE_LEN,fp) == NULL ) goto error0;
      if ( sscanfCtext_ssdd(line, fmt, type, name, &matrix->nrow, &matrix->ncol) != 4 ) goto error22;

   /* Determine matrix type */
      if ( strcmp(type, "int") == 0 ) {
         matrix->type = integerMatrix;
      } else if ( strcmp(type, "float") == 0 ) {
         matrix->type = realMatrix;
      } else if ( strcmp(type, "double") == 0 ) {
         matrix->type = doubleMatrix;
      } else if ( strcmp(type, "char") == 0 ) {
         matrix->type = charMatrix;
      } else {
         goto error22;
      }

   /* Check matrix name */
      if ( matrix->name != NULL ) {
         if ( strcmp(matrix->name, name) != 0 ) goto error24;

      } else {
         len = strlen(name);
         if ( len > 0ul ) {
            len += 1;
            matrix->name = (char *) malloc( len*sizeof(char) );
            if ( matrix->name == NULL ) goto error18;
            strcpy(matrix->name, name);
         } else {
           goto error22;
         }
      }

   return 0;

   /* Error handling */
               
      error0 : if (!file||!matrix)  {amatUnexpectedNull(); return RET_ERR;}
               if ( feof(fp) ) {
                  sprintfC(amatError, amatg_table[0], file->name, "End-Of-File reached.");
                  return RET_EOF;
               } else { 
                  sprintfC(amatError, amatg_table[0], file->name, strerror(errno));
                  return RET_ERR;
               }

      error18: if (!file||!matrix)  {amatUnexpectedNull(); return RET_ERR;}
               sprintfC(amatError, amatg_table[18], file->name, (int) len);
               return RET_ERR;

      error22: if (!file||!matrix)  {amatUnexpectedNull(); return RET_ERR;}
               sprintfC(amatError, amatg_table[22], file->name, line);
               return RET_ERR;

      error24: if (!file||!matrix)  {amatUnexpectedNull(); return RET_ERR;}
               sprintfC(amatError, amatg_table[24], file->name, line, matrix->name);
               return RET_REQ;

} /* End of  amatGetFHeader */



static int amatGetCHeader (AmatGetFile *file, Amatrix *matrix, Amatrix matrixReq) {

   /* Check header of matrix 
      -> file      : File identificator.
      -> matrixReq : Required matrix attributes.
      <- matrix    : Actual matrix attributes.
      <- RETURN    : = 0: no error.
                     = 1: header fully checked; number of matrix elements on  
                          file is zero.
                     = 3: error.
   */

   /* Declarations */
      char    *ch;

      if (!file||!matrix)  {amatUnexpectedNull(); return RET_ERR;}
   /* Check matrix type */
      switch ( matrixReq.type ) {
	     case realRowMatrix:
	     case doubleRowMatrix :
         case integerMatrix:
         case realMatrix   :  
         case doubleMatrix : if ( (matrix->type != integerMatrix) &&
                                  (matrix->type != realMatrix   ) &&
                                  (matrix->type != doubleMatrix )    ) goto error25;
                             break;
         case charMatrix   : if ( matrix->type != charMatrix ) goto error26;
                             break;
         default           : break;   /* do nothing */
      }

   /* Check row and column dimensions */
      if ( matrix->nrow < 0 ) goto error1;
      if ( matrix->ncol < 0 ) goto error2;
      if ( (matrix->nrow)*(matrix->ncol) == 0 ) goto error29;
      if ( matrixReq.nrow > 0 && matrix->nrow != matrixReq.nrow )     goto error27;
      if ( matrixReq.ncol > 0 && matrix->ncol != matrixReq.ncol )     goto error28;
      if ( matrixReq.nrow < 0 && matrix->nrow > Dymola_abs(matrixReq.nrow) ) goto error39;
      if ( matrixReq.ncol < 0 && matrix->ncol > Dymola_abs(matrixReq.ncol) ) goto error40;

   return 0;


   /* Error handling */
      error1 : if (!file||!matrix)  {amatUnexpectedNull(); return RET_ERR;}
               sprintfC(amatError, amatg_table[1], (int) matrix->nrow, matrix->name, file->name);
               return RET_ERR;

      error2 : if (!file||!matrix)  {amatUnexpectedNull(); return RET_ERR;}
               sprintfC(amatError, amatg_table[2], (int) matrix->ncol, matrix->name, file->name);
               return RET_ERR;

      error25: if (!file||!matrix)  {amatUnexpectedNull(); return RET_ERR;}
               switch ( matrix->type ) {
                  case charMatrix: ch = "char";
                                   break;
                  default        : ch = "unknown";
               }
               sprintfC(amatError, amatg_table[25], matrix->name, (int) matrix->nrow,
                       (int) matrix->ncol, file->name, ch); 
               return RET_REQ;

      error26: if (!file||!matrix)  {amatUnexpectedNull(); return RET_ERR;}
               switch ( matrix->type ) {
                  case integerMatrix: ch = "integer"; break;
                  case realRowMatrix: /* Fall thru */
                  case realMatrix   : ch = "real";    break;
				  case doubleRowMatrix: /* Fall thru */
                  case doubleMatrix : ch = "double";  break;
                  default           : ch = "unknown";
               }
               sprintfC(amatError, amatg_table[26], matrix->name, (int) matrix->nrow,
                       (int) matrix->ncol, file->name, ch); 
               return RET_REQ;

      error27: if (!file||!matrix)  {amatUnexpectedNull(); return RET_ERR;}
               sprintfC(amatError, amatg_table[27], matrix->name, (int) matrix->nrow,
                       (int) matrix->ncol, file->name, (int) matrixReq.nrow);
               return RET_REQ;

      error28: if (!file||!matrix)  {amatUnexpectedNull(); return RET_ERR;}
               sprintfC(amatError, amatg_table[28], matrix->name, (int) matrix->nrow,
                       (int) matrix->ncol, file->name, (int) matrixReq.ncol);
               return RET_REQ;

      error29: if (!file||!matrix)  {amatUnexpectedNull(); return RET_ERR;}
               sprintfC(amatError, amatg_table[29], matrix->name, (int) matrix->nrow,
                       (int) matrix->ncol, file->name);
               return RET_MZERO;

      error39: 
               if (!file||!matrix)  {amatUnexpectedNull(); return RET_ERR;}
			   sprintfC(amatError, amatg_table[39], matrix->name, (int) matrix->nrow,
                       (int) matrix->ncol, file->name, (int) Dymola_abs(matrixReq.nrow));
               return RET_REQ;

      error40: 
               if (!file||!matrix)  {amatUnexpectedNull(); return RET_ERR;}
               sprintfC(amatError, amatg_table[40], matrix->name, (int) matrix->nrow,
                       (int) matrix->ncol, file->name, (int) Dymola_abs(matrixReq.ncol));
               return RET_REQ;

} /* End of  amatGetCHeader */



static int amatGetFNMatrix (AmatGetFile *file, Amatrix *matrix) {

   /* Read formatted matrix-data from file */

   /* Declarations */
      FILE       *fp;
      integer     i, j;
      integer    *ip1, *ip2;
      real       *rp1, *rp2;
      doublereal *dp1, *dp2, dp;

   /* Read numeric data of matrix from file */
	  if (!file||!matrix)  {amatUnexpectedNull(); return RET_ERR;}

      fp  = file->fp;
      if ( matrix->type == integerMatrix ) {
         ip1 = matrix->data.i;
         ip2 = ip1;
         for ( i=0; i<matrix->nrow; i++) {
            for ( j=0; j<matrix->ncol; j++ ) {
               if ( amatGetSkip(file) != 0 ) return RET_ERR;
               if ( fscanfClf(fp,&dp) != 1 ) goto error6;
               *ip2 = (integer) dp;
               ip2 += matrix->nrow;
            }
            ip1++;
            ip2 = ip1;
         }

      } else if ( matrix->type == realMatrix ) {
         rp1 = matrix->data.r;
         rp2 = rp1;
         for ( i=0; i<matrix->nrow; i++) {
            for ( j=0; j<matrix->ncol; j++ ) {
               if ( amatGetSkip(file) != 0 ) return RET_ERR;
               if ( fscanfCg(fp,rp2) != 1 ) goto error6;
               rp2 += matrix->nrow;
            }
            rp1++;
            rp2 = rp1;
         }
	  } else if ( matrix->type == realRowMatrix ||matrix->type == doubleRowMatrix ) {
		integer readto=matrix->nrow;
		integer offset=0;
		if (matrix->nrowallocated<0) {
			readto=-matrix->nrowallocated-1;
			offset=readto-matrix->nrow;
		}
        for ( i=matrix->nrowread; i<readto; i++) {
			integer index=i-offset;if (index<0) index=0;
			if (matrix->type == realRowMatrix) {
				for ( j=0; j<matrix->ncol; j++ ) {
					double dp;
					if ( amatGetSkip(file) != 0 ) return RET_ERR;
					if ( fscanfClf(fp,&dp) != 1 ) goto error6;
					matrix->data.rrow[index][j] = (real) dp;
				}
			} else {
				/* doubleRowMatrix */
				for ( j=0; j<matrix->ncol; j++ ) {
					if ( amatGetSkip(file) != 0 ) return RET_ERR;
					if ( fscanfClf(fp,&(matrix->data.drow[index][j])) != 1 ) goto error6;
				}
			}
		}
		matrix->nrowread = i;
      } else {
         dp1 = matrix->data.d;
         dp2 = dp1;
         for ( i=0; i<matrix->nrow; i++) {
            for ( j=0; j<matrix->ncol; j++ ) {
               if ( amatGetSkip(file) != 0 ) return RET_ERR;
               if ( fscanfClf(fp,dp2) != 1 ) goto error6;
               dp2 += matrix->nrow;
            }
            dp1++;
            dp2 = dp1;
         }
      }

      return 0;


   /* Error handling */
      error6 : if (!file||!matrix)  {amatUnexpectedNull(); return RET_ERR;}
               if ( feof(fp) ) {
                  sprintfC(amatError, amatg_table[6], matrix->name, (int) matrix->nrow, 
                          (int) matrix->ncol, file->name, ": End-Of-File reached."); 
                  return RET_EOF;
               } else { 
                  sprintfC(amatError, amatg_table[6], matrix->name, (int) matrix->nrow, 
                          (int) matrix->ncol, file->name, " ");
                  return RET_ERR;
               }

} /* End of  amatGetFNMatrix */


static void fixQuotes(char*s) {
	int i=0,j=0;
	for(;;) {
		if (s[i]=='\\' && s[i+1]=='n') {s[j]='\n';i+=2;j++;}
		else if (s[i]=='\\' && s[i+1]=='\\') {s[j]='\\';i+=2;j++;}
		else {s[j]=s[i];if (!s[j]) break;i++;j++;}
	}
}
static int amatGetFText (AmatGetFile *file, Amatrix *matrix) {

   /* Read text matrix in formatted form from file */

   /* Declarations */
      FILE       *fp;
      integer     nrow, ncol;
      size_t      lenStr, len, lenmax;
      integer     i, ii;
      char      **textVec, *str, *p;
      int         c;

	  if (!file||!matrix)  {amatUnexpectedNull(); return RET_ERR;}

   /* Initialize output data */
      textVec = matrix->data.c;
      nrow    = matrix->nrow;
      ncol    = matrix->ncol;
      fp      = file->fp;

   /* Allocate storage for string with maximum length*2 + 2 ('\n\0') */
      lenStr = (size_t) ( ncol*2 + 2 );
      str    = (char *) malloc( lenStr*sizeof(char) );
      if ( str == NULL ) {len=lenStr; goto error13;}
 
   /* Initialize text pointer vector with NULL */
      for (i=0; i<nrow; i++) textVec[i] = NULL;

   /* Read text */
      lenmax = 0;
      for ( i=0; i<nrow; i++) {
         /* Read next string into temporary storage "str" */
            if ( fgets(str,(int)lenStr,fp) == NULL ) goto error14;

			fixQuotes(str);
         /* Determine length of string */
            len = strlen(str);
            if ( str[len-1] == '\n' ) {
               str[len-1] = '\0';
               len--;
            } else {
               /* read upto End-of-Line */ 
                  c = getc(fp);
                  while( c != '\n' && c != EOF ) c = getc(fp);
            }

         /* Allocate storage for string and copy string into storage */
            p = (char *) malloc( (size_t) (len+1)*sizeof(char) );
            if ( p == NULL ) {
               for (ii=0; ii<nrow; ii++) free(textVec[ii]);
               free(str);
               len = len+1;
               goto error13;
            }
            textVec[i] = strcpy(p, str);
            lenmax = Dymola_max(len,lenmax);
      }      

      free(str);
      matrix->ncol = (integer) lenmax;

   return 0;


   /* Error handling */
      error13 : if (!file||!matrix)  {amatUnexpectedNull(); return RET_ERR;}
                sprintfC(amatError, amatg_table[13], matrix->name, (int) matrix->nrow, 
                        (int) matrix->ncol, file->name, (int) len);
                return RET_ERR;

      error14 : if (!file||!matrix)  {amatUnexpectedNull(); return RET_ERR;}
                if ( feof(fp) ) {
                   sprintfC(amatError, amatg_table[14], matrix->name, (int) matrix->nrow, 
                           (int) matrix->ncol,  file->name, "End-Of-File reached."); 
                } else { 
                   sprintfC(amatError, amatg_table[14], matrix->name, (int) matrix->nrow, 
                           (int) matrix->ncol, file->name, strerror(errno));
                }
                for (i=0; i<nrow; i++) free(textVec[i]);
                free(str);
                return RET_ERR;

} /* End of  amatGetFText */



static int amatGetBMatrix (AmatGetFile *file, Amatrix *matrix, AmatType typeOnFile, int prec,int removeSpace){

   /* Read numeric or text matrix in binary format from file */
	if (!file||!matrix)  {amatUnexpectedNull(); return RET_ERR;}
   /* Read data */
      switch ( typeOnFile ) {
         case integerMatrix: return amatGetBIMatrix (file, matrix);
         case realMatrix   : return amatGetBRMatrix (file, matrix);
         case doubleMatrix : return amatGetBDMatrix (file, matrix);
         case charMatrix   : return amatGetBText    (file, matrix, prec, removeSpace);
      }

   /* Error handling */
      sprintfC(amatError, amatg_table[11], file->name, (int) typeOnFile);
      return RET_ERR;

} /* End of  amatGetBMatrix */



static int amatGetBHeader (AmatGetFile *file, Amatrix *matrix,
                           char *name, integer lname, int *prec) {

   /* Read matrix header from file in amatBinary format */

   /* Declarations */
      FILE    *fp;
      Fmatrix  x;
      size_t   one=1, len, len2;
      int      type, i;
      int     mach, machReq;
      integer  ii;

	  if (!file||!matrix||!name)  {amatUnexpectedNull(); return RET_ERR;}
   /* Read struct and determine length of name */
      fp = file->fp;
      if ( fread(&x, sizeof(Fmatrix), one, fp) != one ) goto error0a;

     /* Determine data type of matrix */
	  /* Before actually checking name, since it detects incorrectly formatted files */
      machReq = MATTYPE/1000L;
      mach    = x.type / 1000L;
      if ( mach != machReq ) goto error20;  /* machine type   */

   /* Read and check name */
      if ( matrix->name == NULL ) {
         /* Allocate storage and read name into storage */
            len  = (size_t) x.namlen;
            len2 = Dymola_max(1,len);
            matrix->name = (char *) malloc( len2*sizeof(char) ); 
            if ( matrix->name == NULL ) goto error18;
            if ( len > 0 ) {
               if ( fread(matrix->name, sizeof(char), len, fp) != len ) goto error0;
            } else {
               matrix->name[0] = '\0';
            }

      } else {
         /* Read name into "name" storage. If name on file is too long, truncate it */
            if ( x.namlen > 0 ) {
               if ( x.namlen <= lname+1 ) {
                  if ( fread(name, sizeof(char), x.namlen, fp) != (size_t) x.namlen ) goto error0;
               } else {
                  if ( fread(name, sizeof(char), lname, fp) != (size_t) lname ) goto error0;
                  name[lname] = '\0';
                  for (ii=0; ii<x.namlen-lname; ii++) {if ( fgetc(fp) == EOF ) goto error0;}
               }
            } else {
              name[0] = '\0';
            }
         
         /* Compare name with required name */
            if ( strcmp(matrix->name, name) != 0 ) goto error31;
      }

      i       = (int) ((x.type % 1000L) % 100L);
      *prec   = i / 10;
      type    = i % 10;
     if ( type == 0 ) {                    /* numeric matrix */
         if ( *prec == 0 ) {
            matrix->type = doubleMatrix;
         } else if ( *prec == 1 ) {
            matrix->type = realMatrix;
         } else if ( *prec == 2 ) {
            matrix->type = integerMatrix;
         } else {
            goto error3;
         }
      } else if ( type == 1 ) {             /* character matrix */
         matrix->type = charMatrix;
      } else {
         goto error4;
      }

   /* store row and column dimension */
      if ( file->bstruct == binNormal ) {
         matrix->nrow = x.mrows;
         matrix->ncol = x.ncols;
      } else {
         matrix->nrow = x.ncols;
         matrix->ncol = x.mrows;
      }

   return 0;


   /* Error handling */
      error0a: if (!file||!matrix||!name)  {amatUnexpectedNull(); return RET_ERR;}
               if ( feof(fp) ) {
                  sprintfC(amatError, amatg_table[0], file->name, "End-Of-File reached."); 
                  return RET_EOF;
               } else { 
                  sprintfC(amatError, amatg_table[0], file->name, strerror(errno));
                  return RET_ERR;
               }

      error0 : if (!file||!matrix||!name)  {amatUnexpectedNull(); return RET_ERR;}
               if ( feof(fp) ) {
                  sprintfC(amatError, amatg_table[0], file->name, "End-Of-File reached."); 
               } else { 
                  sprintfC(amatError, amatg_table[0], file->name, strerror(errno));
               }
               return RET_ERR;

      error3 : if (!file||!matrix||!name)  {amatUnexpectedNull(); return RET_ERR;}
               sprintfC(amatError, amatg_table[3], *prec, matrix->name, file->name);
               return RET_ERR;

      error4 : if (!file||!matrix||!name)  {amatUnexpectedNull(); return RET_ERR;}
               sprintfC(amatError, amatg_table[4], type, matrix->name, file->name);
               return RET_ERR;

      error18: if (!file||!matrix||!name)  {amatUnexpectedNull(); return RET_ERR;}
               sprintfC(amatError, amatg_table[18], file->name, (int) len2, strerror(errno));
               return RET_ERR;

      error20: if (!file||!matrix||!name)  {amatUnexpectedNull(); return RET_ERR;}
               sprintfC(amatError, amatg_table[20], file->name);
               return RET_ERR;

      error31: if (!file||!matrix||!name)  {amatUnexpectedNull(); return RET_ERR;}
               sprintfC(amatError, amatg_table[31], name, file->name, matrix->name);
               return RET_REQ;

} /* End of  amatGetBHeader */



static int amatGetBIMatrix (AmatGetFile *file, Amatrix *matrix) {

   /* Read matrix-data from file. The data is stored as integer matrix
      in Matlab binary format */

   /* Declarations */
      FILE       *fp;
      size_t      one=1, ii;
      integer     i, j;
      doublereal *dp1, *dp2;
      real       *rp1, *rp2;
      integer    *ip1, *ip2, ip;
      size_t      len;

	  if (!file||!matrix)  {amatUnexpectedNull(); return RET_ERR;}
 
   /* Read numeric data of matrix from file */
      fp  = file->fp;
      len = (size_t) matrix->nrow * matrix->ncol;

      if ( file->bstruct == binNormal ) {  /* normal matrix structure */
         if ( matrix->type == integerMatrix ) {
            /* get integer matrix in column dense format as integer matrix */
               if ( fread(matrix->data.i, sizeof(integer), len, fp) != len ) goto error6;

         } else if ( matrix->type == realMatrix ) {
            /* get integer matrix in column dense format as real matrix */
               rp1 = matrix->data.r;
               for ( ii=0; ii<len; ii++ ) {
                     if ( fread(&ip, sizeof(integer), one, fp) != one ) goto error6;
                     *rp1++ = (real) ip;
               }
		 } else if ( matrix->type == realRowMatrix ) {
			 for (j=0;j<matrix->ncol;j++) {
				 for(i=0;i<matrix->nrow;i++) {
                    if ( fread(&ip, sizeof(integer), one, fp) != one ) goto error6;
					matrix->data.rrow[i][j]= (real) ip;
				 }
			 }
		 } else if ( matrix->type == doubleRowMatrix ) {
			 for (j=0;j<matrix->ncol;j++) {
				 for(i=0;i<matrix->nrow;i++) {
                    if ( fread(&ip, sizeof(integer), one, fp) != one ) goto error6;
					matrix->data.drow[i][j]= (doublereal) ip;
				 }
			 }
         } else {
            /* get integer matrix in column dense format as doublereal matrix */
               dp1 = matrix->data.d;
               for ( ii=0; ii<len; ii++ ) {
                     if ( fread(&ip, sizeof(integer), one, fp) != one ) goto error6;
                     *dp1++ = (doublereal) ip;
               }
         }

      } else {    /* matrix has to be transposed */
         if ( matrix->type == integerMatrix ) {
            /* get integer matrix in row dense format as integer matrix */
               ip1 = matrix->data.i;
               ip2 = ip1;
               for ( i=0; i<matrix->nrow; i++) {
                  for ( j=0; j<matrix->ncol; j++ ) {
                     if ( fread(ip2, sizeof(integer), one, fp) != one ) goto error6;
                     ip2 += matrix->nrow;
                  }
                  ip1++;
                  ip2 = ip1;
               }

		 } else if ( matrix->type == realRowMatrix || matrix->type == doubleRowMatrix) {
			 int readto=matrix->nrow;
			 int offset=0;
			 if (matrix->nrowallocated<0) {
				 readto=-matrix->nrowallocated-1;
				 offset=readto-matrix->nrow;
			 }
			 if ( matrix->type == realRowMatrix ) {
				 for(i=matrix->nrowread;i<readto;i++) {
					 int index=i-offset;if(index<0) index=0;
					 for (j=0;j<matrix->ncol;j++) {
						 if ( fread(&ip, sizeof(integer), one, fp) != one ) goto error6;
						 matrix->data.rrow[index][j]= (real) ip;
					 }
				 }
			 } else {
				 /* doubleRowMatrix */
				 for(i=matrix->nrowread;i<readto;i++) {
					 int index=i-offset;if(index<0) index=0;
					 for (j=0;j<matrix->ncol;j++) {
						 if ( fread(&ip, sizeof(integer), one, fp) != one ) goto error6;
						 matrix->data.drow[index][j]= (doublereal) ip;
					 }
				 }
			 }
			 matrix->nrowread = i;
         } else if ( matrix->type == realMatrix ) {
            /* get integer matrix in row dense format as real matrix */
               rp1 = matrix->data.r;
               rp2 = rp1;
               for ( i=0; i<matrix->nrow; i++) {
                  for ( j=0; j<matrix->ncol; j++ ) {
                     if ( fread(&ip, sizeof(integer), one, fp) != one ) goto error6;
                     *rp2 =  (real) ip;
                      rp2 += matrix->nrow;
                  }
                  rp1++;
                  rp2 = rp1;
               }

         } else {
            /* get integer matrix in row dense format as doublereal matrix */
               dp1 = matrix->data.d;
               dp2 = dp1;
               for ( i=0; i<matrix->nrow; i++) {
                  for ( j=0; j<matrix->ncol; j++ ) {
                     if ( fread(&ip, sizeof(integer), one, fp) != one ) goto error6;
                     *dp2 =  (doublereal) ip;
                      dp2 += matrix->nrow;
                  }
                  dp1++;
                  dp2 = dp1;
               }
         }
      }

      return 0;


   /* Error handling */
      error6 : if (!file||!matrix)  {amatUnexpectedNull(); return RET_ERR;}
               if ( feof(fp) ) {
                  sprintfC(amatError, amatg_table[6], matrix->name, (int) matrix->nrow, 
                          (int) matrix->ncol, file->name, "End-Of-File reached."); 
               } else { 
                  sprintfC(amatError, amatg_table[6], matrix->name, (int) matrix->nrow, 
                          (int) matrix->ncol, file->name, strerror(errno));
               }
               return RET_ERR;

} /* End of  amatGetBIMatrix */



static int amatGetBRMatrix (AmatGetFile *file, Amatrix *matrix) {

   /* Read matrix-data from file. The data is stored as real matrix
      in Matlab binary format */

   /* Declarations */
      FILE       *fp;
      size_t      one=1, ii;
      integer     i, j;
      doublereal *dp1, *dp2;
      real       *rp1, *rp2, rp;
      integer    *ip1, *ip2;
      size_t      len;
 
	  if (!file||!matrix)  {amatUnexpectedNull(); return RET_ERR;}
   /* Read numeric data of matrix from file */
      fp  = file->fp;
      len = (size_t) matrix->nrow * matrix->ncol;

      if ( file->bstruct == binNormal ) {  /* normal matrix structure */
         if ( matrix->type == realMatrix ) {
            /* get real matrix in column dense format as real matrix */
               if ( fread(matrix->data.r, sizeof(real), len, fp) != len ) goto error6;

         } else if ( matrix->type == integerMatrix ) {
            /* get real matrix in column dense format as integer matrix */
               ip1 = matrix->data.i;
               for ( ii=0; ii<len; ii++ ) {
                     if ( fread(&rp, sizeof(real), one, fp) != one ) goto error6;
                     *ip1++ = (integer) rp;
               }
		 } else if ( matrix->type == realRowMatrix) {
           /* get real matrix in column dense format as real row-stored matrix */
 			 for(j=0;j<matrix->ncol;j++) {
				 for(i=0;i<matrix->nrow;i++) {
					 if (fread(&rp,sizeof(real),one,fp) != one) goto error6;
					 matrix->data.rrow[i][j] = rp;
				 }
			 }
		 } else if ( matrix->type == doubleRowMatrix) {
           /* get real matrix in column dense format as doublereal row-stored matrix */
 			 for(j=0;j<matrix->ncol;j++) {
				 for(i=0;i<matrix->nrow;i++) {
					 if (fread(&rp,sizeof(real),one,fp) != one) goto error6;
					 matrix->data.drow[i][j] = rp;
				 }
			 }
         } else {
            /* get real matrix in column dense format as doublereal matrix */
               dp1 = matrix->data.d;
               for ( ii=0; ii<len; ii++ ) {
                     if ( fread(&rp, sizeof(real), one, fp) != one ) goto error6;
                     *dp1++ = (doublereal) rp;
               }
         }

      } else {    /* matrix has to be transposed */
         if ( matrix->type == realMatrix ) {
            /* get real matrix in row dense format as real matrix */
               rp1 = matrix->data.r;
               rp2 = rp1;
               for ( i=0; i<matrix->nrow; i++) {
                  for ( j=0; j<matrix->ncol; j++ ) {
                     if ( fread(rp2, sizeof(real), one, fp) != one ) goto error6;
                     rp2 += matrix->nrow;
                  }
                  rp1++;
                  rp2 = rp1;
               }

         } else if ( matrix->type == integerMatrix ) {
            /* get real matrix in row dense format as integer matrix */
               ip1 = matrix->data.i;
               ip2 = ip1;
               for ( i=0; i<matrix->nrow; i++) {
                  for ( j=0; j<matrix->ncol; j++ ) {
                     if ( fread(&rp, sizeof(real), one, fp) != one ) goto error6;
                     *ip2 =  (integer) rp;
                      ip2 += matrix->nrow;
                  }
                  ip1++;
                  ip2 = ip1;
               }
		 } else if (matrix->type == realRowMatrix || matrix->type == doubleRowMatrix) {
           /* get real matrix in row dense format as doublereal row-stored matrix */
			   real buffer[20];
			   real *mbuff;
			   integer readto=matrix->nrow;
			   integer offset=0;
			   if (matrix->nrowallocated<0) {
				   readto=-matrix->nrowallocated-1;
				   offset=readto-matrix->nrow;
			   }
			   for(i=matrix->nrowread;i<readto;i++) {
				  integer index=i-offset;if(index<0) index=0;
				  if (amatStopRequest!=0 && (*amatStopRequest)()) {
					goto error44;
				  }
				  if (matrix->type == realRowMatrix) {
					  for ( j=0; j<matrix->ncol; ) {
						  size_t todo = Dymola_min(matrix->ncol-j,sizeof(buffer)/sizeof(*buffer));
						  if ( fread(&buffer[0], sizeof(real), todo, fp) != todo ) goto error6;
						  for(mbuff=buffer;mbuff<buffer+todo;j++,mbuff++) {
							  matrix->data.rrow[index][j] = *mbuff;
						  }
					  }
				  } else {
					  /* doubleRowMatrix */
					  for ( j=0; j<matrix->ncol; ) {
						  size_t todo = Dymola_min(matrix->ncol-j,sizeof(buffer)/sizeof(*buffer));
						  if ( fread(&buffer[0], sizeof(real), todo, fp) != todo ) goto error6;
						  for(mbuff=buffer;mbuff<buffer+todo;j++,mbuff++) {
							  matrix->data.drow[index][j] = (doublereal) *mbuff;
						  }
					  }
				  }
               }
			   matrix->nrowread = i;
         } else {
            /* get real matrix in row dense format as doublereal matrix */
			   real buffer[20];
			   real *mbuff;
               dp1 = matrix->data.d;

               for ( i=0; i<matrix->nrow; i++) {
				  dp2 = dp1;
                  for ( j=0; j<matrix->ncol; ) {
					 size_t todo = Dymola_min(matrix->ncol-j,sizeof(buffer)/sizeof(*buffer));
                     if ( fread(&buffer[0], sizeof(real), todo, fp) != todo ) goto error6;
					 for(mbuff=buffer;mbuff<buffer+todo;j++,mbuff++) {
						*dp2 =  (doublereal) *mbuff;
						dp2 += matrix->nrow;
					 }
                  }
                  dp1++;
               }
         }
      }

      return 0;


   /* Error handling */
      error6 : 
               if (!file||!matrix)  {amatUnexpectedNull(); return RET_ERR;}
               if ( feof(fp) ) {
                  sprintfC(amatError, amatg_table[6], matrix->name, (int) matrix->nrow, 
                          (int) matrix->ncol, file->name, "End-Of-File reached."); 
               } else { 
                  sprintfC(amatError, amatg_table[6], matrix->name, (int) matrix->nrow, 
                          (int) matrix->ncol, file->name, strerror(errno));
               }
               return RET_ERR;
	error44:
			   if (!file||!matrix)  {amatUnexpectedNull(); return RET_ERR;}
			   sprintfC(amatError, amatg_table[44], file->name);
			   return RET_ERR;

} /* End of  amatGetBRMatrix */



static int amatGetBDMatrix (AmatGetFile *file, Amatrix *matrix){

   /* Read matrix-data from file. The data is stored as doublereal matrix
      in Matlab binary format */

   /* Declarations */
      FILE       *fp;
      size_t      one=1, ii;
      integer     i, j;
      doublereal *dp1, *dp2, dp;
      real       *rp1, *rp2;
      integer    *ip1, *ip2;
      size_t      len;

	  if (!file||!matrix)  {amatUnexpectedNull(); return RET_ERR;}

   /* Read numeric data of matrix from file */
      fp  = file->fp;
      len = (size_t) matrix->nrow * matrix->ncol;

      if ( file->bstruct == binNormal ) {  /* normal matrix structure */
         if ( matrix->type == doubleMatrix ) {
            /* get doublereal matrix in column dense format as doublereal matrix */
               if ( fread(matrix->data.d, sizeof(doublereal), len, fp) != len ) goto error6;

		 } else if (matrix->type == realRowMatrix) {
			 /* get a doublereal matrix in column dense format as a real row-stored matrix */
			 for(j=0;j<matrix->ncol;j++) {
				 for(i=0;i<matrix->nrow;i++) {
					 if (fread(&dp,sizeof(doublereal),one,fp) != one) goto error6;
					 matrix->data.rrow[i][j] = (real) dp;
				 }
			 }
		 } else if (matrix->type == doubleRowMatrix) {
			 /* get a doublereal matrix in column dense format as a doublereal row-stored matrix */
			 for(j=0;j<matrix->ncol;j++) {
				 for(i=0;i<matrix->nrow;i++) {
					 if (fread(&(matrix->data.drow[i][j]),sizeof(doublereal),one,fp) != one) goto error6;
				 }
			 }
         } else if ( matrix->type == realMatrix ) {
            /* get doublereal matrix in column dense format as real matrix */
               rp1 = matrix->data.r;
               for ( ii=0; ii<len; ii++ ) {
                     if ( fread(&dp, sizeof(doublereal), one, fp) != one ) goto error6;
                     *rp1++ = (real) dp;
               }

         } else {
            /* get doublereal matrix in column dense format as integer matrix */
               ip1 = matrix->data.i;
               for ( ii=0; ii<len; ii++ ) {
                     if ( fread(&dp, sizeof(doublereal), one, fp) != one ) goto error6;
                     *ip1++ = (integer) dp;
               }
         }

      } else {    /* matrix has to be transposed */
         if ( matrix->type == doubleMatrix ) {
            /* get doublereal matrix in row dense format as doublereal matrix */
			   doublereal buffer[20];
			   doublereal *mbuff;
               dp1 = matrix->data.d;

               for ( i=0; i<matrix->nrow; i++) {
				  dp2 = dp1;
                  for ( j=0; j<matrix->ncol;) {
					 size_t todo = Dymola_min(matrix->ncol-j,sizeof(buffer)/sizeof(*buffer));
                     if ( fread(&buffer[0], sizeof(doublereal), todo, fp) != todo ) goto error6;
					 for(mbuff=buffer;mbuff<buffer+todo;mbuff++,j++) {
						 *dp2 = *mbuff;
                         dp2 += matrix->nrow;
					 }
                  }
                  dp1++;
               }
         } else if ( matrix->type == realRowMatrix || matrix->type == doubleRowMatrix ) {
			 /* get a doublereal matrix in row dense format as a doublereal row-stored matrix */
			 /* no need to transpose */
			 integer readto=matrix->nrow;
			 integer offset=0;
			 if (matrix->nrowallocated<0) {
				 readto=-matrix->nrowallocated-1;
				 offset=readto-matrix->nrow;
			 }
			 if (matrix->type == realRowMatrix) {
				 for(i=matrix->nrowread;i<readto;i++) {
					 integer index=i-offset;if(index<0) index=0;
					 for (j=0;j<matrix->ncol;j++) {
						 if (fread(&dp,sizeof(doublereal), one, fp) != one ) goto error6;
						 matrix->data.drow[index][j] = (real) dp;
					 }
				 }
			 } else {
				 /* doubleRowMatrix */ 
				 for(i=matrix->nrowread;i<readto;i++) {
					 integer index=i-offset;if(index<0) index=0;
					 if (fread(matrix->data.drow[index],sizeof(doublereal), matrix->ncol, fp) != (size_t)(matrix->ncol) ) goto error6;
				 }
			 }
			 matrix->nrowread = i;
         } else if ( matrix->type == realMatrix ) {
            /* get doublereal matrix in row dense format as real matrix */
               rp1 = matrix->data.r;
               rp2 = rp1;
               for ( i=0; i<matrix->nrow; i++) {
                  for ( j=0; j<matrix->ncol; j++ ) {
                     if ( fread(&dp, sizeof(doublereal), one, fp) != one ) goto error6;
                     *rp2 =  (real) dp;
                      rp2 += matrix->nrow;
                  }
                  rp1++;
                  rp2 = rp1;
               }
         } else {
            /* get doublereal matrix in row dense format as integer matrix */
               ip1 = matrix->data.i;
               ip2 = ip1;
               for ( i=0; i<matrix->nrow; i++) {
                  for ( j=0; j<matrix->ncol; j++ ) {
                     if ( fread(&dp, sizeof(doublereal), one, fp) != one ) goto error6;
                     *ip2 =  (integer) dp;
                      ip2 += matrix->nrow;
                  }
                  ip1++;
                  ip2 = ip1;
               }
         }
      }

      return 0;


   /* Error handling */
      error6 : 
               if (!file||!matrix)  {amatUnexpectedNull(); return RET_ERR;}
			   if ( feof(fp) ) {
                  sprintfC(amatError, amatg_table[6], matrix->name, (int) matrix->nrow, 
                          (int) matrix->ncol, file->name, "End-Of-File reached."); 
               } else { 
                  sprintfC(amatError, amatg_table[6], matrix->name, (int) matrix->nrow, 
                          (int) matrix->ncol, file->name, strerror(errno));
               }
               return RET_ERR;

} /* End of  BamatGetBDMatrix */



static int amatGetBText (AmatGetFile *file, Amatrix *matrix, int prec,int removeSpace){

   /* Read text matrix (Matlab binary format) from file */

   /* Declarations */
      FILE       *fp;
      size_t      one=1, lenTemp, lenTempByte, lenVec;
      size_t      lenStr, lenStrByte, len, lenmax, nrow, ncol;
      integer     i, j, ii;
      char       *textTemp, **textVec, *str;
      char       *p, **pVec, *pTemp, *pStr;
      doublereal  dp;

	  if (!file||!matrix)  {amatUnexpectedNull(); return RET_ERR;}
   /* Check for type of text matrix */
      if ( !(prec==0 || prec==5) ) goto error8;

   /* Initialize output data */
      textTemp = NULL;
      textVec  = matrix->data.c;
      lenVec   = (size_t) matrix->nrow;
      nrow     = (size_t) matrix->nrow;
      ncol     = (size_t) matrix->ncol;
      fp       = file->fp;

   /* Read matrix from file */
      if ( file->bstruct == binNormal ) {  /* normal matrix structure */
         /* Read full matrix temporarily */

         /* Allocate storage for temporary text matrix */
            lenTemp     = (size_t) nrow * ncol;
            lenTempByte = lenTemp * sizeof( char );
			textTemp    = (char *) malloc( lenTempByte?lenTempByte:1 ); /* Avoid null-size */
            if ( textTemp == NULL ) {len = lenTemp; goto error13;}

         /* Read text matrix from file */
            if ( prec == 5 ) {
               /* text stored as character */
                  if ( fread(textTemp, one, lenTempByte, fp) != lenTempByte ) {
                     free(textTemp);
                     goto error14;
                  }
            } else {
               /* text stored as doublereal */
                  p = textTemp; 
                  for (i=0; i<(integer) lenTemp; i++) {
                     if ( fread(&dp, sizeof(doublereal), one, fp) != one ) {
                        free(textTemp);
                        goto error14;
                     }
                     *p++ = (char) dp;
                  }
            }
      }

      /* Allocate storage for string with maximum length + 1 */
         lenStr     = (size_t) ( ncol + 1 );
         lenStrByte = lenStr * sizeof( char );
         str        = (char *) malloc( lenStrByte );
         if ( str == NULL ) {
            free(textTemp);
            len = lenStr;
            goto error13;
         }
 
      /* Initialize text pointer vector with NULL */
         for (i=0; i<(integer) lenVec; i++) textVec[i] = NULL;

      /* Read text */
         pVec   = textVec;
         pTemp  = textTemp;
         lenmax = 0;
         for ( i=0; i<(integer) nrow; i++) {
            /* Read next string into temporary storage "str" */
               if ( file->bstruct == binNormal ) {
                  pStr = str;
                  p    = pTemp;
                  for ( j=0; j<(integer) ncol; j++ ) {
                      *pStr++ = *p;
                       p += nrow;
                  }
                  pTemp++;

               } else {
                  if ( prec == 5 ) {
                     /* text stored as character */
                        if ( fread(str, sizeof(char), ncol, fp) != ncol ) {
                           for (ii=0; ii<(integer) lenVec; ii++) free(textVec[ii]);
                           free(textTemp);
                           free(str);
                           goto error14;
                        }
                  } else {
                     /* text stored as doublereal */
                        p = str;
                        for (j=0; j<(integer) ncol; j++) {
                           if ( fread(&dp, sizeof(doublereal), one, fp) != one ) {
                              for (ii=0; ii<(integer) lenVec; ii++) free(textVec[ii]);
                              free(textTemp);
                              free(str);
                              goto error14;
                           }
                           *p++ = (char) dp;
                        }
                  }
               }
         
            /* Remove trailing blanks and determine length of string */
               len = 0;
               for ( j=ncol-1; j>=0; j--) {
				  if ( removeSpace ? (isspaceC(((unsigned char*)str)[j]) == 0) : (str[j]!='\0') ) {
                     len = j + 1;               /* length of string without '\0' */
                     break;
                  }
               }

            /* Allocate storage for string and copy string into storage */
              if ( len > 0ul ) {
                  p = (char *) malloc( (size_t) (len+1)*sizeof(char) );
                  if ( p == NULL ) {
                     for (ii=0; ii<(integer) lenVec; ii++) free(textVec[ii]);
                     free(textTemp);
                     free(str);
                     len = len+1;
                     goto error13;
                  }
                  str[len] = '\0';
                  *pVec = strcpy(p, str);
               }
               pVec++;
               lenmax = Dymola_max(len,lenmax);
         }      

         free(str);
         free(textTemp);
         matrix->ncol = (integer) lenmax;

      return 0;


   /* Error handling */
      error8  : 
                if (!file||!matrix)  {amatUnexpectedNull(); return RET_ERR;}
				sprintfC(amatError, amatg_table[8], prec, matrix->name, 
                        (int) matrix->nrow, (int) matrix->ncol, file->name);
                return RET_ERR;

      error13 : 
                if (!file||!matrix)  {amatUnexpectedNull(); return RET_ERR;}
				sprintfC(amatError, amatg_table[13], matrix->name, (int) matrix->nrow, 
                        (int) matrix->ncol, file->name, (int) len);
                return RET_ERR;

      error14 : 
                if (!file||!matrix)  {amatUnexpectedNull(); return RET_ERR;}
				if ( feof(fp) ) {
                   sprintfC(amatError, amatg_table[14], matrix->name, (int) matrix->nrow, 
                           (int) matrix->ncol,  file->name, "End-Of-File reached."); 
                } else { 
                   sprintfC(amatError, amatg_table[14], matrix->name, (int) matrix->nrow, 
                           (int) matrix->ncol, file->name, strerror(errno));
                }
                return RET_ERR;

} /* End of  amatGetBText */



static int amatGetSkip(AmatGetFile *file) {

   /* Skip blanks, tabs and comments
      <- RETURN: = 0: no error occured.
                 = 3: EOF occured.
                 = 4: other error occured.
   */

      FILE *fp;
      int   ch;
	  if (!file)  {amatUnexpectedNull(); return RET_ERR;}
      fp = file->fp;
      while(1) {
         /* Read next character */
            if ( (ch = getc(fp)) == EOF ) goto error0;

         /* If character = blank, tab, newline then read next character */
            if ( isspaceC(ch) ) continue;

         /* If character = '#' then skip characters upto newline,
            otherwise break loop */
            if ( ch == '#' ) {
               do {
                  if ( (ch = getc(fp)) == EOF ) goto error0;
               } while ( ch != '\n' );
            } else {
               if ( ungetc(ch, fp) == EOF ) goto error0;
               break;
            }
      }
      return 0;

   /* Error handling */
      error0 : 
               if (!file)  {amatUnexpectedNull(); return RET_ERR;}
			   if ( feof(fp) != 0 ) { 
                  sprintfC(amatError, amatg_table[0], file->name, "End-of-File reached");
                  return RET_EOF;
               } else {
                  sprintfC(amatError, amatg_table[0], file->name, strerror(errno));
                  return RET_ERR;
               }

} /* End of  amatGetSkip */







static int amatPutFImatrix (AmatPutFile *file, const char *name, integer nrow, integer ncol, const integer     data[], const char *matDescr, const char *const elemDescr[], int nblank);
static int amatPutFRmatrix (AmatPutFile *file, const char *name, integer nrow, integer ncol, const real        data[], const char *matDescr, const char *const elemDescr[], int nblank);
static int amatPutFRmatrix2(AmatPutFile *file, const char *name, integer nrow, integer ncol, const Amatrix     matrix, const char *matDescr, const char *const elemDescr[], int nblank);
static int amatPutFDmatrix (AmatPutFile *file, const char *name, integer nrow, integer ncol, const doublereal  data[], const char *matDescr, const char *const elemDescr[], int nblank);
static int amatPutFDmatrix2(AmatPutFile *file, const char *name, integer nrow, integer ncol, const Amatrix     matrix, const char *matDescr, const char *const elemDescr[], int nblank);
static int amatPutFCmatrix (AmatPutFile *file, const char *name, integer nrow, integer ncol, const char       *const data[], const char *matDescr);
static int amatPutFRinit   (AmatPutFile *file);
static int amatPutFRclose  (AmatPutFile *file);
static int amatPutFRreal   (AmatPutFile *file, const real row[]);
static int amatPutFRdouble (AmatPutFile *file, const doublereal row[]);

static int amatPutBImatrix (AmatPutFile *file, const char *name, integer nrow, integer ncol, const integer     data[]);
static int amatPutBRmatrix (AmatPutFile *file, const char *name, integer nrow, integer ncol, const real        data[]);
static int amatPutBRmatrix2(AmatPutFile *file, const char *name, integer nrow, integer ncol, Amatrix     matrix);
static int amatPutBDmatrix (AmatPutFile *file, const char *name, integer nrow, integer ncol, const doublereal  data[]);
static int amatPutBDmatrix2(AmatPutFile *file, const char *name, integer nrow, integer ncol, Amatrix     matrix);
static int amatPutBCmatrix (AmatPutFile *file, const char *name, integer nrow, integer ncol, const char       *const data[], char padUsing);
static int amatPutBRinit   (AmatPutFile *file);
static int amatPutBRclose  (AmatPutFile *file);


/* Error messages */
   static const char *const amatp_table[] = {
      /* 0 */   "Error opening file \"%.400s\": %.400s\n"
      /* 1 */  ,"Error writing class description of class \"%.400s\" to\n"
                "file \"%.400s\":\n%.400s\n"
      /* 2 */  ,"Error should never occur when writing matrix \"%.400s(%d,%d)\" to\n"
                "file \"%.400s\" (matrix.type = %d)\n"
      /* 3 */  ,"Error when writing matrix \"%.400s(%d,%d)\" to file \"%.400s\":\n%.400s\n"
      /* 4 */  ,"Error should never occur when writing matrix \"%.400s\" to\n"
                "file \"%.400s\" (file.format = %d)\n"
      /* 5 */  ,"Not possible to allocate enough storage in order to open\n"
                "file \"%.400s\" for writing of matrices:\n%.400s\n"
      /* 6 */  ,"Not possible to allocate enough storage in order to open\n"
                "file \"%.400s\" for writing of matrix \"%.400s\":\n%.400s\n"
      /* 7 */  ,"Wrong argument \"type\" (= %d) in function amatPutRinit,\n"
                "when initializing writing of matrix \"%.400s\" to file \"%.400s\".\n"
                "(type = realMatrix or doubleMatrix required).\n"
      /* 8 */  ,"Error when writing matrix \"%.400s\" row-wise to file \"%.400s\":\n%.400s\n"
      /* 9 */  ,"In order to write matrix \"%.400s\" row-wise to file \"%.400s\",\n"
                "the file has to be opened with \"bstruct = binTrans\" and not\n"
                "with \"bstruct = binNormal\" (function amatPutOpen).\n"
     /* 10 */  ,"Error when writing row %d of matrix \"%.400s\" row-wise\n"
                "to file \"%.400s\": %.400s\n"
     /* 11 */  ,"Error when writing matrix \"%.400s(%d,%d)\" to file \"%.400s\":\n"
                "The datatype of the matrix is not specified (matrix.type = voidMatrix).\n"
     /* 12 */  ,"Not possible to use functions amatPutRxxx, in order to write\n"
               ,"matrix \"%.400s\" row-wise, because output stream is no file, but stdout.\n"
	 /* 13 */  ,"%s is not of type double{Row}matrix"
   };



int amatWrite (char *fileName, AmatFileType ftype, Amatrix matrix) {
 
   /* Write matrix to file */
      AmatPutFile file;
      if ( amatPutOpen   (fileName, ftype, binNormal, "", "", "", &file) != 0 ) return 1;
      if ( amatPutMatrix (&file, matrix) != 0 ) { amatPutClose(&file); return 1; }
      amatPutClose(&file);
      return 0;

} /* End of  amatWrite */




int amatWriteAll (const char *fileName, AmatFileType ftype, char *className,
                  char *version, char *descr, AclassWrite c[], integer dim_c) {

   /* Declarations */
      AmatPutFile file;
      integer     i;

   /* Open file */
      if ( amatPutOpen(fileName, ftype, binNormal, className, version, descr, &file) != 0 ) return 1;

   /* Write matrices to file */
      for (i=0; i<dim_c; i++) {
         if ( amatPutMatrixD(&file, *(c[i].matrix), c[i].matDescr, c[i].elemDescr, c[i].nblank) != 0 ) {
            /* Error */
               amatPutClose(&file);
               return 1;
         }
      }

   /* Close file */
      amatPutClose(&file);
      return 0;

} /* End of  amatWriteAll */


static const fpos_t null_fpos;
/* Needed for initialization.*/

int amatPutOpen (const char *fileName, AmatFileType format, AmatBinStruct bstruct,
				 char *className,  char *version, char *descr, AmatPutFile *file){
   
   /* Open file for writing of matrices */
   
   /* Declarations */
      FILE      *fp;
      const char      *mode;
      Amatrix    classmat;
      char      * classdata[4];
      size_t     len;

	  if (!file)  {amatUnexpectedNull(); return RET_ERR;}

   /* Initialize file */
      file->fp      = NULL;
      file->format  = format;
      file->bstruct = binNormal;
      file->matName = NULL;
      file->posHead = null_fpos;
      file->type    = voidMatrix;
      file->ncol    = 0;
      file->nrowAct = 0;
      file->noRow   = 0;

   /* Determine file mode */
      if ( format == amatBinary ) mode = "wb";
      else                        mode = "w";

   /* Open file */
      if ( fileName == NULL || fileName[0] == '\0' ) {
         /* use stdout */
            fp = stdout;
      } else {
         if ( (fp = fopen(fileName, mode)) == NULL ) goto error0;
         file->noRow = 1;
      }
      file->fp = fp;    

   /* Allocate storage for file-name and store it in "file" */
      if ( fp != stdout ) {
         len = strlen(fileName) + 1;
         file->name = (char *) malloc( len*sizeof(char) );
         if ( file->name == NULL ) goto error5;
         strcpy( file->name, fileName );
      } else {
         file->name = NULL;
      }

   /* Define class */

   /* Write file version number to file, in case of ASCII file format */
      if ( format == amatASCII ) {
         if ( fprintfC(fp,"#1\n") < 0 ) goto error1;
      }

   /* Write class description to file (always with bstruct=binNormal) */
      if ( className != NULL && className[0] != '\0' ) {
         classmat.name   = "Aclass";
         classmat.nrow   = (format == amatASCII) ? 3 : 4;
         classmat.ncol   = 0;
         classmat.type   = charMatrix;
         classmat.data.c = classdata;
         classdata[0]    = className;
         classdata[1]    = version ? version : " ";
         classdata[2]    = descr   ? descr   : " ";
         classdata[3]    = (bstruct == binNormal) ? "binNormal" : "binTrans";
         if ( amatPutMatrix(file, classmat) != 0 ) goto error1;
      }

   /* Set bstruct for next writing of matrices */
      file->bstruct = bstruct;
   
   return 0;

   /* Error handling */ 
      error0: sprintfC(amatError, amatp_table[0], fileName, strerror(errno));
              return 1;

      error1: sprintfC(amatError, amatp_table[1], className, fileName, strerror(errno));
              amatPutClose(file);
              return 1;

      error5: sprintfC(amatError, amatp_table[5], fileName, strerror(errno));
              fclose(fp);
              return 1;

} /* End of  amatPutOpen */



void amatPutClose (AmatPutFile *file) {

   /* Close file */

   if (file && file->fp != NULL ) {
      fclose(file->fp);
      free(file->name);
      free(file->matName);
      file->fp      = NULL;
      file->name    = NULL;
      file->matName = NULL;
   }
}


void amatPutConnect(FILE *fp, char *fileName, AmatPutFile *file) {
	if (!file)  {amatUnexpectedNull(); return;}
   /* Connect AmatPut stream to an open ASCII file stream */
      file->name    = fileName;
      file->fp      = fp;
      file->format  = amatASCII;
      file->bstruct = binNormal;
      file->noRow   = 0;
      file->matName = NULL;
}



static int amatPutMatrixD2(AmatPutFile* file, Amatrix matrix, const char*matDescr, const char*const elemDescr[],int nblank,char padUsing) {
   /* Write matrix on file "file" including a description text */
	  if (!file)  {amatUnexpectedNull(); return RET_ERR;}

      if ( file->format == amatASCII ) {
         /* Write matrix in ASCII format */ 
            switch ( matrix.type ) {
               case voidMatrix   : goto error11;
               case integerMatrix:  return amatPutFImatrix (file, matrix.name, matrix.nrow, matrix.ncol, matrix.data.i, matDescr, elemDescr, nblank);
               case realMatrix   :  return amatPutFRmatrix (file, matrix.name, matrix.nrow, matrix.ncol, matrix.data.r, matDescr, elemDescr, nblank);
               case doubleMatrix :  return amatPutFDmatrix (file, matrix.name, matrix.nrow, matrix.ncol, matrix.data.d, matDescr, elemDescr, nblank);
               case charMatrix   :  return amatPutFCmatrix (file, matrix.name, matrix.nrow, matrix.ncol, (const char   *const*)(matrix.data.c), matDescr); /* cast-away? */
			   case realRowMatrix:  return amatPutFRmatrix2(file, matrix.name, matrix.nrow, matrix.ncol, matrix,matDescr,elemDescr,nblank);
			   case doubleRowMatrix:return amatPutFDmatrix2(file, matrix.name, matrix.nrow, matrix.ncol, matrix,matDescr,elemDescr,nblank);
               default           : goto error2;
            }

      } else {
         /* Write matrix in binary format */
            switch ( matrix.type ) {
               case voidMatrix   : goto error11;
               case integerMatrix:  return amatPutBImatrix (file, matrix.name, matrix.nrow, matrix.ncol, matrix.data.i);
               case realMatrix   :  return amatPutBRmatrix (file, matrix.name, matrix.nrow, matrix.ncol, matrix.data.r);
               case doubleMatrix :  return amatPutBDmatrix (file, matrix.name, matrix.nrow, matrix.ncol, matrix.data.d);
               case charMatrix   :  return amatPutBCmatrix (file, matrix.name, matrix.nrow, matrix.ncol, (const char   *const*)(matrix.data.c),padUsing); /* cast-away? */
			   case realRowMatrix:  return amatPutBRmatrix2(file, matrix.name, matrix.nrow, matrix.ncol, matrix);
			   case doubleRowMatrix:return amatPutBDmatrix2(file, matrix.name, matrix.nrow, matrix.ncol, matrix);
               default           : goto error2;
            }
      } 

   /* Error handling */ 
      error2:  
               if (!file)  {amatUnexpectedNull(); return RET_ERR;}
			   sprintfC(amatError, amatp_table[2], matrix.name, (int) matrix.nrow,
                       (int) matrix.ncol, file->name, (int) matrix.type);
               return 1;

      error11: 
               if (!file)  {amatUnexpectedNull(); return RET_ERR;}
			   sprintfC(amatError, amatp_table[11], matrix.name, (int) matrix.nrow,
                       (int) matrix.ncol, file->name);
               return 1;

} /* End of  amatPutMatrixD2 */

int amatPutMatrixD (AmatPutFile *file, Amatrix matrix, const char *matDescr, const char *const elemDescr[], int nblank) {
	return amatPutMatrixD2(file,matrix,matDescr,elemDescr,nblank,' ');
}




int amatPutMatrix (AmatPutFile *file, Amatrix matrix) {

   /* Write matrix to file */
      char *s1;
      const char *const*s2;
      s1 = NULL;
      s2 = NULL;
      return amatPutMatrixD (file, matrix, s1, s2, 0);
}

int amatPutMatrixPadding (AmatPutFile *file, Amatrix matrix,char padForText) {

   /* Write matrix to file */
      char *s1;
      const char *const*s2;
      s1 = NULL;
      s2 = NULL;
      return amatPutMatrixD2 (file, matrix, s1, s2, 0, padForText);
}



int amatPutRinit (AmatPutFile *file, const char *name, integer ncol, AmatType type) {
   
   /* Initialize writing of numeric matrix row-wise to file */

   /* Declarations */
      size_t len;

	  if (!file||!name)  {amatUnexpectedNull(); return RET_ERR;}
   /* Check arguments */
      if ( file->noRow == 0 ) goto error12;
      assumption(ncol > 0);
	  /* Should not store as RowMatrix*/
	  if (type == realRowMatrix) {
		  type = realMatrix;
	  } else if (type == doubleRowMatrix) {
		  type = doubleMatrix;
	  } 
      if ( type != realMatrix && type != doubleMatrix) goto error7;
      if ( file->bstruct != binTrans ) goto error9;

   /* Store matrix name in "file" (for error messages) */
      len = strlen(name) + 1;
      file->matName = (char *) malloc( len*sizeof(char) );
      if ( file->matName == NULL ) goto error6;
      strcpy( file->matName, name );

   /* Initialize common data */
      file->type    = type;
      file->ncol    = ncol;
      file->nrowAct = 0;

   /* Initialize header and write header to file */
      switch ( file->format ) {
         case amatASCII : return amatPutFRinit (file);
         case amatBinary: return amatPutBRinit (file);
         default        : goto error4;
      }

   /* Error handling */
      error4:  
               if (!file||!name)  {amatUnexpectedNull(); return RET_ERR;}
			   sprintfC(amatError, amatp_table[4], name, file->name, (int) file->format);
               return 1;

      error6:  
               if (!file||!name)  {amatUnexpectedNull(); return RET_ERR;}
			   sprintfC(amatError, amatp_table[6], file->name, name, strerror(errno));
               return 1;

      error7:  
               if (!file||!name)  {amatUnexpectedNull(); return RET_ERR;}
			   sprintfC(amatError, amatp_table[7], type, name, file->name);
               return 1;

      error9:  
               if (!file||!name)  {amatUnexpectedNull(); return RET_ERR;}
			   sprintfC(amatError, amatp_table[9], name, file->name);
               return 1;

      error12: 
               if (!file||!name)  {amatUnexpectedNull(); return RET_ERR;}
			   sprintfC(amatError, amatp_table[12], name);
               return 1;

} /* End of  amatPutRinit */



int amatPutRclose (AmatPutFile *file) {

   /* Close row-wise writing of matrix */
	  if (!file)  {return 0;}
      if ( !file->fp ) return 0;

      switch ( file->format ) {
         case amatASCII : if ( putc('\n', file->fp) == EOF ) goto error8;
                          free(file->matName);
                          file->matName = NULL;
                          return amatPutFRclose(file);
         case amatBinary: free(file->matName);
                          file->matName = NULL;
                          return amatPutBRclose(file);
         default        : goto error4;
      }

   /* Error handling */
      error4: if (file) sprintfC(amatError, amatp_table[4], file->matName, file->name, (int) file->format);
              return 1;

      error8: if (file) sprintfC(amatError, amatp_table[8], file->matName, file->name,
                      strerror(errno));
              return 1;

} /* End of  amatPutRclose */


#define ALLOW_READING_OF_DATA_DURING_SIMULATION 1
int amatPutRreal (AmatPutFile *file, const real row[]) {

   /* Write row of real matrix to file */
      size_t nrow;
	  
	  if (!file)  {amatUnexpectedNull(); return RET_ERR;}

      switch ( file->format ) {
	     case amatASCII : {
							int res=amatPutFRreal(file, row);
							if (res==0) {
#if ALLOW_READING_OF_DATA_DURING_SIMULATION
								amatPutFRclose(file);
								fflush(file->fp);
#endif
							}
							return res;
						  }
         case amatBinary: nrow = file->x.mrows;
                          if ( fwrite(row, sizeof(real), nrow, file->fp) != nrow ) goto error10;
                          file->x.ncols++;
#if ALLOW_READING_OF_DATA_DURING_SIMULATION
						  amatPutBRclose(file);
						  fflush(file->fp);
#endif
                          return 0;
         default        : goto error4;
      }

   /* Error handling */
      error4 : 
               if (!file)  {amatUnexpectedNull(); return RET_ERR;}
			   sprintfC(amatError, amatp_table[4], file->matName, file->name, (int) file->format);
               return 1;

      error10: 
               if (!file)  {amatUnexpectedNull(); return RET_ERR;}
			   sprintfC(amatError, amatp_table[10], file->x.ncols+1, file->matName, file->name,
                       strerror(errno));
               return 1;

} /* End of  matPutRreal */



int amatPutRdouble (AmatPutFile *file, const doublereal row[]) {

   /* Write row of doublereal matrix to file */
      size_t nrow;

	  if (!file)  {amatUnexpectedNull(); return RET_ERR;}

      switch ( file->format ) {
	     case amatASCII : {
							int res=amatPutFRdouble(file, row);
							if (res==0)  {
#if ALLOW_READING_OF_DATA_DURING_SIMULATION
								amatPutFRclose(file);
								fflush(file->fp);
#endif
							}
							return res;
						  }
         case amatBinary: nrow = file->x.mrows;
                          if ( fwrite(row, sizeof(doublereal), nrow, file->fp) != nrow ) goto error10;
                          file->x.ncols++;
#if ALLOW_READING_OF_DATA_DURING_SIMULATION
						  amatPutBRclose(file);
						  fflush(file->fp);
#endif
                          return 0;
         default        : goto error4;
      }

   /* Error handling */
      error4 : 
               if (!file)  {amatUnexpectedNull(); return RET_ERR;}
			   sprintfC(amatError, amatp_table[4], file->matName, file->name, (int) file->format);
               return 1;

      error10: 
               if (!file)  {amatUnexpectedNull(); return RET_ERR;}
			   sprintfC(amatError, amatp_table[10], file->x.ncols+1, file->matName, file->name,
                       strerror(errno));
               return 1;

} /* End of  matPutRdouble */



/* Procedures for formatted output ..........................................*/

static int writeStringAndPossiblyEndl(FILE*fp,const char*s,int endLine,int quoteIfNeeded) {
	/* Return <0 if error */
	int anyQuoteNeeded=0;
	int i;
	if (quoteIfNeeded) for(i=0;s[i];++i) {
		if (s[i]=='\\' || s[i]=='\n') anyQuoteNeeded=1;
	}
	if (anyQuoteNeeded) {
		for(i=0;s[i];++i) {
			switch(s[i]) {
			case '\\':fprintf(fp,"\\\\");break;
			case '\n':fprintf(fp,"\\n");break;
			default:fprintf(fp,"%c",s[i]);break;
			}
		}
		if (endLine) return fprintf(fp,"\n");
		return 1;
	}
	return fprintfC(fp,endLine?"%s\n":"%s",s);
}

static int amatPutFImatrix (AmatPutFile *file, const char *name, integer nrow, integer ncol, const integer data[],
                            const char *matDescr, const char *const elemDescr[], int nblank) {

   /* Write integer matrix in formatted form on file */

   /* Declarations */
      FILE    *fp;
      const integer *p, *prow;
      integer  i, j, nb;
      integer  imax;
      int      sign, digits, k;
      char     fmt[20];

	  if (!file||!name)  {amatUnexpectedNull(); return RET_ERR;}

   /* Assertions */
      assumption(nrow >= 0);
      assumption(ncol >= 0);

   /* Write matrix type, name and dimensions to file */
      fp = file->fp;
      if ( matDescr!=NULL && matDescr[0]!='\0' ) { if ( fprintfC(fp, "%s\n", matDescr) < 0 ) goto error3; }
      if ( fprintfC(fp,"int %s(%d,%d)\n", name, (int) nrow, (int) ncol) < 0 ) goto error3;

   /* Determine biggest element of matrix */
      imax = 0;
      sign = 0;
      for (i=0; i<nrow*ncol; i++) {
         if ( Dymola_abs(data[i]) > imax ) imax = Dymola_abs(data[i]);
         if ( data[i] < 0 ) sign = 1;
      }

   /* Determine number of digits of maximal element (+ optionally sign) */
      digits = 1 + sign;
      while( imax > 9) {
         imax = imax/10;
         digits++;
      }

   /* Build format */
      sprintfC(fmt," %%%dd", digits);

   /* Write matrix data to file */
      p    = data;
      prow = data;
      for (i=0; i<nrow; i++) {
         k = 1;
         for (j=0; j<ncol; j++) {
             if ( fprintfC(fp, fmt, (int) *p) < 0 ) goto error3;
             p += nrow;
             k++;
             /* insert a new-line after a certain amount of numbers */
                if ( (k > 7) && (j < ncol-1) ) {
                   k = 1;
                   if ( fprintfC(fp,"\n    ") < 0 ) goto error3;
                }
         }
         if ( elemDescr != NULL ) { 
            if ( nblank > 0 ) {
               nb = nblank;
               while( --nb ) {if ( putc(' ', fp) == EOF ) goto error3;}
               if ( fputs("# ", fp) == EOF ) goto error3; 
            }
            if ( writeStringAndPossiblyEndl(fp, elemDescr[i],0,nblank>0) <0) goto error3;
         }
         if ( fprintfC(fp,"\n") < 0 ) goto error3;
         prow++;
         p = prow;
      }

   if ( putc('\n',fp) == EOF ) goto error3;
   return 0;

   /* Error handling */
      error3: 
              if (!file||!name)  {amatUnexpectedNull(); return RET_ERR;}
			  sprintfC(amatError, amatp_table[3], name, (int) nrow, (int) ncol,
                      file->name, strerror(errno));
              return 1;

} /* End of  amatPutFImatrix */



static int amatPutFRmatrix (AmatPutFile *file, const char *name, integer nrow, integer ncol, const real data[],
                            const char *matDescr, const char *const elemDescr[], int nblank) {

   /* Write real matrix in formatted form on file */

   /* Declarations */
      FILE    *fp;
      const real    *p, *prow;
	  real pr;
      integer  i, j, nb;
      double   pg;
      int      k;

	  if (!file||!name)  {amatUnexpectedNull(); return RET_ERR;}
   /* Assertions */
      assumption(nrow >= 0);
      assumption(ncol >= 0);

   /* Write matrix type, name and dimensions to file */
      fp = file->fp;
      if ( matDescr!=NULL && matDescr[0]!='\0' ) { if ( fprintfC(fp, "%s\n", matDescr) < 0 ) goto error3; }
      if ( fprintfC(fp,"float %s(%d,%d)\n", name, (int) nrow, (int) ncol) < 0 ) goto error3;

   /* Write matrix data to file */
      p    = data;
      prow = data;
      for (i=0; i<nrow; i++) {
         k = 1;
         for (j=0; j<ncol; j++) {
             /* Check whether number is an integer number and in the range -100 000 ...100 000 */
			    double x=*p;
			    if (x>=FLT_MAX) {
					x=FLT_MAX;
				} else if (x<=-FLT_MAX) {
					x=-FLT_MAX;
				}
                pr = (real) modf((double) x, &pg);
	            if ( pr == 0.0 && pg >=-100000.0 && pg <= 100000.0) {
                   if ( fprintfC(fp," %7d       ", (int) x) < 0 ) goto error3;
                } else {
                   if ( fprintfC(fp," %14.7E", (double) x) < 0 ) goto error3;
                }
             p += nrow;
             k++;
             /* insert a new-line after a certain amount of numbers */
                if ( (k > 5) && (j < ncol-1) ) {
                   k = 1;
                   if ( fprintfC(fp,"\n    ") < 0 ) goto error3;
                }
         }
         if ( elemDescr != NULL ) { 
            if ( nblank > 0 ) {
               nb = nblank;
               while( --nb ) {if ( putc(' ', fp) == EOF ) goto error3;}
               if ( fputs("# ", fp) == EOF ) goto error3; 
            }
            if ( writeStringAndPossiblyEndl(fp, elemDescr[i],0,nblank>0) < 0 ) goto error3;
         }
         if ( fprintfC(fp,"\n") < 0 ) goto error3;
         prow++;
         p = prow;
      }

   if ( putc('\n',fp) == EOF ) goto error3;
   return 0;

   /* Error handling */
      error3: 
              if (!file||!name)  {amatUnexpectedNull(); return RET_ERR;}
			  sprintfC(amatError, amatp_table[3], name, (int) nrow, (int) ncol,
                      file->name, strerror(errno));
              return 1;

} /* End of  amatPutFRmatrix */

static int amatPutFRmatrix2 (AmatPutFile *file, const char *name, integer nrow, integer ncol, const Amatrix matrix,
                            const char *matDescr, const char *const elemDescr[], int nblank) {

   /* Write real matrix in formatted form on file */

   /* Declarations */
      FILE    *fp;
	  real pr;
      integer  i, j, nb;
	  int      isRowFormat = (matrix.type == realRowMatrix);
      double   pg;
      int      k;

	  if (!file||!name)  {amatUnexpectedNull(); return RET_ERR;}
   /* Assertions */
      assumption(nrow >= 0);
      assumption(ncol >= 0);

   /* Write matrix type, name and dimensions to file */
      fp = file->fp;
      if ( matDescr!=NULL && matDescr[0]!='\0' ) { if ( fprintfC(fp, "%s\n", matDescr) < 0 ) goto error3; }
      if ( fprintfC(fp,"float %s(%d,%d)\n", name, (int) nrow, (int) ncol) < 0 ) goto error3;

   /* Write matrix data to file */
      for (i=0; i<nrow; i++) {
         k = 1;
         for (j=0; j<ncol; j++) {
			 double x;
             /* Check whether number is an integer number and in the range -100 000 ...100 000 */
			 if (isRowFormat) {
				 x = matrix.data.rrow[i][j];
			 } else {
				 x = matrix.data.r[j*nrow+i];
			 }
			 if (x>=FLT_MAX) {
				 x=FLT_MAX;
			 } else if (x<=-FLT_MAX) {
				 x=-FLT_MAX;
			 }
			 pr = (real) modf(x, &pg);
			 if ( pr == 0.0 && pg >=-100000.0 && pg <= 100000.0) {
				 if ( fprintfC(fp," %7d       ", (int) x) < 0 ) goto error3;
			 } else {
				 if ( fprintfC(fp," %14.7E", x) < 0 ) goto error3;
			 }
			 k++;
			 /* insert a new-line after a certain amount of numbers */
			 if ( (k > 5) && (j < ncol-1) ) {
				 k = 1;
				 if ( fprintfC(fp,"\n    ") < 0 ) goto error3;
			 }
		 }
         if ( elemDescr != NULL ) { 
            if ( nblank > 0 ) {
               nb = nblank;
               while( --nb ) {if ( putc(' ', fp) == EOF ) goto error3;}
               if ( fputs("# ", fp) == EOF ) goto error3; 
            }
            if ( writeStringAndPossiblyEndl(fp, elemDescr[i],0,nblank>0) < 0 ) goto error3;
         }
         if ( fprintfC(fp,"\n") < 0 ) goto error3;
      }

   if ( putc('\n',fp) == EOF ) goto error3;
   return 0;

   /* Error handling */
      error3: 
              if (!file||!name)  {amatUnexpectedNull(); return RET_ERR;}
			  sprintfC(amatError, amatp_table[3], name, (int) nrow, (int) ncol,
                      file->name, strerror(errno));
              return 1;

} /* End of  amatPutFRmatrix2 */


static int amatPutFDmatrix (AmatPutFile *file, const char *name, integer nrow, integer ncol, const doublereal data[],
							 const char *matDescr, const char *const elemDescr[], int nblank) {

   /* Write doublereal matrix in formatted form on file */

   /* Declarations */
      FILE      *fp;
      const doublereal *p, *prow;
      integer     i, j, nb, ndigit, nnum;
      doublereal  pr;
      double      pg;
      int         k;
      int        *type=NULL;
      char        str[10];

	  if (!file||!name)  {amatUnexpectedNull(); return RET_ERR;}

   /* Assertions */
      assumption(nrow >= 0);
      assumption(ncol >= 0);

   /* Write matrix type, name and dimensions to file */
      fp = file->fp;
      if ( matDescr!=NULL && matDescr[0]!='\0' ) { if ( fprintfC(fp, "%s\n", matDescr) < 0 ) goto error3; }
      if ( fprintfC(fp,"double %s(%d,%d)\n", name, (int) nrow, (int) ncol) < 0 ) goto error3;

   /* Allocate storage for an integer vector with ncol elements and
      initialize vector with 0 (= column of double numbers).
      This vector is used to produce nicer output.
   */
      type = (int *) calloc( (size_t) ncol, sizeof(int) );
      if ( type == NULL ) goto error3b;

   /* Determine, whether a column of the matrix contains only integer
      numbers and if this is the case, store the number of digits
      in the type vector
   */
      for (j=0; j<ncol; j++) {
         p = data + j*nrow;
         for (i=0; i<nrow; i++) {
            pr = modf(*p, &pg); 
            if ( pr != 0.0 ) {
               /* number is a double number */
                  type[j] = 0;
                  break;
            } else {
               /* number is an integer number */
                  pg     = fabs(pg);
                  ndigit =  1;
                  nnum   = 10;
                  for (k=0; k<8; k++) {
                     if ( pg < nnum ) {
                        type[j] = type[j] > ndigit ? type[j] : ndigit;
                        break;
                     }
                     ndigit++;
                     nnum *= 10;
                  }
                  if ( ndigit >= 8 ) {
                     type[j] = 0;
                     break;
                  }
            }
            p++;
         }
      }
   
   /* Write matrix data to file */
      p    = data;
      prow = data;
      for (i=0; i<nrow; i++) {
         k = 0;
         for (j=0; j<ncol; j++) {
			 double x=*p;
			 if (x>=DBL_MAX) {
				x=DBL_MAX;
			} else if (x<=-DBL_MAX) {
				x= -DBL_MAX;
			}
             /* If column contains only integers, write an integer */
                if ( type[j] > 0 ) {
                   /* write maximum number of digits of column in string */
                      sprintfC(str," %%%dd", type[j]+1);
                   if ( fprintfC(fp, str, (int) x) < 0 ) goto error3;
                   k = k+type[j]+2;
                } else {
                   /* check, whether number is an integer number */
                      pr = modf(x, &pg);
                      if ( pr == 0.0 && pg >=-100000.0 && pg <= 100000.0) {
                         if ( fprintfC(fp," %7d                ", (int) x) < 0 ) goto error3;
                      } else { 
                         if ( fprintfC(fp," %23.16E", (double) x) < 0 ) goto error3;
                      }
                      k = k+24;
                }
             p += nrow;
             /* insert a new-line after a certain amount of numbers */
                if ( (k > 70) && (j < ncol-1) ) {
                   k = 0;
                   if ( fprintfC(fp,"\n") < 0 ) goto error3;
                }
         }
         if ( elemDescr != NULL ) { 
            if ( nblank > 0 ) {
               nb = nblank;
               while( --nb ) {if ( putc(' ', fp) == EOF ) goto error3;}
               if ( fputs("# ", fp) == EOF ) goto error3; 
            }
            if ( writeStringAndPossiblyEndl(fp, elemDescr[i],0,nblank>0) < 0 ) goto error3;
         }
         if ( fprintfC(fp,"\n") < 0 ) goto error3;
         prow++;
         p = prow;
      }

   if ( putc('\n',fp) == EOF ) goto error3;
   free(type);
   return 0;

   /* Error handling */
      error3 : 
               if (!file||!name)  {amatUnexpectedNull(); return RET_ERR;}
			   sprintfC(amatError, amatp_table[3], name, (int) nrow, (int) ncol,
                       file->name, strerror(errno));
               free(type);
               return 1;

      error3b: 
               if (!file||!name)  {amatUnexpectedNull(); return RET_ERR;}
			   sprintfC(amatError, amatp_table[3], name, (int) nrow, (int) ncol,
                       file->name, "not enough memory\n");
               return 1;

} /* End of  amatPutFDmatrix */

static int amatPutFDmatrix2 (AmatPutFile *file, const char *name, integer nrow, integer ncol, const Amatrix matrix, 
                            const char *matDescr, const char *const elemDescr[], int nblank) {

   /* Write doublereal matrix in formatted form on file */

   /* Declarations */
      FILE      *fp;
	  double pstar;
      integer     i, j, nb, ndigit, nnum;
      doublereal  pr;
	  int         isRowFormat = (matrix.type == doubleRowMatrix);
      double      pg;
      int         k;
      int        *type=NULL;
      char        str[10];

	  if (!file||!name)  {amatUnexpectedNull(); return RET_ERR;}

   /* Assertions */
      assumption(nrow >= 0);
      assumption(ncol >= 0);

	  if (matrix.type != doubleMatrix && matrix.type != doubleRowMatrix) {
		  goto error13;
	  }

   /* Write matrix type, name and dimensions to file */
      fp = file->fp;
      if ( matDescr!=NULL && matDescr[0]!='\0' ) { if ( fprintfC(fp, "%s\n", matDescr) < 0 ) goto error3; }
      if ( fprintfC(fp,"double %s(%d,%d)\n", name, (int) nrow, (int) ncol) < 0 ) goto error3;

   /* Allocate storage for an integer vector with ncol elements and
      initialize vector with 0 (= column of double numbers).
      This vector is used to produce nicer output.
   */
      type = (int *) calloc( (size_t) ncol, sizeof(int) );
      if ( type == NULL ) goto error3b;

   /* Determine, whether a column of the matrix contains only integer
      numbers and if this is the case, store the number of digits
      in the type vector
   */
      for (j=0; j<ncol; j++) {
         for (i=0; i<nrow; i++) {
			if (isRowFormat) 
				pstar = matrix.data.drow[i][j];
			else
				pstar = matrix.data.d[j*nrow+i];
            pr = modf(pstar, &pg); 
            if ( pr != 0.0 ) {
               /* number is a double number */
                  type[j] = 0;
                  break;
            } else {
               /* number is an integer number */
                  pg     = fabs(pg);
                  ndigit =  1;
                  nnum   = 10;
                  for (k=0; k<8; k++) {
                     if ( pg < nnum ) {
                        type[j] = type[j] > ndigit ? type[j] : ndigit;
                        break;
                     }
                     ndigit++;
                     nnum *= 10;
                  }
                  if ( ndigit >= 8 ) {
                     type[j] = 0;
                     break;
                  }
            }
         }
      }
   
   /* Write matrix data to file */
      for (i=0; i<nrow; i++) {
         k = 0;
         for (j=0; j<ncol; j++) {
			 if (isRowFormat) 
				pstar = matrix.data.drow[i][j];
			else
				pstar = matrix.data.d[j*nrow+i];
             /* If column contains only integers, write an integer */
                if ( type[j] > 0 ) {
                   /* write maximum number of digits of column in string */
                      sprintfC(str," %%%dd", type[j]+1);
                   if ( fprintfC(fp, str, (int) pstar) < 0 ) goto error3;
                   k = k+type[j]+2;
                } else {
                   /* check, whether number is an integer number */
                      pr = modf(pstar, &pg);
                      if ( pr == 0.0 && pg >=-100000.0 && pg <= 100000.0) {
                         if ( fprintfC(fp," %7d                ", (int) pstar) < 0 ) goto error3;
                      } else { 
                         if ( fprintfC(fp," %23.16E", (double) pstar) < 0 ) goto error3;
                      }
                      k = k+24;
                }
             /* insert a new-line after a certain amount of numbers */
                if ( (k > 70) && (j < ncol-1) ) {
                   k = 0;
                   if ( fprintfC(fp,"\n") < 0 ) goto error3;
                }
         }
         if ( elemDescr != NULL ) { 
            if ( nblank > 0 ) {
               nb = nblank;
               while( --nb ) {if ( putc(' ', fp) == EOF ) goto error3;}
               if ( fputs("# ", fp) == EOF ) goto error3; 
            }
            if ( writeStringAndPossiblyEndl(fp, elemDescr[i],0,nblank>0)) goto error3;
         }
         if ( fprintfC(fp,"\n") < 0 ) goto error3;
      }

   if ( putc('\n',fp) == EOF ) goto error3;
   free(type);
   return 0;

   /* Error handling */
      error3 : 
               if (!file||!name)  {amatUnexpectedNull(); return RET_ERR;}
			   sprintfC(amatError, amatp_table[3], name, (int) nrow, (int) ncol,
                       file->name, strerror(errno));
               free(type);
               return 1;

      error3b: 
               if (!file||!name)  {amatUnexpectedNull(); return RET_ERR;}
			   sprintfC(amatError, amatp_table[3], name, (int) nrow, (int) ncol,
                       file->name, "not enough memory\n");
               return 1;

      error13: 
               if (!file||!name)  {amatUnexpectedNull(); return RET_ERR;}
			   sprintfC(amatError, amatp_table[13], matrix.name);
			   return 1;

} /* End of  amatPutFDmatrix2 */



static int amatPutFCmatrix (AmatPutFile *file, const char *name, integer nrow, integer ncol, const char *const data[],
                            const char *matDescr){
                          
   /* Write char matrix in formatted form on file */

   /* Declarations */
      FILE     *fp;
      integer   i, maxdim, idum;
      const char     *s;

	  if (!file||!name)  {amatUnexpectedNull(); return RET_ERR;}

   /* Assertions */
      assumption(nrow >= 0);
      assumption(ncol >= 0);
      if ( nrow == 0 ) ncol = 0;

   /* Determine maximum length of the strings, if necessary */
      if ( ncol == 0 ) {
         maxdim = (integer) 0;
         for (i=0; i<nrow; i++) {
            if ( (idum=strlen(data[i])) > maxdim ) maxdim = idum;
         }
      } else {
         maxdim = ncol;
      }

   /* Write matrix name and dimensions to file */
      fp = file->fp;
      if ( matDescr!=NULL && matDescr[0]!='\0' ) { if ( fprintfC(fp, "%s\n", matDescr) < 0 ) goto error3; }
      if ( fprintfC(fp,"char %s(%d,%d)\n", name, (int) nrow, (int) maxdim) < 0 ) goto error3;
  
   /* Write vector of strings onto file */
      if ( data != NULL ) {
         for (i=0; i<nrow; i++) {
              s = data[i];
              if ( s == NULL ) {
                 if ( fprintfC(fp,"\n") < 0 )     goto error3;
              } else {
				 if ( writeStringAndPossiblyEndl(fp,s,1,1) < 0 ) goto error3;
              }
         }
      } 

   if ( putc('\n',fp) == EOF ) goto error3;
   return 0;

   /* Error handling */
      error3: 
              if (!file||!name)  {amatUnexpectedNull(); return RET_ERR;}
			  sprintfC(amatError, amatp_table[3], name, (int) nrow, (int) maxdim,
                      file->name, strerror(errno));
              return 1;

} /* End of  amatPutFCmatrix */

                           

static int amatPutFRinit (AmatPutFile *file) {

   /* Initialize row-wise writing of numeric matrix to ASCII file */

   /* Declarations */
      FILE *fp;
	  if (!file)  {amatUnexpectedNull(); return RET_ERR;}
   /* Write matrix type and matrix name to file */
      fp = file->fp;
      if ( file->type == doubleMatrix ) {
         if ( fprintfC(fp,"double ") < 0 ) goto error8;
      } else {
         if ( fprintfC(fp,"float ") < 0 ) goto error8;
      }
      if ( fprintfC(fp,"%s",file->matName) < 0 ) goto error8;

   /* Inquire position of stream */
      if ( fgetpos(fp, &(file->posHead)) != 0 ) goto error8;

   /* Write initial dimensions and enough blanks for actual ones
      (supply an additional blank before '(', because there is a bug
      on MSDOS/gcc with fgetpos/fsetpos). 
   */
      if ( fprintfC(fp," (0,%d)                             \n",
                      (int) file->ncol) < 0 ) goto error8;

   return 0;
   
   /* Error handling */
      error8: 
              if (!file)  {amatUnexpectedNull(); return RET_ERR;}
			  sprintfC(amatError, amatp_table[8], file->matName, file->name, strerror(errno));
              return 1;

} /* End of  amatPutFRinit */



static int amatPutFRclose (AmatPutFile *file) {

   /* Close row-wise writing of matrix to ASCII file */

   /* Declarations */
      FILE    *fp;
      fpos_t   PosEnd;

	  if (!file)  {amatUnexpectedNull(); return RET_ERR;}
   /* Inquire actual position of file */
      fp = file->fp;
      if ( fgetpos(fp, &PosEnd) != 0 ) goto error8;
    
   /* Position file to start of header */
      if ( fsetpos(fp, &file->posHead) != 0 ) goto error8;

   /* Write actual dimensions */
      if ( fprintfC(fp,"(%d,%d) ", (int) file->nrowAct,
                                    (int) file->ncol) < 0 ) goto error8;

   /* Position file to previous position */
      if ( fsetpos(fp, &PosEnd) != 0 ) goto error8;

   return 0;
   
   /* Error handling */
      error8: 
              if (!file)  {amatUnexpectedNull(); return RET_ERR;}
			  sprintfC(amatError, amatp_table[8], file->matName, file->name, strerror(errno));
              return 1;

} /* End of  amatPutFRclose */



static int amatPutFRreal (AmatPutFile *file, const real row[]) {

   /* Write row of real matrix to ASCII file */

   /* Declarations */
      FILE       *fp;
      integer     i;
      int         k;

	  if (!file)  {amatUnexpectedNull(); return RET_ERR;}
   /* Store vector */
      fp = file->fp;
      k  = 1;
      for (i=0; i<file->ncol; i++) {
         if ( fprintfC(fp," %14.7E", (double) row[i]) < 0 ) goto error10;
         k++;
         /* insert a new-line after a certain amount of numbers */
            if ( (k > 5) && (i < file->ncol-1) ) {
               k = 1;
               if ( fprintfC(fp,"\n    ") < 0 ) goto error10;
            }
      }
      if ( fprintfC(fp,"\n") < 0 ) goto error10;

      (file->nrowAct)++;

      return 0;
   
   /* Error handling */
      error10: 
               if (!file)  {amatUnexpectedNull(); return RET_ERR;}
			   sprintfC(amatError, amatp_table[10], file->nrowAct+1, file->matName, 
                       file->name, strerror(errno));
               return 1;
 
} /* End of  amatPutFRreal */



static int amatPutFRdouble (AmatPutFile *file, const doublereal row[]) {

   /* Write row of doublereal matrix to ASCII file */

   /* Declarations */
      FILE       *fp;
      integer     i;
      int         k;

	  if (!file)  {amatUnexpectedNull(); return RET_ERR;}
   /* Store vector */
      fp = file->fp;
      k  = 1;
      for (i=0; i<file->ncol; i++) {
         if ( fprintfC(fp," %23.16E", row[i]) < 0 ) goto error10;
         k++;
         /* insert a new-line after a certain amount of numbers */
            if ( (k > 3) && (i < file->ncol-1) ) {
               k = 1;
               if ( fprintfC(fp,"\n    ") < 0 ) goto error10;
            }
      }
      if ( fprintfC(fp,"\n") < 0 ) goto error10;

      (file->nrowAct)++;

      return 0;
   
   /* Error handling */
      error10: 
               if (!file)  {amatUnexpectedNull(); return RET_ERR;}
			   sprintfC(amatError, amatp_table[10], file->nrowAct+1, file->matName, 
                       file->name, strerror(errno));
               return 1;
 
} /* End of  amatPutFRdouble */

            

/* Procedures for Matlab binary output ......................................*/

static int amatPutBImatrix (AmatPutFile *file, const char *name, integer nrow, integer ncol, const integer mat[]) {

   /* Write integer matrix in Matlab binary format to file */

   /* Declarations */
      FILE    *fp;
      Fmatrix *x;
      size_t   len, one=1;
      integer  i, j;
      const integer *p1, *p2;

	  if (!file||!name)  {amatUnexpectedNull(); return RET_ERR;}
   /* Assertions */
      assumption(nrow >= 0);
      assumption(ncol >= 0);

   /* Define struct Fmatrix */
      x = &file->x;
      x->type   = MATTYPE + 20;
      x->imagf  = 0L;
      len       = strlen(name) + 1;
      x->namlen = (int) len;
      if ( file->bstruct == binNormal ) {
         x->mrows  = (int) nrow;
         x->ncols  = (int) ncol;
      } else {                         /* Write transposed matrix */
         x->mrows  = (int) ncol;
         x->ncols  = (int) nrow;
      }

   /* Write header to file */
      fp = file->fp;
      if ( fwrite(x, sizeof(Fmatrix), one, fp) != one ) goto error3;
      if ( fwrite(name, sizeof(char), len, fp) != len ) goto error3;

   /* Write data to file */
      if ( file->bstruct == binNormal ) {
         /* Write column-wise */
            len = nrow*ncol;
            if ( len > 0ul ) {
               if ( fwrite(mat, sizeof(integer), len, fp) != len ) goto error3;
            }
      } else {
         /* Write row-wise */
            p1 = mat;
            p2 = mat;
            for (i=0; i<nrow; i++) {
               for (j=0; j<ncol; j++) {
                  if ( fwrite(p1, sizeof(integer), one, fp) != one ) goto error3;
                  p1 += nrow;
               }
               p2++;
               p1 = p2;
            }

      }
      return 0;

   /* Error handling */
      error3: 
              if (!file||!name)  {amatUnexpectedNull(); return RET_ERR;}
			  sprintfC(amatError, amatp_table[3], name, (int) nrow, (int) ncol,
                      file->name, strerror(errno));
              return 1;

} /* End of  amatPutBImatrix */



static int amatPutBRmatrix (AmatPutFile *file, const char *name, integer nrow, integer ncol, const real mat[]) {

   /* Write real matrix in Matlab binary format to file */

   /* Declarations */
      FILE    *fp;
      Fmatrix *x;
      size_t   len, one=1;
      integer  i, j;
      const real    *p1, *p2;

	  if (!file||!name)  {amatUnexpectedNull(); return RET_ERR;}
   /* Assertions */
      assumption(nrow >= 0);
      assumption(ncol >= 0);

   /* Define struct Fmatrix */
      x = &file->x;
      x->type   = MATTYPE + 10;
      x->imagf  = 0L;
      len       = strlen(name) + 1;
      x->namlen = (int) len;
      if ( file->bstruct == binNormal ) {
         x->mrows  = (int) nrow;
         x->ncols  = (int) ncol;
      } else {                         /* Write transposed matrix */
         x->mrows  = (int) ncol;
         x->ncols  = (int) nrow;
      }

   /* Write header to file */
      fp = file->fp;
      if ( fwrite(x, sizeof(Fmatrix), one, fp) != one ) goto error3;
      if ( fwrite(name, sizeof(char), len, fp) != len ) goto error3;

   /* Write data to file */
      if ( file->bstruct == binNormal ) {
         /* Write column-wise */
            len = nrow*ncol;
            if ( len > 0ul ) {
               if ( fwrite(mat, sizeof(real), len, fp) != len ) goto error3;
            }
      } else {
         /* Write row-wise */
            p1 = mat;
            p2 = mat;
            for (i=0; i<nrow; i++) {
               for (j=0; j<ncol; j++) {
                  if ( fwrite(p1, sizeof(real), one, fp) != one ) goto error3;
                  p1 += nrow;
               }
               p2++;
               p1 = p2;
            }
      }
      return 0;

   /* Error handling */
      error3: 
              if (!file||!name)  {amatUnexpectedNull(); return RET_ERR;}
			  sprintfC(amatError, amatp_table[3], name, (int) nrow, (int) ncol,
                      file->name, strerror(errno));
              return 1;

} /* End of  amatPutBRmatrix */


static int amatPutBRmatrix2(AmatPutFile *file, const char *name, integer nrow, integer ncol, Amatrix matrix) {

   /* Write real matrix in Matlab binary format to file */

   /* Declarations */
      FILE       *fp;
      Fmatrix    *x;
      size_t      len, one=1;
      integer     i, j;
      int         isRowMatrix = (matrix.type == realRowMatrix);

	  if (!file||!name)  {amatUnexpectedNull(); return RET_ERR;}

   /* Assertions */
      assumption(nrow >= 0);
      assumption(ncol >= 0);

   /* Define struct Fmatrix */
      x = &file->x;
      x->type   = MATTYPE + 10;
      x->imagf  = 0L;
      len       = strlen(name) + 1;
      x->namlen = (int) len;
      if ( file->bstruct == binNormal ) {
         x->mrows  = (int) nrow;
         x->ncols  = (int) ncol;
      } else {                         /* Write transposed matrix */
         x->mrows  = (int) ncol;
         x->ncols  = (int) nrow;
      }

   /* Write header to file */
      fp = file->fp;
      if ( fwrite(x, sizeof(Fmatrix), one, fp) != one ) goto error3;
      if ( fwrite(name, sizeof(char), len, fp) != len ) goto error3;

   /* Write data to file */
      if ( file->bstruct == binNormal ) {
         /* Write column-wise */
            len = nrow*ncol;
            if ( len > 0ul ) { 
				if (isRowMatrix) {
					for(j=0; j<ncol; j++) {
						for(i=0; i<nrow; i++) {
							if (fwrite(&(matrix.data.rrow[i][j]), sizeof(real), one, fp)!=one) goto error3;
						}
					}
				} else {
				   if ( fwrite(matrix.data.r, sizeof(real), len, fp) != len ) goto error3;
				}
            }
      } else {
         /* Write row-wise */
            for (i=0; i<nrow; i++) {
				if (isRowMatrix) {
					if (fwrite(matrix.data.rrow[i],sizeof(real), ncol, fp ) != (size_t)(ncol)) goto error3;
				} else {
					real*p1 = matrix.data.r+i;
					for (j=0; j<ncol; j++) {
						if ( fwrite(p1, sizeof(real), one, fp) != one ) goto error3;
						p1 += nrow;
					}
				}
            }
      }
      return 0;

   /* Error handling */
      error3: 
              if (!file||!name)  {amatUnexpectedNull(); return RET_ERR;}
			  sprintfC(amatError, amatp_table[3], name, (int) nrow, (int) ncol,
                      file->name, strerror(errno));
              return 1;

} /* End of  amatPutBRmatrix2 */


static int amatPutBDmatrix (AmatPutFile *file, const char *name, integer nrow, integer ncol, const doublereal mat[]) {

   /* Write doublereal matrix in Matlab binary format to file */

   /* Declarations */
      FILE       *fp;
      Fmatrix    *x;
      size_t      len, one=1;
      integer     i, j;
      const doublereal *p1, *p2;

   /* Assertions */
      assumption(nrow >= 0);
      assumption(ncol >= 0);

	  if (!file||!name)  {amatUnexpectedNull(); return RET_ERR;}

   /* Define struct Fmatrix */
      x = &file->x;
      x->type   = MATTYPE;
      x->imagf  = 0L;
      len       = strlen(name) + 1;
      x->namlen = (int) len;
      if ( file->bstruct == binNormal ) {
         x->mrows  = (int) nrow;
         x->ncols  = (int) ncol;
      } else {                         /* Write transposed matrix */
         x->mrows  = (int) ncol;
         x->ncols  = (int) nrow;
      }

   /* Write header to file */
      fp = file->fp;
      if ( fwrite(x, sizeof(Fmatrix), one, fp) != one ) goto error3;
      if ( fwrite(name, sizeof(char), len, fp) != len ) goto error3;

   /* Write data to file */
      if ( file->bstruct == binNormal ) {
         /* Write column-wise */
            len = nrow*ncol;
            if ( len > 0ul ) { 
               if ( fwrite(mat, sizeof(doublereal), len, fp) != len ) goto error3;
            }
      } else {
         /* Write row-wise */
            p1 = mat;
            p2 = mat;
            for (i=0; i<nrow; i++) {
               for (j=0; j<ncol; j++) {
                  if ( fwrite(p1, sizeof(doublereal), one, fp) != one ) goto error3;
                  p1 += nrow;
               }
               p2++;
               p1 = p2;
            }
      }
      return 0;

   /* Error handling */
      error3: 
              if (!file||!name)  {amatUnexpectedNull(); return RET_ERR;}
			  sprintfC(amatError, amatp_table[3], name, (int) nrow, (int) ncol,
                      file->name, strerror(errno));
              return 1;

} /* End of  amatPutBDmatrix */

static int amatPutBDmatrix2(AmatPutFile *file, const char *name, integer nrow, integer ncol, Amatrix matrix) {

   /* Write doublereal matrix in Matlab binary format to file */

   /* Declarations */
      FILE       *fp;
      Fmatrix    *x;
      size_t      len, one=1;
      integer     i, j;
      int         isRowMatrix = (matrix.type == doubleRowMatrix);

   /* Assertions */
	  if (!file||!name)  {amatUnexpectedNull(); return RET_ERR;}

      assumption(nrow >= 0);
      assumption(ncol >= 0);

	  if (matrix.type != doubleMatrix && matrix.type != doubleRowMatrix) {
		  goto error13;
	  }

   /* Define struct Fmatrix */
      x = &file->x;
      x->type   = MATTYPE;
      x->imagf  = 0L;
      len       = strlen(name) + 1;
      x->namlen = (int) len;
      if ( file->bstruct == binNormal ) {
         x->mrows  = (int) nrow;
         x->ncols  = (int) ncol;
      } else {                         /* Write transposed matrix */
         x->mrows  = (int) ncol;
         x->ncols  = (int) nrow;
      }

   /* Write header to file */
      fp = file->fp;
      if ( fwrite(x, sizeof(Fmatrix), one, fp) != one ) goto error3;
      if ( fwrite(name, sizeof(char), len, fp) != len ) goto error3;

   /* Write data to file */
      if ( file->bstruct == binNormal ) {
         /* Write column-wise */
            len = nrow*ncol;
            if ( len > 0ul ) { 
				if (isRowMatrix) {
					for(j=0; j<ncol; j++) {
						for(i=0; i<nrow; i++) {
							if (fwrite(&(matrix.data.drow[i][j]), sizeof(doublereal), one, fp)!=one) goto error3;
						}
					}
				} else {
				   if ( fwrite(matrix.data.d, sizeof(doublereal), len, fp) != len ) goto error3;
				}
            }
      } else {
         /* Write row-wise */
            for (i=0; i<nrow; i++) {
				if (isRowMatrix) {
					if (fwrite(matrix.data.drow[i],sizeof(doublereal), ncol, fp ) != (size_t)(ncol)) goto error3;
				} else {
					doublereal*p1 = matrix.data.d+i;
					for (j=0; j<ncol; j++) {
						if ( fwrite(p1, sizeof(doublereal), one, fp) != one ) goto error3;
						p1 += nrow;
					}
				}
            }
      }
      return 0;

   /* Error handling */
      error3: 
              if (!file||!name)  {amatUnexpectedNull(); return RET_ERR;}
			  sprintfC(amatError, amatp_table[3], name, (int) nrow, (int) ncol,
                      file->name, strerror(errno));
              return 1;

      error13: 
              if (!file||!name)  {amatUnexpectedNull(); return RET_ERR;}
			  sprintfC(amatError, amatp_table[13], matrix.name);
			  return 1;

} /* End of  amatPutBDmatrix2 */


static int amatPutBCmatrix (AmatPutFile *file, const char *name, integer nrow, integer ncol, const char *const mat[], char padUsing) {

   /* Write char matrix in Matlab binary format to file */

   /* Declarations */
      FILE    *fp;
      Fmatrix *x;
      integer  maxdim, idum, i, j, jmax;
      size_t   len, one=1;
      const char    *p;
	  char c;

	  if (!file||!name)  {amatUnexpectedNull(); return RET_ERR;}

   /* Assertions */
      assumption(nrow >= 0);
      assumption(ncol >= 0);
      if ( nrow == 0 ) ncol = 0;

   /* Determine maximum length of the strings */
      maxdim = (integer) 0;
      for (i=0; i<nrow; i++) {
		  if ( (idum=mat[i] ? strlen(mat[i]) : 0) > maxdim ) maxdim = idum;
      } 

   /* Define struct Fmatrix */
      x = &file->x;
      x->type   = MATTYPE + 51;
      x->imagf  = 0L;
      len       = strlen(name) + 1;
      x->namlen = (int) len;
      if ( file->bstruct == binNormal ) {
         x->mrows  = (int) nrow;
         x->ncols  = (int) maxdim;
      } else {                         /* Write transposed matrix */
         x->mrows  = (int) maxdim;
         x->ncols  = (int) nrow;
      }

   /* Write header to file */
      fp = file->fp;
      if ( fwrite(x, sizeof(Fmatrix), one, fp) != one ) goto error3;
      if ( fwrite(name, sizeof(char), len, fp) != len ) goto error3;

   /* Write data to file */
      if ( file->bstruct == binNormal ) {
         /* Write column-wise */
            for (j=0; j<maxdim; j++) {
		  	  for (i=0; i<nrow;i++) {
			    p = mat[i];
				if (p==0) {
				  p="";
				}
                if ( j <= (integer) strlen( p ) ) {
                  c = p[j];
                  if ( putc(c,fp) == EOF ) goto error3;
                } else {
                  if ( putc(padUsing,fp) == EOF ) goto error3;
                }
              }
            }

      } else {
         /* Write row-wise */
            for (i=0;i<nrow;i++) {
               /* Write string */
				  p = mat[i];
				  if (p==0) p="";
                  if ( fputs(p,fp) == EOF ) goto error3;

               /* Write trailing blanks */
                  jmax = maxdim - strlen(p);
                  for (j=0;j<jmax;j++) {
                    if ( putc(padUsing,fp) == EOF ) goto error3;
                  }
            }
      }
      return 0;

   /* Error handling */
      error3: 
              if (!file||!name)  {amatUnexpectedNull(); return RET_ERR;}
			  sprintfC(amatError, amatp_table[3], name, (int) nrow, (int) ncol,
                      file->name, strerror(errno));
              return 1;
 
} /* End of  amatPutBCmatrix */ 



static int amatPutBRinit (AmatPutFile *file) {

   /* Initialize row-wise writing of numeric matrix to ASCII file */

   /* Declarations */
      FILE    *fp;
      Fmatrix *x;
      size_t   len, one=1;

	  if (!file)  {amatUnexpectedNull(); return RET_ERR;}

   /* Inquire position of stream */
      fp = file->fp;
      if ( fgetpos(fp, &(file->posHead)) != 0 ) goto error8;

   /* Define struct Fmatrix */
      x = &file->x;
      x->type   = ( file->type == doubleMatrix) ? MATTYPE : MATTYPE + 10;
      x->mrows  = (int) file->ncol;
      x->ncols  = 0L;
      x->imagf  = 0L;
      len       = strlen(file->matName) + 1;
      x->namlen = (int) len;

   /* Write header to file */
      if ( fwrite(x, sizeof(Fmatrix), one, fp) != one ) goto error8;
      if ( fwrite(file->matName, sizeof(char), len, fp) != len ) goto error8;

   return 0;
   
   /* Error handling */
      error8: 
              if (!file)  {amatUnexpectedNull(); return RET_ERR;}
			  sprintfC(amatError, amatp_table[8], file->matName, file->name, strerror(errno));
              return 1;

} /* End of  amatPutBRinit */



static int amatPutBRclose (AmatPutFile *file) {

   /* Close row-wise writing of matrix to Matlab binary file */

   /* Declarations */
      FILE    *fp;
      Fmatrix *x;
      fpos_t   PosEnd;
      size_t   one=1;

	  if (!file)  {amatUnexpectedNull(); return RET_ERR;}

   /* Inquire actual position of file */
      fp = file->fp;
      if ( fgetpos(fp, &PosEnd) != 0 ) goto error8;
    
   /* Position file to start of header */
      if ( fsetpos(fp, &file->posHead) != 0 ) goto error8;

   /* Write actual header */
      x = &file->x;
      if ( fwrite(x, sizeof(Fmatrix), one, fp) != one ) goto error8;

   /* Position file to previous position */
      if ( fsetpos(fp, &PosEnd) != 0 ) goto error8;

   return 0;
   
   /* Error handling */
      error8: 
              if (!file)  {amatUnexpectedNull(); return RET_ERR;}
			  sprintfC(amatError, amatp_table[8], file->matName, file->name, strerror(errno));
              return 1;

} /* End of  amatPutBRclose */

#endif


/* Error messages */
   static const char *const amatu_table[] = {
     /* 0 */    "Not possible to allocate storage for %d character pointers\n"
                "in order to generate a default text matrix.\n"
     /* 1 */   ,"Not possible to resize text matrix from %d rows to %d rows.\n"
     /* 2 */   ,"Error should never occur when generating default text matrix\n"
                "with procedure amatTextGenElem (number of written characters = %d,"
                "should be %d).\n"
     /* 3 */   ,"Not possible to allocate storage for %d characters, in\n"
                "order to generate a default name for text matrix element %d.\n"
     /* 4 */   ,"Not possible to allocate storage for a matrix.\n"
     /* 5 */   ,"Not possible to allocate storage for %d integer pointers.\n"
     /* 6 */   ,"Name \"%.400s\" in matrix \"%.400s\" is unknown.\n"
   };


DYMOLA_STATIC Amatrix *amatNew(void) {
   /* Allocate storage for a new matrix and initialize it */

   /* Declarations */
      Amatrix *matrix;

   /* Allocate storage */
      matrix = (Amatrix *) malloc( sizeof(Amatrix) );
      if ( matrix == NULL ) goto error4;

   /* Initialize storage */
      amatInit(matrix);
      return matrix;

   /* Error handling */
      error4: sprintfC(amatError, amatu_table[4]);
              return NULL;
}



DYMOLA_STATIC void amatInit(Amatrix *mat) {
   /* Initialize an Amatrix object */
	  if (!mat) return;
      mat->name   = NULL;
      mat->nrow   = 0;
      mat->ncol   = 0;
      mat->type   = voidMatrix;
	  mat->nrowallocated =0;
	  mat->nrowread = 0;
      mat->data.v = NULL;
}


DYMOLA_STATIC void amatInitValue(char *name, integer nrow, integer ncol,
                   AmatType type, AmatData data, Amatrix *mat) {
    /* Initialize an Amatrix object with the specified values. 
       is overwritten. */
	   if (!mat) return;
       mat->name   = name;
       mat->nrow   = nrow;
       mat->ncol   = ncol;
       mat->type   = type;
	   mat->nrowallocated =0;
	   mat->nrowread = 0;
       mat->data.v = data.v;
}


DYMOLA_STATIC void amatDel (Amatrix *mat) {
   
   /* Delete matrix */ 
      if ( mat == NULL ) return;
      free(mat->name);
      if ( mat->type == charMatrix ) {
         amatTextDel(mat->data.c, mat->nrow);
      } else {
		 if ((mat->type == realRowMatrix || mat->type == doubleRowMatrix) && mat->data.v) {
			 int i;
			 if (mat->type == realRowMatrix) {
				 for(i=0;i<mat->nrow;++i) {
					 free(mat->data.rrow[i]);
					 mat->data.rrow[i]=0;
				 }
			 } else {
				 for(i=0;i<mat->nrow;++i) {
					 free(mat->data.drow[i]);
					 mat->data.drow[i]=0;
				 }
			 }
		 }
         free(mat->data.v);
      }
      amatInit(mat);

} /* End of  amatDel */

static char **amatTextGenElem (char *text[], char *name, integer nbeg, integer nend);


DYMOLA_STATIC void amatTextDel (char *text[], integer nrow) {

   /* Delete text matrix */
      integer  i;

      if ( text == NULL ) return;
      for (i=0; i<nrow; i++) {
          free( text[i] );
		  text[i]=NULL;
      }
      free (text);
	  text=NULL;

} /* End of  amatTextDel */



DYMOLA_STATIC char **amatTextGen (char *name, integer nrow) {

   /* Generate default text */

   /* Declarations */
      char **text;

   /* Check nrow */
      assumption(nrow > 0);

   /* Allocate storage for pointer text */
      text = (char **) malloc( (size_t) nrow*sizeof(char *) );
      if ( text == NULL ) goto error0;

   /* Generate default text */
      return amatTextGenElem(text, name, (integer) 0, nrow);

   /* Error handling */
      error0: sprintfC(amatError, amatu_table[0], (int) nrow);
              return NULL;

} /* End of  amatTextGen */



DYMOLA_STATIC char **amatTextResize (char *text[], integer nrowIn, integer nrowOut, char *name) {

   /* Resize text */

   /* Declarations */
      integer   i, iend;
      char    **vec;

   /* Resize text */
      if ( nrowIn == nrowOut ) {
         /* do nothing */
            vec = text;

      } else if ( nrowIn >  nrowOut ) {
         /* remove last rows */
            for (i=nrowOut; i<nrowIn; i++) free( text[i] );
            vec = (char **) realloc( (void *)text, (size_t) nrowOut*sizeof(char *) );
            if ( vec == NULL ) {iend=nrowOut; goto error1;}

      } else {
         /* add additional rows */
            vec = (char **) realloc( (void *)text, (size_t) nrowOut*sizeof(char *) );
            if ( vec == NULL ) {iend=nrowIn; goto error1;}
            vec = amatTextGenElem(vec, name, nrowIn, nrowOut);
      }
      return vec;

   /* Error handling */
      error1: sprintfC(amatError, amatu_table[1], (int) nrowIn, (int) nrowOut);
              for (i=0; i<iend; i++) free( vec[i] );
              free( vec );
              return NULL;

} /* End of  amatTextResize */



static char **amatTextGenElem (char *text[], char *name, integer nbeg, integer nend) {

   /* Generate default text from [nbeg] ... [nend-1]. If an error occurs, 
      free complete text */

   /* Declarations */
      size_t    lenName;
      char     *str;
      integer   i, j, len;
      int       digits, digitLimit;

   /* Check arguments */
      assumption(nbeg >= 0) ;
      assumption(nbeg <= nend);
 
   /* Determine length of name */
      lenName = strlen(name);

   /* Determine starting number of digits of nbeg+1 */
      j = nbeg+1;
      digits     = 0;
      digitLimit = 1;
      while (j > 0) {
         j = j/10;
         digits++; 
         digitLimit *= 10;
      }

   /* Generate default names */
      for (i=nbeg; i<nend; i++) {
         /* determine number of digits of i+1 */
            if ( i+1 >= digitLimit ) {
               digits++;
               digitLimit *= 10;
            }

         /* allocate storage for element name */
            len = (integer) (lenName + digits + 1);
            str = (char *) malloc( (size_t) len );
            if ( str == NULL ) {        /* release allocated storage */
               for (j=0; j<i; j++) free( text[j] );
               free( text );
			   text = NULL;
               goto error3;
            }
            text[i] = str;

         /* store default element name */
            if ( (j=(integer) sprintfC(str, "%s%d", name, (int) i+1)) != (len-1) ) {
               for (j=0; j<=i; j++) free( text[j] );
               free( text );
			   text = NULL;
               goto error2;
            }
      }

   return text;

   /* Error handling */
      error2: sprintfC(amatError, amatu_table[2], (int) j, (int) len); 
              return NULL;

      error3: sprintfC(amatError, amatu_table[3], (int) len, (int) i);
              return NULL;

} /* End of  amatTextGenElem */


DYMOLA_STATIC int amatTextFind (Amatrix tmat, const char *str) {

   /* Find C-index of string "str" within text field of "tmat". */
      char *str2;
      int   i;
      int   info = -1;
	  int last =tmat.nrowread;

   /* Initialize */
      if ( str == NULL ) return -1;

   /* Search from last position upto end */
      if ( last < 0   ) last = 0;
      for (i=last; i<tmat.nrow; i++) {
         str2 = tmat.data.c[i];
         if ( (str2 != NULL)  &&  strcmp(str2, str) == 0) goto found;
      }

   /* Search from first position upto last position */
      if ( last >= tmat.nrow ) last = tmat.nrow;
      for (i=0; i<last; i++) {
         str2 = tmat.data.c[i];
         if ( (str2 != NULL)  &&  strcmp(str2, str) == 0) goto found;
      }

   /* No comparision found */
      return info;

   /* Comparision found */
      found: tmat.nrowread=i;
			 info = i;
             return info;
}

DYMOLA_STATIC integer *amatTextIndex(Amatrix tmat1, Amatrix tmat2, int message) {

   /* Build index from two text objects. */
      integer *index = NULL;
      integer  i, j, ii, jbeg;

   /* Allocate storage for index */
      index = (integer *) malloc( tmat2.nrow*sizeof(integer *) );
      if ( index == NULL ) goto error5;

   /* Build index */
      jbeg = 0;
      for (i=0; i<tmat2.nrow; i++) {
         /* search from last position upto end */
            for (j=jbeg; j<tmat1.nrow; j++) {
               if ( strcmp(tmat1.data.c[j], tmat2.data.c[i]) == 0) {
                  index[i] = j;
                  jbeg = j+1;
                  if ( jbeg >= tmat1.nrow ) jbeg = 0;
                  goto next;
               }
            }

         /* search from first position upto last position */
            for (j=0; j<jbeg-1; j++) {
               if ( strcmp(tmat1.data.c[j], tmat2.data.c[i]) == 0) {
                  index[i] = j;
                  jbeg = j+1;
                  if ( jbeg >= tmat1.nrow ) jbeg = 0;
                  goto next;
               }
            }

         /* no comparision found */
            if ( message ) {
               ii = i;
               goto error6;
            } else {
               index[i] = -1;
            }

         /* comparison successful */
            next: ;
      }
      return index;

   /* Error handling */
      error5: sprintfC(amatError, amatu_table[5], (int) tmat2.nrow*sizeof(integer *));
              return NULL;

      error6: free(index);
              sprintfC(amatError, amatu_table[6], tmat2.data.c[ii], tmat2.name);
              return NULL;

} /* End of  amatTextIndex */


/* usertab.c

   User-defined function to define Dymola interpolation tables.
   The data for the user-defined tables are provided via 
   include-file "usertab.h".

   USUALLY, THIS FILE (usertab.c) NEED NOT TO BE CHANGED.

   Author : Martin Otter, DLR.
   Version: 1.0, 1997-09-30: implemented.
*/
/*
 * Copyright (C) 1997-2001 Dynasim AB.
 * All rights reserved.
 *
 */

#include "usertab.h"
#include "sprwat.h"

static int userTabFindName(UsertabTableElement tableDef[], int ntable, char *name);

DYMOLA_STATIC int usertab(char *tableName, int nipo, int dim[], int *colWise,
            double **table) {

   /* Define a table by statically storing the table in function usertab.
      This function can be adapted by the user to his/her needs.
      A 1D-table is defined as a matrix where
        - the first column is the abszissa data
        - the other columns are the ordinate data to
          be interpolated with respect to the first column.
      A 2D-table is defined as a matrix where
        - the first column (without first element), i.e., table(2:,1),
          is the first abscissa (u(1)),
        - the first row (without first element), i.e., table(1,2:),
          is the second abscissa (u(2)),
        - the other elements, i.e., table(i,j) with i>=2,j>=2,
          are the corresponding ordinate values.

      -> tableName  : Name of table.
      -> nipo       : = 0: time-table required (time interpolation).
                      = 1: 1D-table required.
                      = 2: 2D-table required.
      <- dim        : Actual values of dimensions.
      <- colWise    : = 0: table stored row-wise    (row_1, row_2, ..., row_n).
                    : = 1: table stored column-wise (column_1, column_2, ...).
      <- table      : Pointer to value vector.
      <- RETURN: = 0: No error.
                 = 1: An error occured. An error message is printed
                      from "usertab" with function "DymosimError".
   */

      char mess[500];
      int  ID;

   /* Search names */
      ID = userTabFindName(tableDef, N_TABLEDEF, tableName);
      if ( ID < 0 ) return 1;

   /* Check interpolation type */
      if ( nipo != tableDef[ID].nipo ) {
         sprintf(mess,"%dD-interpolation required for %dD-table \"%s\"\n",
                      nipo, tableDef[ID].nipo, tableName);
         DymosimError(mess);
         return 1;
      }

   /* Return desired values */
      *table   = tableDef[ID].value;
      *colWise = 0;
      dim[0] = tableDef[ID].dim[0];         
      dim[1] = tableDef[ID].dim[1];         

      return 0;
}


static int userTabFindName(UsertabTableElement tableDef[], int ntable, char *name) {

   /* Find table element with the given name. Return the index
      with respect to the found element or -1, in case it is
      not found
   */
      int  i;
      char mess[500];

      for (i=0; i<ntable; i++) {
         if ( strcmp(tableDef[i].name, name) == 0 ) return i;
      }

   /* Error */
      sprintf(mess, "The table matrix \"%s\" was not found in the user\n"
                    "supplied function \"usertab\" on file \"usertab.c\".\n",
                    name);
      DymosimError(mess);
      return (-1);              
}

/*
 * Copyright (C) 1997-2001 Dynasim AB.
 * All rights reserved.
 *
 */
/* dymc.c

   Include all relevant C-Files, in order that the number of
   files is reduced (problems with batch-files under
   Windows'95 and Windows 3.11).

   Author : Martin Otter, DLR.
   Release: 1997-11-26: implemented and tested
            1998-10-05: modified for realtime
            1999-01-20: Uses DYMOLA_DSPACE
*/

#if defined(DYMOLA_DSPACE)
/* Realtime system */
#define NO_FILE
#include <dsdefs.h>
#else
#include <stdio.h>
#endif
#include "sprwat.h"

#include "matutil.c"
#include "delay.c"
#include "dymtable.c"
#include "dymf2c.c"
#include "matrixop.c"

/* Start MSL Table Wrappers */

DYMOLA_STATIC int ModelicaTables_CombiTimeTable_init(const char* tableName, const char* fileName, 
                                        double const *table, size_t nRow, size_t nColumn,
                                        double startTime, int smoothness,
                                        int extrapolation) {
  int tableID = (int) dymTableTimeIni2(0.0, startTime, smoothness-1, extrapolation-1, 
                                       tableName, fileName, table, nRow, nColumn, 0.0);
  return tableID;
}

DYMOLA_STATIC void ModelicaTables_CombiTimeTable_close(int tableID) {
  ;
};

DYMOLA_STATIC double ModelicaTables_CombiTimeTable_interpolate(int tableID, int icol, double u) {
  return dymTableTimeIpo2(tableID, icol, u);
}

DYMOLA_STATIC double ModelicaTables_CombiTimeTable_minimumTime(int tableID) {
  int tableID2 = (int) tableID;
  return dymTableTimeTmin(tableID2);
}

DYMOLA_STATIC double ModelicaTables_CombiTimeTable_maximumTime(int tableID) {
  return dymTableTimeTmax(tableID);
}

DYMOLA_STATIC int ModelicaTables_CombiTable1D_init(const char* tableName, const char* fileName, 
                                       double const *table, size_t nRow, size_t nColumn, 
                                       int smoothness) {
  int tableID = (int) dymTableInit(1.0, smoothness-1, tableName, fileName, table, nRow, nColumn, 0.0);
  return tableID;
}

DYMOLA_STATIC void ModelicaTables_CombiTable1D_close(int tableID) {
  ;
};

DYMOLA_STATIC double ModelicaTables_CombiTable1D_interpolate(int tableID, int icol, double u) {
  return dymTableIpo1((double) tableID, icol, u);
}

DYMOLA_STATIC int ModelicaTables_CombiTable2D_init(const char* tableName, const char* fileName,
                                       double const *table, size_t nRow, size_t nColumn, 
                                       int smoothness) {
  int tableID = (int) dymTableInit(2.0, smoothness-1, tableName, fileName, table, nRow, nColumn, 0.0);
  return tableID;
}

DYMOLA_STATIC void ModelicaTables_CombiTable2D_close(int tableID) {
  ;
};

DYMOLA_STATIC double ModelicaTables_CombiTable2D_interpolate(int tableID, double u1, double u2) {
  return dymTableIpo2((double) tableID, u1, u2);
}

/* End MSL Table Wrappers */


/* MSL 3.2.1 Tables */

#if (defined(RT) || defined(NRT)) && !defined(NO_FILE_SYSTEM) && !defined(LABCAR)
#define NO_FILE_SYSTEM 1
#endif

#if !defined(NO_FILE_SYSTEM)
#include "ModelicaMatIO.c"
#endif
#include "ModelicaStandardTables.c"

#if defined(NO_FILE_SYSTEM)
#undef NO_FILE_SYSTEM
#endif

/* End MSL 3.2.1 Tables */
/* dymf.f -- translated by f2c (version 19961209).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/
/*
 * Copyright (C) 1993-2008 Dynasim AB.
 * All rights reserved.
 *
 * Parts of the code is copyrighted Minpack.
 *
 * Minpack Copyright Notice (1999) University of Chicago.  All rights reserved
 *  
 * Redistribution and use in source and binary forms, with or
 * without modification, are permitted provided that the
 * following conditions are met:
 *  
 * 1. Redistributions of source code must retain the above
 * copyright notice, this list of conditions and the following
 * disclaimer.
 *  
 * 2. Redistributions in binary form must reproduce the above
 * copyright notice, this list of conditions and the following
 * disclaimer in the documentation and/or other materials
 * provided with the distribution.
 *  
 * 3. The end-user documentation included with the
 * redistribution, if any, must include the following
 * acknowledgment:
 *  
 *    "This product includes software developed by the
 *    University of Chicago, as Operator of Argonne National
 *    Laboratory.
 *  
 * Alternately, this acknowledgment may appear in the software
 * itself, if and wherever such third-party acknowledgments
 * normally appear.
 *  
 * 4. WARRANTY DISCLAIMER. THE SOFTWARE IS SUPPLIED "AS IS"
 * WITHOUT WARRANTY OF ANY KIND. THE COPYRIGHT HOLDER, THE
 * UNITED STATES, THE UNITED STATES DEPARTMENT OF ENERGY, AND
 * THEIR EMPLOYEES: (1) DISCLAIM ANY WARRANTIES, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO ANY IMPLIED WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, TITLE
 * OR NON-INFRINGEMENT, (2) DO NOT ASSUME ANY LEGAL LIABILITY
 * OR RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR
 * USEFULNESS OF THE SOFTWARE, (3) DO NOT REPRESENT THAT USE OF
 * THE SOFTWARE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS, (4)
 * DO NOT WARRANT THAT THE SOFTWARE WILL FUNCTION
 * UNINTERRUPTED, THAT IT IS ERROR-FREE OR THAT ANY ERRORS WILL
 * BE CORRECTED.
 *  
 * 5. LIMITATION OF LIABILITY. IN NO EVENT WILL THE COPYRIGHT
 * HOLDER, THE UNITED STATES, THE UNITED STATES DEPARTMENT OF
 * ENERGY, OR THEIR EMPLOYEES: BE LIABLE FOR ANY INDIRECT,
 * INCIDENTAL, CONSEQUENTIAL, SPECIAL OR PUNITIVE DAMAGES OF
 * ANY KIND OR NATURE, INCLUDING BUT NOT LIMITED TO LOSS OF
 * PROFITS OR LOSS OF DATA, FOR ANY REASON WHATSOEVER, WHETHER
 * SUCH LIABILITY IS ASSERTED ON THE BASIS OF CONTRACT, TORT
 * (INCLUDING NEGLIGENCE OR STRICT LIABILITY), OR OTHERWISE,
 * EVEN IF ANY OF SAID PARTIES HAS BEEN WARNED OF THE
 * POSSIBILITY OF SUCH LOSS OR DAMAGES.
 */

#include <math.h>
#include "f2c.h"
#if defined(DYMOLAB) || defined(OPENGL)
#define DEBUG_NL 0
void DymosimMessageDouble(const char*c,double d) {
}
#else
#define DEBUG_NL 0
#endif

#if defined(DYMOLA_DSPACE) || defined(NO_FILE)
#include <dsdefs.h>
#else
#include <stdio.h>
#endif

#ifdef DEBUG_NL
static double LocalTime;
#endif

#include "libdssetup.h"
#include "localeless.h"

LIBDS_API int dymlqr_(const integer *fact,doublereal *a,const integer*lda,const integer*n,doublereal *b,const doublereal *tol,
					  integer *itype,doublereal *work,integer *iwork,integer *info);
LIBDS_API int dymli1_(integer *sysnr,const integer *fact,doublereal *a,const integer*lda,const integer*n,doublereal * b,doublereal *time,
					  integer *event, 
	integer *printpriority,doublereal *dwork,integer *iwork,integer *ierr);
LIBDS_API int dymli2_(integer *sysnr,integer *fact,doublereal *a,const integer*lda,const integer*n,doublereal *b,
					  doublereal *time,integer *event, integer *printpriority,doublereal *dwork,
					  integer *iwork,integer *factor,integer *ierr);
LIBDS_API int dymli3_(integer *sysnr,integer *fact,doublereal *a,const integer*lda,const integer*n,doublereal *b,
					  doublereal *time,integer *event,integer *printpriority,doublereal *dwork,
					  integer *iwork,integer *factor,const char*const*varnames,integer *ierr);
LIBDS_API int dymli4_(integer *sysnr,integer *fact,doublereal *a,const integer*lda,const integer*n,doublereal *b,doublereal *time,
					  integer *event,integer *printpriority,doublereal *dwork,integer *iwork,integer *factor,const char*const*varnames,integer *ierr,int*fEvent);

#if !defined(RT)
#include <time.h>
LIBDS_API double myclock(void) {
	static clock_t oldValue=0;
	static double extra=0;
	clock_t t=clock();
	if (t<oldValue) {
		/* Assume that times are represented by integer type, and that they work */
		/* using mod 2^32 and not mod 2^31 */
		if (t==-1) return 0; /* No timer*/
		extra+=pow(256.0,sizeof(clock_t));
	}
	oldValue=t;
	return extra+t;
}
#endif

#include <float.h>
/* Constants used to implement dlamch_(). DB 1998-10-05. */

#define __dj_ENFORCE_FUNCTION_CALLS
/* The define guards missing files for djgpp */
#include <ctype.h>
#undef  __dj_ENFORCE_FUNCTION_CALLS
#define lsame_(ca,cb,nca,ncb) ((*ca==*cb)||(toupperC(*(const unsigned char*)ca)==toupperC(*(const unsigned char*)cb)))

/* Table of constant values */

static const integer c__1 = 1;
static const doublereal c_b246 = -1.;
static const integer c_n1 = -1;
static const doublereal c_b263 = 1.;
static const doublereal c_b362 = 0.;
static const doublereal c_b438 = .5;
static const integer c__0 = 0;
static const integer c__2 = 2;
static const logical c_false = FALSE_;
static const integer c__3 = 3;
static const integer c__4 = 4;


DYMOLA_STATIC void dymnl2_Fast(const integer m, const integer n, doublereal* a, const integer lda, const doublereal* v, const doublereal* w);
DYMOLA_STATIC void dymnl3_Fast(const integer m, const integer n, doublereal*s, const integer ls, const doublereal*u, doublereal*v, doublereal*w, logical*sing,
				 const integer*ipivots);
DYMOLA_STATIC void dymnl4_Fast(const integer n, const doublereal* r__, const integer lr, const doublereal*diag, const doublereal*qtb, 
				 const doublereal delta, doublereal*x, doublereal*wa1, doublereal*wa2,
				 const integer*ipivots,int inJacobian);
DYMOLA_STATIC void dymnl7_Fast(const integer m, const integer n, doublereal* a, const integer lda, 
				 const logical pivot, 
			 integer*ipvt, const integer lipvt, doublereal*sigma, doublereal*acnorm, 
	doublereal*wa,const char*const*varnames,doublereal*x,int printEvent);
DYMOLA_STATIC doublereal dnrm2_Fast(const integer np, const doublereal *x, const integer incxp);
DYMOLA_STATIC doublereal dnrm2_Fast1(const integer np, const doublereal *x);
DYMOLA_STATIC void dymnl6_Fast(const integer m, const integer n, doublereal* q, const integer ldq, doublereal*wa);
DYMOLA_STATIC void dymnl5_Fast(integer*irev, const integer n, doublereal*x, doublereal*fvec, doublereal*fjac, 
			const integer ldfjac, const integer ml, const integer mu, const doublereal*nominalx, 
	doublereal*wa1, doublereal*wa2, doublereal*dsave, integer*isave,logical useCentral,logical printDetails,
	const char*const*varnames,int check);
DYMOLA_STATIC void dymnl1_Fast(integer*infrev, const integer iopt, const integer n, doublereal*x, 
			 doublereal*fvec, doublereal*fjac, const integer ldfjac, const doublereal xtol, 
			 const doublereal xtoldesired, const integer maxfev, const integer ml, 
			 const integer mu, const doublereal*nominalx, doublereal*diag, doublereal *diag2,
			 const doublereal factor, integer*info, integer*nfev, integer*njev, 
			 doublereal *r__, const integer lr,doublereal* qtf, 
			 doublereal * wa1, doublereal*wa2, doublereal*wa3, doublereal*wa4, 
			 const integer jopt, doublereal*dsave, integer*isave,const char*const*varnames,
			 integer*ipivots,int printEvent,int inJacobian);
DYMOLA_STATIC void dymnlinf2_(integer n,const doublereal*sol,const char*const*varnames);
DYMOLA_STATIC void dymnlinf1_(const integer*sysnr, const integer*n, const doublereal*sol, const doublereal*res, 
				const integer*nfunc, const integer*njac, const integer*nmax,const char*const*varnames);
DYMOLA_STATIC void dymnlerr_(const doublereal*time, const integer*sysnr, const integer*info, const integer*n, 
			   const integer*iopt, const doublereal*tol, const doublereal*sol, const doublereal*res, 
	const integer*nfunc, const integer*njac, const integer*nmax,const char*const*varnames);
/* --- */
LIBDS_API int dgemv_(const char*trans, const integer*m, const integer*n, const double*alpha, double*a, const integer*lda, 
		double*x, const integer*incx, const double*beta, double*y, const integer*incy, ftnlen trans_len);
DYMOLA_STATIC int dgetf2Dummy_(const integer*m, const integer*n, doublereal*a, const integer*lda,integer*ipiv,integer*hysteresis,integer*info);
DYMOLA_STATIC int xerbla_(char*srname,integer*info,ftnlen srname_len);
DYMOLA_STATIC int dswap_(const integer*n,doublereal*dx,const integer*incx,doublereal*dy,const integer*incy);
DYMOLA_STATIC int dscal_(const integer*n,const doublereal*da,doublereal*dx,const integer*incx);
DYMOLA_STATIC int dger_(const integer*m,const integer*n,const doublereal*alpha,doublereal*x,const integer*incx,doublereal*y,const integer*incy, doublereal*a,const integer*lda);
DYMOLA_STATIC int dgetrfDummy_(const integer*m,const integer*n,doublereal *a,const integer*lda,integer * ipiv,integer * hysteresis,integer * info);
DYMOLA_STATIC integer ilaenv_(const integer *ispec,const char *name__,const char *opts,const integer*n1,const integer*n2,const integer*n3,const integer*n4,ftnlen name_len,ftnlen opts_len);
DYMOLA_STATIC int dlaswp_(const integer*n,doublereal*a,const integer*lda,const integer *k1,const integer * k2,const integer * ipiv,const integer*incx);
DYMOLA_STATIC int dtrsm_(const char*side,const char*uplo,const char*transa,const char*diag,const integer*m,const integer*n,const doublereal*alpha,doublereal*a,const integer*lda, 
		   doublereal*b,const integer*ldb,ftnlen side_len,ftnlen uplo_len,ftnlen transa_len,ftnlen diag_len);
DYMOLA_STATIC int dgemm_(const char*transa,const char*transb,const integer*m,const integer*n,const integer * k,const doublereal*alpha,doublereal * a,const integer*lda,doublereal * b,const integer*ldb, 
	const doublereal*beta,doublereal * c__, const integer*ldc, ftnlen transa_len, ftnlen transb_len);

DYMOLA_STATIC doublereal dasum_(const integer*n,doublereal * dx,const integer*incx);
LIBDS_API int daxpy_(const integer*n,doublereal * da,doublereal *dx, const integer*incx,doublereal * dy,const integer*incy);
static void daxpy1_(const integer*np,doublereal * dap,doublereal * dx,doublereal * dy);

DYMOLA_STATIC int dgeequ_(const integer*m, const integer*n, doublereal *a, const integer *
	lda, doublereal *r, doublereal *c, doublereal *rowcnd, doublereal *
	colcnd, doublereal *amax, integer *info);
DYMOLA_STATIC int dlaqge_(const integer*m, const integer*n, doublereal *a, const integer *
	lda, doublereal *r, doublereal *c, doublereal *rowcnd, doublereal *
	colcnd, doublereal *amax, char *equed);
DYMOLA_STATIC int dcopy_(const integer*n,const doublereal *dx,const integer*incx,doublereal *dy,const integer*incy);
static doublereal ddot1_(const integer*np, doublereal *dx,doublereal * dy);
DYMOLA_STATIC doublereal ddot_(const integer*n,doublereal *dx,const integer*incx,doublereal *dy,const integer*incy);
DYMOLA_STATIC int dtrsv_(const char*uplo,const char*trans,const char * diag,const integer*n,doublereal * a,const integer*lda,doublereal *x,const integer*incx, ftnlen uplo_len, 
	ftnlen trans_len, ftnlen diag_len);
DYMOLA_STATIC integer idamax_(const integer*n,doublereal*dx,const integer*incx);
DYMOLA_STATIC int dgecon_(char *norm,const integer*n,doublereal *a,const integer*lda,doublereal *anorm,doublereal *rcond,doublereal *work,integer* iwork,integer* info,
	 ftnlen norm_len);
struct dlacon_clean_struct {
	integer i__,j,iter;
	doublereal temp;
	integer jump;
	integer jlast;
	doublereal altsgn,estold;
};
DYMOLA_STATIC int dlacon_(const integer*n,doublereal * v,doublereal * x,integer*isgn,doublereal *est,integer*kase);
DYMOLA_STATIC int dlacon_clean(const integer*n,doublereal * v,doublereal * x,integer*isgn,doublereal *est,integer*kase,struct dlacon_clean_struct *staticarea);
DYMOLA_STATIC int dlatrs_(const char*uplo,const char*trans,const char *diag,const char *normin,const integer*n,doublereal * a,const integer*lda,doublereal *x,doublereal *scale, 
	doublereal *cnorm,integer*info,ftnlen uplo_len,ftnlen trans_len,ftnlen diag_len,ftnlen normin_len);
DYMOLA_STATIC int drscl_(const integer*n,doublereal*sa,doublereal*sx,const integer*incx);
DYMOLA_STATIC int dgeqpf_(const integer*m,const integer*n,doublereal *a,const integer*lda,integer *jpvt,doublereal *tau,doublereal *work,integer *info);
DYMOLA_STATIC int dgeqr2_(const integer*m,const integer*n,doublereal *a,const integer*lda,doublereal *tau,doublereal *work,integer *info);
DYMOLA_STATIC int dorm2r_(const char*side,const char*trans,const integer*m,const integer*n,const integer *k,doublereal *a,const integer*lda,doublereal *tau, 
			doublereal *c__,const integer*ldc, doublereal *work,integer *info,ftnlen side_len,ftnlen trans_len);
DYMOLA_STATIC int dlarfg_(const integer*n, doublereal*alpha,doublereal * x,const integer*incx,doublereal *tau);
DYMOLA_STATIC int dlarf_(const char*side,const integer*m,const integer*n,doublereal *v,const integer*incv,doublereal *tau,doublereal *c__,const integer*ldc,doublereal *work, 
	ftnlen side_len);
DYMOLA_STATIC int dgetrf_(const integer*m,const integer*n,doublereal *a,const integer*lda,integer *ipiv,integer *info);
DYMOLA_STATIC int dgetf2_(const integer*m,const integer*n,doublereal *a,const integer*lda,integer *ipiv,integer *info);
DYMOLA_STATIC int dgetrs_(const char*trans,const integer*n,const integer*nrhs,doublereal *a,const integer*lda,integer *ipiv,doublereal *b,const integer*ldb,integer *info, 
	ftnlen trans_len);
DYMOLA_STATIC int dlabad_(doublereal *small,doublereal *large);
DYMOLA_STATIC int dlaic1_(const integer *job,integer *j,doublereal *x,doublereal *sest,doublereal *w,doublereal *gamma,doublereal *sestpr,doublereal *s,doublereal *c__);
DYMOLA_STATIC doublereal dlange_(char *norm,const integer*m,const integer*n,doublereal *a,const integer*lda,doublereal *work,ftnlen norm_len);
DYMOLA_STATIC int dlassq_(const integer*n,doublereal *x,const integer*incx,doublereal *scale,doublereal *sumsq);
DYMOLA_STATIC doublereal dlapy2_(doublereal *x,doublereal * y);
DYMOLA_STATIC int dlatzm_(const char*side,const integer*m,const integer*n,doublereal *v,const integer*incv,doublereal *tau,doublereal * c1,doublereal *c2,const integer*ldc,doublereal *work, 
	ftnlen side_len);
DYMOLA_STATIC int dtzrqf_(const integer*m,const integer*n,doublereal *a,const integer*lda,doublereal *tau,integer *info);
DYMOLA_STATIC logical (lsame_)(const char*ca,const char*cb,ftnlen ca_len,ftnlen cb_len);
DYMOLA_STATIC int dgbfa_(doublereal *abd,const integer*lda,const integer*n,const integer*ml,const integer*mu,integer *ipvt,integer *info);
DYMOLA_STATIC int dgbsl_(doublereal *abd,const integer*lda,const integer*n,const integer*ml,const integer*mu,integer *ipvt,doublereal *b,const integer *job);
DYMOLA_STATIC int dgeco_(doublereal *a,const integer*lda,const integer*n,integer *ipvt,doublereal * rcond,doublereal *z__);
DYMOLA_STATIC int dgefa_(doublereal *a,const integer*lda,const integer*n,integer *ipvt,integer *info);
DYMOLA_STATIC int dgesl_(doublereal *a,const integer*lda,const integer*n,const integer *ipvt,doublereal *b,const integer *job);
DYMOLA_STATIC int dymres_(doublereal *a,const integer*lda,const integer*n,doublereal *b,integer *ierr);
DYMOLA_STATIC int dymsol_(doublereal *a,const integer*lda,const integer*n,doublereal *b,doublereal *dwork,integer *iwork,integer *ierr);
DYMOLA_STATIC int dymlin_(integer *infrev,const integer*n,doublereal *sol,doublereal *res,doublereal *jac,
			integer *ljac,doublereal *dwork,integer *iwork,doublereal *time,
			integer *event,integer *printpriority,integer *sysnr,integer *ierr);
DYMOLA_STATIC int dymnl_(integer *infrev,integer *iopt,const integer*n,doublereal *x,doublereal *fvec,doublereal *fjac,const doublereal *tol,
		   const doublereal *toldesired,integer *info,doublereal *dum,integer*ldum,integer *idum,const doublereal *epsfcn);
DYMOLA_STATIC int dymnon_(integer *infrev,integer *qiopt,integer *qnnl,doublereal *qsol,doublereal *qres,doublereal *qjac,
			doublereal *qtol, 
	integer *qinfo,doublereal*qd,integer *qi,integer *printevent,integer *qnlnr,doublereal *time,integer *qnlfunc,integer *qnljac,
	integer *qnlmax,integer *ierr___);
DYMOLA_STATIC int dymnon4_(integer *qinfrev,integer *qiopt,integer *qnnl,doublereal *qsol,doublereal *qres,doublereal *qjac,doublereal *qtol,
			 doublereal *qtoldesired,integer *qinfo,doublereal *qd,integer *qi,integer *idum,
			 integer *printevent,integer *qnlnr,doublereal *time,integer *qnlfunc,integer *qnljac,integer *qnlmax,const char*const*varnames,integer*ierr___);
DYMOLA_STATIC int handleevent_(const char *relation,integer*m,doublereal *cp,doublereal *cn,doublereal *c__,integer *init,
				 integer *printevent,doublereal *eps,doublereal *t,logical *anyevent,ftnlen relation_len);
LIBDS_API int booleanchanged_(const char *varname,logical *var,logical *newvar,integer *printevent, 
	logical *anyevent,ftnlen varname_len);

#ifndef GODESS
DYMOLA_STATIC doublereal dasum_(n, dx, incx)
const integer*n;
doublereal *dx;
const integer*incx;
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal ret_val, d__1, d__2, d__3, d__4, d__5, d__6;

    /* Local variables */
    integer i__, m, mp1;
    doublereal dtemp;
    integer nincx;


/*     takes the sum of the absolute values. */
/*     jack dongarra, linpack, 3/11/78. */
/*     modified 3/93 to return if incx .le. 0. */


    /* Parameter adjustments */
    --dx;

    /* Function Body */
    ret_val = 0.;
    dtemp = 0.;
    if (*n <= 0 || *incx <= 0) {
	return ret_val;
    }
    if (*incx == 1) {
	goto L20;
    }

/*        code for increment not equal to 1 */

    nincx = *n * *incx;
    i__1 = nincx;
    i__2 = *incx;
    for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
	dtemp += (d__1 = dx[i__], Dymola_abs(d__1));
/* L10: */
    }
    ret_val = dtemp;
    return ret_val;

/*        code for increment equal to 1 */


/*        clean-up loop */

L20:
    m = *n % 6;
    if (m == 0) {
	goto L40;
    }
    i__2 = m;
    for (i__ = 1; i__ <= i__2; ++i__) {
	dtemp += (d__1 = dx[i__], Dymola_abs(d__1));
/* L30: */
    }
    if (*n < 6) {
	goto L60;
    }
L40:
    mp1 = m + 1;
    i__2 = *n;
    for (i__ = mp1; i__ <= i__2; i__ += 6) {
	dtemp = dtemp + (d__1 = dx[i__], Dymola_abs(d__1)) + (d__2 = dx[i__ + 1], 
		Dymola_abs(d__2)) + (d__3 = dx[i__ + 2], Dymola_abs(d__3)) + (d__4 = dx[i__ 
		+ 3], Dymola_abs(d__4)) + (d__5 = dx[i__ + 4], Dymola_abs(d__5)) + (d__6 = 
		dx[i__ + 5], Dymola_abs(d__6));
/* L50: */
    }
L60:
    ret_val = dtemp;
    return ret_val;
} /* dasum_ */

/* Subroutine */ 
LIBDS_API int daxpy_(n, da, dx, incx, dy, incy)
const integer*n;
doublereal *da, *dx;
const integer*incx;
doublereal *dy;
const integer*incy;
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    integer i__, m, ix, iy, mp1;


/*     constant times a vector plus a vector. */
/*     uses unrolled loops for increments equal to one. */
/*     jack dongarra, linpack, 3/11/78. */


    /* Parameter adjustments */
    --dy;
    --dx;

    /* Function Body */
    if (*n <= 0) {
	return 0;
    }
    if (*da == 0.) {
	return 0;
    }
    if (*incx == 1 && *incy == 1) {
	goto L20;
    }

/*        code for unequal increments or equal increments */
/*          not equal to 1 */

    ix = 1;
    iy = 1;
    if (*incx < 0) {
	ix = (-(*n) + 1) * *incx + 1;
    }
    if (*incy < 0) {
	iy = (-(*n) + 1) * *incy + 1;
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	dy[iy] += *da * dx[ix];
	ix += *incx;
	iy += *incy;
/* L10: */
    }
    return 0;

/*        code for both increments equal to 1 */


/*        clean-up loop */

L20:
    m = *n % 4;
    if (m == 0) {
	goto L40;
    }
    i__1 = m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	dy[i__] += *da * dx[i__];
/* L30: */
    }
    if (*n < 4) {
	return 0;
    }
L40:
    mp1 = m + 1;
    i__1 = *n;
    for (i__ = mp1; i__ <= i__1; i__ += 4) {
	dy[i__] += *da * dx[i__];
	dy[i__ + 1] += *da * dx[i__ + 1];
	dy[i__ + 2] += *da * dx[i__ + 2];
	dy[i__ + 3] += *da * dx[i__ + 3];
/* L50: */
    }
    return 0;
} /* daxpy_ */
#endif

static
#if defined(_MSC_VER) && _MSC_VER>=1200
__inline
#elif __GNUC__
__inline
#endif
void daxpy1_(np, dap, dx, dy)
const integer*np;
doublereal *dap, *dx, *dy;
{
    const integer  n = *np;
    const doublereal da = *dap;
    integer i;
	
	/* constant times a vector plus a vector. */
	/* special case for increment=1, no loop unrolling. */
	
    if (da == 0.) return;
	
    for (i = 0; i < n; ++i) {
		dy[i] += da * dx[i];
    }
} /* daxpy1_ */

#ifndef GODESS
/* Subroutine */ DYMOLA_STATIC int dcopy_(n, dx, incx, dy, incy)
const integer*n;
const doublereal *dx;
const integer*incx;
doublereal *dy;
const integer*incy;
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    integer i__, m, ix, iy, mp1;


/*     copies a vector, x, to a vector, y. */
/*     uses unrolled loops for increments equal to one. */
/*     jack dongarra, linpack, 3/11/78. */


    /* Parameter adjustments */
    --dy;
    --dx;

    /* Function Body */
    if (*n <= 0) {
	return 0;
    }
    if (*incx == 1 && *incy == 1) {
	goto L20;
    }

/*        code for unequal increments or equal increments */
/*          not equal to 1 */

    ix = 1;
    iy = 1;
    if (*incx < 0) {
	ix = (-(*n) + 1) * *incx + 1;
    }
    if (*incy < 0) {
	iy = (-(*n) + 1) * *incy + 1;
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	dy[iy] = dx[ix];
	ix += *incx;
	iy += *incy;
/* L10: */
    }
    return 0;

/*        code for both increments equal to 1 */


/*        clean-up loop */

L20:
    m = *n % 7;
    if (m == 0) {
	goto L40;
    }
    i__1 = m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	dy[i__] = dx[i__];
/* L30: */
    }
    if (*n < 7) {
	return 0;
    }
L40:
    mp1 = m + 1;
    i__1 = *n;
    for (i__ = mp1; i__ <= i__1; i__ += 7) {
	dy[i__] = dx[i__];
	dy[i__ + 1] = dx[i__ + 1];
	dy[i__ + 2] = dx[i__ + 2];
	dy[i__ + 3] = dx[i__ + 3];
	dy[i__ + 4] = dx[i__ + 4];
	dy[i__ + 5] = dx[i__ + 5];
	dy[i__ + 6] = dx[i__ + 6];
/* L50: */
    }
    return 0;
} /* dcopy_ */
#endif /* GODESS */

static
#if defined(_MSC_VER) && _MSC_VER>=1200
__inline
#elif __GNUC__
__inline
#endif
doublereal ddot1_(np, dx, dy)
const integer*np;
doublereal *dx;
doublereal *dy;
{
    const integer n = *np;
    doublereal dtemp = 0;
    integer i;
	
	/*     forms the dot product of two vectors. */
	/*     special case for increment=1, no loop unrolling. */
	
    for (i = 0; i < n; ++i) {
		dtemp += dx[i] * dy[i];
    }
    return dtemp;
} /* ddot1_ */

#ifndef GODESS
DYMOLA_STATIC doublereal ddot_(n, dx, incx, dy, incy)
const integer*n;
doublereal *dx;
const integer*incx;
doublereal *dy;
const integer*incy;
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    integer i__, ix, iy;
    doublereal dtemp;


/*     forms the dot product of two vectors. */
/*     uses unrolled loops for increments equal to one. */
/*     jack dongarra, linpack, 3/11/78. */


    /* Parameter adjustments */
    --dy;
    --dx;

    /* Function Body */
    dtemp = 0.;
    if (*incx == 1 && *incy == 1) {
		i__1=*n;
		for(i__=1;i__+4<=i__1;i__+=5) {
			dtemp += dx[i__] * dy[i__] + dx[i__ + 1] * dy[i__ + 1] + dx[
				i__ + 2] * dy[i__ + 2] + dx[i__ + 3] * dy[i__ + 3] + dx[i__ + 
				4] * dy[i__ + 4];
		}
		for(;i__<=i__1;i__++)
			dtemp+=dx[i__]*dy[i__];
    } else {
		
		/*        code for unequal increments or equal increments */
		/*          not equal to 1 */
		
		ix = 1;
		iy = 1;
		if (*incx < 0) {
			ix = (-(*n) + 1) * *incx + 1;
		}
		if (*incy < 0) {
			iy = (-(*n) + 1) * *incy + 1;
		}
		i__1 = *n;
		for (i__ = 1; i__ <= i__1; ++i__) {
			dtemp += dx[ix] * dy[iy];
			ix += *incx;
			iy += *incy;
			/* L10: */
		}
	}
    return dtemp;
} /* ddot_ */

/* Subroutine */ DYMOLA_STATIC int dgemm_(transa, transb, m, n, k, alpha, a, lda, b, ldb, 
	beta, c__, ldc, transa_len, transb_len)
const char*transa, *transb;
const integer*m, *n, *k;
const doublereal*alpha;
doublereal*a;
const integer*lda;
doublereal *b;
const integer*ldb;
const doublereal*beta;
doublereal*c__;
const integer*ldc;
ftnlen transa_len;
ftnlen transb_len;
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, i__1, i__2, 
	    i__3;

    /* Local variables */
    integer i__, j, l, info;
    logical nota, notb;
    doublereal temp;
    integer ncola;
    integer nrowa, nrowb;

/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  DGEMM  performs one of the matrix-matrix operations */

/*     C := alpha*op( A )*op( B ) + beta*C, */

/*  where  op( X ) is one of */

/*     op( X ) = X   or   op( X ) = X', */

/*  alpha and beta are scalars, and A, B and C are matrices, with op( A ) */
/*  an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix. */

/*  Parameters */
/*  ========== */

/*  TRANSA - CHARACTER*1. */
/*           On entry, TRANSA specifies the form of op( A ) to be used in */
/*           the matrix multiplication as follows: */

/*              TRANSA = 'N' or 'n',  op( A ) = A. */

/*              TRANSA = 'T' or 't',  op( A ) = A'. */

/*              TRANSA = 'C' or 'c',  op( A ) = A'. */

/*           Unchanged on exit. */

/*  TRANSB - CHARACTER*1. */
/*           On entry, TRANSB specifies the form of op( B ) to be used in */
/*           the matrix multiplication as follows: */

/*              TRANSB = 'N' or 'n',  op( B ) = B. */

/*              TRANSB = 'T' or 't',  op( B ) = B'. */

/*              TRANSB = 'C' or 'c',  op( B ) = B'. */

/*           Unchanged on exit. */

/*  M      - INTEGER. */
/*           On entry,  M  specifies  the number  of rows  of the  matrix */
/*           op( A )  and of the  matrix  C.  M  must  be at least  zero. */
/*           Unchanged on exit. */

/*  N      - INTEGER. */
/*           On entry,  N  specifies the number  of columns of the matrix */
/*           op( B ) and the number of columns of the matrix C. N must be */
/*           at least zero. */
/*           Unchanged on exit. */

/*  K      - INTEGER. */
/*           On entry,  K  specifies  the number of columns of the matrix */
/*           op( A ) and the number of rows of the matrix op( B ). K must */
/*           be at least  zero. */
/*           Unchanged on exit. */

/*  ALPHA  - DOUBLE PRECISION. */
/*           On entry, ALPHA specifies the scalar alpha. */
/*           Unchanged on exit. */

/*  A      - DOUBLE PRECISION array of DIMENSION ( LDA, ka ), where ka is */
/*           k  when  TRANSA = 'N' or 'n',  and is  m  otherwise. */
/*           Before entry with  TRANSA = 'N' or 'n',  the leading  m by k */
/*           part of the array  A  must contain the matrix  A,  otherwise */
/*           the leading  k by m  part of the array  A  must contain  the */
/*           matrix A. */
/*           Unchanged on exit. */

/*  LDA    - INTEGER. */
/*           On entry, LDA specifies the first dimension of A as declared */
/*           in the calling (sub) program. When  TRANSA = 'N' or 'n' then */
/*           LDA must be at least  max( 1, m ), otherwise  LDA must be at */
/*           least  max( 1, k ). */
/*           Unchanged on exit. */

/*  B      - DOUBLE PRECISION array of DIMENSION ( LDB, kb ), where kb is */
/*           n  when  TRANSB = 'N' or 'n',  and is  k  otherwise. */
/*           Before entry with  TRANSB = 'N' or 'n',  the leading  k by n */
/*           part of the array  B  must contain the matrix  B,  otherwise */
/*           the leading  n by k  part of the array  B  must contain  the */
/*           matrix B. */
/*           Unchanged on exit. */

/*  LDB    - INTEGER. */
/*           On entry, LDB specifies the first dimension of B as declared */
/*           in the calling (sub) program. When  TRANSB = 'N' or 'n' then */
/*           LDB must be at least  max( 1, k ), otherwise  LDB must be at */
/*           least  max( 1, n ). */
/*           Unchanged on exit. */

/*  BETA   - DOUBLE PRECISION. */
/*           On entry,  BETA  specifies the scalar  beta.  When  BETA  is */
/*           supplied as zero then C need not be set on input. */
/*           Unchanged on exit. */

/*  C      - DOUBLE PRECISION array of DIMENSION ( LDC, n ). */
/*           Before entry, the leading  m by n  part of the array  C must */
/*           contain the matrix  C,  except when  beta  is zero, in which */
/*           case C need not be set on entry. */
/*           On exit, the array  C  is overwritten by the  m by n  matrix */
/*           ( alpha*op( A )*op( B ) + beta*C ). */

/*  LDC    - INTEGER. */
/*           On entry, LDC specifies the first dimension of C as declared */
/*           in  the  calling  (sub)  program.   LDC  must  be  at  least */
/*           max( 1, m ). */
/*           Unchanged on exit. */


/*  Level 3 Blas routine. */

/*  -- Written on 8-February-1989. */
/*     Jack Dongarra, Argonne National Laboratory. */
/*     Iain Duff, AERE Harwell. */
/*     Jeremy Du Croz, Numerical Algorithms Group Ltd. */
/*     Sven Hammarling, Numerical Algorithms Group Ltd. */


/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. Local Scalars .. */
/*     .. Parameters .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not */
/*     transposed and set  NROWA, NCOLA and  NROWB  as the number of rows */
/*     and  columns of  A  and the  number of  rows  of  B  respectively. */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = a_dim1 + 1;
    a -= a_offset;
    b_dim1 = *ldb;
    b_offset = b_dim1 + 1;
    b -= b_offset;
    c_dim1 = *ldc;
    c_offset = c_dim1 + 1;
    c__ -= c_offset;

    /* Function Body */
    nota = lsame_(transa, "N", 1, 1);
    notb = lsame_(transb, "N", 1, 1);
    if (nota) {
	nrowa = *m;
	ncola = *k;
    } else {
	nrowa = *k;
	ncola = *m;
    }
    if (notb) {
	nrowb = *k;
    } else {
	nrowb = *n;
    }

/*     Test the input parameters. */

    info = 0;
    if (! nota && ! lsame_(transa, "C", 1, 1) && ! lsame_(transa, "T", 1, 
	    1)) {
	info = 1;
    } else if (! notb && ! lsame_(transb, "C", 1, 1) && ! lsame_(transb, 
	    "T", 1, 1)) {
	info = 2;
    } else if (*m < 0) {
	info = 3;
    } else if (*n < 0) {
	info = 4;
    } else if (*k < 0) {
	info = 5;
    } else if (*lda < Dymola_max(1,nrowa)) {
	info = 8;
    } else if (*ldb < Dymola_max(1,nrowb)) {
	info = 10;
    } else if (*ldc < Dymola_max(1,*m)) {
	info = 13;
    }
    if (info != 0) {
	xerbla_("DGEMM ", &info, 6);
	return 0;
    }

/*     Quick return if possible. */

    if (*m == 0 || *n == 0 || ((*alpha == 0. || *k == 0) && *beta == 1.)) {
	return 0;
    }

/*     And if  alpha.eq.zero. */

    if (*alpha == 0.) {
	if (*beta == 0.) {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		i__2 = *m;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    c__[i__ + j * c_dim1] = 0.;
/* L10: */
		}
/* L20: */
	    }
	} else {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		i__2 = *m;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    c__[i__ + j * c_dim1] = *beta * c__[i__ + j * c_dim1];
/* L30: */
		}
/* L40: */
	    }
	}
	return 0;
    }

/*     Start the operations. */

    if (notb) {
	if (nota) {

/*           Form  C := alpha*A*B + beta*C. */

	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		if (*beta == 0.) {
		    i__2 = *m;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			c__[i__ + j * c_dim1] = 0.;
/* L50: */
		    }
		} else if (*beta != 1.) {
		    i__2 = *m;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			c__[i__ + j * c_dim1] = *beta * c__[i__ + j * c_dim1];
/* L60: */
		    }
		}
		i__2 = *k;
		for (l = 1; l <= i__2; ++l) {
		    if (b[l + j * b_dim1] != 0.) {
			temp = *alpha * b[l + j * b_dim1];
			i__3 = *m;
			for (i__ = 1; i__ <= i__3; ++i__) {
			    c__[i__ + j * c_dim1] += temp * a[i__ + l * 
				    a_dim1];
/* L70: */
			}
		    }
/* L80: */
		}
/* L90: */
	    }
	} else {

/*           Form  C := alpha*A'*B + beta*C */

	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		i__2 = *m;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    temp = 0.;
		    i__3 = *k;
		    for (l = 1; l <= i__3; ++l) {
			temp += a[l + i__ * a_dim1] * b[l + j * b_dim1];
/* L100: */
		    }
		    if (*beta == 0.) {
			c__[i__ + j * c_dim1] = *alpha * temp;
		    } else {
			c__[i__ + j * c_dim1] = *alpha * temp + *beta * c__[
				i__ + j * c_dim1];
		    }
/* L110: */
		}
/* L120: */
	    }
	}
    } else {
	if (nota) {

/*           Form  C := alpha*A*B' + beta*C */

	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		if (*beta == 0.) {
		    i__2 = *m;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			c__[i__ + j * c_dim1] = 0.;
/* L130: */
		    }
		} else if (*beta != 1.) {
		    i__2 = *m;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			c__[i__ + j * c_dim1] = *beta * c__[i__ + j * c_dim1];
/* L140: */
		    }
		}
		i__2 = *k;
		for (l = 1; l <= i__2; ++l) {
		    if (b[j + l * b_dim1] != 0.) {
			temp = *alpha * b[j + l * b_dim1];
			i__3 = *m;
			for (i__ = 1; i__ <= i__3; ++i__) {
			    c__[i__ + j * c_dim1] += temp * a[i__ + l * 
				    a_dim1];
/* L150: */
			}
		    }
/* L160: */
		}
/* L170: */
	    }
	} else {

/*           Form  C := alpha*A'*B' + beta*C */

	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		i__2 = *m;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    temp = 0.;
		    i__3 = *k;
		    for (l = 1; l <= i__3; ++l) {
			temp += a[l + i__ * a_dim1] * b[j + l * b_dim1];
/* L180: */
		    }
		    if (*beta == 0.) {
			c__[i__ + j * c_dim1] = *alpha * temp;
		    } else {
			c__[i__ + j * c_dim1] = *alpha * temp + *beta * c__[
				i__ + j * c_dim1];
		    }
/* L190: */
		}
/* L200: */
	    }
	}
    }

    return 0;

/*     End of DGEMM . */

} /* dgemm_ */

/* Subroutine */ 
LIBDS_API int dgemv_(trans, m, n, alpha, a, lda, x, incx, beta, y, 
	incy, trans_len)
const char*trans;
const integer*m, *n;
const doublereal*alpha;
doublereal *a;
const integer*lda;
doublereal *x;
const integer*incx;
const doublereal*beta;
doublereal*y;
const integer*incy;
ftnlen trans_len;
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;

    /* Local variables */
    integer i__, j, ix, iy, jx, jy, kx, ky, info;
    doublereal temp;
    integer lenx, leny;

/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  DGEMV  performs one of the matrix-vector operations */

/*     y := alpha*A*x + beta*y,   or   y := alpha*A'*x + beta*y, */

/*  where alpha and beta are scalars, x and y are vectors and A is an */
/*  m by n matrix. */

/*  Parameters */
/*  ========== */

/*  TRANS  - CHARACTER*1. */
/*           On entry, TRANS specifies the operation to be performed as */
/*           follows: */

/*              TRANS = 'N' or 'n'   y := alpha*A*x + beta*y. */

/*              TRANS = 'T' or 't'   y := alpha*A'*x + beta*y. */

/*              TRANS = 'C' or 'c'   y := alpha*A'*x + beta*y. */

/*           Unchanged on exit. */

/*  M      - INTEGER. */
/*           On entry, M specifies the number of rows of the matrix A. */
/*           M must be at least zero. */
/*           Unchanged on exit. */

/*  N      - INTEGER. */
/*           On entry, N specifies the number of columns of the matrix A. */
/*           N must be at least zero. */
/*           Unchanged on exit. */

/*  ALPHA  - DOUBLE PRECISION. */
/*           On entry, ALPHA specifies the scalar alpha. */
/*           Unchanged on exit. */

/*  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ). */
/*           Before entry, the leading m by n part of the array A must */
/*           contain the matrix of coefficients. */
/*           Unchanged on exit. */

/*  LDA    - INTEGER. */
/*           On entry, LDA specifies the first dimension of A as declared */
/*           in the calling (sub) program. LDA must be at least */
/*           max( 1, m ). */
/*           Unchanged on exit. */

/*  X      - DOUBLE PRECISION array of DIMENSION at least */
/*           ( 1 + ( n - 1 )*Dymola_abs( INCX ) ) when TRANS = 'N' or 'n' */
/*           and at least */
/*           ( 1 + ( m - 1 )*Dymola_abs( INCX ) ) otherwise. */
/*           Before entry, the incremented array X must contain the */
/*           vector x. */
/*           Unchanged on exit. */

/*  INCX   - INTEGER. */
/*           On entry, INCX specifies the increment for the elements of */
/*           X. INCX must not be zero. */
/*           Unchanged on exit. */

/*  BETA   - DOUBLE PRECISION. */
/*           On entry, BETA specifies the scalar beta. When BETA is */
/*           supplied as zero then Y need not be set on input. */
/*           Unchanged on exit. */

/*  Y      - DOUBLE PRECISION array of DIMENSION at least */
/*           ( 1 + ( m - 1 )*Dymola_abs( INCY ) ) when TRANS = 'N' or 'n' */
/*           and at least */
/*           ( 1 + ( n - 1 )*Dymola_abs( INCY ) ) otherwise. */
/*           Before entry with BETA non-zero, the incremented array Y */
/*           must contain the vector y. On exit, Y is overwritten by the */
/*           updated vector y. */

/*  INCY   - INTEGER. */
/*           On entry, INCY specifies the increment for the elements of */
/*           Y. INCY must not be zero. */
/*           Unchanged on exit. */


/*  Level 2 Blas routine. */

/*  -- Written on 22-October-1986. */
/*     Jack Dongarra, Argonne National Lab. */
/*     Jeremy Du Croz, Nag Central Office. */
/*     Sven Hammarling, Nag Central Office. */
/*     Richard Hanson, Sandia National Labs. */


/*     .. Parameters .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters. */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = a_dim1 + 1;
    a -= a_offset;
    --x;
    --y;

    /* Function Body */
    info = 0;
    if (! lsame_(trans, "N", 1, 1) && ! lsame_(trans, "T", 1, 1) && ! 
	    lsame_(trans, "C", 1, 1)) {
	info = 1;
    } else if (*m < 0) {
	info = 2;
    } else if (*n < 0) {
	info = 3;
    } else if (*lda < Dymola_max(1,*m)) {
	info = 6;
    } else if (*incx == 0) {
	info = 8;
    } else if (*incy == 0) {
	info = 11;
    }
    if (info != 0) {
	xerbla_("DGEMV ", &info, 6);
	return 0;
    }

/*     Quick return if possible. */

    if (*m == 0 || *n == 0 || (*alpha == 0. && *beta == 1.)) {
	return 0;
    }

/*     Set  LENX  and  LENY, the lengths of the vectors x and y, and set */
/*     up the start points in  X  and  Y. */

    if (lsame_(trans, "N", 1, 1)) {
	lenx = *n;
	leny = *m;
    } else {
	lenx = *m;
	leny = *n;
    }
    if (*incx > 0) {
	kx = 1;
    } else {
	kx = 1 - (lenx - 1) * *incx;
    }
    if (*incy > 0) {
	ky = 1;
    } else {
	ky = 1 - (leny - 1) * *incy;
    }

/*     Start the operations. In this version the elements of A are */
/*     accessed sequentially with one pass through A. */

/*     First form  y := beta*y. */

    if (*beta != 1.) {
	if (*incy == 1) {
	    if (*beta == 0.) {
		i__1 = leny;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    y[i__] = 0.;
/* L10: */
		}
	    } else {
		i__1 = leny;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    y[i__] = *beta * y[i__];
/* L20: */
		}
	    }
	} else {
	    iy = ky;
	    if (*beta == 0.) {
		i__1 = leny;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    y[iy] = 0.;
		    iy += *incy;
/* L30: */
		}
	    } else {
		i__1 = leny;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    y[iy] = *beta * y[iy];
		    iy += *incy;
/* L40: */
		}
	    }
	}
    }
    if (*alpha == 0.) {
	return 0;
    }
    if (lsame_(trans, "N", 1, 1)) {

/*        Form  y := alpha*A*x + y. */

	jx = kx;
	if (*incy == 1) {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		if (x[jx] != 0.) {
		    temp = *alpha * x[jx];
		    i__2 = *m;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			y[i__] += temp * a[i__ + j * a_dim1];
/* L50: */
		    }
		}
		jx += *incx;
/* L60: */
	    }
	} else {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		if (x[jx] != 0.) {
		    temp = *alpha * x[jx];
		    iy = ky;
		    i__2 = *m;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			y[iy] += temp * a[i__ + j * a_dim1];
			iy += *incy;
/* L70: */
		    }
		}
		jx += *incx;
/* L80: */
	    }
	}
    } else {

/*        Form  y := alpha*A'*x + y. */

	jy = ky;
	if (*incx == 1) {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		temp = 0.;
		i__2 = *m;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    temp += a[i__ + j * a_dim1] * x[i__];
/* L90: */
		}
		y[jy] += *alpha * temp;
		jy += *incy;
/* L100: */
	    }
	} else {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		temp = 0.;
		ix = kx;
		i__2 = *m;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    temp += a[i__ + j * a_dim1] * x[ix];
		    ix += *incx;
/* L110: */
		}
		y[jy] += *alpha * temp;
		jy += *incy;
/* L120: */
	    }
	}
    }

    return 0;

/*     End of DGEMV . */

} /* dgemv_ */

/* Subroutine */ DYMOLA_STATIC int dger_(m, n, alpha, x, incx, y, incy, a, lda)
const integer*m, *n;
const doublereal*alpha;
doublereal *x;
const integer*incx;
doublereal *y;
const integer*incy;
doublereal *a;
const integer*lda;
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;

    /* Local variables */
    integer i__, j, ix, jy, kx, info;
    doublereal temp;
	const integer mp=*m;
	const integer np=*n;
	const integer incxp=*incx;
	const integer incyp=*incy;

/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  DGER   performs the rank 1 operation */

/*     A := alpha*x*y' + A, */

/*  where alpha is a scalar, x is an m element vector, y is an n element */
/*  vector and A is an m by n matrix. */

/*  Parameters */
/*  ========== */

/*  M      - INTEGER. */
/*           On entry, M specifies the number of rows of the matrix A. */
/*           M must be at least zero. */
/*           Unchanged on exit. */

/*  N      - INTEGER. */
/*           On entry, N specifies the number of columns of the matrix A. */
/*           N must be at least zero. */
/*           Unchanged on exit. */

/*  ALPHA  - DOUBLE PRECISION. */
/*           On entry, ALPHA specifies the scalar alpha. */
/*           Unchanged on exit. */

/*  X      - DOUBLE PRECISION array of dimension at least */
/*           ( 1 + ( m - 1 )*Dymola_abs( INCX ) ). */
/*           Before entry, the incremented array X must contain the m */
/*           element vector x. */
/*           Unchanged on exit. */

/*  INCX   - INTEGER. */
/*           On entry, INCX specifies the increment for the elements of */
/*           X. INCX must not be zero. */
/*           Unchanged on exit. */

/*  Y      - DOUBLE PRECISION array of dimension at least */
/*           ( 1 + ( n - 1 )*Dymola_abs( INCY ) ). */
/*           Before entry, the incremented array Y must contain the n */
/*           element vector y. */
/*           Unchanged on exit. */

/*  INCY   - INTEGER. */
/*           On entry, INCY specifies the increment for the elements of */
/*           Y. INCY must not be zero. */
/*           Unchanged on exit. */

/*  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ). */
/*           Before entry, the leading m by n part of the array A must */
/*           contain the matrix of coefficients. On exit, A is */
/*           overwritten by the updated matrix. */

/*  LDA    - INTEGER. */
/*           On entry, LDA specifies the first dimension of A as declared */
/*           in the calling (sub) program. LDA must be at least */
/*           max( 1, m ). */
/*           Unchanged on exit. */


/*  Level 2 Blas routine. */

/*  -- Written on 22-October-1986. */
/*     Jack Dongarra, Argonne National Lab. */
/*     Jeremy Du Croz, Nag Central Office. */
/*     Sven Hammarling, Nag Central Office. */
/*     Richard Hanson, Sandia National Labs. */


/*     .. Parameters .. */
/*     .. Local Scalars .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters. */

    /* Parameter adjustments */
    --x;
    --y;
    a_dim1 = *lda;
    a_offset = a_dim1 + 1;
    a -= a_offset;

    /* Function Body */
    info = 0;
    if (mp < 0) {
	info = 1;
    } else if (np < 0) {
	info = 2;
    } else if (incxp == 0) {
	info = 5;
    } else if (incyp == 0) {
	info = 7;
    } else if (*lda < Dymola_max(1,mp)) {
	info = 9;
    }
    if (info != 0) {
	xerbla_("DGER  ", &info, 6);
	return 0;
    }

/*     Quick return if possible. */

    if (mp == 0 || np == 0 || *alpha == 0.) {
	return 0;
    }

/*     Start the operations. In this version the elements of A are */
/*     accessed sequentially with one pass through A. */

    if (incyp > 0) {
	jy = 1;
    } else {
	jy = 1 - (np - 1) * incyp;
    }
    if (incxp == 1) {
	i__1 = np;
	for (j = 1; j <= i__1; ++j) {
	    if (y[jy] != 0.) {
		temp = *alpha * y[jy];
		i__2 = mp;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    a[i__ + j * a_dim1] += x[i__] * temp;
/* L10: */
		}
	    }
	    jy += incyp;
/* L20: */
	}
    } else {
	if (incxp > 0) {
	    kx = 1;
	} else {
	    kx = 1 - (mp - 1) * incxp;
	}
	i__1 = np;
	for (j = 1; j <= i__1; ++j) {
	    if (y[jy] != 0.) {
		temp = *alpha * y[jy];
		ix = kx;
		i__2 = mp;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    a[i__ + j * a_dim1] += x[ix] * temp;
		    ix += incxp;
/* L30: */
		}
	    }
	    jy += incyp;
/* L40: */
	}
    }

    return 0;

/*     End of DGER  . */

} /* dger_ */

#if 0
DYMOLA_STATIC doublereal dnrm2_(n, dx, incx)
const integer*n;
doublereal *dx;
const integer*incx;
{
    /* Initialized data */

    static const doublereal zero = 0.;
    static const doublereal one = 1.;
    static const doublereal cutlo = 8.232e-11;
    static const doublereal cuthi = 1.304e19;

    /* Format strings */
    static char fmt_30[] = "";
    static char fmt_50[] = "";
    static char fmt_70[] = "";
    static char fmt_110[] = "";

    /* System generated locals */
    integer i__1;
    doublereal ret_val, d__1;

    /* Local variables */
    integer i__, j, ix;
    doublereal sum, xmax;
    integer next;
    doublereal hitest;

    /* Assigned format variables */
    static char *next_fmt;

    /* Parameter adjustments */
    --dx;

    /* Function Body */

/*     euclidean norm of the n-vector stored in dx() with storage */
/*     increment incx . */
/*     if    n .le. 0 return with result = 0. */
/*     if n .ge. 1 then incx must be .ge. 1 */

/*           c.l.lawson, 1978 jan 08 */
/*     modified to correct failure to update ix, 1/25/92. */
/*     modified 3/93 to return if incx .le. 0. */

/*     four phase method     using two built-in constants that are */
/*     hopefully applicable to all machines. */
/*         cutlo = maximum of  dsqrt(u/eps)  over all known machines. */
/*         cuthi = minimum of  dsqrt(v)      over all known machines. */
/*     where */
/*         eps = smallest no. such that eps + 1. .gt. 1. */
/*         u   = smallest positive no.   (underflow limit) */
/*         v   = largest  no.            (overflow  limit) */

/*     brief outline of algorithm.. */

/*     phase 1    scans zero components. */
/*     move to phase 2 when a component is nonzero and .le. cutlo */
/*     move to phase 3 when a component is .gt. cutlo */
/*     move to phase 4 when a component is .ge. cuthi/m */
/*     where m = n for x() real and m = 2*n for complex. */

/*     values for cutlo and cuthi.. */
/*     from the environmental parameters listed in the imsl converter */
/*     document the limiting values are as follows.. */
/*     cutlo, s.p.   u/eps = 2**(-102) for  honeywell.  close seconds are */
/*                   univac and dec at 2**(-103) */
/*                   thus cutlo = 2**(-51) = 4.44089e-16 */
/*     cuthi, s.p.   v = 2**127 for univac, honeywell, and dec. */
/*                   thus cuthi = 2**(63.5) = 1.30438e19 */
/*     cutlo, d.p.   u/eps = 2**(-67) for honeywell and dec. */
/*                   thus cutlo = 2**(-33.5) = 8.23181d-11 */
/*     cuthi, d.p.   same as s.p.  cuthi = 1.30438d19 */
/*     data cutlo, cuthi / 8.232d-11,  1.304d19 / */
/*     data cutlo, cuthi / 4.441e-16,  1.304e19 / */

    if (*n > 0 && *incx > 0) {
	goto L10;
    }
    ret_val = zero;
    goto L300;

L10:
    next = 0;
    next_fmt = fmt_30;
    sum = zero;
    i__ = 1;
    ix = 1;
/*                                                 begin main loop */
L20:
    switch ((int)next) {
	case 0: goto L30;
	case 1: goto L50;
	case 2: goto L70;
	case 3: goto L110;
    }
L30:
    if ((d__1 = dx[i__], Dymola_abs(d__1)) > cutlo) {
	goto L85;
    }
    next = 1;
    next_fmt = fmt_50;
    xmax = zero;

/*                        phase 1.  sum is zero */

L50:
    if (dx[i__] == zero) {
	goto L200;
    }
    if ((d__1 = dx[i__], Dymola_abs(d__1)) > cutlo) {
	goto L85;
    }

/*                                prepare for phase 2. */
    next = 2;
    next_fmt = fmt_70;
    goto L105;

/*                                prepare for phase 4. */

L100:
    ix = j;
    next = 3;
    next_fmt = fmt_110;
    sum = sum / dx[i__] / dx[i__];
L105:
    xmax = (d__1 = dx[i__], Dymola_abs(d__1));
    goto L115;

/*                   phase 2.  sum is small. */
/*                             scale to avoid destructive underflow. */

L70:
    if ((d__1 = dx[i__], Dymola_abs(d__1)) > cutlo) {
	goto L75;
    }

/*                     common code for phases 2 and 4. */
/*                     in phase 4 sum is large.  scale to avoid overflow. */

L110:
    if ((d__1 = dx[i__], Dymola_abs(d__1)) <= xmax) {
	goto L115;
    }
/* Computing 2nd power */
    d__1 = xmax / dx[i__];
    sum = one + sum * (d__1 * d__1);
    xmax = (d__1 = dx[i__], Dymola_abs(d__1));
    goto L200;

L115:
/* Computing 2nd power */
    d__1 = dx[i__] / xmax;
    sum += d__1 * d__1;
    goto L200;


/*                  prepare for phase 3. */

L75:
    sum = sum * xmax * xmax;


/*     for real or d.p. set hitest = cuthi/n */
/*     for complex      set hitest = cuthi/(2*n) */

L85:
    hitest = cuthi / (real) (*n);

/*                   phase 3.  sum is mid-range.  no scaling. */

    i__1 = *n;
    for (j = ix; j <= i__1; ++j) {
	if ((d__1 = dx[i__], Dymola_abs(d__1)) >= hitest) {
	    goto L100;
	}
/* Computing 2nd power */
	d__1 = dx[i__];
	sum += d__1 * d__1;
	i__ += *incx;
/* L95: */
    }
    ret_val = sqrt(sum);
    goto L300;

L200:
    ++ix;
    i__ += *incx;
    if (ix <= *n) {
	goto L20;
    }

/*              end of main loop. */

/*              compute square root and adjust for scaling. */

    ret_val = xmax * sqrt(sum);
L300:
    return ret_val;
} /* dnrm2_ */
#else

/* From blas (accompanying lapack)*/
/* Used instead from 2000-08-07 since it can handle +/-inf */
/* and generate sensible results in contrast to the other version of dnrm2 */
DYMOLA_STATIC doublereal dnrm2_(const integer*n, const doublereal *x, const integer*incx)
{


    /* System generated locals */
    const integer np=*n,incxp=*incx;
    doublereal ret_val, d__1;

    /* Local variables */
     doublereal norm, scale, absxi;
     integer ix;
     doublereal ssq;


/*  DNRM2 returns the euclidean norm of a vector via the function   
    name, so that   

       DNRM2 := sqrt( x'*x )   



    -- This version written on 25-October-1982.   
       Modified on 14-October-1993 to inline the call to DLASSQ.   
       Sven Hammarling, Nag Ltd.   


    
   Parameter adjustments   
       Function Body */
#define X(I) x[(I)-1]


    if (np < 1 || incxp < 1) {
	norm = 0.;
    } else if (np == 1) {
	norm = fabs(X(1));
    } else {
	scale = 0.;
	ssq = 1.;
/*        The following loop is equivalent to this call to the LAPACK 
  
          auxiliary routine:   
          CALL DLASSQ( N, X, INCX, SCALE, SSQ ) */

	for (ix = 1; incxp < 0 ? ix >= (np-1)*incxp+1 : 
	     ix <= (np-1)*incxp+1; ix += incxp) {
			 double Xix;
			 Xix=X(ix);
#if defined(_MSC_VER)
			 if (!_finite(Xix)) return DBL_MAX;
#else
			 if (!(Xix==Xix)) return DBL_MAX;
#endif
	    if (Xix != 0.) {
		absxi = (d__1 = Xix, fabs(d__1));
		if (scale < absxi) {
/* Computing 2nd power */
		    d__1 = scale / absxi;
		    ssq = ssq * (d__1 * d__1) + 1.;
		    scale = absxi;
		} else {
/* Computing 2nd power */
		    d__1 = absxi / scale;
		    ssq += d__1 * d__1;
		}
	    }
/* L10: */
	}
	norm = scale * sqrt(ssq);
    }

    ret_val = norm;
    return ret_val;

/*     End of DNRM2. */
}
#undef X
#endif
#endif

/* From blas (accompanying lapack)*/
/* Used instead from 2000-08-07 since it can handle +/-inf */
/* and generate sensible results in contrast to the other version of dnrm2 */
DYMOLA_STATIC doublereal dnrm2_Fast(const integer np, const doublereal *x, const integer incxp)
{

	doublereal d__1;

    /* Local variables */
     doublereal norm, scale, absxi;
     integer ix;
     doublereal ssq;


/*  DNRM2 returns the euclidean norm of a vector via the function   
    name, so that   

       DNRM2 := sqrt( x'*x )   



    -- This version written on 25-October-1982.   
       Modified on 14-October-1993 to inline the call to DLASSQ.   
       Sven Hammarling, Nag Ltd.   


    
   Parameter adjustments   
       Function Body */
#define X(I) x[(I)-1]


    if (np < 1 || incxp < 1) {
	norm = 0.;
    } else if (np == 1) {
	norm = fabs(X(1));
    } else {
	scale = 0.;
	ssq = 1.;
/*        The following loop is equivalent to this call to the LAPACK 
  
          auxiliary routine:   
          CALL DLASSQ( N, X, INCX, SCALE, SSQ ) */

	for (ix = 1; incxp < 0 ? ix >= (np-1)*incxp+1 : 
	     ix <= (np-1)*incxp+1; ix += incxp) {
			 double Xix;
			 Xix=X(ix);
#if defined(_MSC_VER)
			 if (!_finite(Xix)) return DBL_MAX;
#else
			 if (!(Xix==Xix)) return DBL_MAX;
#endif
	    if (Xix != 0.) {
		absxi = fabs(Xix);
		if (scale < absxi) {
/* Computing 2nd power */
		    d__1 = scale / absxi;
		    ssq = ssq * (d__1 * d__1) + 1.;
		    scale = absxi;
		} else {
/* Computing 2nd power */
		    d__1 = absxi / scale;
		    ssq += d__1 * d__1;
		}
	    }
/* L10: */
	}
	norm = scale * sqrt(ssq);
    }

    return norm;

/*     End of DNRM2. */
}
#undef X
/* From blas (accompanying lapack)*/
/* Used instead from 2000-08-07 since it can handle +/-inf */
/* and generate sensible results in contrast to the other version of dnrm2 */
DYMOLA_STATIC doublereal dnrm2_Fast1(const integer np, const doublereal *x)
{

	doublereal d__1;

    /* Local variables */
     doublereal norm, scale, absxi;
     integer ix;
     doublereal ssq,Xix;

     /*  DNRM2 returns the euclidean norm of a vector via the function  name, so that  */

     /*       DNRM2 := sqrt( x'*x )   */



     /*    -- This version written on 25-October-1982.   */
     /*       Modified on 14-October-1993 to inline the call to DLASSQ.   */
     /*       Sven Hammarling, Nag Ltd.   */


    
     /*   Parameter adjustments   */
     /*       Function Body */
#define X(I) x[(I)-1]

	 if (np<1) {
		 norm=0;
	 } else if (np==1) {
		 norm=fabs(X(1));
	 } else {
		 scale = 0.0;
		 ssq = 1.0;
/*        The following loop is equivalent to this call to the LAPACK 
		 
		   auxiliary routine:   
		 CALL DLASSQ( N, X, INCX, SCALE, SSQ ) */
		 
		 for (ix=1; ix <= np; ix += 1) {
			 Xix=X(ix);
#if defined(_MSC_VER)
			 if (!_finite(Xix)) return DBL_MAX;
#else
			 if (!(Xix==Xix)) return DBL_MAX;
#endif
			 if (Xix != 0.) {
				 absxi = fabs(Xix);
				 if (scale < absxi) {
					 /* Computing 2nd power */
					 d__1 = scale / absxi;
					 ssq = ssq * (d__1 * d__1) + 1.;
					 scale = absxi;
				 } else {
					 /* Computing 2nd power */
					 d__1 = absxi / scale;
					 ssq += d__1 * d__1;
				 }
			 }
			 /* L10: */
		 }
		 norm = scale * sqrt(ssq);
	 }

    return norm;

/*     End of DNRM2. */
}
#undef X

#ifndef GODESS
/* Subroutine */ DYMOLA_STATIC int dscal_(n, da, dx, incx)
const integer*n;
const doublereal *da;
doublereal *dx;
const integer*incx;
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    integer i__, m, mp1, nincx;


/*     scales a vector by a constant. */
/*     uses unrolled loops for increment equal to one. */
/*     jack dongarra, linpack, 3/11/78. */
/*     modified 3/93 to return if incx .le. 0. */


    /* Parameter adjustments */
    --dx;

    /* Function Body */
    if (*n <= 0 || *incx <= 0) {
	return 0;
    }
    if (*incx == 1) {
	goto L20;
    }

/*        code for increment not equal to 1 */

    nincx = *n * *incx;
    i__1 = nincx;
    i__2 = *incx;
    for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
	dx[i__] = *da * dx[i__];
/* L10: */
    }
    return 0;

/*        code for increment equal to 1 */


/*        clean-up loop */

L20:
    m = *n % 5;
    if (m == 0) {
	goto L40;
    }
    i__2 = m;
    for (i__ = 1; i__ <= i__2; ++i__) {
	dx[i__] = *da * dx[i__];
/* L30: */
    }
    if (*n < 5) {
	return 0;
    }
L40:
    mp1 = m + 1;
    i__2 = *n;
    for (i__ = mp1; i__ <= i__2; i__ += 5) {
	dx[i__] = *da * dx[i__];
	dx[i__ + 1] = *da * dx[i__ + 1];
	dx[i__ + 2] = *da * dx[i__ + 2];
	dx[i__ + 3] = *da * dx[i__ + 3];
	dx[i__ + 4] = *da * dx[i__ + 4];
/* L50: */
    }
    return 0;
} /* dscal_ */

/* Subroutine */ DYMOLA_STATIC int dswap_(n, dx, incx, dy, incy)
const integer*n;
doublereal *dx;
const integer*incx;
doublereal *dy;
const integer*incy;
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    integer i__, m, ix, iy, mp1;
    doublereal dtemp;


/*     interchanges two vectors. */
/*     uses unrolled loops for increments equal one. */
/*     jack dongarra, linpack, 3/11/78. */


    /* Parameter adjustments */
    --dy;
    --dx;

    /* Function Body */
    if (*n <= 0) {
	return 0;
    }
    if (*incx == 1 && *incy == 1) {
	goto L20;
    }

/*       code for unequal increments or equal increments not equal */
/*         to 1 */

    ix = 1;
    iy = 1;
    if (*incx < 0) {
	ix = (-(*n) + 1) * *incx + 1;
    }
    if (*incy < 0) {
	iy = (-(*n) + 1) * *incy + 1;
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	dtemp = dx[ix];
	dx[ix] = dy[iy];
	dy[iy] = dtemp;
	ix += *incx;
	iy += *incy;
/* L10: */
    }
    return 0;

/*       code for both increments equal to 1 */


/*       clean-up loop */

L20:
    m = *n % 3;
    if (m == 0) {
	goto L40;
    }
    i__1 = m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	dtemp = dx[i__];
	dx[i__] = dy[i__];
	dy[i__] = dtemp;
/* L30: */
    }
    if (*n < 3) {
	return 0;
    }
L40:
    mp1 = m + 1;
    i__1 = *n;
    for (i__ = mp1; i__ <= i__1; i__ += 3) {
	dtemp = dx[i__];
	dx[i__] = dy[i__];
	dy[i__] = dtemp;
	dtemp = dx[i__ + 1];
	dx[i__ + 1] = dy[i__ + 1];
	dy[i__ + 1] = dtemp;
	dtemp = dx[i__ + 2];
	dx[i__ + 2] = dy[i__ + 2];
	dy[i__ + 2] = dtemp;
/* L50: */
    }
    return 0;
} /* dswap_ */

/* Subroutine */ DYMOLA_STATIC int dtrsm_(side, uplo, transa, diag, m, n, alpha, a, lda, b, 
	ldb, side_len, uplo_len, transa_len, diag_len)
const char*side, *uplo, *transa, *diag;
const integer*m, *n;
const doublereal*alpha;
doublereal*a;
const integer*lda;
doublereal *b;
const integer*ldb;
ftnlen side_len;
ftnlen uplo_len;
ftnlen transa_len;
ftnlen diag_len;
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, i__1, i__2, i__3;

    /* Local variables */
    integer i__, j, k, info;
    doublereal temp;
    logical lside;
    integer nrowa;
    logical upper;
    logical nounit;

/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  DTRSM  solves one of the matrix equations */

/*     op( A )*X = alpha*B,   or   X*op( A ) = alpha*B, */

/*  where alpha is a scalar, X and B are m by n matrices, A is a unit, or */
/*  non-unit,  upper or lower triangular matrix  and  op( A )  is one  of */

/*     op( A ) = A   or   op( A ) = A'. */

/*  The matrix X is overwritten on B. */

/*  Parameters */
/*  ========== */

/*  SIDE   - CHARACTER*1. */
/*           On entry, SIDE specifies whether op( A ) appears on the left */
/*           or right of X as follows: */

/*              SIDE = 'L' or 'l'   op( A )*X = alpha*B. */

/*              SIDE = 'R' or 'r'   X*op( A ) = alpha*B. */

/*           Unchanged on exit. */

/*  UPLO   - CHARACTER*1. */
/*           On entry, UPLO specifies whether the matrix A is an upper or */
/*           lower triangular matrix as follows: */

/*              UPLO = 'U' or 'u'   A is an upper triangular matrix. */

/*              UPLO = 'L' or 'l'   A is a lower triangular matrix. */

/*           Unchanged on exit. */

/*  TRANSA - CHARACTER*1. */
/*           On entry, TRANSA specifies the form of op( A ) to be used in */
/*           the matrix multiplication as follows: */

/*              TRANSA = 'N' or 'n'   op( A ) = A. */

/*              TRANSA = 'T' or 't'   op( A ) = A'. */

/*              TRANSA = 'C' or 'c'   op( A ) = A'. */

/*           Unchanged on exit. */

/*  DIAG   - CHARACTER*1. */
/*           On entry, DIAG specifies whether or not A is unit triangular */
/*           as follows: */

/*              DIAG = 'U' or 'u'   A is assumed to be unit triangular. */

/*              DIAG = 'N' or 'n'   A is not assumed to be unit */
/*                                  triangular. */

/*           Unchanged on exit. */

/*  M      - INTEGER. */
/*           On entry, M specifies the number of rows of B. M must be at */
/*           least zero. */
/*           Unchanged on exit. */

/*  N      - INTEGER. */
/*           On entry, N specifies the number of columns of B.  N must be */
/*           at least zero. */
/*           Unchanged on exit. */

/*  ALPHA  - DOUBLE PRECISION. */
/*           On entry,  ALPHA specifies the scalar  alpha. When  alpha is */
/*           zero then  A is not referenced and  B need not be set before */
/*           entry. */
/*           Unchanged on exit. */

/*  A      - DOUBLE PRECISION array of DIMENSION ( LDA, k ), where k is m */
/*           when  SIDE = 'L' or 'l'  and is  n  when  SIDE = 'R' or 'r'. */
/*           Before entry  with  UPLO = 'U' or 'u',  the  leading  k by k */
/*           upper triangular part of the array  A must contain the upper */
/*           triangular matrix  and the strictly lower triangular part of */
/*           A is not referenced. */
/*           Before entry  with  UPLO = 'L' or 'l',  the  leading  k by k */
/*           lower triangular part of the array  A must contain the lower */
/*           triangular matrix  and the strictly upper triangular part of */
/*           A is not referenced. */
/*           Note that when  DIAG = 'U' or 'u',  the diagonal elements of */
/*           A  are not referenced either,  but are assumed to be  unity. */
/*           Unchanged on exit. */

/*  LDA    - INTEGER. */
/*           On entry, LDA specifies the first dimension of A as declared */
/*           in the calling (sub) program.  When  SIDE = 'L' or 'l'  then */
/*           LDA  must be at least  max( 1, m ),  when  SIDE = 'R' or 'r' */
/*           then LDA must be at least max( 1, n ). */
/*           Unchanged on exit. */

/*  B      - DOUBLE PRECISION array of DIMENSION ( LDB, n ). */
/*           Before entry,  the leading  m by n part of the array  B must */
/*           contain  the  right-hand  side  matrix  B,  and  on exit  is */
/*           overwritten by the solution matrix  X. */

/*  LDB    - INTEGER. */
/*           On entry, LDB specifies the first dimension of B as declared */
/*           in  the  calling  (sub)  program.   LDB  must  be  at  least */
/*           max( 1, m ). */
/*           Unchanged on exit. */


/*  Level 3 Blas routine. */


/*  -- Written on 8-February-1989. */
/*     Jack Dongarra, Argonne National Laboratory. */
/*     Iain Duff, AERE Harwell. */
/*     Jeremy Du Croz, Numerical Algorithms Group Ltd. */
/*     Sven Hammarling, Numerical Algorithms Group Ltd. */


/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. Local Scalars .. */
/*     .. Parameters .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters. */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = a_dim1 + 1;
    a -= a_offset;
    b_dim1 = *ldb;
    b_offset = b_dim1 + 1;
    b -= b_offset;

    /* Function Body */
    lside = lsame_(side, "L", 1, 1);
    if (lside) {
	nrowa = *m;
    } else {
	nrowa = *n;
    }
    nounit = lsame_(diag, "N", 1, 1);
    upper = lsame_(uplo, "U", 1, 1);

    info = 0;
    if (! lside && ! lsame_(side, "R", 1, 1)) {
	info = 1;
    } else if (! upper && ! lsame_(uplo, "L", 1, 1)) {
	info = 2;
    } else if (! lsame_(transa, "N", 1, 1) && ! lsame_(transa, "T", 1, 1) 
	    && ! lsame_(transa, "C", 1, 1)) {
	info = 3;
    } else if (! lsame_(diag, "U", 1, 1) && ! lsame_(diag, "N", 1, 1)) {
	info = 4;
    } else if (*m < 0) {
	info = 5;
    } else if (*n < 0) {
	info = 6;
    } else if (*lda < Dymola_max(1,nrowa)) {
	info = 9;
    } else if (*ldb < Dymola_max(1,*m)) {
	info = 11;
    }
    if (info != 0) {
	xerbla_("DTRSM ", &info, 6);
	return 0;
    }

/*     Quick return if possible. */

    if (*n == 0) {
	return 0;
    }

/*     And when  alpha.eq.zero. */

    if (*alpha == 0.) {
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    i__2 = *m;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		b[i__ + j * b_dim1] = 0.;
/* L10: */
	    }
/* L20: */
	}
	return 0;
    }

/*     Start the operations. */

    if (lside) {
	if (lsame_(transa, "N", 1, 1)) {

/*           Form  B := alpha*inv( A )*B. */

	    if (upper) {
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    if (*alpha != 1.) {
			i__2 = *m;
			for (i__ = 1; i__ <= i__2; ++i__) {
			    b[i__ + j * b_dim1] = *alpha * b[i__ + j * b_dim1]
				    ;
/* L30: */
			}
		    }
		    for (k = *m; k >= 1; --k) {
			if (b[k + j * b_dim1] != 0.) {
			    if (nounit) {
				b[k + j * b_dim1] /= a[k + k * a_dim1];
			    }
			    i__2 = k - 1;
			    for (i__ = 1; i__ <= i__2; ++i__) {
				b[i__ + j * b_dim1] -= b[k + j * b_dim1] * a[
					i__ + k * a_dim1];
/* L40: */
			    }
			}
/* L50: */
		    }
/* L60: */
		}
	    } else {
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    if (*alpha != 1.) {
			i__2 = *m;
			for (i__ = 1; i__ <= i__2; ++i__) {
			    b[i__ + j * b_dim1] = *alpha * b[i__ + j * b_dim1]
				    ;
/* L70: */
			}
		    }
		    i__2 = *m;
		    for (k = 1; k <= i__2; ++k) {
			if (b[k + j * b_dim1] != 0.) {
			    if (nounit) {
				b[k + j * b_dim1] /= a[k + k * a_dim1];
			    }
			    i__3 = *m;
			    for (i__ = k + 1; i__ <= i__3; ++i__) {
				b[i__ + j * b_dim1] -= b[k + j * b_dim1] * a[
					i__ + k * a_dim1];
/* L80: */
			    }
			}
/* L90: */
		    }
/* L100: */
		}
	    }
	} else {

/*           Form  B := alpha*inv( A' )*B. */

	    if (upper) {
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    i__2 = *m;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			temp = *alpha * b[i__ + j * b_dim1];
			i__3 = i__ - 1;
			for (k = 1; k <= i__3; ++k) {
			    temp -= a[k + i__ * a_dim1] * b[k + j * b_dim1];
/* L110: */
			}
			if (nounit) {
			    temp /= a[i__ + i__ * a_dim1];
			}
			b[i__ + j * b_dim1] = temp;
/* L120: */
		    }
/* L130: */
		}
	    } else {
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    for (i__ = *m; i__ >= 1; --i__) {
			temp = *alpha * b[i__ + j * b_dim1];
			i__2 = *m;
			for (k = i__ + 1; k <= i__2; ++k) {
			    temp -= a[k + i__ * a_dim1] * b[k + j * b_dim1];
/* L140: */
			}
			if (nounit) {
			    temp /= a[i__ + i__ * a_dim1];
			}
			b[i__ + j * b_dim1] = temp;
/* L150: */
		    }
/* L160: */
		}
	    }
	}
    } else {
	if (lsame_(transa, "N", 1, 1)) {

/*           Form  B := alpha*B*inv( A ). */

	    if (upper) {
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    if (*alpha != 1.) {
			i__2 = *m;
			for (i__ = 1; i__ <= i__2; ++i__) {
			    b[i__ + j * b_dim1] = *alpha * b[i__ + j * b_dim1]
				    ;
/* L170: */
			}
		    }
		    i__2 = j - 1;
		    for (k = 1; k <= i__2; ++k) {
			if (a[k + j * a_dim1] != 0.) {
			    i__3 = *m;
			    for (i__ = 1; i__ <= i__3; ++i__) {
				b[i__ + j * b_dim1] -= a[k + j * a_dim1] * b[
					i__ + k * b_dim1];
/* L180: */
			    }
			}
/* L190: */
		    }
		    if (nounit) {
			temp = 1. / a[j + j * a_dim1];
			i__2 = *m;
			for (i__ = 1; i__ <= i__2; ++i__) {
			    b[i__ + j * b_dim1] = temp * b[i__ + j * b_dim1];
/* L200: */
			}
		    }
/* L210: */
		}
	    } else {
		for (j = *n; j >= 1; --j) {
		    if (*alpha != 1.) {
			i__1 = *m;
			for (i__ = 1; i__ <= i__1; ++i__) {
			    b[i__ + j * b_dim1] = *alpha * b[i__ + j * b_dim1]
				    ;
/* L220: */
			}
		    }
		    i__1 = *n;
		    for (k = j + 1; k <= i__1; ++k) {
			if (a[k + j * a_dim1] != 0.) {
			    i__2 = *m;
			    for (i__ = 1; i__ <= i__2; ++i__) {
				b[i__ + j * b_dim1] -= a[k + j * a_dim1] * b[
					i__ + k * b_dim1];
/* L230: */
			    }
			}
/* L240: */
		    }
		    if (nounit) {
			temp = 1. / a[j + j * a_dim1];
			i__1 = *m;
			for (i__ = 1; i__ <= i__1; ++i__) {
			    b[i__ + j * b_dim1] = temp * b[i__ + j * b_dim1];
/* L250: */
			}
		    }
/* L260: */
		}
	    }
	} else {

/*           Form  B := alpha*B*inv( A' ). */

	    if (upper) {
		for (k = *n; k >= 1; --k) {
		    if (nounit) {
			temp = 1. / a[k + k * a_dim1];
			i__1 = *m;
			for (i__ = 1; i__ <= i__1; ++i__) {
			    b[i__ + k * b_dim1] = temp * b[i__ + k * b_dim1];
/* L270: */
			}
		    }
		    i__1 = k - 1;
		    for (j = 1; j <= i__1; ++j) {
			if (a[j + k * a_dim1] != 0.) {
			    temp = a[j + k * a_dim1];
			    i__2 = *m;
			    for (i__ = 1; i__ <= i__2; ++i__) {
				b[i__ + j * b_dim1] -= temp * b[i__ + k * 
					b_dim1];
/* L280: */
			    }
			}
/* L290: */
		    }
		    if (*alpha != 1.) {
			i__1 = *m;
			for (i__ = 1; i__ <= i__1; ++i__) {
			    b[i__ + k * b_dim1] = *alpha * b[i__ + k * b_dim1]
				    ;
/* L300: */
			}
		    }
/* L310: */
		}
	    } else {
		i__1 = *n;
		for (k = 1; k <= i__1; ++k) {
		    if (nounit) {
			temp = 1. / a[k + k * a_dim1];
			i__2 = *m;
			for (i__ = 1; i__ <= i__2; ++i__) {
			    b[i__ + k * b_dim1] = temp * b[i__ + k * b_dim1];
/* L320: */
			}
		    }
		    i__2 = *n;
		    for (j = k + 1; j <= i__2; ++j) {
			if (a[j + k * a_dim1] != 0.) {
			    temp = a[j + k * a_dim1];
			    i__3 = *m;
			    for (i__ = 1; i__ <= i__3; ++i__) {
				b[i__ + j * b_dim1] -= temp * b[i__ + k * 
					b_dim1];
/* L330: */
			    }
			}
/* L340: */
		    }
		    if (*alpha != 1.) {
			i__2 = *m;
			for (i__ = 1; i__ <= i__2; ++i__) {
			    b[i__ + k * b_dim1] = *alpha * b[i__ + k * b_dim1]
				    ;
/* L350: */
			}
		    }
/* L360: */
		}
	    }
	}
    }

    return 0;

/*     End of DTRSM . */

} /* dtrsm_ */

/* Subroutine */ DYMOLA_STATIC int dtrsv_(uplo, trans, diag, n, a, lda, x, incx, uplo_len, 
	trans_len, diag_len)
const char*uplo, *trans, *diag;
const integer*n;
doublereal *a;
const integer*lda;
doublereal *x;
const integer*incx;
ftnlen uplo_len;
ftnlen trans_len;
ftnlen diag_len;
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;

    /* Local variables */
    integer i__, j, ix, jx, kx, info;
    doublereal temp;
    logical nounit;

/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  DTRSV  solves one of the systems of equations */

/*     A*x = b,   or   A'*x = b, */

/*  where b and x are n element vectors and A is an n by n unit, or */
/*  non-unit, upper or lower triangular matrix. */

/*  No test for singularity or near-singularity is included in this */
/*  routine. Such tests must be performed before calling this routine. */

/*  Parameters */
/*  ========== */

/*  UPLO   - CHARACTER*1. */
/*           On entry, UPLO specifies whether the matrix is an upper or */
/*           lower triangular matrix as follows: */

/*              UPLO = 'U' or 'u'   A is an upper triangular matrix. */

/*              UPLO = 'L' or 'l'   A is a lower triangular matrix. */

/*           Unchanged on exit. */

/*  TRANS  - CHARACTER*1. */
/*           On entry, TRANS specifies the equations to be solved as */
/*           follows: */

/*              TRANS = 'N' or 'n'   A*x = b. */

/*              TRANS = 'T' or 't'   A'*x = b. */

/*              TRANS = 'C' or 'c'   A'*x = b. */

/*           Unchanged on exit. */

/*  DIAG   - CHARACTER*1. */
/*           On entry, DIAG specifies whether or not A is unit */
/*           triangular as follows: */

/*              DIAG = 'U' or 'u'   A is assumed to be unit triangular. */

/*              DIAG = 'N' or 'n'   A is not assumed to be unit */
/*                                  triangular. */

/*           Unchanged on exit. */

/*  N      - INTEGER. */
/*           On entry, N specifies the order of the matrix A. */
/*           N must be at least zero. */
/*           Unchanged on exit. */

/*  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ). */
/*           Before entry with  UPLO = 'U' or 'u', the leading n by n */
/*           upper triangular part of the array A must contain the upper */
/*           triangular matrix and the strictly lower triangular part of */
/*           A is not referenced. */
/*           Before entry with UPLO = 'L' or 'l', the leading n by n */
/*           lower triangular part of the array A must contain the lower */
/*           triangular matrix and the strictly upper triangular part of */
/*           A is not referenced. */
/*           Note that when  DIAG = 'U' or 'u', the diagonal elements of */
/*           A are not referenced either, but are assumed to be unity. */
/*           Unchanged on exit. */

/*  LDA    - INTEGER. */
/*           On entry, LDA specifies the first dimension of A as declared */
/*           in the calling (sub) program. LDA must be at least */
/*           max( 1, n ). */
/*           Unchanged on exit. */

/*  X      - DOUBLE PRECISION array of dimension at least */
/*           ( 1 + ( n - 1 )*Dymola_abs( INCX ) ). */
/*           Before entry, the incremented array X must contain the n */
/*           element right-hand side vector b. On exit, X is overwritten */
/*           with the solution vector x. */

/*  INCX   - INTEGER. */
/*           On entry, INCX specifies the increment for the elements of */
/*           X. INCX must not be zero. */
/*           Unchanged on exit. */


/*  Level 2 Blas routine. */

/*  -- Written on 22-October-1986. */
/*     Jack Dongarra, Argonne National Lab. */
/*     Jeremy Du Croz, Nag Central Office. */
/*     Sven Hammarling, Nag Central Office. */
/*     Richard Hanson, Sandia National Labs. */


/*     .. Parameters .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters. */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = a_dim1 + 1;
    a -= a_offset;
    --x;

    /* Function Body */
    kx = 0;
    info = 0;
    if (! lsame_(uplo, "U", 1, 1) && ! lsame_(uplo, "L", 1, 1)) {
	info = 1;
    } else if (! lsame_(trans, "N", 1, 1) && ! lsame_(trans, "T", 1, 1) &&
	     ! lsame_(trans, "C", 1, 1)) {
	info = 2;
    } else if (! lsame_(diag, "U", 1, 1) && ! lsame_(diag, "N", 1, 1)) {
	info = 3;
    } else if (*n < 0) {
	info = 4;
    } else if (*lda < Dymola_max(1,*n)) {
	info = 6;
    } else if (*incx == 0) {
	info = 8;
    }
    if (info != 0) {
	xerbla_("DTRSV ", &info, 6);
	return 0;
    }

/*     Quick return if possible. */

    if (*n == 0) {
	return 0;
    }

    nounit = lsame_(diag, "N", 1, 1);

/*     Set up the start point in X if the increment is not unity. This */
/*     will be  ( N - 1 )*INCX  too small for descending loops. */

    if (*incx <= 0) {
	kx = 1 - (*n - 1) * *incx;
    } else if (*incx != 1) {
	kx = 1;
    }

/*     Start the operations. In this version the elements of A are */
/*     accessed sequentially with one pass through A. */

    if (lsame_(trans, "N", 1, 1)) {

/*        Form  x := inv( A )*x. */

	if (lsame_(uplo, "U", 1, 1)) {
	    if (*incx == 1) {
		for (j = *n; j >= 1; --j) {
		    if (x[j] != 0.) {
			if (nounit) {
			    x[j] /= a[j + j * a_dim1];
			}
			temp = x[j];
			for (i__ = j - 1; i__ >= 1; --i__) {
			    x[i__] -= temp * a[i__ + j * a_dim1];
/* L10: */
			}
		    }
/* L20: */
		}
	    } else {
		jx = kx + (*n - 1) * *incx;
		for (j = *n; j >= 1; --j) {
		    if (x[jx] != 0.) {
			if (nounit) {
			    x[jx] /= a[j + j * a_dim1];
			}
			temp = x[jx];
			ix = jx;
			for (i__ = j - 1; i__ >= 1; --i__) {
			    ix -= *incx;
			    x[ix] -= temp * a[i__ + j * a_dim1];
/* L30: */
			}
		    }
		    jx -= *incx;
/* L40: */
		}
	    }
	} else {
	    if (*incx == 1) {
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    if (x[j] != 0.) {
			if (nounit) {
			    x[j] /= a[j + j * a_dim1];
			}
			temp = x[j];
			i__2 = *n;
			for (i__ = j + 1; i__ <= i__2; ++i__) {
			    x[i__] -= temp * a[i__ + j * a_dim1];
/* L50: */
			}
		    }
/* L60: */
		}
	    } else {
		jx = kx;
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    if (x[jx] != 0.) {
			if (nounit) {
			    x[jx] /= a[j + j * a_dim1];
			}
			temp = x[jx];
			ix = jx;
			i__2 = *n;
			for (i__ = j + 1; i__ <= i__2; ++i__) {
			    ix += *incx;
			    x[ix] -= temp * a[i__ + j * a_dim1];
/* L70: */
			}
		    }
		    jx += *incx;
/* L80: */
		}
	    }
	}
    } else {

/*        Form  x := inv( A' )*x. */

	if (lsame_(uplo, "U", 1, 1)) {
	    if (*incx == 1) {
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    temp = x[j];
		    i__2 = j - 1;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			temp -= a[i__ + j * a_dim1] * x[i__];
/* L90: */
		    }
		    if (nounit) {
			temp /= a[j + j * a_dim1];
		    }
		    x[j] = temp;
/* L100: */
		}
	    } else {
		jx = kx;
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    temp = x[jx];
		    ix = kx;
		    i__2 = j - 1;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			temp -= a[i__ + j * a_dim1] * x[ix];
			ix += *incx;
/* L110: */
		    }
		    if (nounit) {
			temp /= a[j + j * a_dim1];
		    }
		    x[jx] = temp;
		    jx += *incx;
/* L120: */
		}
	    }
	} else {
	    if (*incx == 1) {
		for (j = *n; j >= 1; --j) {
		    temp = x[j];
		    i__1 = j + 1;
		    for (i__ = *n; i__ >= i__1; --i__) {
			temp -= a[i__ + j * a_dim1] * x[i__];
/* L130: */
		    }
		    if (nounit) {
			temp /= a[j + j * a_dim1];
		    }
		    x[j] = temp;
/* L140: */
		}
	    } else {
		kx += (*n - 1) * *incx;
		jx = kx;
		for (j = *n; j >= 1; --j) {
		    temp = x[jx];
		    ix = kx;
		    i__1 = j + 1;
		    for (i__ = *n; i__ >= i__1; --i__) {
			temp -= a[i__ + j * a_dim1] * x[ix];
			ix -= *incx;
/* L150: */
		    }
		    if (nounit) {
			temp /= a[j + j * a_dim1];
		    }
		    x[jx] = temp;
		    jx -= *incx;
/* L160: */
		}
	    }
	}
    }

    return 0;

/*     End of DTRSV . */

} /* dtrsv_ */

DYMOLA_STATIC integer idamax_(n, dx, incx)
const integer*n;
doublereal *dx;
const integer*incx;
{
    /* System generated locals */
    integer ret_val, i__1;
    doublereal d__1;

    /* Local variables */
    integer i__, ix;
    doublereal dmax__;


/*     finds the index of element having max. absolute value. */
/*     jack dongarra, linpack, 3/11/78. */
/*     modified 3/93 to return if incx .le. 0. */


    /* Parameter adjustments */
    --dx;

    /* Function Body */
    ret_val = 0;
    if (*n < 1 || *incx <= 0) {
	return ret_val;
    }
    ret_val = 1;
    if (*n == 1) {
	return ret_val;
    }
    if (*incx == 1) {
	goto L20;
    }

/*        code for increment not equal to 1 */

    ix = 1;
    dmax__ = Dymola_abs(dx[1]);
    ix += *incx;
    i__1 = *n;
    for (i__ = 2; i__ <= i__1; ++i__) {
	if ((d__1 = dx[ix], Dymola_abs(d__1)) <= dmax__) {
	    goto L5;
	}
	ret_val = i__;
	dmax__ = (d__1 = dx[ix], Dymola_abs(d__1));
L5:
	ix += *incx;
/* L10: */
    }
    return ret_val;

/*        code for increment equal to 1 */

L20:
    dmax__ = Dymola_abs(dx[1]);
    i__1 = *n;
    for (i__ = 2; i__ <= i__1; ++i__) {
	if ((d__1 = dx[i__], Dymola_abs(d__1)) <= dmax__) {
	    goto L30;
	}
	ret_val = i__;
	dmax__ = (d__1 = dx[i__], Dymola_abs(d__1));
L30:
	;
    }
    return ret_val;
} /* idamax_ */





/* Subroutine */ DYMOLA_STATIC int dgecon_(norm, n, a, lda, anorm, rcond, work, iwork, info,
	 norm_len)
char *norm;
const integer*n;
doublereal *a;
const integer*lda;
doublereal *anorm, *rcond, *work;
integer *iwork, *info;
ftnlen norm_len;
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1;
    doublereal d__1;

    /* Local variables */
    doublereal sl;
    integer ix;
    doublereal su;
    integer kase, kase1;
    doublereal scale;
    doublereal ainvnm;
    logical onenrm;
    char normin[1];
    doublereal smlnum;
	struct dlacon_clean_struct dummy;


/*  -- LAPACK routine (version 1.1) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
/*     Courant Institute, Argonne National Lab, and Rice University */
/*     February 29, 1992 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  DGECON estimates the reciprocal of the condition number of a general */
/*  real matrix A, in either the 1-norm or the infinity-norm, using */
/*  the LU factorization computed by DGETRF. */

/*  An estimate is obtained for norm(inv(A)), and the reciprocal of the */
/*  condition number is computed as */
/*     RCOND = 1 / ( norm(A) * norm(inv(A)) ). */

/*  Arguments */
/*  ========= */

/*  NORM    (input) CHARACTER*1 */
/*          Specifies whether the 1-norm condition number or the */
/*          infinity-norm condition number is required: */
/*          = '1' or 'O':  1-norm; */
/*          = 'I':         Infinity-norm. */

/*  N       (input) INTEGER */
/*          The order of the matrix A.  N >= 0. */

/*  A       (input) DOUBLE PRECISION array, dimension (LDA,N) */
/*          The factors L and U from the factorization A = P*L*U */
/*          as computed by DGETRF. */

/*  LDA     (input) INTEGER */
/*          The leading dimension of the array A.  LDA >= max(1,N). */

/*  ANORM   (input) DOUBLE PRECISION */
/*          If NORM = '1' or 'O', the 1-norm of the original matrix A. */
/*          If NORM = 'I', the infinity-norm of the original matrix A. */

/*  RCOND   (output) DOUBLE PRECISION */
/*          The reciprocal of the condition number of the matrix A, */
/*          computed as RCOND = 1/(norm(A) * norm(inv(A))). */

/*  WORK    (workspace) DOUBLE PRECISION array, dimension (4*N) */

/*  IWORK   (workspace) INTEGER array, dimension (N) */

/*  INFO    (output) INTEGER */
/*          = 0:  successful exit */
/*          < 0:  if INFO = -i, the i-th argument had an illegal value */

/*  ===================================================================== */

/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters. */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = a_dim1 + 1;
    a -= a_offset;
    --work;
    --iwork;

    /* Function Body */
    *info = 0;
    onenrm = *(unsigned char *)norm == '1' || lsame_(norm, "O", 1, 1);
    if (! onenrm && ! lsame_(norm, "I", 1, 1)) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*lda < Dymola_max(1,*n)) {
	*info = -4;
    } else if (*anorm < 0.) {
	*info = -5;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DGECON", &i__1, 6);
	return 0;
    }

/*     Quick return if possible */

    *rcond = 0.;
    if (*n == 0) {
	*rcond = 1.;
	return 0;
    } else if (*anorm == 0.) {
	return 0;
    }

    smlnum = DBL_MIN;

/*     Estimate the norm of inv(A). */

    ainvnm = 0.;
    *(unsigned char *)normin = 'N';
    if (onenrm) {
	kase1 = 1;
    } else {
	kase1 = 2;
    }
    kase = 0;
L10:
    dlacon_clean(n, &work[*n + 1], &work[1], &iwork[1], &ainvnm, &kase,&dummy);
    if (kase != 0) {
	if (kase == kase1) {

/*           Multiply by inv(L). */

	    dlatrs_("Lower", "No transpose", "Unit", normin, n, &a[a_offset], 
		    lda, &work[1], &sl, &work[(*n << 1) + 1], info, 5, 12, 
		    4, 1);

/*           Multiply by inv(U). */

	    dlatrs_("Upper", "No transpose", "Non-unit", normin, n, &a[
		    a_offset], lda, &work[1], &su, &work[*n * 3 + 1], info, 
		    5, 12, 8, 1);
	} else {

/*           Multiply by inv(U'). */

	    dlatrs_("Upper", "Transpose", "Non-unit", normin, n, &a[a_offset],
		     lda, &work[1], &su, &work[*n * 3 + 1], info, 5, 9, 8, 
		    1);

/*           Multiply by inv(L'). */

	    dlatrs_("Lower", "Transpose", "Unit", normin, n, &a[a_offset], 
		    lda, &work[1], &sl, &work[(*n << 1) + 1], info, 5, 9, 
		    4, 1);
	}

/*        Divide X by 1/(SL*SU) if doing so will not cause overflow. */

	scale = sl * su;
	*(unsigned char *)normin = 'Y';
	if (scale != 1.) {
	    ix = idamax_(n, &work[1], &c__1);
	    if (scale < (d__1 = work[ix], Dymola_abs(d__1)) * smlnum || scale == 0.) 
		    {
		goto L20;
	    }
	    drscl_(n, &scale, &work[1], &c__1);
	}
	goto L10;
    }

/*     Compute the estimate of the reciprocal condition number. */

    if (ainvnm != 0.) {
	*rcond = 1. / ainvnm / *anorm;
    }

L20:
    return 0;

/*     End of DGECON */

} /* dgecon_ */

/* Subroutine */ DYMOLA_STATIC int dgeqpf_(m, n, a, lda, jpvt, tau, work, info)
const integer*m, *n;
doublereal *a;
const integer*lda;
integer *jpvt;
doublereal *tau, *work;
integer *info;
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;
    doublereal d__1, d__2;

    /* Local variables */
    integer i__, j, ma, mn;
    doublereal aii;
    integer pvt;
    doublereal temp;
    doublereal temp2;
    integer itemp;


/*  -- LAPACK test routine (version 1.1) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
/*     Courant Institute, Argonne National Lab, and Rice University */
/*     March 31, 1993 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  DGEQPF computes a QR factorization with column pivoting of a */
/*  real M-by-N matrix A: A*P = Q*R. */

/*  Arguments */
/*  ========= */

/*  M       (input) INTEGER */
/*          The number of rows of the matrix A. M >= 0. */

/*  N       (input) INTEGER */
/*          The number of columns of the matrix A. N >= 0 */

/*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*          On entry, the M-by-N matrix A. */
/*          On exit, the upper triangle of the array contains the */
/*          min(M,N)-by-N upper triangular matrix R; the elements */
/*          below the diagonal, together with the array TAU, */
/*          represent the orthogonal matrix Q as a product of */
/*          min(m,n) elementary reflectors. */

/*  LDA     (input) INTEGER */
/*          The leading dimension of the array A. LDA >= max(1,M). */

/*  JPVT    (input/output) INTEGER array, dimension (N) */
/*          On entry, if JPVT(i) .ne. 0, the i-th column of A is permuted */
/*          to the front of A*P (a leading column); if JPVT(i) = 0, */
/*          the i-th column of A is a free column. */
/*          On exit, if JPVT(i) = k, then the i-th column of A*P */
/*          was the k-th column of A. */

/*  TAU     (output) DOUBLE PRECISION array, dimension (min(M,N)) */
/*          The scalar factors of the elementary reflectors. */

/*  WORK    (workspace) DOUBLE PRECISION array, dimension (3*N) */

/*  INFO    (output) INTEGER */
/*          = 0:  successful exit */
/*          < 0:  if INFO = -i, the i-th argument had an illegal value */

/*  Further Details */
/*  =============== */

/*  The matrix Q is represented as a product of elementary reflectors */

/*     Q = H(1) H(2) . . . H(n) */

/*  Each H(i) has the form */

/*     H = I - tau * v * v' */

/*  where tau is a real scalar, and v is a real vector with */
/*  v(1:i-1) = 0 and v(i) = 1; v(i+1:m) is stored on exit in A(i+1:m,i). */

/*  The matrix P is represented in jpvt as follows: If */
/*     jpvt(j) = i */
/*  then the jth column of P is the ith canonical unit vector. */

/*  ===================================================================== */

/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input arguments */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = a_dim1 + 1;
    a -= a_offset;
    --jpvt;
    --tau;
    --work;

    /* Function Body */
    *info = 0;
    if (*m < 0) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*lda < Dymola_max(1,*m)) {
	*info = -4;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DGEQPF", &i__1, 6);
	return 0;
    }

    mn = Dymola_min(*m,*n);

/*     Move initial columns up front */

    itemp = 1;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (jpvt[i__] != 0) {
	    if (i__ != itemp) {
		dswap_(m, &a[i__ * a_dim1 + 1], &c__1, &a[itemp * a_dim1 + 1],
			 &c__1);
		jpvt[i__] = jpvt[itemp];
		jpvt[itemp] = i__;
	    } else {
		jpvt[i__] = i__;
	    }
	    ++itemp;
	} else {
	    jpvt[i__] = i__;
	}
/* L10: */
    }
    --itemp;

/*     Compute the QR factorization and update remaining columns */

    if (itemp > 0) {
	ma = Dymola_min(itemp,*m);
	dgeqr2_(m, &ma, &a[a_offset], lda, &tau[1], &work[1], info);
	if (ma < *n) {
	    i__1 = *n - ma;
	    dorm2r_("Left", "Transpose", m, &i__1, &ma, &a[a_offset], lda, &
		    tau[1], &a[(ma + 1) * a_dim1 + 1], lda, &work[1], info, 
		    4, 9);
	}
    }

    if (itemp < mn) {

/*        Initialize partial column norms. The first n entries of */
/*        work store the exact column norms. */

	i__1 = *n;
	for (i__ = itemp + 1; i__ <= i__1; ++i__) {
	    i__2 = *m - itemp;
	    work[i__] = dnrm2_Fast1(i__2, &a[itemp + 1 + i__ * a_dim1]);
	    work[*n + i__] = work[i__];
/* L20: */
	}

/*        Compute factorization */

	i__1 = mn;
	for (i__ = itemp + 1; i__ <= i__1; ++i__) {

/*           Determine ith pivot column and swap if necessary */

	    i__2 = *n - i__ + 1;
	    pvt = i__ - 1 + idamax_(&i__2, &work[i__], &c__1);

	    if (pvt != i__) {
		dswap_(m, &a[pvt * a_dim1 + 1], &c__1, &a[i__ * a_dim1 + 1], &
			c__1);
		itemp = jpvt[pvt];
		jpvt[pvt] = jpvt[i__];
		jpvt[i__] = itemp;
		work[pvt] = work[i__];
		work[*n + pvt] = work[*n + i__];
	    }

/*           Generate elementary reflector H(i) */

	    if (i__ < *m) {
		i__2 = *m - i__ + 1;
		dlarfg_(&i__2, &a[i__ + i__ * a_dim1], &a[i__ + 1 + i__ * 
			a_dim1], &c__1, &tau[i__]);
	    } else {
		dlarfg_(&c__1, &a[*m + *m * a_dim1], &a[*m + *m * a_dim1], &
			c__1, &tau[*m]);
	    }

	    if (i__ < *n) {

/*              Apply H(i) to A(i:m,i+1:n) from the left */

		aii = a[i__ + i__ * a_dim1];
		a[i__ + i__ * a_dim1] = 1.;
		i__2 = *m - i__ + 1;
		i__3 = *n - i__;
		dlarf_("LEFT", &i__2, &i__3, &a[i__ + i__ * a_dim1], &c__1, &
			tau[i__], &a[i__ + (i__ + 1) * a_dim1], lda, &work[(*
			n << 1) + 1], 4);
		a[i__ + i__ * a_dim1] = aii;
	    }

/*           Update partial column norms */

	    i__2 = *n;
	    for (j = i__ + 1; j <= i__2; ++j) {
		if (work[j] != 0.) {
/* Computing 2nd power */
		    d__2 = (d__1 = a[i__ + j * a_dim1], Dymola_abs(d__1)) / work[j];
		    temp = 1. - d__2 * d__2;
		    temp = Dymola_max(temp,0.);
/* Computing 2nd power */
		    d__1 = work[j] / work[*n + j];
		    temp2 = temp * .05 * (d__1 * d__1) + 1.;
		    if (temp2 == 1.) {
			if (*m - i__ > 0) {
			    i__3 = *m - i__;
			    work[j] = dnrm2_Fast1(i__3, &a[i__ + 1 + j * a_dim1]);
			    work[*n + j] = work[j];
			} else {
			    work[j] = 0.;
			    work[*n + j] = 0.;
			}
		    } else {
			work[j] *= sqrt(temp);
		    }
		}
/* L30: */
	    }

/* L40: */
	}
    }
    return 0;

/*     End of DGEQPF */

} /* dgeqpf_ */

/* Subroutine */ DYMOLA_STATIC int dgeqr2_(m, n, a, lda, tau, work, info)
const integer*m, *n;
doublereal *a;
const integer*lda;
doublereal *tau, *work;
integer *info;
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;

    /* Local variables */
    integer i__, k;
    doublereal aii;

/*  -- LAPACK routine (version 1.1) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
/*     Courant Institute, Argonne National Lab, and Rice University */
/*     February 29, 1992 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  DGEQR2 computes a QR factorization of a real m by n matrix A: */
/*  A = Q * R. */

/*  Arguments */
/*  ========= */

/*  M       (input) INTEGER */
/*          The number of rows of the matrix A.  M >= 0. */

/*  N       (input) INTEGER */
/*          The number of columns of the matrix A.  N >= 0. */

/*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*          On entry, the m by n matrix A. */
/*          On exit, the elements on and above the diagonal of the array */
/*          contain the min(m,n) by n upper trapezoidal matrix R (R is */
/*          upper triangular if m >= n); the elements below the diagonal, */
/*          with the array TAU, represent the orthogonal matrix Q as a */
/*          product of elementary reflectors (see Further Details). */

/*  LDA     (input) INTEGER */
/*          The leading dimension of the array A.  LDA >= max(1,M). */

/*  TAU     (output) DOUBLE PRECISION array, dimension (min(M,N)) */
/*          The scalar factors of the elementary reflectors (see Further */
/*          Details). */

/*  WORK    (workspace) DOUBLE PRECISION array, dimension (N) */

/*  INFO    (output) INTEGER */
/*          = 0: successful exit */
/*          < 0: if INFO = -i, the i-th argument had an illegal value */

/*  Further Details */
/*  =============== */

/*  The matrix Q is represented as a product of elementary reflectors */

/*     Q = H(1) H(2) . . . H(k), where k = min(m,n). */

/*  Each H(i) has the form */

/*     H(i) = I - tau * v * v' */

/*  where tau is a real scalar, and v is a real vector with */
/*  v(1:i-1) = 0 and v(i) = 1; v(i+1:m) is stored on exit in A(i+1:m,i), */
/*  and tau in TAU(i). */

/*  ===================================================================== */

/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input arguments */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = a_dim1 + 1;
    a -= a_offset;
    --tau;
    --work;

    /* Function Body */
    *info = 0;
    if (*m < 0) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*lda < Dymola_max(1,*m)) {
	*info = -4;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DGEQR2", &i__1, 6);
	return 0;
    }

    k = Dymola_min(*m,*n);

    i__1 = k;
    for (i__ = 1; i__ <= i__1; ++i__) {

/*        Generate elementary reflector H(i) to annihilate A(i+1:m,i) */

	i__2 = *m - i__ + 1;
/* Computing MIN */
	i__3 = i__ + 1;
	dlarfg_(&i__2, &a[i__ + i__ * a_dim1], &a[Dymola_min(i__3,*m) + i__ * a_dim1]
		, &c__1, &tau[i__]);
	if (i__ < *n) {

/*           Apply H(i) to A(i:m,i+1:n) from the left */

	    aii = a[i__ + i__ * a_dim1];
	    a[i__ + i__ * a_dim1] = 1.;
	    i__2 = *m - i__ + 1;
	    i__3 = *n - i__;
	    dlarf_("Left", &i__2, &i__3, &a[i__ + i__ * a_dim1], &c__1, &tau[
		    i__], &a[i__ + (i__ + 1) * a_dim1], lda, &work[1], 4);
	    a[i__ + i__ * a_dim1] = aii;
	}
/* L10: */
    }
    return 0;

/*     End of DGEQR2 */

} /* dgeqr2_ */

/* Subroutine */ DYMOLA_STATIC int dgetf2_(m, n, a, lda, ipiv, info)
const integer*m, *n;
doublereal *a;
const integer*lda;
integer *ipiv, *info;
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;
    doublereal d__1;

    /* Local variables */
    integer j, jp;


/*  -- LAPACK routine (version 1.1) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
/*     Courant Institute, Argonne National Lab, and Rice University */
/*     June 30, 1992 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  DGETF2 computes an LU factorization of a general m-by-n matrix A */
/*  using partial pivoting with row interchanges. */

/*  The factorization has the form */
/*     A = P * L * U */
/*  where P is a permutation matrix, L is lower triangular with unit */
/*  diagonal elements (lower trapezoidal if m > n), and U is upper */
/*  triangular (upper trapezoidal if m < n). */

/*  This is the right-looking Level 2 BLAS version of the algorithm. */

/*  Arguments */
/*  ========= */

/*  M       (input) INTEGER */
/*          The number of rows of the matrix A.  M >= 0. */

/*  N       (input) INTEGER */
/*          The number of columns of the matrix A.  N >= 0. */

/*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*          On entry, the m by n matrix to be factored. */
/*          On exit, the factors L and U from the factorization */
/*          A = P*L*U; the unit diagonal elements of L are not stored. */

/*  LDA     (input) INTEGER */
/*          The leading dimension of the array A.  LDA >= max(1,M). */

/*  IPIV    (output) INTEGER array, dimension (Dymola_min(M,N)) */
/*          The pivot indices; for 1 <= i <= min(M,N), row i of the */
/*          matrix was interchanged with row IPIV(i). */

/*  INFO    (output) INTEGER */
/*          = 0: successful exit */
/*          < 0: if INFO = -k, the k-th argument had an illegal value */
/*          > 0: if INFO = k, U(k,k) is exactly zero. The factorization */
/*               has been completed, but the factor U is exactly */
/*               singular, and division by zero will occur if it is used */
/*               to solve a system of equations. */

/*  ===================================================================== */

/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters. */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = a_dim1 + 1;
    a -= a_offset;
    --ipiv;

    /* Function Body */
    *info = 0;
    if (*m < 0) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*lda < Dymola_max(1,*m)) {
	*info = -4;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DGETF2", &i__1, 6);
	return 0;
    }

/*     Quick return if possible */

    if (*m == 0 || *n == 0) {
	return 0;
    }

    i__1 = Dymola_min(*m,*n);
    for (j = 1; j <= i__1; ++j) {

/*        Find pivot and test for singularity. */

	i__2 = *m - j + 1;
	jp = j - 1 + idamax_(&i__2, &a[j + j * a_dim1], &c__1);
	ipiv[j] = jp;
	if (a[jp + j * a_dim1] != 0.) {

/*           Apply the interchange to columns 1:N. */

	    if (jp != j) {
		dswap_(n, &a[j + a_dim1], lda, &a[jp + a_dim1], lda);
	    }

/*           Compute elements J+1:M of J-th column. */

	    if (j < *m) {
		i__2 = *m - j;
		d__1 = 1. / a[j + j * a_dim1];
		dscal_(&i__2, &d__1, &a[j + 1 + j * a_dim1], &c__1);
	    }

	} else if (*info == 0) {

	    *info = j;
	}

	if (j < Dymola_min(*m,*n)) {

/*           Update trailing submatrix. */

	    i__2 = *m - j;
	    i__3 = *n - j;
	    dger_(&i__2, &i__3, &c_b246, &a[j + 1 + j * a_dim1], &c__1, &a[j 
		    + (j + 1) * a_dim1], lda, &a[j + 1 + (j + 1) * a_dim1], 
		    lda);
	}
/* L10: */
    }
    return 0;

/*     End of DGETF2 */

} /* dgetf2_ */


/* Subroutine */ DYMOLA_STATIC int dgetrf_(m, n, a, lda, ipiv, info)
const integer*m, *n;
doublereal *a;
const integer*lda;
integer *ipiv, *info;
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4, i__5;

    /* Local variables */
    integer i__, j, jb, nb;
    integer iinfo;


/*  -- LAPACK routine (version 1.1) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
/*     Courant Institute, Argonne National Lab, and Rice University */
/*     March 31, 1993 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  DGETRF computes an LU factorization of a general M-by-N matrix A */
/*  using partial pivoting with row interchanges. */

/*  The factorization has the form */
/*     A = P * L * U */
/*  where P is a permutation matrix, L is lower triangular with unit */
/*  diagonal elements (lower trapezoidal if m > n), and U is upper */
/*  triangular (upper trapezoidal if m < n). */

/*  This is the right-looking Level 3 BLAS version of the algorithm. */

/*  Arguments */
/*  ========= */

/*  M       (input) INTEGER */
/*          The number of rows of the matrix A.  M >= 0. */

/*  N       (input) INTEGER */
/*          The number of columns of the matrix A.  N >= 0. */

/*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*          On entry, the M-by-N matrix to be factored. */
/*          On exit, the factors L and U from the factorization */
/*          A = P*L*U; the unit diagonal elements of L are not stored. */

/*  LDA     (input) INTEGER */
/*          The leading dimension of the array A.  LDA >= max(1,M). */

/*  IPIV    (output) INTEGER array, dimension (min(M,N)) */
/*          The pivot indices; for 1 <= i <= min(M,N), row i of the */
/*          matrix was interchanged with row IPIV(i). */

/*  INFO    (output) INTEGER */
/*          = 0:  successful exit */
/*          < 0:  if INFO = -i, the i-th argument had an illegal value */
/*          > 0:  if INFO = i, U(i,i) is exactly zero. The factorization */
/*                has been completed, but the factor U is exactly */
/*                singular, and division by zero will occur if it is used */
/*                to solve a system of equations. */

/*  ===================================================================== */

/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters. */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = a_dim1 + 1;
    a -= a_offset;
    --ipiv;

    /* Function Body */
    *info = 0;
    if (*m < 0) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*lda < Dymola_max(1,*m)) {
	*info = -4;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DGETRF", &i__1, 6);
	return 0;
    }

/*     Quick return if possible */

    if (*m == 0 || *n == 0) {
	return 0;
    }

/*     Determine the block size for this environment. */

    nb = ilaenv_(&c__1, "DGETRF", " ", m, n, &c_n1, &c_n1, 6, 1);
    if (nb <= 1 || nb >= Dymola_min(*m,*n)) {

/*        Use unblocked code. */

	dgetf2_(m, n, &a[a_offset], lda, &ipiv[1], info);
    } else {

/*        Use blocked code. */

	i__1 = Dymola_min(*m,*n);
	i__2 = nb;
	for (j = 1; i__2 < 0 ? j >= i__1 : j <= i__1; j += i__2) {
/* Computing MIN */
	    i__3 = Dymola_min(*m,*n) - j + 1;
	    jb = Dymola_min(i__3,nb);

/*           Factor diagonal and subdiagonal blocks and test for e
xact */
/*           singularity. */

	    i__3 = *m - j + 1;
	    dgetf2_(&i__3, &jb, &a[j + j * a_dim1], lda, &ipiv[j], &iinfo);

/*           Adjust INFO and the pivot indices. */

	    if (*info == 0 && iinfo > 0) {
		*info = iinfo + j - 1;
	    }
/* Computing MIN */
	    i__4 = *m, i__5 = j + jb - 1;
	    i__3 = Dymola_min(i__4,i__5);
	    for (i__ = j; i__ <= i__3; ++i__) {
		ipiv[i__] = j - 1 + ipiv[i__];
/* L10: */
	    }

/*           Apply interchanges to columns 1:J-1. */

	    i__3 = j - 1;
	    i__4 = j + jb - 1;
	    dlaswp_(&i__3, &a[a_offset], lda, &j, &i__4, &ipiv[1], &c__1);

	    if (j + jb <= *n) {

/*              Apply interchanges to columns J+JB:N. */

		i__3 = *n - j - jb + 1;
		i__4 = j + jb - 1;
		dlaswp_(&i__3, &a[(j + jb) * a_dim1 + 1], lda, &j, &i__4, &
			ipiv[1], &c__1);

/*              Compute block row of U. */

		i__3 = *n - j - jb + 1;
		dtrsm_("Left", "Lower", "No transpose", "Unit", &jb, &i__3, &
			c_b263, &a[j + j * a_dim1], lda, &a[j + (j + jb) * 
			a_dim1], lda, 4, 5, 12, 4);
		if (j + jb <= *m) {

/*                 Update trailing submatrix. */

		    i__3 = *m - j - jb + 1;
		    i__4 = *n - j - jb + 1;
		    dgemm_("No transpose", "No transpose", &i__3, &i__4, &jb, 
			    &c_b246, &a[j + jb + j * a_dim1], lda, &a[j + (j 
			    + jb) * a_dim1], lda, &c_b263, &a[j + jb + (j + 
			    jb) * a_dim1], lda, 12, 12);
		}
	    }
/* L20: */
	}
    }
    return 0;

/*     End of DGETRF */

} /* dgetrf_ */

/* Subroutine */ DYMOLA_STATIC int dgetrs_(trans, n, nrhs, a, lda, ipiv, b, ldb, info, 
	trans_len)
const char*trans;
const integer*n, *nrhs;
doublereal *a;
const integer*lda;
integer *ipiv;
doublereal *b;
const integer*ldb;
integer *info;
ftnlen trans_len;
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, i__1;

    /* Local variables */
    logical notran;


/*  -- LAPACK routine (version 1.1) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
/*     Courant Institute, Argonne National Lab, and Rice University */
/*     March 31, 1993 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  DGETRS solves a system of linear equations */
/*     A * X = B  or  A' * X = B */
/*  with a general N-by-N matrix A using the LU factorization computed */
/*  by DGETRF. */

/*  Arguments */
/*  ========= */

/*  TRANS   (input) CHARACTER*1 */
/*          Specifies the form of the system of equations: */
/*          = 'N':  A * X = B  (No transpose) */
/*          = 'T':  A'* X = B  (Transpose) */
/*          = 'C':  A'* X = B  (Conjugate transpose = Transpose) */

/*  N       (input) INTEGER */
/*          The order of the matrix A.  N >= 0. */

/*  NRHS    (input) INTEGER */
/*          The number of right hand sides, i.e., the number of columns */
/*          of the matrix B.  NRHS >= 0. */

/*  A       (input) DOUBLE PRECISION array, dimension (LDA,N) */
/*          The factors L and U from the factorization A = P*L*U */
/*          as computed by DGETRF. */

/*  LDA     (input) INTEGER */
/*          The leading dimension of the array A.  LDA >= max(1,N). */

/*  IPIV    (input) INTEGER array, dimension (N) */
/*          The pivot indices from DGETRF; for 1<=i<=N, row i of the */
/*          matrix was interchanged with row IPIV(i). */

/*  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS) */
/*          On entry, the right hand side matrix B. */
/*          On exit, the solution matrix X. */

/*  LDB     (input) INTEGER */
/*          The leading dimension of the array B.  LDB >= max(1,N). */

/*  INFO    (output) INTEGER */
/*          = 0:  successful exit */
/*          < 0:  if INFO = -i, the i-th argument had an illegal value */

/*  ===================================================================== */

/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters. */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = a_dim1 + 1;
    a -= a_offset;
    --ipiv;
    b_dim1 = *ldb;
    b_offset = b_dim1 + 1;
    b -= b_offset;

    /* Function Body */
    *info = 0;
    notran = lsame_(trans, "N", 1, 1);
    if (! notran && ! lsame_(trans, "T", 1, 1) && ! lsame_(trans, "C", 1, 
	    1)) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*nrhs < 0) {
	*info = -3;
    } else if (*lda < Dymola_max(1,*n)) {
	*info = -5;
    } else if (*ldb < Dymola_max(1,*n)) {
	*info = -8;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DGETRS", &i__1, 6);
	return 0;
    }

/*     Quick return if possible */

    if (*n == 0 || *nrhs == 0) {
	return 0;
    }

    if (notran) {

/*        Solve A * X = B. */

/*        Apply row interchanges to the right hand sides. */

	dlaswp_(nrhs, &b[b_offset], ldb, &c__1, n, &ipiv[1], &c__1);

/*        Solve L*X = B, overwriting B with X. */

	dtrsm_("Left", "Lower", "No transpose", "Unit", n, nrhs, &c_b263, &a[
		a_offset], lda, &b[b_offset], ldb, 4, 5, 12, 4);

/*        Solve U*X = B, overwriting B with X. */

	dtrsm_("Left", "Upper", "No transpose", "Non-unit", n, nrhs, &c_b263, 
		&a[a_offset], lda, &b[b_offset], ldb, 4, 5, 12, 8);
    } else {

/*        Solve A' * X = B. */

/*        Solve U'*X = B, overwriting B with X. */

	dtrsm_("Left", "Upper", "Transpose", "Non-unit", n, nrhs, &c_b263, &a[
		a_offset], lda, &b[b_offset], ldb, 4, 5, 9, 8);

/*        Solve L'*X = B, overwriting B with X. */

	dtrsm_("Left", "Lower", "Transpose", "Unit", n, nrhs, &c_b263, &a[
		a_offset], lda, &b[b_offset], ldb, 4, 5, 9, 4);

/*        Apply row interchanges to the solution vectors. */

	dlaswp_(nrhs, &b[b_offset], ldb, &c__1, n, &ipiv[1], &c_n1);
    }

    return 0;

/*     End of DGETRS */

} /* dgetrs_ */
#endif

#ifndef GODESS
/* Subroutine */ DYMOLA_STATIC int dlabad_(small, large)
doublereal *small, *large;
{


/*  -- LAPACK auxiliary routine (version 1.0b) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
/*     Courant Institute, Argonne National Lab, and Rice University */
/*     October 31, 1992 */

/*     .. Scalar Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  DLABAD takes as input the values computed by SLAMCH for underflow and */
/*  overflow, and returns the square root of each of these values if the */
/*  log of LARGE is sufficiently large.  This subroutine is intended to */
/*  identify machines with a large exponent range, such as the Crays, and */
/*  redefine the underflow and overflow limits to be the square roots of */
/*  the values computed by DLAMCH.  This subroutine is needed because */
/*  DLAMCH does not compensate for poor arithmetic in the upper half of */
/*  the exponent range, as is found on a Cray. */

/*  Arguments */
/*  ========= */

/*  SMALL   (input/output) DOUBLE PRECISION */
/*          On entry, the underflow threshold as computed by DLAMCH. */
/*          On exit, if LOG10(LARGE) is sufficiently large, the square */
/*          root of SMALL, otherwise unchanged. */

/*  LARGE   (input/output) DOUBLE PRECISION */
/*          On entry, the overflow threshold as computed by DLAMCH. */
/*          On exit, if LOG10(LARGE) is sufficiently large, the square */
/*          root of LARGE, otherwise unchanged. */

/*  ===================================================================== */

/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     If it looks like we're on a Cray, take the square root of */
/*     SMALL and LARGE to avoid overflow and underflow problems. */

/*      IF( LOG10( LARGE ).GT.2000.D0 ) THEN */
    *small = sqrt(*small);
    *large = sqrt(*large);
/*      END IF */

    return 0;

/*     End of DLABAD */

} /* dlabad_ */

/* Subroutine */ DYMOLA_STATIC int dlacon_(n, v, x, isgn, est, kase)
const integer*n;
doublereal *v, *x;
integer *isgn;
doublereal *est;
integer *kase;
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Builtin functions */

    /* Local variables */
    static integer i__, j, iter;
    static doublereal temp;
    static integer jump;
    static integer jlast;
    static doublereal altsgn, estold;


/*  -- LAPACK auxiliary routine (version 1.1) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
/*     Courant Institute, Argonne National Lab, and Rice University */
/*     February 29, 1992 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  DLACON estimates the 1-norm of a square, real matrix A. */
/*  Reverse communication is used for evaluating matrix-vector products. */

/*  Arguments */
/*  ========= */

/*  N      (input) INTEGER */
/*         The order of the matrix.  N >= 1. */

/*  V      (workspace) DOUBLE PRECISION array, dimension (N) */
/*         On the final return, V = A*W,  where  EST = norm(V)/norm(W) */
/*         (W is not returned). */

/*  X      (input/output) DOUBLE PRECISION array, dimension (N) */
/*         On an intermediate return, X should be overwritten by */
/*               A * X,   if KASE=1, */
/*               A' * X,  if KASE=2, */
/*         and DLACON must be re-called with all the other parameters */
/*         unchanged. */

/*  ISGN   (workspace) INTEGER array, dimension (N) */

/*  EST    (output) DOUBLE PRECISION */
/*         An estimate (a lower bound) for norm(A). */

/*  KASE   (input/output) INTEGER */
/*         On the initial call to DLACON, KASE should be 0. */
/*         On an intermediate return, KASE will be 1 or 2, indicating */
/*         whether X should be overwritten by A * X  or A' * X. */
/*         On the final return from DLACON, KASE will again be 0. */

/*  Further Details */
/*  ======= ======= */

/*  Contributed by Nick Higham, University of Manchester. */
/*  Originally named SONEST, dated March 16, 1988. */

/*  Reference: N.J. Higham, "FORTRAN codes for estimating the one-norm of */
/*  a real or complex matrix, with applications to condition estimation", */
/*  ACM Trans. Math. Soft., vol. 14, no. 4, pp. 381-396, December 1988. */

/*  ===================================================================== */

/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Save statement .. */
/*     .. */
/*     .. Executable Statements .. */

    /* Parameter adjustments */
    --isgn;
    --x;
    --v;

    /* Function Body */
    if (*kase == 0) {
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    x[i__] = 1. / (doublereal) (*n);
/* L10: */
	}
	*kase = 1;
	jump = 1;
	return 0;
    }

    switch ((int)jump) {
	case 1:  goto L20;
	case 2:  goto L40;
	case 3:  goto L70;
	case 4:  goto L110;
	case 5:  goto L140;
    }

/*     ................ ENTRY   (JUMP = 1) */
/*     FIRST ITERATION.  X HAS BEEN OVERWRITTEN BY A*X. */

L20:
    if (*n == 1) {
	v[1] = x[1];
	*est = Dymola_abs(v[1]);
/*        ... QUIT */
	goto L150;
    }
    *est = dasum_(n, &x[1], &c__1);

    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
		if (x[i__]>=0.0) {
			x[i__]=1;
			isgn[i__]=1;
		} else {
			x[i__]=-1;
			isgn[i__]=-1;
		}
/* L30: */
    }
    *kase = 2;
    jump = 2;
    return 0;

/*     ................ ENTRY   (JUMP = 2) */
/*     FIRST ITERATION.  X HAS BEEN OVERWRITTEN BY TRANDPOSE(A)*X. */

L40:
    j = idamax_(n, &x[1], &c__1);
    iter = 2;

/*     MAIN LOOP - ITERATIONS 2,3,...,ITMAX. */

L50:
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	x[i__] = 0.;
/* L60: */
    }
    x[j] = 1.;
    *kase = 1;
    jump = 3;
    return 0;

/*     ................ ENTRY   (JUMP = 3) */
/*     X HAS BEEN OVERWRITTEN BY A*X. */

L70:
    dcopy_(n, &x[1], &c__1, &v[1], &c__1);
    estold = *est;
    *est = dasum_(n, &v[1], &c__1);
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
		if (((x[i__]>=0.0) ? 1 : -1)!=isgn[i__]) {
	/*d__1 = d_sign(&c_b263, &x[i__]);
	if (i_dnnt(&d__1) != isgn[i__]){ */
	    goto L90;
	}
/* L80: */
    }
/*     REPEATED SIGN VECTOR DETECTED, HENCE ALGORITHM HAS CONVERGED. */
    goto L120;

L90:
/*     TEST FOR CYCLING. */
    if (*est <= estold) {
	goto L120;
    }

    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
		if (x[i__]>=0.0) {
			x[i__]=1;
			isgn[i__]=1;
		} else {
			x[i__]=-1;
			isgn[i__]=-1;
		}
				/*
	x[i__] = d_sign(&c_b263, &x[i__]);
	isgn[i__] = i_dnnt(&x[i__]);*/
/* L100: */
    }
    *kase = 2;
    jump = 4;
    return 0;

/*     ................ ENTRY   (JUMP = 4) */
/*     X HAS BEEN OVERWRITTEN BY TRANDPOSE(A)*X. */

L110:
    jlast = j;
    j = idamax_(n, &x[1], &c__1);
    if (x[jlast] != (d__1 = x[j], Dymola_abs(d__1)) && iter < 5) {
	++iter;
	goto L50;
    }

/*     ITERATION COMPLETE.  FINAL STAGE. */

L120:
    altsgn = 1.;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	x[i__] = altsgn * ((doublereal) (i__ - 1) / (doublereal) (*n - 1) + 
		1.);
	altsgn = -altsgn;
/* L130: */
    }
    *kase = 1;
    jump = 5;
    return 0;

/*     ................ ENTRY   (JUMP = 5) */
/*     X HAS BEEN OVERWRITTEN BY A*X. */

L140:
    temp = dasum_(n, &x[1], &c__1) / (doublereal) (*n * 3) * 2.;
    if (temp > *est) {
	dcopy_(n, &x[1], &c__1, &v[1], &c__1);
	*est = temp;
    }

L150:
    *kase = 0;
    return 0;

/*     End of DLACON */

} /* dlacon_ */
#endif

/* Subroutine */ DYMOLA_STATIC int dlacon_clean(n, v, x, isgn, est, kase, staticarea)
const integer*n;
doublereal *v, *x;
integer *isgn;
doublereal *est;
integer *kase;
struct dlacon_clean_struct *staticarea;
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Builtin functions */

    /* Local variables */
#define i__ (staticarea->i__)
#define j (staticarea->j)
#define iter (staticarea->iter)
#define temp (staticarea->temp)
#define jump (staticarea->jump)
#define jlast (staticarea->jlast)
#define altsgn (staticarea->altsgn)
#define estold (staticarea->estold)


/*  -- LAPACK auxiliary routine (version 1.1) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
/*     Courant Institute, Argonne National Lab, and Rice University */
/*     February 29, 1992 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  DLACON estimates the 1-norm of a square, real matrix A. */
/*  Reverse communication is used for evaluating matrix-vector products. */

	/* Cleaned up to avoid internal memory. Requires extra input argument */
/*  Arguments */
/*  ========= */

/*  N      (input) INTEGER */
/*         The order of the matrix.  N >= 1. */

/*  V      (workspace) DOUBLE PRECISION array, dimension (N) */
/*         On the final return, V = A*W,  where  EST = norm(V)/norm(W) */
/*         (W is not returned). */

/*  X      (input/output) DOUBLE PRECISION array, dimension (N) */
/*         On an intermediate return, X should be overwritten by */
/*               A * X,   if KASE=1, */
/*               A' * X,  if KASE=2, */
/*         and DLACON must be re-called with all the other parameters */
/*         unchanged. */

/*  ISGN   (workspace) INTEGER array, dimension (N) */

/*  EST    (output) DOUBLE PRECISION */
/*         An estimate (a lower bound) for norm(A). */

/*  KASE   (input/output) INTEGER */
/*         On the initial call to DLACON, KASE should be 0. */
/*         On an intermediate return, KASE will be 1 or 2, indicating */
/*         whether X should be overwritten by A * X  or A' * X. */
/*         On the final return from DLACON, KASE will again be 0. */

/*  Further Details */
/*  ======= ======= */

/*  Contributed by Nick Higham, University of Manchester. */
/*  Originally named SONEST, dated March 16, 1988. */

/*  Reference: N.J. Higham, "FORTRAN codes for estimating the one-norm of */
/*  a real or complex matrix, with applications to condition estimation", */
/*  ACM Trans. Math. Soft., vol. 14, no. 4, pp. 381-396, December 1988. */

/*  ===================================================================== */

/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Save statement .. */
/*     .. */
/*     .. Executable Statements .. */

    /* Parameter adjustments */
    --isgn;
    --x;
    --v;

    /* Function Body */
    if (*kase == 0) {
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    x[i__] = 1. / (doublereal) (*n);
/* L10: */
	}
	*kase = 1;
	jump = 1;
	return 0;
    }

    switch ((int)jump) {
	case 1:  goto L20;
	case 2:  goto L40;
	case 3:  goto L70;
	case 4:  goto L110;
	case 5:  goto L140;
    }

/*     ................ ENTRY   (JUMP = 1) */
/*     FIRST ITERATION.  X HAS BEEN OVERWRITTEN BY A*X. */

L20:
    if (*n == 1) {
	v[1] = x[1];
	*est = Dymola_abs(v[1]);
/*        ... QUIT */
	goto L150;
    }
    *est = dasum_(n, &x[1], &c__1);

    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
		if (x[i__]>=0.0) {
			x[i__]=1;
			isgn[i__]=1;
		} else {
			x[i__]=-1;
			isgn[i__]=-1;
		}
/* L30: */
    }
    *kase = 2;
    jump = 2;
    return 0;

/*     ................ ENTRY   (JUMP = 2) */
/*     FIRST ITERATION.  X HAS BEEN OVERWRITTEN BY TRANDPOSE(A)*X. */

L40:
    j = idamax_(n, &x[1], &c__1);
    iter = 2;

/*     MAIN LOOP - ITERATIONS 2,3,...,ITMAX. */

L50:
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	x[i__] = 0.;
/* L60: */
    }
    x[j] = 1.;
    *kase = 1;
    jump = 3;
    return 0;

/*     ................ ENTRY   (JUMP = 3) */
/*     X HAS BEEN OVERWRITTEN BY A*X. */

L70:
    dcopy_(n, &x[1], &c__1, &v[1], &c__1);
    estold = *est;
    *est = dasum_(n, &v[1], &c__1);
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
		if (((x[i__]>=0.0) ? 1 : -1)!=isgn[i__]) {
	/*d__1 = d_sign(&c_b263, &x[i__]);
	if (i_dnnt(&d__1) != isgn[i__]){ */
	    goto L90;
	}
/* L80: */
    }
/*     REPEATED SIGN VECTOR DETECTED, HENCE ALGORITHM HAS CONVERGED. */
    goto L120;

L90:
/*     TEST FOR CYCLING. */
    if (*est <= estold) {
	goto L120;
    }

    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
		if (x[i__]>=0.0) {
			x[i__]=1;
			isgn[i__]=1;
		} else {
			x[i__]=-1;
			isgn[i__]=-1;
		}
				/*
	x[i__] = d_sign(&c_b263, &x[i__]);
	isgn[i__] = i_dnnt(&x[i__]);*/
/* L100: */
    }
    *kase = 2;
    jump = 4;
    return 0;

/*     ................ ENTRY   (JUMP = 4) */
/*     X HAS BEEN OVERWRITTEN BY TRANDPOSE(A)*X. */

L110:
    jlast = j;
    j = idamax_(n, &x[1], &c__1);
    if (x[jlast] != (d__1 = x[j], Dymola_abs(d__1)) && iter < 5) {
	++iter;
	goto L50;
    }

/*     ITERATION COMPLETE.  FINAL STAGE. */

L120:
    altsgn = 1.;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	x[i__] = altsgn * ((doublereal) (i__ - 1) / (doublereal) (*n - 1) + 
		1.);
	altsgn = -altsgn;
/* L130: */
    }
    *kase = 1;
    jump = 5;
    return 0;

/*     ................ ENTRY   (JUMP = 5) */
/*     X HAS BEEN OVERWRITTEN BY A*X. */

L140:
    temp = dasum_(n, &x[1], &c__1) / (doublereal) (*n * 3) * 2.;
    if (temp > *est) {
	dcopy_(n, &x[1], &c__1, &v[1], &c__1);
	*est = temp;
    }

L150:
    *kase = 0;
    return 0;

/*     End of DLACON_CLEAN */

} /* dlacon_clean */

#undef i__
#undef j
#undef iter
#undef temp
#undef jump
#undef jlast
#undef altsgn
#undef estold

/* Subroutine */ DYMOLA_STATIC int dlaic1_(job, j, x, sest, w, gamma, sestpr, s, c__)
const integer *job;
integer *j;
doublereal *x, *sest, *w, *gamma, *sestpr, *s, *c__;
{
    /* System generated locals */
    doublereal d__1, d__2, d__3, d__4;

    /* Local variables */
    doublereal b, t, s1, s2, eps, tmp;
    doublereal sine, test, zeta1, zeta2, alpha, norma;
    doublereal absgam, absalp, cosine, absest;
	doublereal sp,cp;


/*  -- LAPACK auxiliary routine (version 1.1) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
/*     Courant Institute, Argonne National Lab, and Rice University */
/*     October 31, 1992 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  DLAIC1 applies one step of incremental condition estimation in */
/*  its simplest version: */

/*  Let x, twonorm(x) = 1, be an approximate singular vector of an j-by-j */
/*  lower triangular matrix L, such that */
/*           twonorm(L*x) = sest */
/*  Then DLAIC1 computes sestpr, s, c such that */
/*  the vector */
/*                  [ s*x ] */
/*           xhat = [  c  ] */
/*  is an approximate singular vector of */
/*                  [ L     0  ] */
/*           Lhat = [ w' gamma ] */
/*  in the sense that */
/*           twonorm(Lhat*xhat) = sestpr. */

/*  Depending on JOB, an estimate for the largest or smallest singular */
/*  value is computed. */

/*  Note that [s c]' and sestpr**2 is an eigenpair of the system */

/*      diag(sest*sest, 0) + [alpha  gamma] * [ alpha ] */
/*                                            [ gamma ] */

/*  where  alpha =  x'*w. */

/*  Arguments */
/*  ========= */

/*  JOB     (input) INTEGER */
/*          = 1: an estimate for the largest singular value is computed. */
/*          = 2: an estimate for the smallest singular value is computed. */

/*  J       (input) INTEGER */
/*          Length of X and W */

/*  X       (input) DOUBLE PRECISION array, dimension (J) */
/*          The j-vector x. */

/*  SEST    (input) DOUBLE PRECISION */
/*          Estimated singular value of j by j matrix L */

/*  W       (input) DOUBLE PRECISION array, dimension (J) */
/*          The j-vector w. */

/*  GAMMA   (input) DOUBLE PRECISION */
/*          The diagonal element gamma. */

/*  SEDTPR  (output) DOUBLE PRECISION */
/*          Estimated singular value of (j+1) by (j+1) matrix Lhat. */

/*  S       (output) DOUBLE PRECISION */
/*          Sine needed in forming xhat. */

/*  C       (output) DOUBLE PRECISION */
/*          Cosine needed in forming xhat. */

/*  ===================================================================== */

/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. Executable Statements .. */

    /* Parameter adjustments */
    --w;
    --x;

    /* Function Body */
    eps = DBL_EPSILON;
    /* alpha = ddot_(j, &x[1], &c__1, &w[1], &c__1); */
    alpha = ddot1_(j, &x[1], &w[1]);

    absalp = Dymola_abs(alpha);
    absgam = Dymola_abs(*gamma);
    absest = Dymola_abs(*sest);

    if (*job == 1) {

/*        Estimating largest singular value */

/*        special cases */

	if (absest == 0.) {
	    s1 = Dymola_max(absgam,absalp);
	    if (s1 == 0.) {
		sp = 0.;
		cp = 1.;
		*sestpr = 0.;
	    } else {
		sp = alpha / s1;
		cp = *gamma / s1;
		tmp = sqrt(sp * sp + cp * cp);
		sp /= tmp;
		cp /= tmp;
		*sestpr = s1 * tmp;
	    }
	} else if (absgam <= eps * absest) {
	    sp = 1.;
	    cp = 0.;
	    tmp = Dymola_max(absest,absalp);
	    s1 = absest / tmp;
	    s2 = absalp / tmp;
	    *sestpr = tmp * sqrt(s1 * s1 + s2 * s2);
	} else if (absalp <= eps * absest) {
	    s1 = absgam;
	    s2 = absest;
	    if (s1 <= s2) {
		sp = 1.;
		cp = 0.;
		*sestpr = s2;
	    } else {
		sp = 0.;
		cp = 1.;
		*sestpr = s1;
	    }
	} else if (absest <= eps * absalp || absest <= eps * absgam) {
	    s1 = absgam;
	    s2 = absalp;
	    if (s1 <= s2) {
		tmp = s1 / s2;
		sp = sqrt(tmp * tmp + 1.);
		*sestpr = s2 * sp;
		cp = *gamma / s2 / sp;
		sp=(alpha>=0.0? 1.0 : -1.0)/ sp;
		/* *s = d_sign(&c_b263, &alpha) / *s; */
	    } else {
		tmp = s2 / s1;
		cp = sqrt(tmp * tmp + 1.);
		*sestpr = s1 * cp;
		sp = alpha / s1 / cp;
		cp = (*gamma>=0.0 ? 1.0 : -1.0) / cp;
		/*	*c__ = d_sign(&c_b263, gamma) / *c__; */
	    }
	} else {

/*           normal case */

	    zeta1 = alpha / absest;
	    zeta2 = *gamma / absest;

	    b = (1. - zeta1 * zeta1 - zeta2 * zeta2) * .5;
	    cp = zeta1 * zeta1;
	    if (b > 0.) {
		t = cp / (b + sqrt(b * b + cp));
	    } else {
		t = sqrt(b * b + cp) - b;
	    }

	    sine = -zeta1 / t;
	    cosine = -zeta2 / (t + 1.);
	    tmp = sqrt(sine * sine + cosine * cosine);
	    sp = sine / tmp;
	    cp = cosine / tmp;
	    *sestpr = sqrt(t + 1.) * absest;
	}

    } else if (*job == 2) {

/*        Estimating smallest singular value */

/*        special cases */

	if (absest == 0.) {
	    *sestpr = 0.;
	    if (Dymola_max(absgam,absalp) == 0.) {
		sine = 1.;
		cosine = 0.;
	    } else {
		sine = -(*gamma);
		cosine = alpha;
	    }
/* Computing MAX */
	    d__1 = Dymola_abs(sine), d__2 = Dymola_abs(cosine);
	    s1 = Dymola_max(d__1,d__2);
	    sp = sine / s1;
	    cp = cosine / s1;
	    tmp = sqrt(sp * sp + cp * cp);
	    sp /= tmp;
	    cp /= tmp;
	} else if (absgam <= eps * absest) {
	    sp = 0.;
	    cp = 1.;
	    *sestpr = absgam;
	} else if (absalp <= eps * absest) {
	    s1 = absgam;
	    s2 = absest;
	    if (s1 <= s2) {
		sp = 0.;
		cp = 1.;
		*sestpr = s1;
	    } else {
		sp = 1.;
		cp = 0.;
		*sestpr = s2;
	    }
	} else if (absest <= eps * absalp || absest <= eps * absgam) {
	    s1 = absgam;
	    s2 = absalp;
	    if (s1 <= s2) {
		tmp = s1 / s2;
		cp = sqrt(tmp * tmp + 1.);
		*sestpr = absest * (tmp / cp);
		sp = -(*gamma / s2) / cp;
		cp = (alpha>=0.0 ? 1.0 : -1.0) / cp;
		/* *c__ = d_sign(&c_b263, &alpha) / *c__; */
	    } else {
		tmp = s2 / s1;
		sp = sqrt(tmp * tmp + 1.);
		*sestpr = absest / sp;
		cp = alpha / s1 / sp;
		sp = -(*gamma >= 0.0 ? 1.0 : -1.0) / sp;
		/* *s = -d_sign(&c_b263, gamma) / *s; */
	    }
	} else {

/*           normal case */

	    zeta1 = alpha / absest;
	    zeta2 = *gamma / absest;

/* Computing MAX */
	    d__3 = zeta1 * zeta1 + 1. + (d__1 = zeta1 * zeta2, Dymola_abs(d__1)), 
		    d__4 = (d__2 = zeta1 * zeta2, Dymola_abs(d__2)) + zeta2 * zeta2;
	    norma = Dymola_max(d__3,d__4);

/*           See if root is closer to zero or to ONE */

	    test = (zeta1 - zeta2) * 2. * (zeta1 + zeta2) + 1.;
	    if (test >= 0.) {

/*              root is close to zero, compute directly */

		b = (zeta1 * zeta1 + zeta2 * zeta2 + 1.) * .5;
		cp = zeta2 * zeta2;
		t = cp / (b + sqrt((d__1 = b * b - cp, Dymola_abs(d__1))));
		sine = zeta1 / (1. - t);
		cosine = -zeta2 / t;
		*sestpr = sqrt(t + eps * 4. * eps * norma) * absest;
	    } else {

/*              root is closer to ONE, shift by that amount */

		b = (zeta2 * zeta2 + zeta1 * zeta1 - 1.) * .5;
		cp = zeta1 * zeta1;
		if (b >= 0.) {
		    t = -(cp) / (b + sqrt(b * b + cp));
		} else {
		    t = b - sqrt(b * b + cp);
		}
		sine = -zeta1 / t;
		cosine = -zeta2 / (t + 1.);
		*sestpr = sqrt(t + 1. + eps * 4. * eps * norma) * absest;
	    }
	    tmp = sqrt(sine * sine + cosine * cosine);
	    sp = sine / tmp;
	    cp = cosine / tmp;

	}
    } else return 0;
	*s = sp;
	*c__ = cp;
    return 0;

/*     End of DLAIC1 */

} /* dlaic1_ */

#ifndef GODESS
DYMOLA_STATIC doublereal dlamch_(char*cmach, ftnlen cmach_len);
DYMOLA_STATIC doublereal dlamch_(cmach, cmach_len)
char *cmach;
ftnlen cmach_len;
{


#if 1
	/* Fast variant, no extra calls */
    switch (toupperC(*cmach)) {
    case 'E':
	return DBL_EPSILON;
    case 'S':
	return DBL_MIN;
    case 'B':
	return FLT_RADIX;
    case 'P':
	return DBL_EPSILON * FLT_RADIX;
    case 'N':
	return DBL_MANT_DIG;
    case 'R':
	return FLT_ROUNDS == 1;
    case 'M':
	return DBL_MIN_10_EXP;
    case 'U':
	return DBL_MIN;
    case 'L':
	return DBL_MAX_10_EXP;
    case 'O':
	return DBL_MAX;
    default:
	return 0;	/* error */
    }
#else
    /* Initialized data */

    static logical first = TRUE_;

    /* System generated locals */
    integer i__1;
    doublereal ret_val;

    /* Builtin functions */
    double pow_di();

    /* Local variables */
    static doublereal t;
    integer it;
    static doublereal rnd, eps, base;
    integer beta;
    static doublereal emin, prec, emax;
    integer imin, imax;
    logical lrnd;
    static doublereal rmin, rmax;
    doublereal rmach;
    doublereal small;
    static doublereal sfmin;
    extern /* Subroutine */ int dlamc2_();


/*  -- LAPACK auxiliary routine (version 1.1) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
/*     Courant Institute, Argonne National Lab, and Rice University */
/*     October 31, 1992 */

/*     .. Scalar Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  DLAMCH determines double precision machine parameters. */

/*  Arguments */
/*  ========= */

/*  CMACH   (input) CHARACTER*1 */
/*          Specifies the value to be returned by DLAMCH: */
/*          = 'E' or 'e',   DLAMCH := eps */
/*          = 'S' or 's ,   DLAMCH := sfmin */
/*          = 'B' or 'b',   DLAMCH := base */
/*          = 'P' or 'p',   DLAMCH := eps*base */
/*          = 'N' or 'n',   DLAMCH := t */
/*          = 'R' or 'r',   DLAMCH := rnd */
/*          = 'M' or 'm',   DLAMCH := emin */
/*          = 'U' or 'u',   DLAMCH := rmin */
/*          = 'L' or 'l',   DLAMCH := emax */
/*          = 'O' or 'o',   DLAMCH := rmax */

/*          where */

/*          eps   = relative machine precision */
/*          sfmin = safe minimum, such that 1/sfmin does not overflow */
/*          base  = base of the machine */
/*          prec  = eps*base */
/*          t     = number of (base) digits in the mantissa */
/*          rnd   = 1.0 when rounding occurs in addition, 0.0 otherwise */
/*          emin  = minimum exponent before (gradual) underflow */
/*          rmin  = underflow threshold - base**(emin-1) */
/*          emax  = largest exponent before overflow */
/*          rmax  = overflow threshold  - (base**emax)*(1-eps) */

/* ===================================================================== */

/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Save statement .. */
/*     .. */
/*     .. Data statements .. */
/*     .. */
/*     .. Executable Statements .. */
    if (first) {
	first = FALSE_;
	dlamc2_(&beta, &it, &lrnd, &eps, &imin, &rmin, &imax, &rmax);
	base = (doublereal) beta;
	t = (doublereal) it;
	if (lrnd) {
	    rnd = 1.;
	    i__1 = 1 - it;
	    eps = pow_di(&base, &i__1) / 2;
	} else {
	    rnd = 0.;
	    i__1 = 1 - it;
	    eps = pow_di(&base, &i__1);
	}
	prec = eps * base;
	emin = (doublereal) imin;
	emax = (doublereal) imax;
	sfmin = rmin;
	small = 1. / rmax;
	if (small >= sfmin) {

/*           Use SMALL plus a bit, to avoid the possibility of rou
nding */
/*           causing overflow when computing  1/sfmin. */

	    sfmin = small * (eps + 1.);
	}
    }
#if 1
    switch (toupperC(*cmach)) {
    case 'E':return eps;
    case 'S':return sfmin;
    case 'B':return base;
    case 'N':return t;
    case 'R':return rnd;
    case 'M':return emin;
    case 'U':return rmin;
    case 'L':return emax;
    case 'O':return rmax;
    default:break;
    };
#else
    if (lsame_(cmach, "E", 1, 1)) {
	rmach = eps;
    } else if (lsame_(cmach, "S", 1, 1)) {
	rmach = sfmin;
    } else if (lsame_(cmach, "B", 1, 1)) {
	rmach = base;
    } else if (lsame_(cmach, "P", 1, 1)) {
	rmach = prec;
    } else if (lsame_(cmach, "N", 1, 1)) {
	rmach = t;
    } else if (lsame_(cmach, "R", 1, 1)) {
	rmach = rnd;
    } else if (lsame_(cmach, "M", 1, 1)) {
	rmach = emin;
    } else if (lsame_(cmach, "U", 1, 1)) {
	rmach = rmin;
    } else if (lsame_(cmach, "L", 1, 1)) {
	rmach = emax;
    } else if (lsame_(cmach, "O", 1, 1)) {
	rmach = rmax;
    }
#endif
    ret_val = rmach;
    return ret_val;

/*     End of DLAMCH */
#endif
} /* dlamch_ */


/* *********************************************************************** */
#if 0
/* Subroutine */ DYMOLA_STATIC int dlamc1_(beta, t, rnd, ieee1)
integer *beta, *t;
logical *rnd, *ieee1;
{
    /* Initialized data */

    static logical first = TRUE_;

    /* System generated locals */
    doublereal d__1, d__2;

    /* Local variables */
    doublereal a, b, c__, f, t1, t2;
    static integer lt;
    doublereal one, qtr;
    static logical lrnd;
    static integer lbeta;
    doublereal savec;
    extern doublereal dlamc3_();
    static logical lieee1;


/*  -- LAPACK auxiliary routine (version 1.1) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
/*     Courant Institute, Argonne National Lab, and Rice University */
/*     October 31, 1992 */

/*     .. Scalar Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  DLAMC1 determines the machine parameters given by BETA, T, RND, and */
/*  IEEE1. */

/*  Arguments */
/*  ========= */

/*  BETA    (output) INTEGER */
/*          The base of the machine. */

/*  T       (output) INTEGER */
/*          The number of ( BETA ) digits in the mantissa. */

/*  RND     (output) LOGICAL */
/*          Specifies whether proper rounding  ( RND = .TRUE. )  or */
/*          chopping  ( RND = .FALSE. )  occurs in addition. This may not */
/*          be a reliable guide to the way in which the machine performs */
/*          its arithmetic. */

/*  IEEE1   (output) LOGICAL */
/*          Specifies whether rounding appears to be done in the IEEE */
/*          'round to nearest' style. */

/*  Further Details */
/*  =============== */

/*  The routine is based on the routine  ENVRON  by Malcolm and */
/*  incorporates suggestions by Gentleman and Marovich. See */

/*     Malcolm M. A. (1972) Algorithms to reveal properties of */
/*        floating-point arithmetic. Comms. of the ACM, 15, 949-951. */

/*     Gentleman W. M. and Marovich S. B. (1974) More on algorithms */
/*        that reveal properties of floating point arithmetic units. */
/*        Comms. of the ACM, 17, 276-277. */

/* ===================================================================== */

/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. Save statement .. */
/*     .. */
/*     .. Data statements .. */
/*     .. */
/*     .. Executable Statements .. */

    if (first) {
	first = FALSE_;
	one = 1.;

/*        LBETA,  LIEEE1,  LT and  LRND  are the  local values  of  BETA, */
/*        IEEE1, T and RND. */

/*        Throughout this routine  we use the function  DLAMC3  to ensure */
/*        that relevant values are  stored and not held in registers,  or */
/*        are not affected by optimizers. */

/*        Compute  a = 2.0**m  with the  smallest positive integer m such */
/*        that */

/*           fl( a + 1.0 ) = a. */

	a = 1.;
	c__ = 1.;

/* +       WHILE( C.EQ.ONE )LOOP */
L10:
	if (c__ == one) {
	    a *= 2;
	    c__ = dlamc3_(&a, &one);
	    d__1 = -a;
	    c__ = dlamc3_(&c__, &d__1);
	    goto L10;
	}
/* +       END WHILE */

/*        Now compute  b = 2.0**m  with the smallest positive integer m */
/*        such that */

/*           fl( a + b ) .gt. a. */

	b = 1.;
	c__ = dlamc3_(&a, &b);

/* +       WHILE( C.EQ.A )LOOP */
L20:
	if (c__ == a) {
	    b *= 2;
	    c__ = dlamc3_(&a, &b);
	    goto L20;
	}
/* +       END WHILE */

/*        Now compute the base.  a and c  are neighbouring floating point */
/*        numbers  in the  interval  ( beta**t, beta**( t + 1 ) )  and so */
/*        their difference is beta. Adding 0.25 to c is to ensure that it */
/*        is truncated to beta and not ( beta - 1 ). */

	qtr = one / 4;
	savec = c__;
	d__1 = -a;
	c__ = dlamc3_(&c__, &d__1);
	lbeta = (integer) (c__ + qtr);

/*        Now determine whether rounding or chopping occurs,  by adding a */
/*        bit  less  than  beta/2  and a  bit  more  than  beta/2  to  a. */

	b = (doublereal) lbeta;
	d__1 = b / 2;
	d__2 = -b / 100;
	f = dlamc3_(&d__1, &d__2);
	c__ = dlamc3_(&f, &a);
	if (c__ == a) {
	    lrnd = TRUE_;
	} else {
	    lrnd = FALSE_;
	}
	d__1 = b / 2;
	d__2 = b / 100;
	f = dlamc3_(&d__1, &d__2);
	c__ = dlamc3_(&f, &a);
	if (lrnd && c__ == a) {
	    lrnd = FALSE_;
	}

/*        Try and decide whether rounding is done in the  IEEE  'round to */
/*        nearest' style. B/2 is half a unit in the last place of the two */
/*        numbers A and SAVEC. Furthermore, A is even, i.e. has last  bit */
/*        zero, and SAVEC is odd. Thus adding B/2 to A should not  change */
/*        A, but adding B/2 to SAVEC should change SAVEC. */

	d__1 = b / 2;
	t1 = dlamc3_(&d__1, &a);
	d__1 = b / 2;
	t2 = dlamc3_(&d__1, &savec);
	lieee1 = t1 == a && t2 > savec && lrnd;

/*        Now find  the  mantissa, t.  It should  be the  integer part of */
/*        log to the base beta of a,  however it is safer to determine  t */
/*        by powering.  So we find t as the smallest positive integer for */
/*        which */

/*           fl( beta**t + 1.0 ) = 1.0. */

	lt = 0;
	a = 1.;
	c__ = 1.;

/* +       WHILE( C.EQ.ONE )LOOP */
L30:
	if (c__ == one) {
	    ++lt;
	    a *= lbeta;
	    c__ = dlamc3_(&a, &one);
	    d__1 = -a;
	    c__ = dlamc3_(&c__, &d__1);
	    goto L30;
	}
/* +       END WHILE */

    }

    *beta = lbeta;
    *t = lt;
    *rnd = lrnd;
    *ieee1 = lieee1;
    return 0;

/*     End of DLAMC1 */

} /* dlamc1_ */


/* *********************************************************************** */

/* Subroutine */ DYMOLA_STATIC int dlamc2_(beta, t, rnd, eps, emin, rmin, emax, rmax)
integer *beta, *t;
logical *rnd;
doublereal *eps;
integer *emin;
doublereal *rmin;
integer *emax;
doublereal *rmax;
{
    /* Initialized data */

    static logical first = TRUE_;
    static logical iwarn = FALSE_;

    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2, d__3, d__4, d__5;

    /* Builtin functions */
    double pow_di();

    /* Local variables */
    doublereal a, b, c__;
    integer i__;
    static integer lt;
    doublereal one, two;
    logical ieee;
    doublereal half;
    logical lrnd;
    static doublereal leps;
    doublereal zero;
    static integer lbeta;
    doublereal rbase;
    static integer lemin, lemax;
    integer gnmin;
    doublereal small;
    integer gpmin;
    doublereal third;
    static doublereal lrmin, lrmax;
    doublereal sixth;
    extern /* Subroutine */ int dlamc1_();
    extern doublereal dlamc3_();
    logical lieee1;
    extern /* Subroutine */ int dlamc4_(), dlamc5_();
    integer ngnmin, ngpmin;


/*  -- LAPACK auxiliary routine (version 1.1) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
/*     Courant Institute, Argonne National Lab, and Rice University */
/*     October 31, 1992 */

/*     .. Scalar Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  DLAMC2 determines the machine parameters specified in its argument */
/*  list. */

/*  Arguments */
/*  ========= */

/*  BETA    (output) INTEGER */
/*          The base of the machine. */

/*  T       (output) INTEGER */
/*          The number of ( BETA ) digits in the mantissa. */

/*  RND     (output) LOGICAL */
/*          Specifies whether proper rounding  ( RND = .TRUE. )  or */
/*          chopping  ( RND = .FALSE. )  occurs in addition. This may not */
/*          be a reliable guide to the way in which the machine performs */
/*          its arithmetic. */

/*  EPS     (output) DOUBLE PRECISION */
/*          The smallest positive number such that */

/*             fl( 1.0 - EPS ) .LT. 1.0, */

/*          where fl denotes the computed value. */

/*  EMIN    (output) INTEGER */
/*          The minimum exponent before (gradual) underflow occurs. */

/*  RMIN    (output) DOUBLE PRECISION */
/*          The smallest normalized number for the machine, given by */
/*          BASE**( EMIN - 1 ), where  BASE  is the floating point value */
/*          of BETA. */

/*  EMAX    (output) INTEGER */
/*          The maximum exponent before overflow occurs. */

/*  RMAX    (output) DOUBLE PRECISION */
/*          The largest positive number for the machine, given by */
/*          BASE**EMAX * ( 1 - EPS ), where  BASE  is the floating point */
/*          value of BETA. */

/*  Further Details */
/*  =============== */

/*  The computation of  EPS  is based on a routine PARANOIA by */
/*  W. Kahan of the University of California at Berkeley. */

/* ===================================================================== */

/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Save statement .. */
/*     .. */
/*     .. Data statements .. */
/*     .. */
/*     .. Executable Statements .. */

    if (first) {
	first = FALSE_;
	zero = 0.;
	one = 1.;
	two = 2.;

/*        LBETA, LT, LRND, LEPS, LEMIN and LRMIN  are the local values of */
/*        BETA, T, RND, EPS, EMIN and RMIN. */

/*        Throughout this routine  we use the function  DLAMC3  to ensure */
/*        that relevant values are stored  and not held in registers,  or */
/*        are not affected by optimizers. */

/*        DLAMC1 returns the parameters  LBETA, LT, LRND and LIEEE1. */

	dlamc1_(&lbeta, &lt, &lrnd, &lieee1);

/*        Start to find EPS. */

	b = (doublereal) lbeta;
	i__1 = -lt;
	a = pow_di(&b, &i__1);
	leps = a;

/*        Try some tricks to see whether or not this is the correct  EPS. */

	b = two / 3;
	half = one / 2;
	d__1 = -half;
	sixth = dlamc3_(&b, &d__1);
	third = dlamc3_(&sixth, &sixth);
	d__1 = -half;
	b = dlamc3_(&third, &d__1);
	b = dlamc3_(&b, &sixth);
	b = Dymola_abs(b);
	if (b < leps) {
	    b = leps;
	}

	leps = 1.;

/* +       WHILE( ( LEPS.GT.B ).AND.( B.GT.ZERO ) )LOOP */
L10:
	if (leps > b && b > zero) {
	    leps = b;
	    d__1 = half * leps;
/* Computing 5th power */
	    d__3 = two, d__4 = d__3, d__3 *= d__3;
/* Computing 2nd power */
	    d__5 = leps;
	    d__2 = d__4 * (d__3 * d__3) * (d__5 * d__5);
	    c__ = dlamc3_(&d__1, &d__2);
	    d__1 = -c__;
	    c__ = dlamc3_(&half, &d__1);
	    b = dlamc3_(&half, &c__);
	    d__1 = -b;
	    c__ = dlamc3_(&half, &d__1);
	    b = dlamc3_(&half, &c__);
	    goto L10;
	}
/* +       END WHILE */

	if (a < leps) {
	    leps = a;
	}

/*        Computation of EPS complete. */

/*        Now find  EMIN.  Let A = + or - 1, and + or - (1 + BASE**(-3
)). */
/*        Keep dividing  A by BETA until (gradual) underflow occurs. T
his */
/*        is detected when we cannot recover the previous A. */

	rbase = one / lbeta;
	small = one;
	for (i__ = 1; i__ <= 3; ++i__) {
	    d__1 = small * rbase;
	    small = dlamc3_(&d__1, &zero);
/* L20: */
	}
	a = dlamc3_(&one, &small);
	dlamc4_(&ngpmin, &one, &lbeta);
	d__1 = -one;
	dlamc4_(&ngnmin, &d__1, &lbeta);
	dlamc4_(&gpmin, &a, &lbeta);
	d__1 = -a;
	dlamc4_(&gnmin, &d__1, &lbeta);
	ieee = FALSE_;

	if (ngpmin == ngnmin && gpmin == gnmin) {
	    if (ngpmin == gpmin) {
		lemin = ngpmin;
/*            ( Non twos-complement machines, no gradual under
flow; */
/*              e.g.,  VAX ) */
	    } else if (gpmin - ngpmin == 3) {
		lemin = ngpmin - 1 + lt;
		ieee = TRUE_;
/*            ( Non twos-complement machines, with gradual und
erflow; */
/*              e.g., IEEE standard followers ) */
	    } else {
		lemin = Dymola_min(ngpmin,gpmin);
/*            ( A guess; no known machine ) */
		iwarn = TRUE_;
	    }

	} else if (ngpmin == gpmin && ngnmin == gnmin) {
	    if ((i__1 = ngpmin - ngnmin, Dymola_abs(i__1)) == 1) {
		lemin = Dymola_max(ngpmin,ngnmin);
/*            ( Twos-complement machines, no gradual underflow
; */
/*              e.g., CYBER 205 ) */
	    } else {
		lemin = Dymola_min(ngpmin,ngnmin);
/*            ( A guess; no known machine ) */
		iwarn = TRUE_;
	    }

	} else if ((i__1 = ngpmin - ngnmin, Dymola_abs(i__1)) == 1 && gpmin == gnmin)
		 {
	    if (gpmin - Dymola_min(ngpmin,ngnmin) == 3) {
		lemin = Dymola_max(ngpmin,ngnmin) - 1 + lt;
/*            ( Twos-complement machines with gradual underflo
w; */
/*              no known machine ) */
	    } else {
		lemin = Dymola_min(ngpmin,ngnmin);
/*            ( A guess; no known machine ) */
		iwarn = TRUE_;
	    }

	} else {
/* Computing MIN */
	    i__1 = Dymola_min(ngpmin,ngnmin), i__1 = Dymola_min(i__1,gpmin);
	    lemin = Dymola_min(i__1,gnmin);
/*         ( A guess; no known machine ) */
	    iwarn = TRUE_;
	}
/* ** */
/* Comment out this if block if EMIN is ok */
	if (iwarn) {
	    first = TRUE_;
/* CC Otter            WRITE( 6, FMT = 9999 )LEMIN */
	}
/* ** */

/*        Assume IEEE arithmetic if we found denormalised  numbers abo
ve, */
/*        or if arithmetic seems to round in the  IEEE style,  determi
ned */
/*        in routine DLAMC1. A true IEEE machine should have both  thi
ngs */
/*        true; however, faulty machines may have one or the other. */

	ieee = ieee || lieee1;

/*        Compute  RMIN by successive division by  BETA. We could comp
ute */
/*        RMIN as BASE**( EMIN - 1 ),  but some machines underflow dur
ing */
/*        this computation. */

	lrmin = 1.;
	i__1 = 1 - lemin;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    d__1 = lrmin * rbase;
	    lrmin = dlamc3_(&d__1, &zero);
/* L30: */
	}

/*        Finally, call DLAMC5 to compute EMAX and RMAX. */

	dlamc5_(&lbeta, &lt, &lemin, &ieee, &lemax, &lrmax);
    }

    *beta = lbeta;
    *t = lt;
    *rnd = lrnd;
    *eps = leps;
    *emin = lemin;
    *rmin = lrmin;
    *emax = lemax;
    *rmax = lrmax;

    return 0;

/*CC Otter 9999 FORMAT( / / ' WARNING. The value EMIN may be incorrect:-',
*/
/* CC Otter     $      '  EMIN = ', I8, / */
/* CC Otter     $      ' If, after inspection, the value EMIN looks', */
/* CC Otter     $      ' acceptable please comment out ', */
/*CC Otter     $      / ' the IF block as marked within the code of routin
e',*/
/*CC Otter     $      ' DLAMC2,', / ' otherwise supply EMIN explicitly.', 
/ )*/

/*     End of DLAMC2 */

} /* dlamc2_ */

/* *********************************************************************** */

DYMOLA_STATIC doublereal dlamc3_(a, b)
doublereal *a, *b;
{
    /* System generated locals */
    doublereal ret_val;


/*  -- LAPACK auxiliary routine (version 1.1) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
/*     Courant Institute, Argonne National Lab, and Rice University */
/*     October 31, 1992 */

/*     .. Scalar Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  DLAMC3  is intended to force  A  and  B  to be stored prior to doing 
*/
/*  the addition of  A  and  B ,  for use in situations where optimizers 
*/
/*  might hold one of these in a register. */

/*  Arguments */
/*  ========= */

/*  A, B    (input) DOUBLE PRECISION */
/*          The values A and B. */

/* ===================================================================== 
*/

/*     .. Executable Statements .. */

    ret_val = *a + *b;

    return ret_val;

/*     End of DLAMC3 */

} /* dlamc3_ */

/* *********************************************************************** */

/* Subroutine */ DYMOLA_STATIC int dlamc4_(emin, start, base)
integer *emin;
doublereal *start;
integer *base;
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Local variables */
    doublereal a;
    integer i__;
    doublereal b1, b2, c1, c2, d1, d2, one, zero, rbase;
    extern doublereal dlamc3_();


/*  -- LAPACK auxiliary routine (version 1.1) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
/*     Courant Institute, Argonne National Lab, and Rice University */
/*     October 31, 1992 */

/*     .. Scalar Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  DLAMC4 is a service routine for DLAMC2. */

/*  Arguments */
/*  ========= */

/*  EMIN    (output) EMIN */
/*          The minimum exponent before (gradual) underflow, computed by 
*/
/*          setting A = START and dividing by BASE until the previous A */
/*          can not be recovered. */

/*  START   (input) DOUBLE PRECISION */
/*          The starting point for determining EMIN. */

/*  BASE    (input) INTEGER */
/*          The base of the machine. */

/* ===================================================================== 
*/

/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. Executable Statements .. */

    a = *start;
    one = 1.;
    rbase = one / *base;
    zero = 0.;
    *emin = 1;
    d__1 = a * rbase;
    b1 = dlamc3_(&d__1, &zero);
    c1 = a;
    c2 = a;
    d1 = a;
    d2 = a;
/* +    WHILE( ( C1.EQ.A ).AND.( C2.EQ.A ).AND. */
/*    $       ( D1.EQ.A ).AND.( D2.EQ.A )      )LOOP */
L10:
    if (c1 == a && c2 == a && d1 == a && d2 == a) {
	--(*emin);
	a = b1;
	d__1 = a / *base;
	b1 = dlamc3_(&d__1, &zero);
	d__1 = b1 * *base;
	c1 = dlamc3_(&d__1, &zero);
	d1 = zero;
	i__1 = *base;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    d1 += b1;
/* L20: */
	}
	d__1 = a * rbase;
	b2 = dlamc3_(&d__1, &zero);
	d__1 = b2 / rbase;
	c2 = dlamc3_(&d__1, &zero);
	d2 = zero;
	i__1 = *base;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    d2 += b2;
/* L30: */
	}
	goto L10;
    }
/* +    END WHILE */

    return 0;

/*     End of DLAMC4 */

} /* dlamc4_ */


/* *********************************************************************** */

/* Subroutine */ DYMOLA_STATIC int dlamc5_(beta, p, emin, ieee, emax, rmax)
integer *beta, *p, *emin;
logical *ieee;
integer *emax;
doublereal *rmax;
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Local variables */
    integer i__;
    doublereal y, z__;
    integer try__, lexp;
    doublereal oldy;
    integer uexp, nbits;
    extern doublereal dlamc3_();
    doublereal recbas;
    integer exbits, expsum;


/*  -- LAPACK auxiliary routine (version 1.1) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
/*     Courant Institute, Argonne National Lab, and Rice University */
/*     October 31, 1992 */

/*     .. Scalar Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  DLAMC5 attempts to compute RMAX, the largest machine floating-point */
/*  number, without overflow.  It assumes that EMAX + Dymola_abs(EMIN) sum */
/*  approximately to a power of 2.  It will fail on machines where this */
/*  assumption does not hold, for example, the Cyber 205 (EMIN = -28625, 
*/
/*  EMAX = 28718).  It will also fail if the value supplied for EMIN is */
/*  too large (i.e. too close to zero), probably with overflow. */

/*  Arguments */
/*  ========= */

/*  BETA    (input) INTEGER */
/*          The base of floating-point arithmetic. */

/*  P       (input) INTEGER */
/*          The number of base BETA digits in the mantissa of a */
/*          floating-point value. */

/*  EMIN    (input) INTEGER */
/*          The minimum exponent before (gradual) underflow. */

/*  IEEE    (input) LOGICAL */
/*          A logical flag specifying whether or not the arithmetic */
/*          system is thought to comply with the IEEE standard. */

/*  EMAX    (output) INTEGER */
/*          The largest exponent before overflow */

/*  RMAX    (output) DOUBLE PRECISION */
/*          The largest machine floating-point number. */

/* ===================================================================== 
*/

/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     First compute LEXP and UEXP, two powers of 2 that bound */
/*     Dymola_abs(EMIN). We then assume that EMAX + Dymola_abs(EMIN) will sum */
/*     approximately to the bound that is closest to Dymola_abs(EMIN). */
/*     (EMAX is the exponent of the required number RMAX). */

    oldy = 0;
    lexp = 1;
    exbits = 1;
L10:
    try__ = lexp << 1;
    if (try__ <= -(*emin)) {
	lexp = try__;
	++exbits;
	goto L10;
    }
    if (lexp == -(*emin)) {
	uexp = lexp;
    } else {
	uexp = try__;
	++exbits;
    }

/*     Now -LEXP is less than or equal to EMIN, and -UEXP is greater */
/*     than or equal to EMIN. EXBITS is the number of bits needed to */
/*     store the exponent. */

    if (uexp + *emin > -lexp - *emin) {
	expsum = lexp << 1;
    } else {
	expsum = uexp << 1;
    }

/*     EXPSUM is the exponent range, approximately equal to */
/*     EMAX - EMIN + 1 . */

    *emax = expsum + *emin - 1;
    nbits = exbits + 1 + *p;

/*     NBITS is the total number of bits needed to store a */
/*     floating-point number. */

    if (nbits % 2 == 1 && *beta == 2) {

/*        Either there are an odd number of bits used to store a */
/*        floating-point number, which is unlikely, or some bits are 
*/
/*        not used in the representation of numbers, which is possible
, */
/*        (e.g. Cray machines) or the mantissa has an implicit bit, */
/*        (e.g. IEEE machines, Dec Vax machines), which is perhaps the
 */
/*        most likely. We have to assume the last alternative. */
/*        If this is true, then we need to reduce EMAX by one because 
*/
/*        there must be some way of representing zero in an implicit-b
it */
/*        system. On machines like Cray, we are reducing EMAX by one 
*/
/*        unnecessarily. */

	--(*emax);
    }

    if (*ieee) {

/*        Assume we are on an IEEE machine which reserves one exponent
 */
/*        for infinity and NaN. */

	--(*emax);
    }

/*     Now create RMAX, the largest machine number, which should */
/*     be equal to (1.0 - BETA**(-P)) * BETA**EMAX . */

/*     First compute 1.0 - BETA**(-P), being careful that the */
/*     result is less than 1.0 . */

    recbas = 1. / *beta;
    z__ = *beta - 1.;
    y = 0.;
    i__1 = *p;
    for (i__ = 1; i__ <= i__1; ++i__) {
	z__ *= recbas;
	if (y < 1.) {
	    oldy = y;
	}
	y = dlamc3_(&y, &z__);
/* L20: */
    }
    if (y >= 1.) {
	y = oldy;
    }

/*     Now multiply by BETA**EMAX to get RMAX. */

    i__1 = *emax;
    for (i__ = 1; i__ <= i__1; ++i__) {
	d__1 = y * *beta;
	y = dlamc3_(&d__1, &c_b362);
/* L30: */
    }

    *rmax = y;
    return 0;

/*     End of DLAMC5 */

} /* dlamc5_ */
#endif

DYMOLA_STATIC doublereal dlange_(norm, m, n, a, lda, work, norm_len)
char *norm;
const integer*m, *n;
doublereal *a;
const integer*lda;
doublereal *work;
ftnlen norm_len;
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;
    doublereal ret_val, d__1, d__2, d__3;


    /* Local variables */
    integer i__, j;
    doublereal sum, scale;
    doublereal value;


/*  -- LAPACK auxiliary routine (version 1.1) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
/*     Courant Institute, Argonne National Lab, and Rice University */
/*     October 31, 1992 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  DLANGE  returns the value of the one norm,  or the Frobenius norm, or 
*/
/*  the  infinity norm,  or the  element of  largest absolute value  of a 
*/
/*  real matrix A. */

/*  Description */
/*  =========== */

/*  DLANGE returns the value */

/*     DLANGE = ( max(Dymola_abs(A(i,j))), NORM = 'M' or 'm' */
/*              ( */
/*              ( norm1(A),         NORM = '1', 'O' or 'o' */
/*              ( */
/*              ( normI(A),         NORM = 'I' or 'i' */
/*              ( */
/*              ( normF(A),         NORM = 'F', 'f', 'E' or 'e' */

/*  where  norm1  denotes the  one norm of a matrix (maximum column sum), 
*/
/*  normI  denotes the  infinity norm  of a matrix  (maximum row sum) and 
*/
/*  normF  denotes the  Frobenius norm of a matrix (square root of sum of 
*/
/*  squares).  Note that  max(Dymola_abs(A(i,j)))  is not a  matrix norm. */

/*  Arguments */
/*  ========= */

/*  NORM    (input) CHARACTER*1 */
/*          Specifies the value to be returned in DLANGE as described */
/*          above. */

/*  M       (input) INTEGER */
/*          The number of rows of the matrix A.  M >= 0.  When M = 0, */
/*          DLANGE is set to zero. */

/*  N       (input) INTEGER */
/*          The number of columns of the matrix A.  N >= 0.  When N = 0, 
*/
/*          DLANGE is set to zero. */

/*  A       (input) DOUBLE PRECISION array, dimension (LDA,N) */
/*          The m by n matrix A. */

/*  LDA     (input) INTEGER */
/*          The leading dimension of the array A.  LDA >= max(M,1). */

/*  WORK    (workspace) DOUBLE PRECISION array, dimension (LWORK), */
/*          where LWORK >= M when NORM = 'I'; otherwise, WORK is not */
/*          referenced. */

/* ===================================================================== 
*/

/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = a_dim1 + 1;
    a -= a_offset;
    --work;

    /* Function Body */
    value = 0;
    if (Dymola_min(*m,*n) == 0) {
	value = 0.;
    } else if (lsame_(norm, "M", 1, 1)) {

/*        Find max(Dymola_abs(A(i,j))). */

	value = 0.;
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    i__2 = *m;
	    for (i__ = 1; i__ <= i__2; ++i__) {
/* Computing MAX */
		d__2 = value, d__3 = (d__1 = a[i__ + j * a_dim1], Dymola_abs(d__1));
		value = Dymola_max(d__2,d__3);
/* L10: */
	    }
/* L20: */
	}
    } else if (lsame_(norm, "O", 1, 1) || *(unsigned char *)norm == '1') {

/*        Find norm1(A). */

	value = 0.;
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    sum = 0.;
	    i__2 = *m;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		sum += (d__1 = a[i__ + j * a_dim1], Dymola_abs(d__1));
/* L30: */
	    }
	    value = Dymola_max(value,sum);
/* L40: */
	}
    } else if (lsame_(norm, "I", 1, 1)) {

/*        Find normI(A). */

	i__1 = *m;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    work[i__] = 0.;
/* L50: */
	}
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    i__2 = *m;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		work[i__] += (d__1 = a[i__ + j * a_dim1], Dymola_abs(d__1));
/* L60: */
	    }
/* L70: */
	}
	value = 0.;
	i__1 = *m;
	for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MAX */
	    d__1 = value, d__2 = work[i__];
	    value = Dymola_max(d__1,d__2);
/* L80: */
	}
    } else if (lsame_(norm, "F", 1, 1) || lsame_(norm, "E", 1, 1)) {

/*        Find normF(A). */

	scale = 0.;
	sum = 1.;
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    dlassq_(m, &a[j * a_dim1 + 1], &c__1, &scale, &sum);
/* L90: */
	}
	value = scale * sqrt(sum);
    }

    ret_val = value;
    return ret_val;

/*     End of DLANGE */

} /* dlange_ */

DYMOLA_STATIC doublereal dlapy2_(x, y)
doublereal *x, *y;
{
    /* System generated locals */
    doublereal ret_val, d__1;


    /* Local variables */
    doublereal w, z__, xabs, yabs;


/*  -- LAPACK auxiliary routine (version 1.1) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
/*     Courant Institute, Argonne National Lab, and Rice University */
/*     October 31, 1992 */

/*     .. Scalar Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  DLAPY2 returns sqrt(x**2+y**2), taking care not to cause unnecessary 
*/
/*  overflow. */

/*  Arguments */
/*  ========= */

/*  X       (input) DOUBLE PRECISION */
/*  Y       (input) DOUBLE PRECISION */
/*          X and Y specify the values x and y. */

/*  ===================================================================== 
*/

/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

    xabs = Dymola_abs(*x);
    yabs = Dymola_abs(*y);
    w = Dymola_max(xabs,yabs);
    z__ = Dymola_min(xabs,yabs);
    if (z__ == 0.) {
	ret_val = w;
    } else {
/* Computing 2nd power */
	d__1 = z__ / w;
	ret_val = w * sqrt(d__1 * d__1 + 1.);
    }
    return ret_val;

/*     End of DLAPY2 */

} /* dlapy2_ */

/* Subroutine */ DYMOLA_STATIC int dlarf_(side, m, n, v, incv, tau, c__, ldc, work, 
	side_len)
const char*side;
const integer*m, *n;
doublereal *v;
const integer*incv;
doublereal *tau, *c__;
const integer*ldc;
doublereal *work;
ftnlen side_len;
{
    /* System generated locals */
    doublereal d__1;


/*  -- LAPACK auxiliary routine (version 1.1) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
/*     Courant Institute, Argonne National Lab, and Rice University */
/*     February 29, 1992 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  DLARF applies a real elementary reflector H to a real m by n matrix */
/*  C, from either the left or the right. H is represented in the form */

/*        H = I - tau * v * v' */

/*  where tau is a real scalar and v is a real vector. */

/*  If tau = 0, then H is taken to be the unit matrix. */

/*  Arguments */
/*  ========= */

/*  SIDE    (input) CHARACTER*1 */
/*          = 'L': form  H * C */
/*          = 'R': form  C * H */

/*  M       (input) INTEGER */
/*          The number of rows of the matrix C. */

/*  N       (input) INTEGER */
/*          The number of columns of the matrix C. */

/*  V       (input) DOUBLE PRECISION array, dimension */
/*                     (1 + (M-1)*Dymola_abs(INCV)) if SIDE = 'L' */
/*                  or (1 + (N-1)*Dymola_abs(INCV)) if SIDE = 'R' */
/*          The vector v in the representation of H. V is not used if */
/*          TAU = 0. */

/*  INCV    (input) INTEGER */
/*          The increment between elements of v. INCV <> 0. */

/*  TAU     (input) DOUBLE PRECISION */
/*          The value tau in the representation of H. */

/*  C       (input/output) DOUBLE PRECISION array, dimension (LDC,N) */
/*          On entry, the m by n matrix C. */
/*          On exit, C is overwritten by the matrix H * C if SIDE = 'L', 
*/
/*          or C * H if SIDE = 'R'. */

/*  LDC     (input) INTEGER */
/*          The leading dimension of the array C. LDC >= max(1,M). */

/*  WORK    (workspace) DOUBLE PRECISION array, dimension */
/*                         (N) if SIDE = 'L' */
/*                      or (M) if SIDE = 'R' */

/*  ===================================================================== 
*/

/*     .. Parameters .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Function Body */
	if (*tau != 0.) {
		if (lsame_(side, "L", 1, 1)) {
		
		/*        Form  H * C */
		
			/*           w := C' * v */
			
			dgemv_("Transpose", m, n, &c_b263, c__, ldc, v, 
				incv, &c_b362, work, &c__1, 9);
			
			/*           C := C - v * w' */
			
			d__1 = -(*tau);
			dger_(m, n, &d__1, v, incv, work, &c__1, c__, 
				ldc);
		} else {
		
			/*        Form  C * H */
			
			/*           w := C * v */
			
			dgemv_("No transpose", m, n, &c_b263, c__, ldc, v, 
				incv, &c_b362, work, &c__1, 12);
			
			/*           C := C - w * v' */
			
			d__1 = -(*tau);
			dger_(m, n, &d__1, work, &c__1, v, incv, c__, 
				ldc);
		}
	}
    return 0;

/*     End of DLARF */

} /* dlarf_ */

/* Subroutine */ DYMOLA_STATIC int dlarfg_(n, alpha, x, incx, tau)
const integer*n;
doublereal*alpha;
doublereal *x;
const integer*incx;
doublereal *tau;
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;


    /* Local variables */
    integer j, knt;
    doublereal beta;
    doublereal xnorm;
    doublereal safmin, rsafmn;


/*  -- LAPACK auxiliary routine (version 1.1) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
/*     Courant Institute, Argonne National Lab, and Rice University */
/*     February 29, 1992 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  DLARFG generates a real elementary reflector H of order n, such */
/*  that */

/*        H * ( alpha ) = ( beta ),   H' * H = I. */
/*            (   x   )   (   0  ) */

/*  where alpha and beta are scalars, and x is an (n-1)-element real */
/*  vector. H is represented in the form */

/*        H = I - tau * ( 1 ) * ( 1 v' ) , */
/*                      ( v ) */

/*  where tau is a real scalar and v is a real (n-1)-element */
/*  vector. */

/*  If the elements of x are all zero, then tau = 0 and H is taken to be 
*/
/*  the unit matrix. */

/*  Otherwise  1 <= tau <= 2. */

/*  Arguments */
/*  ========= */

/*  N       (input) INTEGER */
/*          The order of the elementary reflector. */

/*  ALPHA   (input/output) DOUBLE PRECISION */
/*          On entry, the value alpha. */
/*          On exit, it is overwritten with the value beta. */

/*  X       (input/output) DOUBLE PRECISION array, dimension */
/*                         (1+(N-2)*Dymola_abs(INCX)) */
/*          On entry, the vector x. */
/*          On exit, it is overwritten with the vector v. */

/*  INCX    (input) INTEGER */
/*          The increment between elements of X. INCX <> 0. */

/*  TAU     (output) DOUBLE PRECISION */
/*          The value tau. */

/*  ===================================================================== 
*/

/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Executable Statements .. */

    /* Parameter adjustments */
    --x;

    /* Function Body */
    if (*n <= 1) {
	*tau = 0.;
	return 0;
    }

    i__1 = *n - 1;
    xnorm = dnrm2_Fast(i__1, &x[1], *incx);

    if (xnorm == 0.) {

/*        H  =  I */

	*tau = 0.;
    } else {

/*        general case */

	d__1 = dlapy2_(alpha, &xnorm); /* Cannot be negative */
	beta = -(*alpha>=0.0 ? d__1 : -d__1);
	/* beta = -d_sign(&d__1, alpha); */
	safmin = DBL_MIN;
	if (Dymola_abs(beta) < safmin) {

/*           XNORM, BETA may be inaccurate; scale X and recompute 
them */

	    rsafmn = 1. / safmin;
	    knt = 0;
L10:
	    ++knt;
	    i__1 = *n - 1;
	    dscal_(&i__1, &rsafmn, &x[1], incx);
	    beta *= rsafmn;
	    *alpha *= rsafmn;
	    if (Dymola_abs(beta) < safmin) {
		goto L10;
	    }

/*           New BETA is at most 1, at least SAFMIN */

	    i__1 = *n - 1;
	    xnorm = dnrm2_Fast(i__1, &x[1], *incx);
	    d__1 = dlapy2_(alpha, &xnorm); /* cannot be negative */
		beta = - (*alpha>=0.0 ? d__1 : -d__1);
	    /* beta = -d_sign(&d__1, alpha); */
	    *tau = (beta - *alpha) / beta;
	    i__1 = *n - 1;
	    d__1 = 1. / (*alpha - beta);
	    dscal_(&i__1, &d__1, &x[1], incx);

/*           If ALPHA is subnormal, it may lose relative accuracy 
*/

	    *alpha = beta;
	    i__1 = knt;
	    for (j = 1; j <= i__1; ++j) {
		*alpha *= safmin;
/* L20: */
	    }
	} else {
	    *tau = (beta - *alpha) / beta;
	    i__1 = *n - 1;
	    d__1 = 1. / (*alpha - beta);
	    dscal_(&i__1, &d__1, &x[1], incx);
	    *alpha = beta;
	}
    }

    return 0;

/*     End of DLARFG */

} /* dlarfg_ */

/* Subroutine */ DYMOLA_STATIC int dlassq_(n, x, incx, scale, sumsq)
const integer*n;
doublereal *x;
const integer*incx;
doublereal *scale, *sumsq;
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1;

    /* Local variables */
    integer ix;
    doublereal absxi;


/*  -- LAPACK auxiliary routine (version 1.1) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
/*     Courant Institute, Argonne National Lab, and Rice University */
/*     October 31, 1992 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  DLASSQ  returns the values  scl  and  smsq  such that */

/*     ( scl**2 )*smsq = x( 1 )**2 +...+ x( n )**2 + ( scale**2 )*sumsq, 
*/

/*  where  x( i ) = X( 1 + ( i - 1 )*INCX ). The value of  sumsq  is */
/*  assumed to be non-negative and  scl  returns the value */

/*     scl = max( scale, Dymola_abs( x( i ) ) ). */

/*  scale and sumsq must be supplied in SCALE and SUMSQ and */
/*  scl and smsq are overwritten on SCALE and SUMSQ respectively. */

/*  The routine makes only one pass through the vector x. */

/*  Arguments */
/*  ========= */

/*  N       (input) INTEGER */
/*          The number of elements to be used from the vector X. */

/*  X       (input) DOUBLE PRECISION */
/*          The vector for which a scaled sum of squares is computed. */
/*             x( i )  = X( 1 + ( i - 1 )*INCX ), 1 <= i <= n. */

/*  INCX    (input) INTEGER */
/*          The increment between successive values of the vector X. */
/*          INCX > 0. */

/*  SCALE   (input/output) DOUBLE PRECISION */
/*          On entry, the value  scale  in the equation above. */
/*          On exit, SCALE is overwritten with  scl , the scaling factor 
*/
/*          for the sum of squares. */

/*  SUMSQ   (input/output) DOUBLE PRECISION */
/*          On entry, the value  sumsq  in the equation above. */
/*          On exit, SUMSQ is overwritten with  smsq , the basic sum of */
/*          squares from which  scl  has been factored out. */

/* ===================================================================== 
*/

/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

    /* Parameter adjustments */
    --x;

    /* Function Body */
    if (*n > 0) {
	i__1 = (*n - 1) * *incx + 1;
	i__2 = *incx;
	for (ix = 1; i__2 < 0 ? ix >= i__1 : ix <= i__1; ix += i__2) {
	    if (x[ix] != 0.) {
		absxi = (d__1 = x[ix], Dymola_abs(d__1));
		if (*scale < absxi) {
/* Computing 2nd power */
		    d__1 = *scale / absxi;
		    *sumsq = *sumsq * (d__1 * d__1) + 1;
		    *scale = absxi;
		} else {
/* Computing 2nd power */
		    d__1 = absxi / *scale;
		    *sumsq += d__1 * d__1;
		}
	    }
/* L10: */
	}
    }
    return 0;

/*     End of DLASSQ */

} /* dlassq_ */

/* Subroutine */ DYMOLA_STATIC int dlaswp_(n, a, lda, k1, k2, ipiv, incx)
const integer*n;
doublereal *a;
const integer*lda, *k1, *k2;
const integer *ipiv, *incx;
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1;

    /* Local variables */
    integer i__, ip, ix;


/*  -- LAPACK auxiliary routine (version 1.1) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
/*     Courant Institute, Argonne National Lab, and Rice University */
/*     October 31, 1992 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  DLASWP performs a series of row interchanges on the matrix A. */
/*  One row interchange is initiated for each of rows K1 through K2 of A. 
*/

/*  Arguments */
/*  ========= */

/*  N       (input) INTEGER */
/*          The number of columns of the matrix A. */

/*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*          On entry, the matrix of column dimension N to which the row */
/*          interchanges will be applied. */
/*          On exit, the permuted matrix. */

/*  LDA     (input) INTEGER */
/*          The leading dimension of the array A. */

/*  K1      (input) INTEGER */
/*          The first element of IPIV for which a row interchange will */
/*          be done. */

/*  K2      (input) INTEGER */
/*          The last element of IPIV for which a row interchange will */
/*          be done. */

/*  IPIV    (input) INTEGER array, dimension (M*Dymola_abs(INCX)) */
/*          The vector of pivot indices.  Only the elements in positions 
*/
/*          K1 through K2 of IPIV are accessed. */
/*          IPIV(K) = L implies rows K and L are to be interchanged. */

/*  INCX    (input) INTEGER */
/*          The increment between successive values of IPIV.  If IPIV */
/*          is negative, the pivots are applied in reverse order. */

/* ===================================================================== 
*/

/*     .. Local Scalars .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Interchange row I with row IPIV(I) for each of rows K1 through K2. 
*/

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = a_dim1 + 1;
    a -= a_offset;
    --ipiv;

    /* Function Body */
    if (*incx == 0) {
	return 0;
    }
    if (*incx > 0) {
	ix = *k1;
    } else {
	ix = (1 - *k2) * *incx + 1;
    }
    if (*incx == 1) {
	i__1 = *k2;
	for (i__ = *k1; i__ <= i__1; ++i__) {
	    ip = ipiv[i__];
	    if (ip != i__) {
		dswap_(n, &a[i__ + a_dim1], lda, &a[ip + a_dim1], lda);
	    }
/* L10: */
	}
    } else if (*incx > 1) {
	i__1 = *k2;
	for (i__ = *k1; i__ <= i__1; ++i__) {
	    ip = ipiv[ix];
	    if (ip != i__) {
		dswap_(n, &a[i__ + a_dim1], lda, &a[ip + a_dim1], lda);
	    }
	    ix += *incx;
/* L20: */
	}
    } else if (*incx < 0) {
	i__1 = *k1;
	for (i__ = *k2; i__ >= i__1; --i__) {
	    ip = ipiv[ix];
	    if (ip != i__) {
		dswap_(n, &a[i__ + a_dim1], lda, &a[ip + a_dim1], lda);
	    }
	    ix += *incx;
/* L30: */
	}
    }

    return 0;

/*     End of DLASWP */

} /* dlaswp_ */


/* Subroutine */ DYMOLA_STATIC int dlatrs_(uplo, trans, diag, normin, n, a, lda, x, scale, 
	cnorm, info, uplo_len, trans_len, diag_len, normin_len)
const char*uplo, *trans, *diag, *normin;
const integer*n;
doublereal *a;
const integer*lda;
doublereal *x, *scale, *cnorm;
integer *info;
ftnlen uplo_len;
ftnlen trans_len;
ftnlen diag_len;
ftnlen normin_len;
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;
    doublereal d__1, d__2, d__3;

    /* Local variables */
    integer i__, j;
    doublereal xj, rec, tjj;
    integer jinc;
    doublereal xbnd;
    integer imax;
    doublereal tmax, tjjs, xmax, grow, sumj;
    doublereal tscal, uscal;
    integer jlast;
    logical upper;
    doublereal bignum;
    logical notran;
    integer jfirst;
    doublereal smlnum;
    logical nounit;


/*  -- LAPACK auxiliary routine (version 1.1) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
/*     Courant Institute, Argonne National Lab, and Rice University */
/*     June 30, 1992 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  DLATRS solves one of the triangular systems */
/*  Initially, let M(0) = 0 and G(0) = max{x(i), i=1,...,n}. */

/*  Then for iteration j+1 we have */
/*     M(j+1) <= G(j) / | A(j+1,j+1) | */
/*     G(j+1) <= G(j) + M(j+1) * | A[j+2:n,j+1] | */
/*            <= G(j) ( 1 + CNORM(j+1) / | A(j+1,j+1) | ) */

/*  where CNORM(j+1) is greater than or equal to the infinity-norm of */
/*  column j+1 of A, not counting the diagonal.  Hence */

/*     G(j) <= G(0) product ( 1 + CNORM(i) / | A(i,i) | ) */
/*                  1<=i<=j */
/*  and */

/*     |x(j)| <= ( G(0) / |A(j,j)| ) product ( 1 + CNORM(i) / |A(i,i)| ) 
*/
/*                                   1<=i< j */

/*  Since |x(j)| <= M(j), we use the Level 2 BLAS routine DTRSV if the */
/*  reciprocal of the largest M(j), j=1,..,n, is larger than */
/*  max(underflow, 1/overflow). */

/*  The bound on x(j) is also used to determine when a step in the */
/*  columnwise method can be performed without fear of overflow.  If */
/*  the computed bound is greater than a large constant, x is scaled to */
/*  prevent overflow, but if the bound overflows, x is set to 0, x(j) to 
*/
/*  1, and scale to 0, and a non-trivial solution to A*x = 0 is found. */

/*  Similarly, a row-wise scheme is used to solve A'*x = b.  The basic */
/*  algorithm for A upper triangular is */

/*       for j = 1, ..., n */
/*            x(j) := ( b(j) - A[1:j-1,j]' * x[1:j-1] ) / A(j,j) */
/*       end */

/*  We simultaneously compute two bounds */
/*       G(j) = bound on ( b(i) - A[1:i-1,i]' * x[1:i-1] ), 1<=i<=j */
/*       M(j) = bound on x(i), 1<=i<=j */

/*  The initial values are G(0) = 0, M(0) = max{b(i), i=1,..,n}, and we */
/*  add the constraint G(j) >= G(j-1) and M(j) >= M(j-1) for j >= 1. */
/*  Then the bound on x(j) is */

/*       M(j) <= M(j-1) * ( 1 + CNORM(j) ) / | A(j,j) | */

/*            <= M(0) * product ( ( 1 + CNORM(i) ) / |A(i,i)| ) */
/*                      1<=i<=j */

/*  and we can safely call DTRSV if 1/M(n) and 1/G(n) are both greater */
/*  than max(underflow, 1/overflow). */

/*  ===================================================================== 
*/

/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = a_dim1 + 1;
    a -= a_offset;
    --x;
    --cnorm;

    /* Function Body */
    tjjs = 0;
    *info = 0;
    upper = lsame_(uplo, "U", 1, 1);
    notran = lsame_(trans, "N", 1, 1);
    nounit = lsame_(diag, "N", 1, 1);

/*     Test the input parameters. */

    if (! upper && ! lsame_(uplo, "L", 1, 1)) {
	*info = -1;
    } else if (! notran && ! lsame_(trans, "T", 1, 1) && ! lsame_(trans, 
	    "C", 1, 1)) {
	*info = -2;
    } else if (! nounit && ! lsame_(diag, "U", 1, 1)) {
	*info = -3;
    } else if (! lsame_(normin, "Y", 1, 1) && ! lsame_(normin, "N", 1, 1))
	     {
	*info = -4;
    } else if (*n < 0) {
	*info = -5;
    } else if (*lda < Dymola_max(1,*n)) {
	*info = -7;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DLATRS", &i__1, 6);
	return 0;
    }

/*     Quick return if possible */

    if (*n == 0) {
	return 0;
    }

/*     Determine machine dependent parameters to control overflow. */

    smlnum = DBL_MIN / (DBL_EPSILON*FLT_RADIX);
    bignum = 1. / smlnum;
    *scale = 1.;

    if (lsame_(normin, "N", 1, 1)) {

/*        Compute the 1-norm of each column, not including the diagona
l. */

	if (upper) {

/*           A is upper triangular. */

	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		i__2 = j - 1;
		cnorm[j] = dasum_(&i__2, &a[j * a_dim1 + 1], &c__1);
/* L10: */
	    }
	} else {

/*           A is lower triangular. */

	    i__1 = *n - 1;
	    for (j = 1; j <= i__1; ++j) {
		i__2 = *n - j;
		cnorm[j] = dasum_(&i__2, &a[j + 1 + j * a_dim1], &c__1);
/* L20: */
	    }
	    cnorm[*n] = 0.;
	}
    }

/*     Scale the column norms by TSCAL if the maximum entry in CNORM is */
/*     greater than BIGNUM. */

    imax = idamax_(n, &cnorm[1], &c__1);
    tmax = cnorm[imax];
    if (tmax <= bignum) {
	tscal = 1.;
    } else {
	tscal = 1. / (smlnum * tmax);
	dscal_(n, &tscal, &cnorm[1], &c__1);
    }

/*     Compute a bound on the computed solution vector to see if the */
/*     Level 2 BLAS routine DTRSV can be used. */

    j = idamax_(n, &x[1], &c__1);
    xmax = (d__1 = x[j], Dymola_abs(d__1));
    xbnd = xmax;
    if (notran) {

/*        Compute the growth in A * x = b. */

	if (upper) {
	    jfirst = *n;
	    jlast = 1;
	    jinc = -1;
	} else {
	    jfirst = 1;
	    jlast = *n;
	    jinc = 1;
	}

	if (tscal != 1.) {
	    grow = 0.;
	    goto L50;
	}

	if (nounit) {

/*           A is non-unit triangular. */

/*           Compute GROW = 1/G(j) and XBND = 1/M(j). */
/*           Initially, G(0) = max{x(i), i=1,...,n}. */

	    grow = 1. / Dymola_max(xbnd,smlnum);
	    xbnd = grow;
	    i__1 = jlast;
	    i__2 = jinc;
	    for (j = jfirst; i__2 < 0 ? j >= i__1 : j <= i__1; j += i__2) {

/*              Exit the loop if the growth factor is too smal
l. */

		if (grow <= smlnum) {
		    goto L50;
		}

/*              M(j) = G(j-1) / Dymola_abs(A(j,j)) */

		tjj = (d__1 = a[j + j * a_dim1], Dymola_abs(d__1));
/* Computing MIN */
		d__1 = xbnd, d__2 = Dymola_min(1.,tjj) * grow;
		xbnd = Dymola_min(d__1,d__2);
		if (tjj + cnorm[j] >= smlnum) {

/*                 G(j) = G(j-1)*( 1 + CNORM(j) / Dymola_abs(A(j,
j)) ) */

		    grow *= tjj / (tjj + cnorm[j]);
		} else {

/*                 G(j) could overflow, set GROW to 0. */

		    grow = 0.;
		}
/* L30: */
	    }
	    grow = xbnd;
	} else {

/*           A is unit triangular. */

/*           Compute GROW = 1/G(j), where G(0) = max{x(i), i=1,...
,n}. */

/* Computing MIN */
	    d__1 = 1., d__2 = 1. / Dymola_max(xbnd,smlnum);
	    grow = Dymola_min(d__1,d__2);
	    i__2 = jlast;
	    i__1 = jinc;
	    for (j = jfirst; i__1 < 0 ? j >= i__2 : j <= i__2; j += i__1) {

/*              Exit the loop if the growth factor is too smal
l. */

		if (grow <= smlnum) {
		    goto L50;
		}

/*              G(j) = G(j-1)*( 1 + CNORM(j) ) */

		grow *= 1. / (cnorm[j] + 1.);
/* L40: */
	    }
	}
L50:

	;
    } else {

/*        Compute the growth in A' * x = b. */

	if (upper) {
	    jfirst = 1;
	    jlast = *n;
	    jinc = 1;
	} else {
	    jfirst = *n;
	    jlast = 1;
	    jinc = -1;
	}

	if (tscal != 1.) {
	    grow = 0.;
	    goto L80;
	}

	if (nounit) {

/*           A is non-unit triangular. */

/*           Compute GROW = 1/G(j) and XBND = 1/M(j). */
/*           Initially, M(0) = max{x(i), i=1,...,n}. */

	    grow = 1. / Dymola_max(xbnd,smlnum);
	    xbnd = grow;
	    i__1 = jlast;
	    i__2 = jinc;
	    for (j = jfirst; i__2 < 0 ? j >= i__1 : j <= i__1; j += i__2) {

/*              Exit the loop if the growth factor is too smal
l. */

		if (grow <= smlnum) {
		    goto L80;
		}

/*              G(j) = max( G(j-1), M(j-1)*( 1 + CNORM(j) ) ) 
*/

		xj = cnorm[j] + 1.;
/* Computing MIN */
		d__1 = grow, d__2 = xbnd / xj;
		grow = Dymola_min(d__1,d__2);

/*              M(j) = M(j-1)*( 1 + CNORM(j) ) / Dymola_abs(A(j,j)) 
*/

		tjj = (d__1 = a[j + j * a_dim1], Dymola_abs(d__1));
		if (xj > tjj) {
		    xbnd *= tjj / xj;
		}
/* L60: */
	    }
	    grow = Dymola_min(grow,xbnd);
	} else {

/*           A is unit triangular. */

/*           Compute GROW = 1/G(j), where G(0) = max{x(i), i=1,...
,n}. */

/* Computing MIN */
	    d__1 = 1., d__2 = 1. / Dymola_max(xbnd,smlnum);
	    grow = Dymola_min(d__1,d__2);
	    i__2 = jlast;
	    i__1 = jinc;
	    for (j = jfirst; i__1 < 0 ? j >= i__2 : j <= i__2; j += i__1) {

/*              Exit the loop if the growth factor is too smal
l. */

		if (grow <= smlnum) {
		    goto L80;
		}

/*              G(j) = ( 1 + CNORM(j) )*G(j-1) */

		xj = cnorm[j] + 1.;
		grow /= xj;
/* L70: */
	    }
	}
L80:
	;
    }

    if (grow * tscal > smlnum) {

/*        Use the Level 2 BLAS solve if the reciprocal of the bound on
 */
/*        elements of X is not too small. */

	dtrsv_(uplo, trans, diag, n, &a[a_offset], lda, &x[1], &c__1, 1, 1, 
		1);
    } else {

/*        Use a Level 1 BLAS solve, scaling intermediate results. */

	if (xmax > bignum) {

/*           Scale X so that its components are less than or equal
 to */
/*           BIGNUM in absolute value. */

	    *scale = bignum / xmax;
	    dscal_(n, scale, &x[1], &c__1);
	    xmax = bignum;
	}

	if (notran) {

/*           Solve A * x = b */

	    i__1 = jlast;
	    i__2 = jinc;
	    for (j = jfirst; i__2 < 0 ? j >= i__1 : j <= i__1; j += i__2) {

/*              Compute x(j) = b(j) / A(j,j), scaling x if nec
essary. */

		xj = (d__1 = x[j], Dymola_abs(d__1));
		if (nounit) {
		    tjjs = a[j + j * a_dim1] * tscal;
		} else {
		    tjjs = tscal;
		    if (tscal == 1.) {
			goto L100;
		    }
		}
		tjj = Dymola_abs(tjjs);
		if (tjj > smlnum) {

/*                    Dymola_abs(A(j,j)) > SMLNUM: */

		    if (tjj < 1.) {
			if (xj > tjj * bignum) {

/*                          Scale x by 1/b(j). */

			    rec = 1. / xj;
			    dscal_(n, &rec, &x[1], &c__1);
			    *scale *= rec;
			    xmax *= rec;
			}
		    }
		    x[j] /= tjjs;
		    xj = (d__1 = x[j], Dymola_abs(d__1));
		} else if (tjj > 0.) {

/*                    0 < Dymola_abs(A(j,j)) <= SMLNUM: */

		    if (xj > tjj * bignum) {

/*                       Scale x by (1/Dymola_abs(x(j)))*Dymola_abs(
A(j,j))*BIGNUM */
/*                       to avoid overflow when dividi
ng by A(j,j). */

			rec = tjj * bignum / xj;
			if (cnorm[j] > 1.) {

/*                          Scale by 1/CNORM(j) to
 avoid overflow when */
/*                          multiplying x(j) times
 column j. */

			    rec /= cnorm[j];
			}
			dscal_(n, &rec, &x[1], &c__1);
			*scale *= rec;
			xmax *= rec;
		    }
		    x[j] /= tjjs;
		    xj = (d__1 = x[j], Dymola_abs(d__1));
		} else {

/*                    A(j,j) = 0:  Set x(1:n) = 0, x(j) = 
1, and */
/*                    scale = 0, and compute a solution to
 A*x = 0. */

		    i__3 = *n;
		    for (i__ = 1; i__ <= i__3; ++i__) {
			x[i__] = 0.;
/* L90: */
		    }
		    x[j] = 1.;
		    xj = 1.;
		    *scale = 0.;
		    xmax = 0.;
		}
L100:

/*              Scale x if necessary to avoid overflow when ad
ding a */
/*              multiple of column j of A. */

		if (xj > 1.) {
		    rec = 1. / xj;
		    if (cnorm[j] > (bignum - xmax) * rec) {

/*                    Scale x by 1/(2*Dymola_abs(x(j))). */

			rec *= .5;
			dscal_(n, &rec, &x[1], &c__1);
			*scale *= rec;
		    }
		} else if (xj * cnorm[j] > bignum - xmax) {

/*                 Scale x by 1/2. */

		    dscal_(n, &c_b438, &x[1], &c__1);
		    *scale *= .5;
		}

		if (upper) {
		    if (j > 1) {

/*                    Compute the update */
/*                       x(1:j-1) := x(1:j-1) - x(j) *
 A(1:j-1,j) */

			i__3 = j - 1;
			d__1 = -x[j] * tscal;
			daxpy_(&i__3, &d__1, &a[j * a_dim1 + 1], &c__1, &x[1],
				 &c__1);
			i__3 = j - 1;
			i__ = idamax_(&i__3, &x[1], &c__1);
			xmax = (d__1 = x[i__], Dymola_abs(d__1));
		    }
		} else {
		    if (j < *n) {

/*                    Compute the update */
/*                       x(j+1:n) := x(j+1:n) - x(j) *
 A(j+1:n,j) */

			i__3 = *n - j;
			d__1 = -x[j] * tscal;
			daxpy_(&i__3, &d__1, &a[j + 1 + j * a_dim1], &c__1, &
				x[j + 1], &c__1);
			i__3 = *n - j;
			i__ = j + idamax_(&i__3, &x[j + 1], &c__1);
			xmax = (d__1 = x[i__], Dymola_abs(d__1));
		    }
		}
/* L110: */
	    }

	} else {

/*           Solve A' * x = b */

	    i__2 = jlast;
	    i__1 = jinc;
	    for (j = jfirst; i__1 < 0 ? j >= i__2 : j <= i__2; j += i__1) {

/*              Compute x(j) = b(j) - sum A(k,j)*x(k). */
/*                                    k<>j */

		xj = (d__1 = x[j], Dymola_abs(d__1));
		uscal = tscal;
		rec = 1. / Dymola_max(xmax,1.);
		if (cnorm[j] > (bignum - xj) * rec) {

/*                 If x(j) could overflow, scale x by 1/(2
*XMAX). */

		    rec *= .5;
		    if (nounit) {
			tjjs = a[j + j * a_dim1] * tscal;
		    } else {
			tjjs = tscal;
		    }
		    tjj = Dymola_abs(tjjs);
		    if (tjj > 1.) {

/*                       Divide by A(j,j) when scaling
 x if A(j,j) > 1. */

/* Computing MIN */
			d__1 = 1., d__2 = rec * tjj;
			rec = Dymola_min(d__1,d__2);
			uscal /= tjjs;
		    }
		    if (rec < 1.) {
			dscal_(n, &rec, &x[1], &c__1);
			*scale *= rec;
			xmax *= rec;
		    }
		}

		sumj = 0.;
		if (uscal == 1.) {

/*                 If the scaling needed for A in the dot 
product is 1, */
/*                 call DDOT to perform the dot product. 
*/

		    if (upper) {
			i__3 = j - 1;
			/* sumj = ddot_(&i__3, &a[j * a_dim1 + 1], &c__1, &x[1], &c__1); */
			 sumj = ddot1_(&i__3, &a[j * a_dim1 + 1], &x[1]);
		    } else if (j < *n) {
			i__3 = *n - j;
			/* sumj = ddot_(&i__3, &a[j + 1 + j * a_dim1], &c__1, &x[j + 1], &c__1); */
			sumj = ddot1_(&i__3, &a[j + 1 + j * a_dim1], &x[j + 1]);
		    }
		} else {

/*                 Otherwise, use in-line code for the dot
 product. */

		    if (upper) {
			i__3 = j - 1;
			for (i__ = 1; i__ <= i__3; ++i__) {
			    sumj += a[i__ + j * a_dim1] * uscal * x[i__];
/* L120: */
			}
		    } else if (j < *n) {
			i__3 = *n;
			for (i__ = j + 1; i__ <= i__3; ++i__) {
			    sumj += a[i__ + j * a_dim1] * uscal * x[i__];
/* L130: */
			}
		    }
		}

		if (uscal == tscal) {

/*                 Compute x(j) := ( x(j) - sumj ) / A(j,j
) if 1/A(j,j) */
/*                 was not used to scale the dotproduct. 
*/

		    x[j] -= sumj;
		    xj = (d__1 = x[j], Dymola_abs(d__1));
		    if (nounit) {
			tjjs = a[j + j * a_dim1] * tscal;
		    } else {
			tjjs = tscal;
			if (tscal == 1.) {
			    goto L150;
			}
		    }

/*                    Compute x(j) = x(j) / A(j,j), scalin
g if necessary. */

		    tjj = Dymola_abs(tjjs);
		    if (tjj > smlnum) {

/*                       Dymola_abs(A(j,j)) > SMLNUM: */

			if (tjj < 1.) {
			    if (xj > tjj * bignum) {

/*                             Scale X by 1/ab
s(x(j)). */

				rec = 1. / xj;
				dscal_(n, &rec, &x[1], &c__1);
				*scale *= rec;
				xmax *= rec;
			    }
			}
			x[j] /= tjjs;
		    } else if (tjj > 0.) {

/*                       0 < Dymola_abs(A(j,j)) <= SMLNUM: */

			if (xj > tjj * bignum) {

/*                          Scale x by (1/Dymola_abs(x(j)
))*Dymola_abs(A(j,j))*BIGNUM. */

			    rec = tjj * bignum / xj;
			    dscal_(n, &rec, &x[1], &c__1);
			    *scale *= rec;
			    xmax *= rec;
			}
			x[j] /= tjjs;
		    } else {

/*                       A(j,j) = 0:  Set x(1:n) = 0, 
x(j) = 1, and */
/*                       scale = 0, and compute a solu
tion to A'*x = 0. */

			i__3 = *n;
			for (i__ = 1; i__ <= i__3; ++i__) {
			    x[i__] = 0.;
/* L140: */
			}
			x[j] = 1.;
			*scale = 0.;
			xmax = 0.;
		    }
L150:
		    ;
		} else {

/*                 Compute x(j) := x(j) / A(j,j)  - sumj i
f the dot */
/*                 product has already been divided by 1/A
(j,j). */

		    x[j] = x[j] / tjjs - sumj;
		}
/* Computing MAX */
		d__2 = xmax, d__3 = (d__1 = x[j], Dymola_abs(d__1));
		xmax = Dymola_max(d__2,d__3);
/* L160: */
	    }
	}
	*scale /= tscal;
    }

/*     Scale the column norms by 1/TSCAL for return. */

    if (tscal != 1.) {
	d__1 = 1. / tscal;
	dscal_(n, &d__1, &cnorm[1], &c__1);
    }

    return 0;

/*     End of DLATRS */

} /* dlatrs_ */

/* Subroutine */ DYMOLA_STATIC int dlatzm_(side, m, n, v, incv, tau, c1, c2, ldc, work, 
	side_len)
const char*side;
const integer*m, *n;
doublereal *v;
const integer*incv;
doublereal *tau, *c1, *c2;
const integer*ldc;
doublereal *work;
ftnlen side_len;
{
    /* System generated locals */
    integer c1_dim1, c1_offset, c2_dim1, c2_offset, i__1;
    doublereal d__1;

    /* Local variables */

/*  -- LAPACK routine (version 1.1) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
/*     Courant Institute, Argonne National Lab, and Rice University */
/*     February 29, 1992 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  DLATZM applies a Householder matrix generated by DTZRQF to a matrix. 
*/

/*  Let P = I - tau*u*u',   u = ( 1 ), */
/*                              ( v ) */
/*  where v is an (m-1) vector if SIDE = 'L', or a (n-1) vector if */
/*  SIDE = 'R'. */

/*  If SIDE equals 'L', let */
/*         C = [ C1 ] 1 */
/*             [ C2 ] m-1 */
/*               n */
/*  Then C is overwritten by P*C. */

/*  If SIDE equals 'R', let */
/*         C = [ C1, C2 ] m */
/*                1  n-1 */
/*  Then C is overwritten by C*P. */

/*  Arguments */
/*  ========= */

/*  SIDE    (input) CHARACTER*1 */
/*          = 'L': form P * C */
/*          = 'R': form C * P */

/*  M       (input) INTEGER */
/*          The number of rows of the matrix C. */

/*  N       (input) INTEGER */
/*          The number of columns of the matrix C. */

/*  V       (input) DOUBLE PRECISION array, dimension */
/*                  (1 + (M-1)*Dymola_abs(INCV)) if SIDE = 'L' */
/*                  (1 + (N-1)*Dymola_abs(INCV)) if SIDE = 'R' */
/*          The vector v in the representation of P. V is not used */
/*          if TAU = 0. */

/*  INCV    (input) INTEGER */
/*          The increment between elements of v. INCV <> 0 */

/*  TAU     (input) DOUBLE PRECISION */
/*          The value tau in the representation of P. */

/*  C1      (input/output) DOUBLE PRECISION array, dimension */
/*                         (LDC,N) if SIDE = 'L' */
/*                         (M,1)   if SIDE = 'R' */
/*          On entry, the n-vector C1 if SIDE = 'L', or the m-vector C1 */
/*          if SIDE = 'R'. */

/*          On exit, the first row of P*C if SIDE = 'L', or the first */
/*          column of C*P if SIDE = 'R'. */

/*  C2      (input/output) DOUBLE PRECISION array, dimension */
/*                         (LDC, N)   if SIDE = 'L' */
/*                         (LDC, N-1) if SIDE = 'R' */
/*          On entry, the (m - 1) x n matrix C2 if SIDE = 'L', or the */
/*          m x (n - 1) matrix C2 if SIDE = 'R'. */

/*          On exit, rows 2:m of P*C if SIDE = 'L', or columns 2:m of C*P 
*/
/*          if SIDE = 'R'. */

/*  LDC     (input) INTEGER */
/*          The leading dimension of the arrays C1 and C2. LDC >= (1,M). 
*/

/*  WORK    (workspace) DOUBLE PRECISION array, dimension */
/*                      (N) if SIDE = 'L' */
/*                      (M) if SIDE = 'R' */

/*  ===================================================================== 
*/

/*     .. Parameters .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

    /* Parameter adjustments */
    --v;
    c2_dim1 = *ldc;
    c2_offset = c2_dim1 + 1;
    c2 -= c2_offset;
    c1_dim1 = *ldc;
    c1_offset = c1_dim1 + 1;
    c1 -= c1_offset;
    --work;

    /* Function Body */
    if (Dymola_min(*m,*n) == 0 || *tau == 0.) {
	return 0;
    }

    if (lsame_(side, "L", 1, 1)) {

/*        w := C1 + v' * C2 */

	dcopy_(n, &c1[c1_offset], ldc, &work[1], &c__1);
	i__1 = *m - 1;
	dgemv_("Transpose", &i__1, n, &c_b263, &c2[c2_offset], ldc, &v[1], 
		incv, &c_b263, &work[1], &c__1, 9);

/*        [ C1 ] := [ C1 ] - tau* [ 1 ] * w' */
/*        [ C2 ]    [ C2 ]        [ v ] */

	d__1 = -(*tau);
	daxpy_(n, &d__1, &work[1], &c__1, &c1[c1_offset], ldc);
	i__1 = *m - 1;
	d__1 = -(*tau);
	dger_(&i__1, n, &d__1, &v[1], incv, &work[1], &c__1, &c2[c2_offset], 
		ldc);

    } else if (lsame_(side, "R", 1, 1)) {

/*        w := C1 + C2 * v */

	dcopy_(m, &c1[c1_offset], &c__1, &work[1], &c__1);
	i__1 = *n - 1;
	dgemv_("No transpose", m, &i__1, &c_b263, &c2[c2_offset], ldc, &v[1], 
		incv, &c_b263, &work[1], &c__1, 12);

/*        [ C1, C2 ] := [ C1, C2 ] - tau* w * [ 1 , v'] */

	d__1 = -(*tau);
	daxpy_(m, &d__1, &work[1], &c__1, &c1[c1_offset], &c__1);
	i__1 = *n - 1;
	d__1 = -(*tau);
	dger_(m, &i__1, &d__1, &work[1], &c__1, &v[1], incv, &c2[c2_offset], 
		ldc);
    }

    return 0;

/*     End of DLATZM */

} /* dlatzm_ */

/* Subroutine */ DYMOLA_STATIC int dorm2r_(side, trans, m, n, k, a, lda, tau, c__, ldc, 
	work, info, side_len, trans_len)
const char*side, *trans;
const integer*m, *n, *k;
doublereal *a;
const integer*lda;
doublereal *tau, *c__;
const integer*ldc;
doublereal *work;
integer *info;
ftnlen side_len;
ftnlen trans_len;
{
    /* System generated locals */
    integer a_dim1, a_offset, c_dim1, c_offset, i__1, i__2;

    /* Local variables */
    integer i__, i1, i2, i3, ic, jc, mi, ni, nq;
    doublereal aii;
    logical left;
    logical notran;


/*  -- LAPACK routine (version 1.1) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
/*     Courant Institute, Argonne National Lab, and Rice University */
/*     February 29, 1992 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  DORM2R overwrites the general real m by n matrix C with */

/*        Q * C  if SIDE = 'L' and TRANS = 'N', or */

/*        Q'* C  if SIDE = 'L' and TRANS = 'T', or */

/*        C * Q  if SIDE = 'R' and TRANS = 'N', or */

/*        C * Q' if SIDE = 'R' and TRANS = 'T', */

/*  where Q is a real orthogonal matrix defined as the product of k */
/*  elementary reflectors */

/*        Q = H(1) H(2) . . . H(k) */

/*  as returned by DGEQRF. Q is of order m if SIDE = 'L' and of order n */
/*  if SIDE = 'R'. */

/*  Arguments */
/*  ========= */

/*  SIDE    (input) CHARACTER*1 */
/*          = 'L': apply Q or Q' from the Left */
/*          = 'R': apply Q or Q' from the Right */

/*  TRANS   (input) CHARACTER*1 */
/*          = 'N': apply Q  (No transpose) */
/*          = 'T': apply Q' (Transpose) */

/*  M       (input) INTEGER */
/*          The number of rows of the matrix C. M >= 0. */

/*  N       (input) INTEGER */
/*          The number of columns of the matrix C. N >= 0. */

/*  K       (input) INTEGER */
/*          The number of elementary reflectors whose product defines */
/*          the matrix Q. */
/*          If SIDE = 'L', M >= K >= 0; */
/*          if SIDE = 'R', N >= K >= 0. */

/*  A       (input) DOUBLE PRECISION array, dimension (LDA,K) */
/*          The i-th column must contain the vector which defines the */
/*          elementary reflector H(i), for i = 1,2,...,k, as returned by 
*/
/*          DGEQRF in the first k columns of its array argument A. */
/*          A is modified by the routine but restored on exit. */

/*  LDA     (input) INTEGER */
/*          The leading dimension of the array A. */
/*          If SIDE = 'L', LDA >= max(1,M); */
/*          if SIDE = 'R', LDA >= max(1,N). */

/*  TAU     (input) DOUBLE PRECISION array, dimension (K) */
/*          TAU(i) must contain the scalar factor of the elementary */
/*          reflector H(i), as returned by DGEQRF. */

/*  C       (input/output) DOUBLE PRECISION array, dimension (LDC,N) */
/*          On entry, the m by n matrix C. */
/*          On exit, C is overwritten by Q*C or Q'*C or C*Q' or C*Q. */

/*  LDC     (input) INTEGER */
/*          The leading dimension of the array C. LDC >= max(1,M). */

/*  WORK    (workspace) DOUBLE PRECISION array, dimension */
/*                                   (N) if SIDE = 'L', */
/*                                   (M) if SIDE = 'R' */

/*  INFO    (output) INTEGER */
/*          = 0: successful exit */
/*          < 0: if INFO = -i, the i-th argument had an illegal value */

/*  ===================================================================== 
*/

/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input arguments */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = a_dim1 + 1;
    a -= a_offset;
    --tau;
    c_dim1 = *ldc;
    c_offset = c_dim1 + 1;
    c__ -= c_offset;
    --work;

    /* Function Body */
    ic = jc = 0;
    *info = 0;
    left = lsame_(side, "L", 1, 1);
    notran = lsame_(trans, "N", 1, 1);

/*     NQ is the order of Q */

    if (left) {
	nq = *m;
    } else {
	nq = *n;
    }
    if (! left && ! lsame_(side, "R", 1, 1)) {
	*info = -1;
    } else if (! notran && ! lsame_(trans, "T", 1, 1)) {
	*info = -2;
    } else if (*m < 0) {
	*info = -3;
    } else if (*n < 0) {
	*info = -4;
    } else if (*k < 0 || *k > nq) {
	*info = -5;
    } else if (*lda < Dymola_max(1,nq)) {
	*info = -7;
    } else if (*ldc < Dymola_max(1,*m)) {
	*info = -10;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DORM2R", &i__1, 6);
	return 0;
    }

/*     Quick return if possible */

    if (*m == 0 || *n == 0 || *k == 0) {
	return 0;
    }

    if ((left && ! notran) || (! left && notran)) {
	i1 = 1;
	i2 = *k;
	i3 = 1;
    } else {
	i1 = *k;
	i2 = 1;
	i3 = -1;
    }

    if (left) {
	ni = *n;
	jc = 1;
    } else {
	mi = *m;
	ic = 1;
    }

    i__1 = i2;
    i__2 = i3;
    for (i__ = i1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
	if (left) {

/*           H(i) is applied to C(i:m,1:n) */

	    mi = *m - i__ + 1;
	    ic = i__;
	} else {

/*           H(i) is applied to C(1:m,i:n) */

	    ni = *n - i__ + 1;
	    jc = i__;
	}

/*        Apply H(i) */

	aii = a[i__ + i__ * a_dim1];
	a[i__ + i__ * a_dim1] = 1.;
	dlarf_(side, &mi, &ni, &a[i__ + i__ * a_dim1], &c__1, &tau[i__], &c__[
		ic + jc * c_dim1], ldc, &work[1], 1);
	a[i__ + i__ * a_dim1] = aii;
/* L10: */
    }
    return 0;

/*     End of DORM2R */

} /* dorm2r_ */

/* Subroutine */ DYMOLA_STATIC int drscl_(n, sa, sx, incx)
const integer*n;
doublereal *sa, *sx;
const integer*incx;
{
    doublereal mul, cden;
    logical done;
    doublereal cnum, cden1, cnum1;
    doublereal bignum, smlnum;


/*  -- LAPACK auxiliary routine (version 1.1) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
/*     Courant Institute, Argonne National Lab, and Rice University */
/*     October 31, 1992 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  DRSCL multiplies an n-element real vector x by the real scalar 1/a. */
/*  This is done without overflow or underflow as long as */
/*  the final result x/a does not overflow or underflow. */

/*  Arguments */
/*  ========= */

/*  N       (input) INTEGER */
/*          The number of components of the vector x. */

/*  SA      (input) DOUBLE PRECISION */
/*          The scalar a which is used to divide each component of x. */
/*          SA must be >= 0, or the subroutine will divide by zero. */

/*  SX      (input/output) DOUBLE PRECISION array, dimension */
/*                         (1+(N-1)*Dymola_abs(INCX)) */
/*          The n-element vector x. */

/*  INCX    (input) INTEGER */
/*          The increment between successive values of the vector SX. */
/*          > 0:  SX(1) = X(1) and SX(1+(i-1)*INCX) = x(i),     1< i<= n 
*/
/*          < 0:  SX(1) = X(n) and SX(1+(i-1)*INCX) = x(n-i+1), 1< i<= n 
*/

/* ===================================================================== 
*/

/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Quick return if possible */

    /* Parameter adjustments */
    --sx;

    /* Function Body */
    if (*n <= 0) {
	return 0;
    }

/*     Get machine parameters */

    smlnum = DBL_MIN;
    bignum = 1. / smlnum;
    dlabad_(&smlnum, &bignum);

/*     Initialize the denominator to SA and the numerator to 1. */

    cden = *sa;
    cnum = 1.;

L10:
    cden1 = cden * smlnum;
    cnum1 = cnum / bignum;
    if (Dymola_abs(cden1) > Dymola_abs(cnum) && cnum != 0.) {

/*        Pre-multiply X by SMLNUM if CDEN is large compared to CNUM. 
*/

	mul = smlnum;
	done = FALSE_;
	cden = cden1;
    } else if (Dymola_abs(cnum1) > Dymola_abs(cden)) {

/*        Pre-multiply X by BIGNUM if CDEN is small compared to CNUM. 
*/

	mul = bignum;
	done = FALSE_;
	cnum = cnum1;
    } else {

/*        Multiply X by CNUM / CDEN and return. */

	mul = cnum / cden;
	done = TRUE_;
    }

/*     Scale the vector X by MUL */

    dscal_(n, &mul, &sx[1], incx);

    if (! done) {
	goto L10;
    }

    return 0;

/*     End of DRSCL */

} /* drscl_ */

/* Subroutine */ DYMOLA_STATIC int dtzrqf_(m, n, a, lda, tau, info)
const integer*m, *n;
doublereal *a;
const integer*lda;
doublereal *tau;
integer *info;
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;
    doublereal d__1;

    /* Local variables */
    integer i__, k, m1;


/*  -- LAPACK routine (version 1.1) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
/*     Courant Institute, Argonne National Lab, and Rice University */
/*     March 31, 1993 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  DTZRQF reduces the M-by-N ( M<=N ) real upper trapezoidal matrix A */
/*  to upper triangular form by means of orthogonal transformations. */

/*  The upper trapezoidal matrix A is factored as */

/*     A = ( R  0 ) * Z, */

/*  where Z is an N-by-N orthogonal matrix and R is an M-by-M upper */
/*  triangular matrix. */

/*  Arguments */
/*  ========= */

/*  M       (input) INTEGER */
/*          The number of rows of the matrix A.  M >= 0. */

/*  N       (input) INTEGER */
/*          The number of columns of the matrix A.  N >= M. */

/*  A       (input/output) DOUBLE PRECISION array, dimension */
/*                         (LDA,max(1,N)) */
/*          On entry, the leading M-by-N upper trapezoidal part of the */
/*          array A must contain the matrix to be factorized. */
/*          On exit, the leading M-by-M upper triangular part of A */
/*          contains the upper triangular matrix R, and elements M+1 to */
/*          N of the first M rows of A, with the array TAU, represent the 
*/
/*          orthogonal matrix Z as a product of M elementary reflectors. 
*/

/*  LDA     (input) INTEGER */
/*          The leading dimension of the array A.  LDA >= max(1,M). */

/*  TAU     (output) DOUBLE PRECISION array, dimension (max(1,M)) */
/*          The scalar factors of the elementary reflectors. */

/*  INFO    (output) INTEGER */
/*          = 0:  successful exit */
/*          < 0:  if INFO = -i, the i-th argument had an illegal value */

/*  Further Details */
/*  =============== */

/*  The factorization is obtained by Householder's method.  The kth */
/*  transformation matrix, Z( k ), which is used to introduce zeros into 
*/
/*  the ( m - k + 1 )th row of A, is given in the form */

/*     Z( k ) = ( I     0   ), */
/*              ( 0  T( k ) ) */

/*  where */

/*     T( k ) = I - tau*u( k )*u( k )',   u( k ) = (   1    ), */
/*                                                 (   0    ) */
/*                                                 ( z( k ) ) */

/*  tau is a scalar and z( k ) is an ( n - m ) element vector. */
/*  tau and z( k ) are chosen to annihilate the elements of the kth row */
/*  of X. */

/*  The scalar tau is returned in the kth element of TAU and the vector */
/*  u( k ) in the kth row of A, such that the elements of z( k ) are */
/*  in  a( k, m + 1 ), ..., a( k, n ). The elements of R are returned in 
*/
/*  the upper triangular part of A. */

/*  Z is given by */

/*     Z =  Z( 1 ) * Z( 2 ) * ... * Z( m ). */

/*  ===================================================================== 
*/

/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters. */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = a_dim1 + 1;
    a -= a_offset;
    --tau;

    /* Function Body */
    *info = 0;
    if (*m < 0) {
	*info = -1;
    } else if (*n < *m) {
	*info = -2;
    } else if (*lda < Dymola_max(1,*m)) {
	*info = -4;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DTZRQF", &i__1, 6);
	return 0;
    }

/*     Perform the factorization. */

    if (*m == 0) {
	return 0;
    }
    if (*m == *n) {
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    tau[i__] = 0.;
/* L10: */
	}
    } else {
/* Computing MIN */
	i__1 = *m + 1;
	m1 = Dymola_min(i__1,*n);
	for (k = *m; k >= 1; --k) {

/*           Use a Householder reflection to zero the kth row of A
. */
/*           First set up the reflection. */

	    i__1 = *n - *m + 1;
	    dlarfg_(&i__1, &a[k + k * a_dim1], &a[k + m1 * a_dim1], lda, &tau[
		    k]);

	    if (tau[k] != 0. && k > 1) {

/*              We now perform the operation  A := A*P( k ). 
*/

/*              Use the first ( k - 1 ) elements of TAU to sto
re  a( k ), */
/*              where  a( k ) consists of the first ( k - 1 ) 
elements of */
/*              the  kth column  of  A.  Also  let  B  denote 
 the  first */
/*              ( k - 1 ) rows of the last ( n - m ) columns o
f A. */

		i__1 = k - 1;
		dcopy_(&i__1, &a[k * a_dim1 + 1], &c__1, &tau[1], &c__1);

/*              Form   w = a( k ) + B*z( k )  in TAU. */

		i__1 = k - 1;
		i__2 = *n - *m;
		dgemv_("No transpose", &i__1, &i__2, &c_b263, &a[m1 * a_dim1 
			+ 1], lda, &a[k + m1 * a_dim1], lda, &c_b263, &tau[1],
			 &c__1, 12);

/*              Now form  a( k ) := a( k ) - tau*w */
/*              and       B      := B      - tau*w*z( k )'. */

		i__1 = k - 1;
		d__1 = -tau[k];
		daxpy_(&i__1, &d__1, &tau[1], &c__1, &a[k * a_dim1 + 1], &
			c__1);
		i__1 = k - 1;
		i__2 = *n - *m;
		d__1 = -tau[k];
		dger_(&i__1, &i__2, &d__1, &tau[1], &c__1, &a[k + m1 * a_dim1]
			, lda, &a[m1 * a_dim1 + 1], lda);
	    }
/* L20: */
	}
    }

    return 0;

/*     End of DTZRQF */

} /* dtzrqf_ */

DYMOLA_STATIC integer ilaenv_(ispec, name__, opts, n1, n2, n3, n4, name_len, opts_len)
const integer *ispec;
const char *name__, *opts;
const integer*n1, *n2, *n3, *n4;
ftnlen name_len;
ftnlen opts_len;
{
    /* System generated locals */
    integer ret_val;


    /* Local variables */
    integer i__;
    char c1[1], c2[2], c3[3], c4[2];
    integer ic, nb, iz, nx;
    logical cname, sname;
    integer nbmin;
    char subnam[6];


/*  -- LAPACK auxiliary routine (preliminary version) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
/*     Courant Institute, Argonne National Lab, and Rice University */
/*     February 20, 1992 */

/*     .. Scalar Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  ILAENV is called from the LAPACK routines to choose problem-dependent 
*/
/*  parameters for the local environment.  See ISPEC for a description of 
*/
/*  the parameters. */

/*  This version provides a set of parameters which should give good, */
/*  but not optimal, performance on many of the currently available */
/*  computers.  Users are encouraged to modify this subroutine to set */
/*  the tuning parameters for their particular machine using the option */
/*  and problem size information in the arguments. */

/*  This routine will not function correctly if it is converted to all */
/*  lower case.  Converting it to all upper case is allowed. */

/*  Arguments */
/*  ========= */

/*  ISPEC   (input) INTEGER */
/*          Specifies the parameter to be returned as the value of */
/*          ILAENV. */
/*          = 1: the optimal blocksize; if this value is 1, an unblocked 
*/
/*               algorithm will give the best performance. */
/*          = 2: the minimum block size for which the block routine */
/*               should be used; if the usable block size is less than */
/*               this value, an unblocked routine should be used. */
/*          = 3: the crossover point (in a block routine, for N less */
/*               than this value, an unblocked routine should be used) */
/*          = 4: the number of shifts, used in the nonsymmetric */
/*               eigenvalue routines */
/*          = 5: the minimum column dimension for blocking to be used; */
/*               rectangular blocks must have dimension at least k by m, 
*/
/*               where k is given by ILAENV(2,...) and m by ILAENV(5,...) 
*/
/*          = 6: the crossover point for the SVD (when reducing an m by n 
*/
/*               matrix to bidiagonal form, if max(m,n)/min(m,n) exceeds 
*/
/*               this value, a QR factorization is used first to reduce */
/*               the matrix to a triangular form.) */
/*          = 7: the number of processors */
/*          = 8: the crossover point for the multishift QR and QZ methods 
*/
/*               for nonsymmetric eigenvalue problems. */

/*  NAME    (input) CHARACTER*(*) */
/*          The name of the calling subroutine, in either upper case or */
/*          lower case. */

/*  OPTS    (input) CHARACTER*(*) */
/*          The character options to the subroutine NAME, concatenated */
/*          into a single character string.  For example, UPLO = 'U', */
/*          TRANS = 'T', and DIAG = 'N' for a triangular routine would */
/*          be specified as OPTS = 'UTN'. */

/*  N1      (input) INTEGER */
/*  N2      (input) INTEGER */
/*  N3      (input) INTEGER */
/*  N4      (input) INTEGER */
/*          Problem dimensions for the subroutine NAME; these may not all 
*/
/*          be required. */

/* (ILAENV) (output) INTEGER */
/*          >= 0: the value of the parameter specified by ISPEC */
/*          < 0:  if ILAENV = -k, the k-th argument had an illegal value. 
*/

/*  Further Details */
/*  =============== */

/*  The following conventions have been used when calling ILAENV from the 
*/
/*  LAPACK routines: */
/*  1)  OPTS is a concatenation of all of the character options to */
/*      subroutine NAME, in the same order that they appear in the */
/*      argument list for NAME, even if they are not used in determining 
*/
/*      the value of the parameter specified by ISPEC. */
/*  2)  The problem dimensions N1, N2, N3, N4 are specified in the order 
*/
/*      that they appear in the argument list for NAME.  N1 is used */
/*      first, N2 second, and so on, and unused problem dimensions are */
/*      passed a value of -1. */
/*  3)  The parameter value returned by ILAENV is checked for validity in 
*/
/*      the calling subroutine.  For example, ILAENV is used to retrieve 
*/
/*      the optimal blocksize for STRTRI as follows: */

/*      NB = ILAENV( 1, 'STRTRI', UPLO // DIAG, N, -1, -1, -1 ) */
/*      IF( NB.LE.1 ) NB = MAX( 1, N ) */

/*  ===================================================================== 
*/

/*     .. Local Scalars .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

    switch ((int)*ispec) {
	case 1:  goto L100;
	case 2:  goto L100;
	case 3:  goto L100;
	case 4:  goto L400;
	case 5:  goto L500;
	case 6:  goto L600;
	case 7:  goto L700;
	case 8:  goto L800;
    }

/*     Invalid value for ISPEC */

    ret_val = -1;
    return ret_val;

L100:

/*     Convert NAME to upper case if the first character is lower case. */

    ret_val = 1;
    s_copy(subnam, name__, 6, name_len);
    ic = *(unsigned char *)subnam;
    iz = 'Z';
    if (iz == 90 || iz == 122) {

/*        ASCII character set */

	if (ic >= 97 && ic <= 122) {
	    *(unsigned char *)subnam = (char) (ic - 32);
	    for (i__ = 2; i__ <= 6; ++i__) {
		ic = *(unsigned char *)&subnam[i__ - 1];
		if (ic >= 97 && ic <= 122) {
		    *(unsigned char *)&subnam[i__ - 1] = (char) (ic - 32);
		}
/* L10: */
	    }
	}

    } else if (iz == 233 || iz == 169) {

/*        EBCDIC character set */

	if ((ic >= 129 && ic <= 137) || (ic >= 145 && ic <= 153) || (ic >= 162 && ic <= 169)) {
	    *(unsigned char *)subnam = (char) (ic + 64);
	    for (i__ = 2; i__ <= 6; ++i__) {
		ic = *(unsigned char *)&subnam[i__ - 1];
		if ((ic >= 129 && ic <= 137) || (ic >= 145 && ic <= 153) || (ic >= 162 && ic <= 169)) {
		    *(unsigned char *)&subnam[i__ - 1] = (char) (ic + 64);
		}
/* L20: */
	    }
	}

    } else if (iz == 218 || iz == 250) {

/*        Prime machines:  ASCII+128 */

	if (ic >= 225 && ic <= 250) {
	    *(unsigned char *)subnam = (char) (ic - 32);
	    for (i__ = 2; i__ <= 6; ++i__) {
		ic = *(unsigned char *)&subnam[i__ - 1];
		if (ic >= 225 && ic <= 250) {
		    *(unsigned char *)&subnam[i__ - 1] = (char) (ic - 32);
		}
/* L30: */
	    }
	}
    }

    *(unsigned char *)c1 = *(unsigned char *)subnam;
    sname = *(unsigned char *)c1 == 'S' || *(unsigned char *)c1 == 'D';
    cname = *(unsigned char *)c1 == 'C' || *(unsigned char *)c1 == 'Z';
    if (! (cname || sname)) {
	return ret_val;
    }
    s_copy(c2, subnam + 1, 2, 2);
    s_copy(c3, subnam + 3, 3, 3);
    s_copy(c4, c3 + 1, 2, 2);

    switch ((int)*ispec) {
	case 1:  goto L110;
	case 2:  goto L200;
	case 3:  goto L300;
    }

L110:

/*     ISPEC = 1:  block size */

/*     In these examples, separate code is provided for setting NB for */
/*     real and complex.  We assume that NB will take the same value in */
/*     single or double precision. */

    nb = 1;

    if (s_cmp(c2, "GE", 2, 2) == 0) {
	if (s_cmp(c3, "TRF", 3, 3) == 0) {
	    if (sname) {
		nb = 64;
	    } else {
		nb = 64;
	    }
	} else if (s_cmp(c3, "QRF", 3, 3) == 0 || s_cmp(c3, "RQF", 3, 3) 
		== 0 || s_cmp(c3, "LQF", 3, 3) == 0 || s_cmp(c3, "QLF", 3, 
		3) == 0) {
	    if (sname) {
		nb = 32;
	    } else {
		nb = 32;
	    }
	} else if (s_cmp(c3, "HRD", 3, 3) == 0) {
	    if (sname) {
		nb = 32;
	    } else {
		nb = 32;
	    }
	} else if (s_cmp(c3, "BRD", 3, 3) == 0) {
	    if (sname) {
		nb = 32;
	    } else {
		nb = 32;
	    }
	} else if (s_cmp(c3, "TRI", 3, 3) == 0) {
	    if (sname) {
		nb = 64;
	    } else {
		nb = 64;
	    }
	}
    } else if (s_cmp(c2, "PO", 2, 2) == 0) {
	if (s_cmp(c3, "TRF", 3, 3) == 0) {
	    if (sname) {
		nb = 64;
	    } else {
		nb = 64;
	    }
	}
    } else if (s_cmp(c2, "SY", 2, 2) == 0) {
	if (s_cmp(c3, "TRF", 3, 3) == 0) {
	    if (sname) {
		nb = 64;
	    } else {
		nb = 64;
	    }
	} else if (sname && s_cmp(c3, "TRD", 3, 3) == 0) {
	    nb = 1;
	} else if (sname && s_cmp(c3, "GST", 3, 3) == 0) {
	    nb = 64;
	}
    } else if (cname && s_cmp(c2, "HE", 2, 2) == 0) {
	if (s_cmp(c3, "TRF", 3, 3) == 0) {
	    nb = 64;
	} else if (s_cmp(c3, "TRD", 3, 3) == 0) {
	    nb = 1;
	} else if (s_cmp(c3, "GST", 3, 3) == 0) {
	    nb = 64;
	}
    } else if (sname && s_cmp(c2, "OR", 2, 2) == 0) {
	if (*(unsigned char *)c3 == 'G') {
	    if (s_cmp(c4, "QR", 2, 2) == 0 || s_cmp(c4, "RQ", 2, 2) == 0 
		    || s_cmp(c4, "LQ", 2, 2) == 0 || s_cmp(c4, "QL", 2, 2)
		     == 0 || s_cmp(c4, "HR", 2, 2) == 0 || s_cmp(c4, "TR", 
		    2, 2) == 0 || s_cmp(c4, "BR", 2, 2) == 0) {
		nb = 32;
	    }
	} else if (*(unsigned char *)c3 == 'M') {
	    if (s_cmp(c4, "QR", 2, 2) == 0 || s_cmp(c4, "RQ", 2, 2) == 0 
		    || s_cmp(c4, "LQ", 2, 2) == 0 || s_cmp(c4, "QL", 2, 2)
		     == 0 || s_cmp(c4, "HR", 2, 2) == 0 || s_cmp(c4, "TR", 
		    2, 2) == 0 || s_cmp(c4, "BR", 2, 2) == 0) {
		nb = 32;
	    }
	}
    } else if (cname && s_cmp(c2, "UN", 2, 2) == 0) {
	if (*(unsigned char *)c3 == 'G') {
	    if (s_cmp(c4, "QR", 2, 2) == 0 || s_cmp(c4, "RQ", 2, 2) == 0 
		    || s_cmp(c4, "LQ", 2, 2) == 0 || s_cmp(c4, "QL", 2, 2)
		     == 0 || s_cmp(c4, "HR", 2, 2) == 0 || s_cmp(c4, "TR", 
		    2, 2) == 0 || s_cmp(c4, "BR", 2, 2) == 0) {
		nb = 32;
	    }
	} else if (*(unsigned char *)c3 == 'M') {
	    if (s_cmp(c4, "QR", 2, 2) == 0 || s_cmp(c4, "RQ", 2, 2) == 0 
		    || s_cmp(c4, "LQ", 2, 2) == 0 || s_cmp(c4, "QL", 2, 2)
		     == 0 || s_cmp(c4, "HR", 2, 2) == 0 || s_cmp(c4, "TR", 
		    2, 2) == 0 || s_cmp(c4, "BR", 2, 2) == 0) {
		nb = 32;
	    }
	}
    } else if (s_cmp(c2, "GB", 2, 2) == 0) {
	if (s_cmp(c3, "TRF", 3, 3) == 0) {
	    if (sname) {
		if (*n4 <= 64) {
		    nb = 1;
		} else {
		    nb = 32;
		}
	    } else {
		if (*n4 <= 64) {
		    nb = 1;
		} else {
		    nb = 32;
		}
	    }
	}
    } else if (s_cmp(c2, "PB", 2, 2) == 0) {
	if (s_cmp(c3, "TRF", 3, 3) == 0) {
	    if (sname) {
		if (*n2 <= 64) {
		    nb = 1;
		} else {
		    nb = 32;
		}
	    } else {
		if (*n2 <= 64) {
		    nb = 1;
		} else {
		    nb = 32;
		}
	    }
	}
    } else if (s_cmp(c2, "TR", 2, 2) == 0) {
	if (s_cmp(c3, "TRI", 3, 3) == 0) {
	    if (sname) {
		nb = 64;
	    } else {
		nb = 64;
	    }
	}
    } else if (s_cmp(c2, "LA", 2, 2) == 0) {
	if (s_cmp(c3, "UUM", 3, 3) == 0) {
	    if (sname) {
		nb = 64;
	    } else {
		nb = 64;
	    }
	}
    } else if (sname && s_cmp(c2, "ST", 2, 2) == 0) {
	if (s_cmp(c3, "EBZ", 3, 3) == 0) {
	    nb = 1;
	}
    }
    ret_val = nb;
    return ret_val;

L200:

/*     ISPEC = 2:  minimum block size */

    nbmin = 2;
    if (s_cmp(c2, "GE", 2, 2) == 0) {
	if (s_cmp(c3, "QRF", 3, 3) == 0 || s_cmp(c3, "RQF", 3, 3) == 0 || 
		s_cmp(c3, "LQF", 3, 3) == 0 || s_cmp(c3, "QLF", 3, 3) == 
		0) {
	    if (sname) {
		nbmin = 2;
	    } else {
		nbmin = 2;
	    }
	} else if (s_cmp(c3, "HRD", 3, 3) == 0) {
	    if (sname) {
		nbmin = 2;
	    } else {
		nbmin = 2;
	    }
	} else if (s_cmp(c3, "BRD", 3, 3) == 0) {
	    if (sname) {
		nbmin = 2;
	    } else {
		nbmin = 2;
	    }
	} else if (s_cmp(c3, "TRI", 3, 3) == 0) {
	    if (sname) {
		nbmin = 2;
	    } else {
		nbmin = 2;
	    }
	}
    } else if (s_cmp(c2, "SY", 2, 2) == 0) {
	if (s_cmp(c3, "TRF", 3, 3) == 0) {
	    if (sname) {
		nbmin = 2;
	    } else {
		nbmin = 2;
	    }
	} else if (sname && s_cmp(c3, "TRD", 3, 3) == 0) {
	    nbmin = 2;
	}
    } else if (cname && s_cmp(c2, "HE", 2, 2) == 0) {
	if (s_cmp(c3, "TRD", 3, 3) == 0) {
	    nbmin = 2;
	}
    } else if (sname && s_cmp(c2, "OR", 2, 2) == 0) {
	if (*(unsigned char *)c3 == 'G') {
	    if (s_cmp(c4, "QR", 2, 2) == 0 || s_cmp(c4, "RQ", 2, 2) == 0 
		    || s_cmp(c4, "LQ", 2, 2) == 0 || s_cmp(c4, "QL", 2, 2)
		     == 0 || s_cmp(c4, "HR", 2, 2) == 0 || s_cmp(c4, "TR", 
		    2, 2) == 0 || s_cmp(c4, "BR", 2, 2) == 0) {
		nbmin = 2;
	    }
	} else if (*(unsigned char *)c3 == 'M') {
	    if (s_cmp(c4, "QR", 2, 2) == 0 || s_cmp(c4, "RQ", 2, 2) == 0 
		    || s_cmp(c4, "LQ", 2, 2) == 0 || s_cmp(c4, "QL", 2, 2)
		     == 0 || s_cmp(c4, "HR", 2, 2) == 0 || s_cmp(c4, "TR", 
		    2, 2) == 0 || s_cmp(c4, "BR", 2, 2) == 0) {
		nbmin = 2;
	    }
	}
    } else if (cname && s_cmp(c2, "UN", 2, 2) == 0) {
	if (*(unsigned char *)c3 == 'G') {
	    if (s_cmp(c4, "QR", 2, 2) == 0 || s_cmp(c4, "RQ", 2, 2) == 0 
		    || s_cmp(c4, "LQ", 2, 2) == 0 || s_cmp(c4, "QL", 2, 2)
		     == 0 || s_cmp(c4, "HR", 2, 2) == 0 || s_cmp(c4, "TR", 
		    2, 2) == 0 || s_cmp(c4, "BR", 2, 2) == 0) {
		nbmin = 2;
	    }
	} else if (*(unsigned char *)c3 == 'M') {
	    if (s_cmp(c4, "QR", 2, 2) == 0 || s_cmp(c4, "RQ", 2, 2) == 0 
		    || s_cmp(c4, "LQ", 2, 2) == 0 || s_cmp(c4, "QL", 2, 2)
		     == 0 || s_cmp(c4, "HR", 2, 2) == 0 || s_cmp(c4, "TR", 
		    2, 2) == 0 || s_cmp(c4, "BR", 2, 2) == 0) {
		nbmin = 2;
	    }
	}
    }
    ret_val = nbmin;
    return ret_val;

L300:

/*     ISPEC = 3:  crossover point */

    nx = 0;
    if (s_cmp(c2, "GE", 2, 2) == 0) {
	if (s_cmp(c3, "QRF", 3, 3) == 0 || s_cmp(c3, "RQF", 3, 3) == 0 || 
		s_cmp(c3, "LQF", 3, 3) == 0 || s_cmp(c3, "QLF", 3, 3) == 
		0) {
	    if (sname) {
		nx = 128;
	    } else {
		nx = 128;
	    }
	} else if (s_cmp(c3, "HRD", 3, 3) == 0) {
	    if (sname) {
		nx = 128;
	    } else {
		nx = 128;
	    }
	} else if (s_cmp(c3, "BRD", 3, 3) == 0) {
	    if (sname) {
		nx = 128;
	    } else {
		nx = 128;
	    }
	}
    } else if (s_cmp(c2, "SY", 2, 2) == 0) {
	if (sname && s_cmp(c3, "TRD", 3, 3) == 0) {
	    nx = 1;
	}
    } else if (cname && s_cmp(c2, "HE", 2, 2) == 0) {
	if (s_cmp(c3, "TRD", 3, 3) == 0) {
	    nx = 1;
	}
    } else if (sname && s_cmp(c2, "OR", 2, 2) == 0) {
	if (*(unsigned char *)c3 == 'G') {
	    if (s_cmp(c4, "QR", 2, 2) == 0 || s_cmp(c4, "RQ", 2, 2) == 0 
		    || s_cmp(c4, "LQ", 2, 2) == 0 || s_cmp(c4, "QL", 2, 2)
		     == 0 || s_cmp(c4, "HR", 2, 2) == 0 || s_cmp(c4, "TR", 
		    2, 2) == 0 || s_cmp(c4, "BR", 2, 2) == 0) {
		nx = 128;
	    }
	}
    } else if (cname && s_cmp(c2, "UN", 2, 2) == 0) {
	if (*(unsigned char *)c3 == 'G') {
	    if (s_cmp(c4, "QR", 2, 2) == 0 || s_cmp(c4, "RQ", 2, 2) == 0 
		    || s_cmp(c4, "LQ", 2, 2) == 0 || s_cmp(c4, "QL", 2, 2)
		     == 0 || s_cmp(c4, "HR", 2, 2) == 0 || s_cmp(c4, "TR", 
		    2, 2) == 0 || s_cmp(c4, "BR", 2, 2) == 0) {
		nx = 128;
	    }
	}
    }
    ret_val = nx;
    return ret_val;

L400:

/*     ISPEC = 4:  number of shifts (used by xHSEQR) */

    ret_val = 6;
    return ret_val;

L500:

/*     ISPEC = 5:  minimum column dimension (not used) */

    ret_val = 2;
    return ret_val;

L600:

/*     ISPEC = 6:  crossover point for SVD (used by xGELSS and xGESVD) */

    ret_val = (integer) ((real) Dymola_min(*n1,*n2) * (float)1.6);
    return ret_val;

L700:

/*     ISPEC = 7:  number of processors (not used) */

    ret_val = 1;
    return ret_val;

L800:

/*     ISPEC = 8:  crossover point for multishift (used by xHSEQR) */

    ret_val = 50;
    return ret_val;

/*     End of ILAENV */

} /* ilaenv_ */

#undef lsame_
DYMOLA_STATIC logical lsame_(ca, cb, ca_len, cb_len)
const char *ca, *cb;
ftnlen ca_len;
ftnlen cb_len;
{
    /* System generated locals */
    logical ret_val;

    /* Local variables */
    integer inta, intb, zcode;


/*  -- LAPACK auxiliary routine (version 1.1) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
/*     Courant Institute, Argonne National Lab, and Rice University */
/*     February 29, 1992 */

/*     .. Scalar Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  LSAME returns .TRUE. if CA is the same letter as CB regardless of */
/*  case. */

/*  Arguments */
/*  ========= */

/*  CA      (input) CHARACTER*1 */
/*  CB      (input) CHARACTER*1 */
/*          CA and CB specify the single characters to be compared. */

/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test if the characters are equal */
#if 0
    return (*ca==*cb)||(toupperC(*(unsigned char*)ca)==toupperC(*(unsigned char*)cb)); 
#endif
    ret_val = *(unsigned char *)ca == *(unsigned char *)cb;
    if (ret_val) {
	return ret_val;
    }

/*     Now test for equivalence if both characters are alphabetic. */

    zcode = 'Z';

/*     Use 'Z' rather than 'A' so that ASCII can be detected on Prime */
/*     machines, on which ICHAR returns a value with bit 8 set. */
/*     ICHAR('A') on Prime machines returns 193 which is the same as */
/*     ICHAR('A') on an EBCDIC machine. */

    inta = *(unsigned char *)ca;
    intb = *(unsigned char *)cb;

    if (zcode == 90 || zcode == 122) {

/*        ASCII is assumed - ZCODE is the ASCII code of either lower o
r */
/*        upper case 'Z'. */

	if (inta >= 97 && inta <= 122) {
	    inta += -32;
	}
	if (intb >= 97 && intb <= 122) {
	    intb += -32;
	}

    } else if (zcode == 233 || zcode == 169) {

/*        EBCDIC is assumed - ZCODE is the EBCDIC code of either lower
 or */
/*        upper case 'Z'. */

	if ((inta >= 129 && inta <= 137) || (inta >= 145 && inta <= 153) || (inta >= 162 && inta <= 169)) {
	    inta += 64;
	}
	if ((intb >= 129 && intb <= 137) || (intb >= 145 && intb <= 153) || (intb >= 162 && intb <= 169)) {
	    intb += 64;
	}

    } else if (zcode == 218 || zcode == 250) {

/*        ASCII is assumed, on Prime machines - ZCODE is the ASCII cod
e */
/*        plus 128 of either lower or upper case 'Z'. */

	if (inta >= 225 && inta <= 250) {
	    inta += -32;
	}
	if (intb >= 225 && intb <= 250) {
	    intb += -32;
	}
    }
    ret_val = inta == intb;

/*     RETURN */

/*     End of LSAME */

    return ret_val;
} /* lsame_ */

/* Subroutine */ DYMOLA_STATIC int xerbla_(srname, info, srname_len)
char *srname;
integer *info;
ftnlen srname_len;
{

/*  -- LAPACK auxiliary routine (version 1.1) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
/*     Courant Institute, Argonne National Lab, and Rice University */
/*     February 29, 1992 */

/*     .. Scalar Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  XERBLA  is an error handler for the LAPACK routines. */
/*  It is called by an LAPACK routine if an input parameter has an */
/*  invalid value.  A message is printed and execution stops. */

/*  Installers may consider modifying the STOP statement in order to */
/*  call system-specific exception-handling facilities. */

/*  Arguments */
/*  ========= */

/*  SRNAME  (input) CHARACTER*6 */
/*          The name of the routine which called XERBLA. */

/*  INFO    (input) INTEGER */
/*          The position of the invalid parameter in the parameter list */
/*          of the calling routine. */

/*     .. Executable Statements .. */

/* CC      WRITE( *, FMT = 9999 )SRNAME, INFO */

/* CC      STOP */
    return 0;

/*CC 9999 FORMAT( ' ** On entry to ', A6, ' parameter number ', I2, ' had 
',*/
/* CC     $      'an illegal value' ) */

/*     End of XERBLA */

} /* xerbla_ */
#endif

/* Subroutine */ DYMOLA_STATIC int dgbfa_(abd, lda, n, ml, mu, ipvt, info)
doublereal *abd;
const integer*lda, *n, *ml, *mu;
integer *ipvt, *info;
{
    /* System generated locals */
    integer abd_dim1, abd_offset, i__1, i__2, i__3, i__4;

    /* Local variables */
    integer i__, j, k, l, m;
    doublereal t;
    integer i0, j0, j1, lm, mm, ju, jz, kp1, nm1;


/*     DGBFA FACTORS A DOUBLE PRECISION BAND MATRIX BY ELIMINATION. */

/*     DGBFA IS USUALLY CALLED BY DGBCO, BUT IT CAN BE CALLED */
/*     DIRECTLY WITH A SAVING IN TIME IF  RCOND  IS NOT NEEDED. */

/*     ON ENTRY */

/*        ABD     DOUBLE PRECISION(LDA, N) */
/*                CONTAINS THE MATRIX IN BAND STORAGE.  THE COLUMNS */
/*                OF THE MATRIX ARE STORED IN THE COLUMNS OF  ABD  AND */
/*                THE DIAGONALS OF THE MATRIX ARE STORED IN ROWS */
/*                ML+1 THROUGH 2*ML+MU+1 OF  ABD . */
/*                SEE THE COMMENTS BELOW FOR DETAILS. */

/*        LDA     INTEGER */
/*                THE LEADING DIMENSION OF THE ARRAY  ABD . */
/*                LDA MUST BE .GE. 2*ML + MU + 1 . */

/*        N       INTEGER */
/*                THE ORDER OF THE ORIGINAL MATRIX. */

/*        ML      INTEGER */
/*                NUMBER OF DIAGONALS BELOW THE MAIN DIAGONAL. */
/*                0 .LE. ML .LT. N . */

/*        MU      INTEGER */
/*                NUMBER OF DIAGONALS ABOVE THE MAIN DIAGONAL. */
/*                0 .LE. MU .LT. N . */
/*                MORE EFFICIENT IF  ML .LE. MU . */
/*     ON RETURN */

/*        ABD     AN UPPER TRIANGULAR MATRIX IN BAND STORAGE AND */
/*                THE MULTIPLIERS WHICH WERE USED TO OBTAIN IT. */
/*                THE FACTORIZATION CAN BE WRITTEN  A = L*U  WHERE */
/*                L  IS A PRODUCT OF PERMUTATION AND UNIT LOWER */
/*                TRIANGULAR MATRICES AND  U  IS UPPER TRIANGULAR. */

/*        IPVT    INTEGER(N) */
/*                AN INTEGER VECTOR OF PIVOT INDICES. */

/*        INFO    INTEGER */
/*                = 0  NORMAL VALUE. */
/*                = K  IF  U(K,K) .EQ. 0.0 .  THIS IS NOT AN ERROR */
/*                     CONDITION FOR THIS SUBROUTINE, BUT IT DOES */
/*                     INDICATE THAT DGBSL WILL DIVIDE BY ZERO IF */
/*                     CALLED.  USE  RCOND  IN DGBCO FOR A RELIABLE */
/*                     INDICATION OF SINGULARITY. */

/*     BAND STORAGE */

/*           IF  A  IS A BAND MATRIX, THE FOLLOWING PROGRAM SEGMENT */
/*           WILL SET UP THE INPUT. */

/*                   ML = (BAND WIDTH BELOW THE DIAGONAL) */
/*                   MU = (BAND WIDTH ABOVE THE DIAGONAL) */
/*                   M = ML + MU + 1 */
/*                   DO 20 J = 1, N */
/*                      I1 = MAX0(1, J-MU) */
/*                      I2 = MIN0(N, J+ML) */
/*                      DO 10 I = I1, I2 */
/*                         K = I - J + M */
/*                         ABD(K,J) = A(I,J) */
/*                10    CONTINUE */
/*                20 CONTINUE */

/*           THIS USES ROWS  ML+1  THROUGH  2*ML+MU+1  OF  ABD . */
/*           IN ADDITION, THE FIRST  ML  ROWS IN  ABD  ARE USED FOR */
/*           ELEMENTS GENERATED DURING THE TRIANGULARIZATION. */
/*           THE TOTAL NUMBER OF ROWS NEEDED IN  ABD  IS  2*ML+MU+1 . */
/*           THE  ML+MU BY ML+MU  UPPER LEFT TRIANGLE AND THE */
/*           ML BY ML  LOWER RIGHT TRIANGLE ARE NOT REFERENCED. */

/*     LINPACK. THIS VERSION DATED 08/14/78 . */
/*     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB. */

/*     SUBROUTINES AND FUNCTIONS */

/*     BLAS DAXPY,DSCAL,IDAMAX */
/*     FORTRAN MAX0,MIN0 */

/*     INTERNAL VARIABLES */



    /* Parameter adjustments */
    abd_dim1 = *lda;
    abd_offset = abd_dim1 + 1;
    abd -= abd_offset;
    --ipvt;

    /* Function Body */
    m = *ml + *mu + 1;
    *info = 0;

/*     ZERO INITIAL FILL-IN COLUMNS */

    j0 = *mu + 2;
    j1 = Dymola_min(*n,m) - 1;
    if (j1 < j0) {
	goto L30;
    }
    i__1 = j1;
    for (jz = j0; jz <= i__1; ++jz) {
	i0 = m + 1 - jz;
	i__2 = *ml;
	for (i__ = i0; i__ <= i__2; ++i__) {
	    abd[i__ + jz * abd_dim1] = 0.;
/* L10: */
	}
/* L20: */
    }
L30:
    jz = j1;
    ju = 0;

/*     GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING */

    nm1 = *n - 1;
    if (nm1 < 1) {
	goto L130;
    }
    i__1 = nm1;
    for (k = 1; k <= i__1; ++k) {
	kp1 = k + 1;

/*        ZERO NEXT FILL-IN COLUMN */

	++jz;
	if (jz > *n) {
	    goto L50;
	}
	if (*ml < 1) {
	    goto L50;
	}
	i__2 = *ml;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    abd[i__ + jz * abd_dim1] = 0.;
/* L40: */
	}
L50:

/*        FIND L = PIVOT INDEX */

/* Computing MIN */
	i__2 = *ml, i__3 = *n - k;
	lm = Dymola_min(i__2,i__3);
	i__2 = lm + 1;
	l = idamax_(&i__2, &abd[m + k * abd_dim1], &c__1) + m - 1;
	ipvt[k] = l + k - m;

/*        ZERO PIVOT IMPLIES THIS COLUMN ALREADY TRIANGULARIZED */

	if (abd[l + k * abd_dim1] == 0.) {
	    goto L100;
	}

/*           INTERCHANGE IF NECESSARY */

	if (l == m) {
	    goto L60;
	}
	t = abd[l + k * abd_dim1];
	abd[l + k * abd_dim1] = abd[m + k * abd_dim1];
	abd[m + k * abd_dim1] = t;
L60:

/*           COMPUTE MULTIPLIERS */

	t = -1. / abd[m + k * abd_dim1];
	dscal_(&lm, &t, &abd[m + 1 + k * abd_dim1], &c__1);

/*           ROW ELIMINATION WITH COLUMN INDEXING */

/* Computing MIN */
/* Computing MAX */
	i__3 = ju, i__4 = *mu + ipvt[k];
	i__2 = Dymola_max(i__3,i__4);
	ju = Dymola_min(i__2,*n);
	mm = m;
	if (ju < kp1) {
	    goto L90;
	}
	i__2 = ju;
	for (j = kp1; j <= i__2; ++j) {
	    --l;
	    --mm;
	    t = abd[l + j * abd_dim1];
	    if (l == mm) {
		goto L70;
	    }
	    abd[l + j * abd_dim1] = abd[mm + j * abd_dim1];
	    abd[mm + j * abd_dim1] = t;
L70:
	    daxpy_(&lm, &t, &abd[m + 1 + k * abd_dim1], &c__1, &abd[mm + 1 + 
		    j * abd_dim1], &c__1);
/* L80: */
	}
L90:
	goto L110;
L100:
	*info = k;
L110:
/* L120: */
	;
    }
L130:
    ipvt[*n] = *n;
    if (abd[m + *n * abd_dim1] == 0.) {
	*info = *n;
    }
    return 0;
} /* dgbfa_ */

/* Subroutine */ DYMOLA_STATIC int dgbsl_(abd, lda, n, ml, mu, ipvt, b, job)
doublereal *abd;
const integer*lda, *n, *ml, *mu;
integer *ipvt;
doublereal *b;
const integer *job;
{
    /* System generated locals */
    integer abd_dim1, abd_offset, i__1, i__2, i__3;

    /* Local variables */
    integer k, l, m;
    doublereal t;
    integer kb, la, lb, lm, nm1;


/*     DGBSL SOLVES THE DOUBLE PRECISION BAND SYSTEM */
/*     A * X = B  OR  TRANS(A) * X = B */
/*     USING THE FACTORS COMPUTED BY DGBCO OR DGBFA. */

/*     ON ENTRY */

/*        ABD     DOUBLE PRECISION(LDA, N) */
/*                THE OUTPUT FROM DGBCO OR DGBFA. */

/*        LDA     INTEGER */
/*                THE LEADING DIMENSION OF THE ARRAY  ABD . */

/*        N       INTEGER */
/*                THE ORDER OF THE ORIGINAL MATRIX. */

/*        ML      INTEGER */
/*                NUMBER OF DIAGONALS BELOW THE MAIN DIAGONAL. */

/*        MU      INTEGER */
/*                NUMBER OF DIAGONALS ABOVE THE MAIN DIAGONAL. */

/*        IPVT    INTEGER(N) */
/*                THE PIVOT VECTOR FROM DGBCO OR DGBFA. */

/*        B       DOUBLE PRECISION(N) */
/*                THE RIGHT HAND SIDE VECTOR. */

/*        JOB     INTEGER */
/*                = 0         TO SOLVE  A*X = B , */
/*                = NONZERO   TO SOLVE  TRANS(A)*X = B , WHERE */
/*                            TRANS(A)  IS THE TRANSPOSE. */

/*     ON RETURN */

/*        B       THE SOLUTION VECTOR  X . */

/*     ERROR CONDITION */

/*        A DIVISION BY ZERO WILL OCCUR IF THE INPUT FACTOR CONTAINS A */
/*        ZERO ON THE DIAGONAL.  TECHNICALLY THIS INDICATES SINGULARITY */
/*        BUT IT IS OFTEN CAUSED BY IMPROPER ARGUMENTS OR IMPROPER */
/*        SETTING OF LDA .  IT WILL NOT OCCUR IF THE SUBROUTINES ARE */
/*        CALLED CORRECTLY AND IF DGBCO HAS SET RCOND .GT. 0.0 */
/*        OR DGBFA HAS SET INFO .EQ. 0 . */

/*     TO COMPUTE  INVERSE(A) * C  WHERE  C  IS A MATRIX */
/*     WITH  P  COLUMNS */
/*           CALL DGBCO(ABD,LDA,N,ML,MU,IPVT,RCOND,Z) */
/*           IF (RCOND IS TOO SMALL) GO TO ... */
/*           DO 10 J = 1, P */
/*              CALL DGBSL(ABD,LDA,N,ML,MU,IPVT,C(1,J),0) */
/*        10 CONTINUE */

/*     LINPACK. THIS VERSION DATED 08/14/78 . */
/*     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB. */

/*     SUBROUTINES AND FUNCTIONS */

/*     BLAS DAXPY,DDOT */
/*     FORTRAN MIN0 */

/*     INTERNAL VARIABLES */


    /* Parameter adjustments */
    abd_dim1 = *lda;
    abd_offset = abd_dim1 + 1;
    abd -= abd_offset;
    --ipvt;
    --b;

    /* Function Body */
    m = *mu + *ml + 1;
    nm1 = *n - 1;
    if (*job != 0) {
	goto L50;
    }

/*        JOB = 0 , SOLVE  A * X = B */
/*        FIRST SOLVE L*Y = B */

    if (*ml == 0) {
	goto L30;
    }
    if (nm1 < 1) {
	goto L30;
    }
    i__1 = nm1;
    for (k = 1; k <= i__1; ++k) {
/* Computing MIN */
	i__2 = *ml, i__3 = *n - k;
	lm = Dymola_min(i__2,i__3);
	l = ipvt[k];
	t = b[l];
	if (l == k) {
	    goto L10;
	}
	b[l] = b[k];
	b[k] = t;
L10:
	daxpy_(&lm, &t, &abd[m + 1 + k * abd_dim1], &c__1, &b[k + 1], &c__1);
/* L20: */
    }
L30:

/*        NOW SOLVE  U*X = Y */

    i__1 = *n;
    for (kb = 1; kb <= i__1; ++kb) {
	k = *n + 1 - kb;
	b[k] /= abd[m + k * abd_dim1];
	lm = Dymola_min(k,m) - 1;
	la = m - lm;
	lb = k - lm;
	t = -b[k];
	daxpy_(&lm, &t, &abd[la + k * abd_dim1], &c__1, &b[lb], &c__1);
/* L40: */
    }
    goto L100;
L50:

/*        JOB = NONZERO, SOLVE  TRANS(A) * X = B */
/*        FIRST SOLVE  TRANS(U)*Y = B */

    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	lm = Dymola_min(k,m) - 1;
	la = m - lm;
	lb = k - lm;
	/* t = ddot_(&lm, &abd[la + k * abd_dim1], &c__1, &b[lb], &c__1); */
	t = ddot1_(&lm, &abd[la + k * abd_dim1], &b[lb]);
	b[k] = (b[k] - t) / abd[m + k * abd_dim1];
/* L60: */
    }

/*        NOW SOLVE TRANS(L)*X = Y */

    if (*ml == 0) {
	goto L90;
    }
    if (nm1 < 1) {
	goto L90;
    }
    i__1 = nm1;
    for (kb = 1; kb <= i__1; ++kb) {
	k = *n - kb;
/* Computing MIN */
	i__2 = *ml, i__3 = *n - k;
	lm = Dymola_min(i__2,i__3);
	/* b[k] += ddot_(&lm, &abd[m + 1 + k * abd_dim1], &c__1, &b[k + 1], &c__1); */
	b[k] += ddot1_(&lm, &abd[m + 1 + k * abd_dim1], &b[k + 1]);
	l = ipvt[k];
	if (l == k) {
	    goto L70;
	}
	t = b[l];
	b[l] = b[k];
	b[k] = t;
L70:
/* L80: */
	;
    }
L90:
L100:
    return 0;
} /* dgbsl_ */

/* Subroutine */ DYMOLA_STATIC int dgeco_(a, lda, n, ipvt, rcond, z__)
doublereal *a;
const integer*lda, *n;
integer *ipvt;
doublereal *rcond, *z__;
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;
    doublereal d__1, d__2;

    /* Builtin functions */

    /* Local variables */
    integer j, k, l;
    doublereal s, t;
    integer kb;
    doublereal ek, sm, wk;
    integer kp1;
    doublereal wkm;
    integer info;
    doublereal anorm;
    doublereal ynorm;


/*     DGECO FACTORS A DOUBLE PRECISION MATRIX BY GAUSSIAN ELIMINATION */
/*     AND ESTIMATES THE CONDITION OF THE MATRIX. */

/*     IF  RCOND  IS NOT NEEDED, DGEFA IS SLIGHTLY FASTER. */
/*     TO SOLVE  A*X = B , FOLLOW DGECO BY DGESL. */
/*     TO COMPUTE  INVERSE(A)*C , FOLLOW DGECO BY DGESL. */
/*     TO COMPUTE  DETERMINANT(A) , FOLLOW DGECO BY DGEDI. */
/*     TO COMPUTE  INVERSE(A) , FOLLOW DGECO BY DGEDI. */

/*     ON ENTRY */

/*        A       DOUBLE PRECISION(LDA, N) */
/*                THE MATRIX TO BE FACTORED. */

/*        LDA     INTEGER */
/*                THE LEADING DIMENSION OF THE ARRAY  A . */

/*        N       INTEGER */
/*                THE ORDER OF THE MATRIX  A . */

/*     ON RETURN */

/*        A       AN UPPER TRIANGULAR MATRIX AND THE MULTIPLIERS */
/*                WHICH WERE USED TO OBTAIN IT. */
/*                THE FACTORIZATION CAN BE WRITTEN  A = L*U  WHERE */
/*                L  IS A PRODUCT OF PERMUTATION AND UNIT LOWER */
/*                TRIANGULAR MATRICES AND  U  IS UPPER TRIANGULAR. */

/*        IPVT    INTEGER(N) */
/*                AN INTEGER VECTOR OF PIVOT INDICES. */

/*        RCOND   DOUBLE PRECISION */
/*                AN ESTIMATE OF THE RECIPROCAL CONDITION OF  A . */
/*                FOR THE SYSTEM  A*X = B , RELATIVE PERTURBATIONS */
/*                IN  A  AND  B  OF SIZE  EPSILON  MAY CAUSE */
/*                RELATIVE PERTURBATIONS IN  X  OF SIZE  EPSILON/RCOND . 
*/
/*                IF  RCOND  IS SO SMALL THAT THE LOGICAL EXPRESSION */
/*                           1.0 + RCOND .EQ. 1.0 */
/*                IS TRUE, THEN  A  MAY BE SINGULAR TO WORKING */
/*                PRECISION.  IN PARTICULAR,  RCOND  IS ZERO  IF */
/*                EXACT SINGULARITY IS DETECTED OR THE ESTIMATE */
/*                UNDERFLOWS. */

/*        Z       DOUBLE PRECISION(N) */
/*                A WORK VECTOR WHOSE CONTENTS ARE USUALLY UNIMPORTANT. */
/*                IF  A  IS CLOSE TO A SINGULAR MATRIX, THEN  Z  IS */
/*                AN APPROXIMATE NULL VECTOR IN THE SENSE THAT */
/*                NORM(A*Z) = RCOND*NORM(A)*NORM(Z) . */

/*     LINPACK. THIS VERSION DATED 08/14/78 . */
/*     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB. */

/*     SUBROUTINES AND FUNCTIONS */

/*     LINPACK DGEFA */
/*     BLAS DAXPY,DDOT,DSCAL,DASUM */
/*     FORTRAN DABS,DMAX1,DSIGN */

/*     INTERNAL VARIABLES */



/*     COMPUTE 1-NORM OF A */

    /* Parameter adjustments */
    --z__;
    --ipvt;
    a_dim1 = *lda;
    a_offset = a_dim1 + 1;
    a -= a_offset;

    /* Function Body */
    anorm = 0.;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
/* Computing MAX */
	d__1 = anorm, d__2 = dasum_(n, &a[j * a_dim1 + 1], &c__1);
	anorm = Dymola_max(d__1,d__2);
/* L10: */
    }

/*     FACTOR */

    dgefa_(&a[a_offset], lda, n, &ipvt[1], &info);

/*     RCOND = 1/(NORM(A)*(ESTIMATE OF NORM(INVERSE(A)))) . */
/*     ESTIMATE = NORM(Z)/NORM(Y) WHERE  A*Z = Y  AND  TRANS(A)*Y = E . */
/*     TRANS(A)  IS THE TRANSPOSE OF A .  THE COMPONENTS OF  E  ARE */
/*     CHOSEN TO CAUSE MAXIMUM LOCAL GROWTH IN THE ELEMENTS OF W  WHERE */
/*     TRANS(U)*W = E .  THE VECTORS ARE FREQUENTLY RESCALED TO AVOID */
/*     OVERFLOW. */

/*     SOLVE TRANS(U)*W = E */

    ek = 1.;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	z__[j] = 0.;
/* L20: */
    }
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	if (z__[k] != 0.) {
	    d__1 = -z__[k];
	    ek = d_sign(&ek, &d__1);
	}
	if ((d__1 = ek - z__[k], Dymola_abs(d__1)) <= (d__2 = a[k + k * a_dim1], Dymola_abs(
		d__2))) {
	    goto L30;
	}
	s = (d__1 = a[k + k * a_dim1], Dymola_abs(d__1)) / (d__2 = ek - z__[k], Dymola_abs(
		d__2));
	dscal_(n, &s, &z__[1], &c__1);
	ek = s * ek;
L30:
	wk = ek - z__[k];
	wkm = -ek - z__[k];
	s = Dymola_abs(wk);
	sm = Dymola_abs(wkm);
	if (a[k + k * a_dim1] == 0.) {
	    goto L40;
	}
	wk /= a[k + k * a_dim1];
	wkm /= a[k + k * a_dim1];
	goto L50;
L40:
	wk = 1.;
	wkm = 1.;
L50:
	kp1 = k + 1;
	if (kp1 > *n) {
	    goto L90;
	}
	i__2 = *n;
	for (j = kp1; j <= i__2; ++j) {
	    sm += (d__1 = z__[j] + wkm * a[k + j * a_dim1], Dymola_abs(d__1));
	    z__[j] += wk * a[k + j * a_dim1];
	    s += (d__1 = z__[j], Dymola_abs(d__1));
/* L60: */
	}
	if (s >= sm) {
	    goto L80;
	}
	t = wkm - wk;
	wk = wkm;
	i__2 = *n;
	for (j = kp1; j <= i__2; ++j) {
	    z__[j] += t * a[k + j * a_dim1];
/* L70: */
	}
L80:
L90:
	z__[k] = wk;
/* L100: */
    }
    s = 1. / dasum_(n, &z__[1], &c__1);
    dscal_(n, &s, &z__[1], &c__1);

/*     SOLVE TRANS(L)*Y = W */

    i__1 = *n;
    for (kb = 1; kb <= i__1; ++kb) {
	k = *n + 1 - kb;
	if (k < *n) {
	    i__2 = *n - k;
	    /* z__[k] += ddot_(&i__2, &a[k + 1 + k * a_dim1], &c__1, &z__[k + 1], &c__1); */
	    z__[k] += ddot1_(&i__2, &a[k + 1 + k * a_dim1], &z__[k + 1]);
	}
	if ((d__1 = z__[k], Dymola_abs(d__1)) <= 1.) {
	    goto L110;
	}
	s = 1. / (d__1 = z__[k], Dymola_abs(d__1));
	dscal_(n, &s, &z__[1], &c__1);
L110:
	l = ipvt[k];
	t = z__[l];
	z__[l] = z__[k];
	z__[k] = t;
/* L120: */
    }
    s = 1. / dasum_(n, &z__[1], &c__1);
    dscal_(n, &s, &z__[1], &c__1);

    ynorm = 1.;

/*     SOLVE L*V = Y */

    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	l = ipvt[k];
	t = z__[l];
	z__[l] = z__[k];
	z__[k] = t;
	if (k < *n) {
	    i__2 = *n - k;
	    daxpy_(&i__2, &t, &a[k + 1 + k * a_dim1], &c__1, &z__[k + 1], &
		    c__1);
	}
	if ((d__1 = z__[k], Dymola_abs(d__1)) <= 1.) {
	    goto L130;
	}
	s = 1. / (d__1 = z__[k], Dymola_abs(d__1));
	dscal_(n, &s, &z__[1], &c__1);
	ynorm = s * ynorm;
L130:
/* L140: */
	;
    }
    s = 1. / dasum_(n, &z__[1], &c__1);
    dscal_(n, &s, &z__[1], &c__1);
    ynorm = s * ynorm;

/*     SOLVE  U*Z = V */

    i__1 = *n;
    for (kb = 1; kb <= i__1; ++kb) {
	k = *n + 1 - kb;
	if ((d__1 = z__[k], Dymola_abs(d__1)) <= (d__2 = a[k + k * a_dim1], Dymola_abs(d__2)
		)) {
	    goto L150;
	}
	s = (d__1 = a[k + k * a_dim1], Dymola_abs(d__1)) / (d__2 = z__[k], Dymola_abs(d__2))
		;
	dscal_(n, &s, &z__[1], &c__1);
	ynorm = s * ynorm;
L150:
	if (a[k + k * a_dim1] != 0.) {
	    z__[k] /= a[k + k * a_dim1];
	}
	if (a[k + k * a_dim1] == 0.) {
	    z__[k] = 1.;
	}
	t = -z__[k];
	i__2 = k - 1;
	daxpy_(&i__2, &t, &a[k * a_dim1 + 1], &c__1, &z__[1], &c__1);
/* L160: */
    }
/*     MAKE ZNORM = 1.0 */
    s = 1. / dasum_(n, &z__[1], &c__1);
    dscal_(n, &s, &z__[1], &c__1);
    ynorm = s * ynorm;

    if (anorm != 0.) {
	*rcond = ynorm / anorm;
    }
    if (anorm == 0.) {
	*rcond = 0.;
    }
    return 0;
} /* dgeco_ */

/* Subroutine */ DYMOLA_STATIC int dgefa_(a, lda, n, ipvt, info)
doublereal *a;
const integer*lda, *n;
integer *ipvt, *info;
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;

    /* Local variables */
    integer j, k, l;
    doublereal t;
    integer kp1, nm1;


/*     DGEFA FACTORS A DOUBLE PRECISION MATRIX BY GAUSSIAN ELIMINATION. */

/*     DGEFA IS USUALLY CALLED BY DGECO, BUT IT CAN BE CALLED */
/*     DIRECTLY WITH A SAVING IN TIME IF  RCOND  IS NOT NEEDED. */
/*     (TIME FOR DGECO) = (1 + 9/N)*(TIME FOR DGEFA) . */

/*     ON ENTRY */

/*        A       DOUBLE PRECISION(LDA, N) */
/*                THE MATRIX TO BE FACTORED. */

/*        LDA     INTEGER */
/*                THE LEADING DIMENSION OF THE ARRAY  A . */

/*        N       INTEGER */
/*                THE ORDER OF THE MATRIX  A . */

/*     ON RETURN */

/*        A       AN UPPER TRIANGULAR MATRIX AND THE MULTIPLIERS */
/*                WHICH WERE USED TO OBTAIN IT. */
/*                THE FACTORIZATION CAN BE WRITTEN  A = L*U  WHERE */
/*                L  IS A PRODUCT OF PERMUTATION AND UNIT LOWER */
/*                TRIANGULAR MATRICES AND  U  IS UPPER TRIANGULAR. */

/*        IPVT    INTEGER(N) */
/*                AN INTEGER VECTOR OF PIVOT INDICES. */

/*        INFO    INTEGER */
/*                = 0  NORMAL VALUE. */
/*                = K  IF  U(K,K) .EQ. 0.0 .  THIS IS NOT AN ERROR */
/*                     CONDITION FOR THIS SUBROUTINE, BUT IT DOES */
/*                     INDICATE THAT DGESL OR DGEDI WILL DIVIDE BY ZERO */
/*                     IF CALLED.  USE  RCOND  IN DGECO FOR A RELIABLE */
/*                     INDICATION OF SINGULARITY. */

/*     LINPACK. THIS VERSION DATED 08/14/78 . */
/*     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB. */

/*     SUBROUTINES AND FUNCTIONS */

/*     BLAS DAXPY,DSCAL,IDAMAX */

/*     INTERNAL VARIABLES */



/*     GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING */

    /* Parameter adjustments */
    --ipvt;
    a_dim1 = *lda;
    a_offset = a_dim1 + 1;
    a -= a_offset;

    /* Function Body */
    *info = 0;
    nm1 = *n - 1;
    if (nm1 < 1) {
	goto L70;
    }
    i__1 = nm1;
    for (k = 1; k <= i__1; ++k) {
	kp1 = k + 1;

/*        FIND L = PIVOT INDEX */

	i__2 = *n - k + 1;
	l = idamax_(&i__2, &a[k + k * a_dim1], &c__1) + k - 1;
	ipvt[k] = l;

/*        ZERO PIVOT IMPLIES THIS COLUMN ALREADY TRIANGULARIZED */

	if (a[l + k * a_dim1] == 0.) {
	    goto L40;
	}

/*           INTERCHANGE IF NECESSARY */

	if (l == k) {
	    goto L10;
	}
	t = a[l + k * a_dim1];
	a[l + k * a_dim1] = a[k + k * a_dim1];
	a[k + k * a_dim1] = t;
L10:

/*           COMPUTE MULTIPLIERS */

	t = -1. / a[k + k * a_dim1];
	i__2 = *n - k;
	dscal_(&i__2, &t, &a[k + 1 + k * a_dim1], &c__1);

/*           ROW ELIMINATION WITH COLUMN INDEXING */

	i__2 = *n;
	for (j = kp1; j <= i__2; ++j) {
	    t = a[l + j * a_dim1];
	    if (l == k) {
		goto L20;
	    }
	    a[l + j * a_dim1] = a[k + j * a_dim1];
	    a[k + j * a_dim1] = t;
L20:
	    i__3 = *n - k;
	    /* daxpy_(&i__3, &t, &a[k + 1 + k * a_dim1], &c__1, &a[k + 1 + j * a_dim1], &c__1); */
	    daxpy1_(&i__3, &t, &a[k + 1 + k * a_dim1], &a[k + 1 + j * a_dim1]);
/* L30: */
	}
	goto L50;
L40:
	*info = k;
L50:
/* L60: */
	;
    }
L70:
    ipvt[*n] = *n;
    if (a[*n + *n * a_dim1] == 0.) {
	*info = *n;
    }
    return 0;
} /* dgefa_ */

/* Subroutine */ DYMOLA_STATIC int dgesl_(a, lda, n, ipvt, b, job)
doublereal *a;
const integer*lda, *n, *ipvt;
doublereal *b;
const integer *job;
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;

    /* Local variables */
    integer k, l;
    doublereal t;
    integer kb, nm1;


/*     DGESL SOLVES THE DOUBLE PRECISION SYSTEM */
/*     A * X = B  OR  TRANS(A) * X = B */
/*     USING THE FACTORS COMPUTED BY DGECO OR DGEFA. */

/*     ON ENTRY */

/*        A       DOUBLE PRECISION(LDA, N) */
/*                THE OUTPUT FROM DGECO OR DGEFA. */

/*        LDA     INTEGER */
/*                THE LEADING DIMENSION OF THE ARRAY  A . */

/*        N       INTEGER */
/*                THE ORDER OF THE MATRIX  A . */

/*        IPVT    INTEGER(N) */
/*                THE PIVOT VECTOR FROM DGECO OR DGEFA. */

/*        B       DOUBLE PRECISION(N) */
/*                THE RIGHT HAND SIDE VECTOR. */

/*        JOB     INTEGER */
/*                = 0         TO SOLVE  A*X = B , */
/*                = NONZERO   TO SOLVE  TRANS(A)*X = B  WHERE */
/*                            TRANS(A)  IS THE TRANSPOSE. */

/*     ON RETURN */

/*        B       THE SOLUTION VECTOR  X . */

/*     ERROR CONDITION */

/*        A DIVISION BY ZERO WILL OCCUR IF THE INPUT FACTOR CONTAINS A */
/*        ZERO ON THE DIAGONAL.  TECHNICALLY THIS INDICATES SINGULARITY */
/*        BUT IT IS OFTEN CAUSED BY IMPROPER ARGUMENTS OR IMPROPER */
/*        SETTING OF LDA .  IT WILL NOT OCCUR IF THE SUBROUTINES ARE */
/*        CALLED CORRECTLY AND IF DGECO HAS SET RCOND .GT. 0.0 */
/*        OR DGEFA HAS SET INFO .EQ. 0 . */

/*     TO COMPUTE  INVERSE(A) * C  WHERE  C  IS A MATRIX */
/*     WITH  P  COLUMNS */
/*           CALL DGECO(A,LDA,N,IPVT,RCOND,Z) */
/*           IF (RCOND IS TOO SMALL) GO TO ... */
/*           DO 10 J = 1, P */
/*              CALL DGESL(A,LDA,N,IPVT,C(1,J),0) */
/*        10 CONTINUE */

/*     LINPACK. THIS VERSION DATED 08/14/78 . */
/*     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB. */

/*     SUBROUTINES AND FUNCTIONS */

/*     BLAS DAXPY,DDOT */

/*     INTERNAL VARIABLES */


    /* Parameter adjustments */
    --b;
    a_dim1 = *lda;
    a_offset = a_dim1 + 1;
    a -= a_offset;
    --ipvt;

    /* Function Body */
    nm1 = *n - 1;
    if (*job != 0) {
	goto L50;
    }

/*        JOB = 0 , SOLVE  A * X = B */
/*        FIRST SOLVE  L*Y = B */

    if (nm1 < 1) {
	goto L30;
    }
    i__1 = nm1;
    for (k = 1; k <= i__1; ++k) {
	l = ipvt[k];
	t = b[l];
	if (l == k) {
	    goto L10;
	}
	b[l] = b[k];
	b[k] = t;
L10:
	i__2 = *n - k;
	/* daxpy_(&i__2, &t, &a[k + 1 + k * a_dim1], &c__1, &b[k + 1], &c__1); */
	daxpy1_(&i__2, &t, &a[k + 1 + k * a_dim1], &b[k + 1]);
/* L20: */
    }
L30:

/*        NOW SOLVE  U*X = Y */

    i__1 = *n;
    for (kb = 1; kb <= i__1; ++kb) {
	k = *n + 1 - kb;
	b[k] /= a[k + k * a_dim1];
	t = -b[k];
	i__2 = k - 1;
	/* daxpy_(&i__2, &t, &a[k * a_dim1 + 1], &c__1, &b[1], &c__1); */
	daxpy1_(&i__2, &t, &a[k * a_dim1 + 1], &b[1]);
/* L40: */
    }
    goto L100;
L50:

/*        JOB = NONZERO, SOLVE  TRANS(A) * X = B */
/*        FIRST SOLVE  TRANS(U)*Y = B */

    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	i__2 = k - 1;
	/* t = ddot_(&i__2, &a[k * a_dim1 + 1], &c__1, &b[1], &c__1); */
	t = ddot1_(&i__2, &a[k * a_dim1 + 1], &b[1]);
	b[k] = (b[k] - t) / a[k + k * a_dim1];
/* L60: */
    }

/*        NOW SOLVE TRANS(L)*X = Y */

    if (nm1 < 1) {
	goto L90;
    }
    i__1 = nm1;
    for (kb = 1; kb <= i__1; ++kb) {
	k = *n - kb;
	i__2 = *n - k;
	/* b[k] += ddot_(&i__2, &a[k + 1 + k * a_dim1], &c__1, &b[k + 1], &c__1); */
	b[k] += ddot1_(&i__2, &a[k + 1 + k * a_dim1], &b[k + 1]);
	l = ipvt[k];
	if (l == k) {
	    goto L70;
	}
	t = b[l];
	b[l] = b[k];
	b[k] = t;
L70:
/* L80: */
	;
    }
L90:
L100:
    return 0;
} /* dgesl_ */

/* Subroutine */ DYMOLA_STATIC int dymres_(a, lda, n, b, ierr)
doublereal *a;
const integer*lda, *n;
doublereal *b;
integer *ierr;
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;

    /* Local variables */
    integer i__, j;


/*     Reset matrix A and vector b (columnwise) */


    /* Parameter adjustments */
    --b;
    a_dim1 = *lda;
    a_offset = a_dim1 + 1;
    a -= a_offset;

    /* Function Body */
    if (*n > *lda) {
	dymosimmessageSev_(1, "--- Too large matrix. MaxMat in dsmodel.c needs to be increased.", 64);
	*ierr = 1;
    } else {
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    i__2 = *n;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		a[i__ + j * a_dim1] = 0.;
/* L10: */
	    }
	    b[j] = 0.;
/* L20: */
	}
	*ierr = 0;
    }
    return 0;
} /* dymres_ */



/* Subroutine */ DYMOLA_STATIC int dymsol_(a, lda, n, b, dwork, iwork, ierr)
doublereal *a;
const integer*lda, *n;
doublereal *b, *dwork;
integer *iwork, *ierr;
{
    /* System generated locals */
    integer a_dim1, a_offset;

    /* Local variables */
    doublereal rcond;



/*     Solves the equation */
/*       A*x = b */

    /* Parameter adjustments */
    --iwork;
    --dwork;
    --b;
    a_dim1 = *lda;
    a_offset = a_dim1 + 1;
    a -= a_offset;

    /* Function Body */
    *ierr = 0;
    dgeco_(&a[a_offset], lda, n, &iwork[1], &rcond, &dwork[1]);
    if (rcond * .01 + 1. <= 1.) {
	dymosimmessageSev_(1, "--- Error in Dymola model:", 26);
	dymosimmessageSev_(1, "---    System of equations is singular.", 39);
	dymosimmessagedoubleSev_(1, "---   Condition number = ", &rcond, 25);
	*ierr = 1;
    } else {
	dgesl_(&a[a_offset], lda, n, &iwork[1], &b[1], &c__0);
    }
    return 0;
} /* dymsol_ */


#if 0
/* Subroutine */ int dymsos_(a, lda, n, b, dwork, iwork, ierr)
doublereal *a;
const integer*lda, *n;
doublereal *b, *dwork;
integer *iwork, *ierr;
{
    /* System generated locals */
    integer a_dim1, a_offset;

    /* Local variables */
    doublereal rcond;


/*     Solves the equation */
/*       A*x = b */
/*     where A is symmetric. */

/*      CALL DSICO(A, LDA, N, IWORK, RCOND, DWORK) */
    /* Parameter adjustments */
    --iwork;
    --dwork;
    --b;
    a_dim1 = *lda;
    a_offset = a_dim1 + 1;
    a -= a_offset;

    /* Function Body */
    if (rcond * .01 + 1. == 1.) {
	dymosimmessage_("--- Error in Dymola model:", 26);
	dymosimmessage_("---    System of equations is singular.", 39);
	dymosimmessagedouble_("---   Condition number = ", &rcond, 25);
	*ierr = 1;
    } else {
/*        CALL DSISL(A, LDA, N, IWORK, B) */
    }
    return 0;
} /* dymsos_ */
#endif



/* Subroutine */ DYMOLA_STATIC int dymlin_(infrev, n, sol, res, jac, ljac, dwork, iwork, 
	time, event, printpriority, sysnr, ierr)
integer *infrev;
const integer*n;
doublereal *sol, *res, *jac;
integer *ljac;
doublereal *dwork;
integer *iwork;
doublereal *time;
integer *event, *printpriority, *sysnr, *ierr;
{
    /* System generated locals */
    integer jac_dim1, jac_offset, i__1;

    /* Local variables */
    integer i__, icol;


/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */

/* Solve linear system "A*x = b", where code to compute */
/* "Res = A*x - b" is provided (x=Sol) */

/* InfRev: IN, OUT, INTEGER, reverse communication parameter */
/*         InfRev is an integer variable that must be set by the user to 
*/
/*         a non-positive value initially, in order to indicate */
/*         the start of the calculation. The value of InfRev on each */
/*         return indicates the reason for the return. If it is 1,3 or 4, 
*/
/*         then the user's program has to provide new values of Res, */
/*         and then recall the subroutine with InfRev still set to this */
/*         number. */
/*         InfRev on return */
/*            = 0  : Solver terminated */
/*            = 1,3: Calculate Res */

/* N     : IN, INTEGER, N.gt.0 */
/*         Number of residues and solution variables. */

/* Sol   : OUT, DOUBLE(N) */
/*         On output QSol contains the values for which Res has to be */
/*         computed in the next iteration step (InfRev=1,3) or the */
/*         final solution (InfRev=0). */

/* Res   : IN, DOUBLE(N) */
/*         Res has to contain the residues evaluated at Sol. */
/*         (INFREV=1,3) */

/* Jac   : OUT, DOUBLE(N,N) */
/*         Jac contains the A-matrix after successful termination. */

/* dwork : OUT, DOUBLE(ldwork) */
/*         work array of length ldwork. */
/*         ldwork must not be less than N*N+5*N */

/* iwork : IN, INTEGER (liwork) */
/*         work array of length liwork. */
/*         liwork must not be less than 2*N. */

/*  Time    (input) double */
/*          Actualt time instant (used in error message) */

/* Event : IN, INTEGER */
/*         = 0: DYMSOL called during continuous integration. */
/*         = 1: DYMSOL called at an event instant. */

/* Printpriority: IN, INTEGER */
/*         If PrintPriority >= 80 and Event=1, a message is printed */
/*         if A is singular but with consistent equations. */

/* SysNr : IN, INTEGER */
/*         The number of the linear equation system (used to identify the 
*/
/*         system of equations in the Dymola model */

/* ierr  : OUT, INTEGER */
/*         = 0: successful exit */
/*        = 1: an error occured; an error message has been already printed
.*/

/* Author: M. Otter, DLR Oberpfaffenhofen, March 1998. */

/* ====================================================================== 
*/


    /* Parameter adjustments */
    --iwork;
    --dwork;
    --res;
    --sol;
    jac_dim1 = *ljac;
    jac_offset = jac_dim1 + 1;
    jac -= jac_offset;

    /* Function Body */
    *ierr = 0;
    if (*infrev < 0) {
/*        -- First call. Compute right hand side vector b */
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    sol[i__] = 0.;
/* L10: */
	}
	*infrev = 1;
	return 0;

    } else if (*infrev == 1) {
/*        -- Right hand side vector computed in last step. */
/*        -- Store it in double work array. */
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    dwork[i__] = -res[i__];
/* L20: */
	}

/*        -- Initialize computation of all the columns of the A-Matrix
 */
	iwork[1] = 1;
	sol[1] = 1.;
	i__1 = *n;
	for (i__ = 2; i__ <= i__1; ++i__) {
	    sol[i__] = 0.;
/* L30: */
	}
	*infrev = 3;
	return 0;

    } else if (*infrev == 3) {
/*        -- Store computed column in A-Matrix */
	icol = iwork[1];
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    jac[i__ + icol * jac_dim1] = res[i__] + dwork[i__];
/* L40: */
	}

/*        -- Compute next column */
	++icol;
	iwork[1] = icol;
	if (icol <= *n) {
	    i__1 = *n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		sol[i__] = 0.;
/* L50: */
	    }
	    sol[icol] = 1.;
	    *infrev = 3;
	    return 0;
	}

/*        -- Solve linear system of equations */
	dymli1_(sysnr, &c__0, &jac[jac_offset], ljac, n, &dwork[1], time, 
		event, printpriority, &dwork[*n + 1], &iwork[1], ierr);
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    sol[i__] = dwork[i__];
/* L60: */
	}
	*infrev = 0;
	return 0;
    } else {
/*        -- Wrong input argument */
	DymosimMessageIntSev_(1, "... Wrong input argument of DYMLIN, InvRev = ", 
		infrev, 45);
	*infrev = 0;
	*ierr = 1;
	return 0;
    }
    return 0;
} /* dymlin_ */




/* Subroutine */ LIBDS_API int dymli1_(sysnr, fact, a, lda, n, b, time, event, 
	printpriority, dwork, iwork, ierr)
integer *sysnr;
const integer *fact;
doublereal *a;
const integer*lda, *n;
doublereal *b, *time;
integer *event, *printpriority;
doublereal *dwork;
integer *iwork, *ierr;
{
    /* Initialized data */

    /* System generated locals */
    integer a_dim1, a_offset, i__1;

    /* Local variables */
    integer i__;
    doublereal eps;
    integer info;

/*     .. */

/*  Purpose */
/*  ======= */

/*  DYMSOL computes the solution to a real system of linear equations */
/*     A * X = B, */
/*  where A is an N-by-N matrix and X and B are N-by-1 matrices. */

/*  If FACT = 0, the LU decomposition with partial pivoting and row */
/*  interchanges is used to factor A as */
/*     A = P * L * U, */
/*  where P is a permutation matrix, L is unit lower triangular, and U is 
*/
/*  upper triangular.  If A is nonsingular, then the factored form of A */
/*  is used to solve the system of equations A * X = B. */
/*  If FACT = 1 or if A is singular, the QR decomposition with */
/*  column pivoting is used to factor A as */
/*     A*P = Q * R */
/*  where P is a permutation matrix, Q is an orthogonal matrix, and R is 
*/
/*  upper triangular. The rank of A is determined and the compatibility */
/*  of the linear system is checked. If the system is compatible a */
/*  minimum norm solution is computed using the QR decomposition of A. */

/*  If A is singular but with consistent equations, a message is */
/*  printed at event instants and if the print priority is >= 80. */


/*  Arguments */
/*  ========= */

/*  SysNr   (input) INTEGER */
/*          The number of the linear equation system (used to identify the
 */
/*          system of equations in the Dymola model) */

/*  FACT    (input) INTEGER */
/*          Option to specify the solution method. */
/*          FACT = 0, use the LU-decomposition to solve AX = B or */
/*                    the QR-decomposition if A is singular. */
/*          FACT = 1, use the QR-decomposition to solve AX = B. */

/*  A       (input/output) DOUBLE PRECISION array, dimension (N,N) */
/*          On entry, the N-by-N coefficient matrix A. */
/*          On exit, information on LU or QR decomposition of A. */

/*  LDA     (input) INTEGER */
/*          The declared leading dimension of A. LDA >= 1. */

/*  N       (input) INTEGER */
/*          The number of linear equations, i.e., the order of the */
/*          matrix A.  1 <= N <= LDA. */

/*  B       (input/output) DOUBLE PRECISION array, dimension (N) */
/*          On entry, the N-by-1 matrix of right hand side matrix B. */
/*          On exit, if INFO = 0, the N-by-1 solution matrix X. */

/*  Time    (input) double */
/*          Actual time instant (used in error message) */

/*  Event   (input) INTEGER */
/*          = 0: DYMSOL called during continuous integration. */
/*          = 1: DYMSOL called at an event instant. */

/*  PrintPriority (input) INTEGER */
/*          If PrintPriority >= 80 and Event=1, a message is printed */
/*          if A is singular but with consistent equations. */

/*  WORK    (output) REAL work array, dimension (LWORK), where */
/*           LWORK = N*N+4*N if FACT = 0 */
/*           LWORK = 4*N     if FACT = 1 */

/*  IWORK   (output) INTEGER work array, dimension (2*N) */

/*  IERR    (output) INTEGER */
/*          = 0:  successful exit */
/*                (A regular or singular with consistent equations) */
/*          = 1:  error exit (A singular, AX=B not compatible). */
/*                a message has already been computed. */

/* Author: M. Otter, DLR Oberpfaffenhofen, Feb. 1998. */

/* ====================================================================== 
*/

		eps=DBL_EPSILON*1e7;




    /* Parameter adjustments */
    --iwork;
    --dwork;
    --b;
    a_dim1 = *lda;
    a_offset = a_dim1 + 1;
    a -= a_offset;

    /* Function Body */

    *ierr = 0;
/* Inquire machine epsilon */

    dymlqr_(fact, &a[a_offset], lda, n, &b[1], &eps, &iwork[1], &dwork[1], &
	    iwork[*n + 1], &info);
    if (info == 1 && *event == 1 && (*printpriority & (1<<12))) {
	dymosimmessagedoubleSev_(1, "... System of equations is singular at Time = "
		, time, 46);
	DymosimMessageIntSev_(1, "...    Linear system of equations number = ", 
		sysnr, 43);
	dymosimmessageSev_(1, "...    Variables X(i) which cannot be uniquely computed:", 56);
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    if (iwork[i__] == 1) {
		DymosimMessageIntSev_(1, "...     ", &i__, 8);
	    }
/* L10: */
	}
	dymosimmessageSev_(1, "... Simulation continues with minimum norm solution."
		, 52);
    } else if (info == 2) {
	dymosimmessageSev_(1, " ", 1);
	dymosimmessagedoubleSev_(1, "LINEAR SYSTEM OF EQUATIONS IS SINGULAR AT TIME = ", time, 49);
	DymosimMessageIntSev_(1, "   linear system of equations number = ", sysnr, 
		39);
	dymosimmessageSev_(1, " ", 1);
	*ierr = 1;
    }

    return 0;
} /* dymli1_ */
LIBDS_API int dymlqr3_(integer *fact, doublereal*a, const integer*lda, const integer*n, doublereal*b, doublereal*tol, 
                    integer*pivots, doublereal*work, integer*iwork, integer*factor,integer*info,integer*sysnr,
					const char*const*varnames,integer*fEvent);

/* Subroutine */ LIBDS_API int dymli4_(sysnr, fact, a, lda, n, b, time, event, 
	printpriority, dwork, iwork, factor,varnames,ierr,fEvent)
integer *sysnr, *fact;
doublereal *a;
const integer*lda, *n;
doublereal *b, *time;
integer *event, *printpriority;
doublereal *dwork;
integer *iwork, *factor;
const char*const*varnames;
integer *ierr;
int *fEvent;
{
    /* Initialized data */

    /* System generated locals */
    /* integer a_dim1, a_offset;*/
	integer i__1;

    /* Local variables */
    integer i__;
    doublereal eps;
    integer info;

/*     .. */

/*  Purpose */
/*  ======= */

/*  DYMSOL computes the solution to a real system of linear equations */
/*     A * X = B, */
/*  where A is an N-by-N matrix and X and B are N-by-1 matrices. */

/*  If FACT = 0, the LU decomposition with partial pivoting and row */
/*  interchanges is used to factor A as */
/*     A = P * L * U, */
/*  where P is a permutation matrix, L is unit lower triangular, and U is 
*/
/*  upper triangular.  If A is nonsingular, then the factored form of A */
/*  is used to solve the system of equations A * X = B. */
/*  If FACT = 1 or if A is singular, the QR decomposition with */
/*  column pivoting is used to factor A as */
/*     A*P = Q * R */
/*  where P is a permutation matrix, Q is an orthogonal matrix, and R is 
*/
/*  upper triangular. The rank of A is determined and the compatibility */
/*  of the linear system is checked. If the system is compatible a */
/*  minimum norm solution is computed using the QR decomposition of A. */

/*  If A is singular but with consistent equations, a message is */
/*  printed at event instants and if the print priority is >= 80. */


/*  Arguments */
/*  ========= */

/*  SysNr   (input) INTEGER */
/*          The number of the linear equation system (used to identify the
 */
/*          system of equations in the Dymola model) */

/*  FACT    (input) INTEGER */
/*          Option to specify the solution method. */
/*          FACT&1 = 0, use the LU-decomposition to solve AX = B or */
/*                    the QR-decomposition if A is singular. */
/*          FACT&1 = 1, use the QR-decomposition to solve AX = B. */
/*          FACT&2 = 0, do not perform row and column scaling */
/*          FACT&2 = 2, perform row and column before LU/QR-decomposition */

/*  A       (input) DOUBLE PRECISION array, dimension (N,N) */
/*          On entry, the N-by-N coefficient matrix A. */

/*  LDA     (input) INTEGER */
/*          The declared leading dimension of A. LDA >= 1. */

/*  N       (input) INTEGER */
/*          The number of linear equations, i.e., the order of the */
/*          matrix A.  1 <= N <= LDA. */

/*  B       (input/output) DOUBLE PRECISION array, dimension (N) */
/*          On entry, the N-by-1 matrix of right hand side matrix B. */
/*          On exit, if INFO = 0, the N-by-1 solution matrix X. */

/*  Time    (input) double */
/*          Actual time instant (used in error message) */

/*  Event   (input) INTEGER */
/*          = 0: DYMSOL called during continuous integration. */
/*          = 1: DYMSOL called at an event instant. */

/*  PrintPriority (input) INTEGER */
/*          If PrintPriority >= 80 and Event=1, a message is printed */
/*          if A is singular but with consistent equations. */

/*  WORK    (output) REAL work array, dimension (LWORK), where */
/*           LWORK = N*N+6*N */

/*  IWORK   (output) INTEGER work array, dimension (2*N+1) */

  /* Factor (input/output)
     Factor&3     = 0: Factorize matrix
                  = 1: Use existing LU-decomposition
                  = 2: Use existing QR-decomposition
     Factor&256   =256: Use row-equlibrization
           &512   =512: Use column-equilibrization
             Updated after each call.
             The LU-factorization is stored in the first N*N elements 
             of WORK and the first N elements of IWORK.
             The QR-factorization is stored in the first N*N+N elements
             of WORK and the first N+1 elements of IWORK.
	     The equilib information is stored in the 2*N last elements
	     of WORK */

/*  IERR    (output) INTEGER */
/*          = 0:  successful exit */
/*                (A regular or singular with consistent equations) */
/*          = 1:  error exit (A singular, AX=B not compatible). */
/*                a message has already been computed. */

/* Author: M. Otter, DLR Oberpfaffenhofen, Feb. 1998. */
/* Modified to reuse factorization: H. Olsson, Dynasim, Feb 1999. */
/* Modified to include variable names: H. Olsson, Dynasim, Sep 2000 */

/* ====================================================================== 
*/





    /* Shortened parameter adjustments */
    --iwork;
 
    /* Function Body */

    *ierr = 0;
/* Inquire machine epsilon */
    if (1) {
	/* eps = dlamch_("E", 1) * 1e3; */
        eps = DBL_EPSILON * 1e7;
		if (fEvent && *fEvent) eps = eps*1e4;
#ifdef DYMOSIM
		{
			extern int Check5(char*key);
			if (!Check5("")) {
				DymosimMessage("");
				*ierr = 1;
				return 0;
			}
		}
#endif
    }
    dymlqr3_(fact, a, lda, n, b, &eps, &iwork[1], dwork, &
	    iwork[*n + 2], factor,&info,sysnr,varnames,fEvent);
    if ((info == 1 && *event == 1 && *printpriority >= 80)||info==2) {
		if (info==2)
		dymosimmessagedoubleSev_(1, "LINEAR SYSTEM OF EQUATIONS IS SINGULAR AT TIME = ", time, 49);
		else
	dymosimmessagedoubleSev_(1, "... System of equations is singular at Time = "
		, time, 46);
	DymosimMessageIntSev_(1, "...    Linear system of equations number = ", 
		sysnr, 43);
		if (varnames==0) {
			dymosimmessageSev_(1, "...    Variables X(i) which cannot be uniquely computed:", 56);
			i__1 = *n;
			for (i__ = 1; i__ <= i__1; ++i__) {
				if (iwork[(*n)+1+i__] == 1) {
					DymosimMessageInt_("...     ", &i__, 8);
				}
/* L10: */
			}
		} else {
			dymosimmessageSev_(1, "... Variables which cannot be uniquely computed:", 48);
			i__1 = *n;
			for (i__ = 1; i__ <= i__1; ++i__) {
				if (iwork[(*n)+1+i__] == 1) {
					char str[200];
					int i;
					for(i=0;varnames[i__-1][i]!='\0';i++) {
						str[i]=varnames[i__-1][i];
						if (i>=140) {
							strcpy(str+i,"...");i+=3;
							break;
						}
					}
					strcpy(str+i," = ");i+=3;
					dymosimmessagedoubleSev_(1, str,b+(i__-1),i);
				}
			}
		}
		if (info==1) 
	dymosimmessageSev_(1, "... Simulation continues with minimum norm solution."
		, 52);
		else {
			dymosimmessageSev_(1, "... NOT ACCEPTING SINCE TOO LARGE RESDIUAL", 42);
			*ierr = 1;
		}
	}

    return 0;
} /* dymli4_ */

LIBDS_API int dymli3_(sysnr, fact, a, lda, n, b, time, event, 
	printpriority, dwork, iwork, factor,varnames,ierr)
integer *sysnr, *fact;
doublereal *a;
const integer*lda, *n;
doublereal *b, *time;
integer *event, *printpriority;
doublereal *dwork;
integer *iwork, *factor;
const char*const*varnames;
integer *ierr;
{
	return dymli4_(sysnr,fact,a,lda,n,b,time,event,printpriority,dwork,iwork,factor,varnames,ierr,0);
}
LIBDS_API int dymli2_(sysnr, fact, a, lda, n, b, time, event, 
	printpriority, dwork, iwork, factor,ierr)
integer *sysnr, *fact;
doublereal *a;
const integer*lda, *n;
doublereal *b, *time;
integer *event, *printpriority;
doublereal *dwork;
integer *iwork, *factor, *ierr;
{
	return dymli3_(sysnr,fact,a,lda,n,b,time,event,printpriority,dwork,iwork,factor,0,ierr);
}

/* Subroutine */ LIBDS_API int dymlqr_(fact, a, lda, n, b, tol, itype, work, iwork, 
	info)
const integer *fact;
doublereal *a;
const integer*lda, *n;
doublereal *b;
const doublereal *tol;
integer *itype;
doublereal *work;
integer *iwork, *info;
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;
    doublereal d__1;

    /* Local variables */
    integer i__, j, k;
    doublereal c1, c2;
    integer n1, n2;
    doublereal s1, s2, t1, t2;
    integer iw, it1, it2, job, rank, ierr;
    doublereal smin, smax;
    integer info2;
    doublereal rcond, anorm;
    integer ismin, ismax;
    doublereal sminpr, smaxpr;


/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/* C      INTEGER            ITYPE( * ), IWORK( * ) */
/* C      DOUBLE PRECISION   A( LDA, * ), B( * ), WORK( * ) */
/*     .. */

/*  Purpose */
/*  ======= */

/*  DYMLQR computes the solution to a real system of linear equations */
/*     A * X = B, */
/*  where A is an N-by-N matrix and X and B are N-by-1 matrices. */

/*  If FACT = 0, the LU decomposition with partial pivoting and row */
/*  interchanges is used to factor A as */
/*     A = P * L * U, */
/*  where P is a permutation matrix, L is unit lower triangular, and U is 
*/
/*  upper triangular.  If A is nonsingular, then the factored form of A */
/*  is used to solve the system of equations A * X = B. */
/*  If FACT = 1 or if A is singular, the QR decomposition with */
/*  column pivoting is used to factor A as */
/*     A*P = Q * R */
/*  where P is a permutation matrix, Q is an orthogonal matrix, and R is 
*/
/*  upper triangular. The rank of A is determined and the compatibility */
/*  of the linear system is checked. If the system is compatible a */
/*  minimum norm solution is computed using the QR decomposition of A. */


/*  Arguments */
/*  ========= */

/*  FACT    (input) INTEGER */
/*          Option to specify the solution method. */
/*          FACT = 0, use the LU-decomposition to solve AX = B or */
/*                    the QR-decomposition if A is singular. */
/*          FACT = 1, use the QR-decomposition to solve AX = B. */

/*  A       (input/output) DOUBLE PRECISION array, dimension (N,N) */
/*          On entry, the N-by-N coefficient matrix A. */
/*          On exit, information on LU or QR decomposition of A. */

/*  LDA     (input) INTEGER */
/*          The declared leading dimension of A. LDA >= 1. */

/*  N       (input) INTEGER */
/*          The number of linear equations, i.e., the order of the */
/*          matrix A.  1 <= N <= LDA. */

/*  B       (input/output) DOUBLE PRECISION array, dimension (N) */
/*          On entry, the N-by-1 matrix of right hand side matrix B. */
/*          On exit, if INFO = 0, the N-by-1 solution matrix X. */

/*  TOL     (input) DOUBLE PRECISION */
/*           Tolerance for system compatibility test. */

/*  ITYPE   (output) INTEGER array, dimension (N) */
/*          Specifies the redundant/nonredundant components of X: */
/*          ITYPE(k) = 0, if the k-th component is nonredundant */
/*          ITYPE(k) = 1, if the k-th component is redundant */

/*  WORK    (output) REAL work array, dimension (LWORK), where */
/*           LWORK = N*N+4*N if FACT = 0 */
/*           LWORK = 4*N     if FACT = 1 */

/*  IWORK   (output) INTEGER work array, dimension (2*N) */

/*  INFO    (output) INTEGER */
/*          = 0:  successful exit (A regular) */
/*          = 1:  successful exit (A singular, AX=B compatible) */
/*          = 2:  error exit (A singular, AX=B not compatible). */


/* Author: A. Varga, M. Otter, DLR Oberpfaffenhofen, Feb. 1998. */

/* ====================================================================== 
*/

/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Executable Statements .. */

    /* Parameter adjustments */
    --iwork;
    --work;
    --itype;
    --b;
    a_dim1 = *lda;
    a_offset = a_dim1 + 1;
    a -= a_offset;

    /* Function Body */
    *info = 0;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	itype[i__] = 0;
/* L10: */
    }
    n1 = *n + 1;
    n2 = (*n << 1) + 1;
    iw = *n * *n + 1;

/*     Compute the LU factorization of A. */

    job = *fact;
    if (job == 0) {

/*         Save A. */

/* CC       CALL DCOPY( N*N, A, 1, WORK, 1 ) */
	k = 0;
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    i__2 = *n;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		++k;
		work[k] = a[i__ + j * a_dim1];
/* L11: */
	    }
/* L12: */
	}

/*         Compute the LU factorization of A. */

	dgetrf_(n, n, &work[1], n, &iwork[1], &ierr);

/*         Estimate reciprocal condition number and check for singular
ity. */

	if (ierr == 0) {
	    anorm = dlange_("1", n, n, &work[1], n, &work[iw], 1);
	    dgecon_("1", n, &work[1], n, &anorm, &rcond, &work[iw], &iwork[n1]
		    , &info2, 1);

/*            real workspace 4*N; integer workspace N */

	    if (rcond <= *tol) {
		ierr = 1;
	    }
	}
	if (ierr == 0) {

/*            Solve the system A*X = B, overwriting B with X. */

	    dgetrs_("No transpose", n, &c__1, &work[1], n, &iwork[1], &b[1], 
		    n, &info2, 12);
	} else {

/*            Singular A, use QR decomposition. */

	    job = 1;
	}
    }

    if (job == 1) {
	*info = 1;

/*        Compute the QR factorization with column pivoting of A: */
/*          A * P = Q * R */

	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    iwork[i__] = 0;
/* L15: */
	}
	dgeqpf_(n, n, &a[a_offset], lda, &iwork[1], &work[1], &work[n1], &
		ierr);

/*        workspace 3*N. Details of Householder rotations stored */
/*        in WORK(1:N). */


/*        B := Q' * B */

	dorm2r_("Left", "Transpose", n, &c__1, n, &a[a_offset], lda, &work[1],
		 &b[1], n, &c1, &info2, 4, 9);

/*        Determine RANK using incremental condition estimation */

	ismin = n1;
	ismax = n2;
	work[ismin] = 1.;
	work[ismax] = 1.;
	smax = (d__1 = a[a_dim1 + 1], Dymola_abs(d__1));
	smin = smax;
	if ((d__1 = a[a_dim1 + 1], Dymola_abs(d__1)) == 0.) {
	    rank = 0;
	} else {
	    rank = 1;
	}

L20:
	if (rank < *n) {
	    i__ = rank + 1;
	    dlaic1_(&c__2, &rank, &work[ismin], &smin, &a[i__ * a_dim1 + 1], &
		    a[i__ + i__ * a_dim1], &sminpr, &s1, &c1);
	    dlaic1_(&c__1, &rank, &work[ismax], &smax, &a[i__ * a_dim1 + 1], &
		    a[i__ + i__ * a_dim1], &smaxpr, &s2, &c2);

	    if (smaxpr * *tol <= sminpr) {
		i__1 = rank;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    work[ismin + i__ - 1] = s1 * work[ismin + i__ - 1];
		    work[ismax + i__ - 1] = s2 * work[ismax + i__ - 1];
/* L30: */
		}
		work[ismin + rank] = c1;
		work[ismax + rank] = c2;
		smin = sminpr;
		smax = smaxpr;
		++rank;
		goto L20;
	    }
	}

/*        Check compatibility of the system */

	if (rank >= *n) {
	    *info = 0;
	}
	i__1 = *n;
	for (i__ = rank + 1; i__ <= i__1; ++i__) {
	    if ((d__1 = b[i__], Dymola_abs(d__1)) > *tol) {
		*info = 2;
		return 0;
	    }
	    itype[i__] = 1;
	    b[i__] = 0.;
/* L40: */
	}
	if (rank == 0) {
	    return 0;
	}

/*        Logically partition R = [ R11 R12 ] */
/*                                [  0  R22 ] */
/*        where R11 = R(1:RANK,1:RANK) */

/*        [R11,R12] = [ T11, 0 ] * Y */
	if (rank < *n) {
	    dtzrqf_(&rank, n, &a[a_offset], lda, &work[n1], &ierr);
	}

/*        Details of Householder rotations stored in WORK(N1:2*N) */

/*        B(1:RANK) := inv(T11) * B(1:RANK) */

	dtrsm_("Left", "Upper", "No transpose", "Non-unit", &rank, &c__1, &
		c_b263, &a[a_offset], lda, &b[1], n, 4, 5, 12, 8);

/*        B := Y' * B */

	if (rank < *n) {
	    i__1 = rank;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		i__2 = *n - rank + 1;
		dlatzm_("Left", &i__2, &c__1, &a[i__ + (rank + 1) * a_dim1], 
			n, &work[*n + i__], &b[i__], &b[rank + 1], n, &c1, 4)
			;
/* L50: */
	    }
	}

/*        B := P * B */

	--n2;
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    work[(*n << 1) + i__] = 1.;
/* L60: */
	}
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    if (work[n2 + i__] == 1.) {
		if (iwork[i__] != i__) {
		    k = i__;
		    t1 = b[k];
		    t2 = b[iwork[k]];
		    it1 = itype[k];
		    it2 = itype[iwork[k]];
L70:
		    b[iwork[k]] = t1;
		    itype[iwork[k]] = it1;
		    work[n2 + k] = 0.;
		    t1 = t2;
		    it1 = it2;
		    k = iwork[k];
		    t2 = b[iwork[k]];
		    it2 = itype[iwork[k]];
		    if (iwork[k] != i__) {
			goto L70;
		    }
		    b[i__] = t1;
		    itype[i__] = it1;
		    work[n2 + k] = 0.;
		}
	    }
/* L80: */
	}
    }
    return 0;

/*     End of DYMLQR */

} /* dymlqr_ */

/* Subroutine */ LIBDS_API int dymlqr3_(fact, a, lda, n, b, tol, ipivot, work, iwork, factor,
	info,sysnr,varnames,fEvent)
integer *fact;
doublereal *a;
const integer*lda, *n;
doublereal *b, *tol;
integer *ipivot;
doublereal *work;
integer *iwork, *factor, *info;
integer *sysnr;
const char*const*varnames;
integer *fEvent;
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;
    doublereal d__1;

    /* Local variables */
    integer i__, j, k;
    doublereal c1, c2;
    integer n1, n2;
	integer workTau,workTau2,workWork,workB;
    doublereal s1, s2, t1, t2;
    integer iw, it1, it2, factjob, eqjob, rank, ierr;
    doublereal smin, smax;
    integer info2;
    doublereal rcond, anorm;
    integer ismin, ismax;
    integer iscale,iequilib;
    doublereal eqrow,eqcol,eqmax;
	integer rank_def;
	integer rank_max;
    char equed;
    doublereal sminpr, smaxpr;


/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/* C      INTEGER            ITYPE( * ), IWORK( * ) */
/* C      DOUBLE PRECISION   A( LDA, * ), B( * ), WORK( * ) */
/*     .. */

/*  Purpose */
/*  ======= */

/*  DYMLQR computes the solution to a real system of linear equations */
/*     A * X = B, */
/*  where A is an N-by-N matrix and X and B are N-by-1 matrices. */

/*  If FACT = 0, the LU decomposition with partial pivoting and row */
/*  interchanges is used to factor A as */
/*     A = P * L * U, */
/*  where P is a permutation matrix, L is unit lower triangular, and U is 
*/
/*  upper triangular.  If A is nonsingular, then the factored form of A */
/*  is used to solve the system of equations A * X = B. */
/*  If FACT = 1 or if A is singular, the QR decomposition with */
/*  column pivoting is used to factor A as */
/*     A*P = Q * R */
/*  where P is a permutation matrix, Q is an orthogonal matrix, and R is 
*/
/*  upper triangular. The rank of A is determined and the compatibility */
/*  of the linear system is checked. If the system is compatible a */
/*  minimum norm solution is computed using the QR decomposition of A. */


/*  Arguments */
/*  ========= */

/*  FACT    (input) INTEGER */
/*          Option to specify the solution method. (See also Factor). */
/*          FACT&1 = 0, use the LU-decomposition to solve AX = B or */
/*                    the QR-decomposition if A is singular. */
/*          FACT&1 = 1, use the QR-decomposition to solve AX = B. */
/*          FACT&2 = 2, automatic row and column scaling */

/*  A       (input) DOUBLE PRECISION array, dimension (N,N) */
/*          On entry, the N-by-N coefficient matrix A. */

/*  LDA     (input) INTEGER */
/*          The declared leading dimension of A. LDA >= 1. */

/*  N       (input) INTEGER */
/*          The number of linear equations, i.e., the order of the */
/*          matrix A.  1 <= N <= LDA. */

/*  B       (input/output) DOUBLE PRECISION array, dimension (N) */
/*          On entry, the N-by-1 matrix of right hand side matrix B. */
/*          On exit, if INFO = 0, the N-by-1 solution matrix X. */

/*  TOL     (input) DOUBLE PRECISION */
/*           Tolerance for system compatibility test. */

/*  IPIVOT  (input/output) INTEGER array, dimension (N+1) */
/*          Pivot information for LU and QR. The last elements */
/*          is the rank for the QR-factorization */

/*  WORK    (input/output) REAL work array, dimension (LWORK), where */
/*           LWORK = N*N+6*N */

/*  IWORK   (output) INTEGER work array, dimension (N) */
/*          Specifies the redundant/nonredundant components of X: */
/*          ITYPE(k) = 0, if the k-th component is nonredundant */
/*          ITYPE(k) = 1, if the k-th component is redundant */

/*  FACTOR  (input) INTEGER */
/*    Factor&3     = 0: Factorize matrix
                  = 1: Use existing LU-decomposition
                  = 2: Use existing QR-decomposition
     Factor&256   =256: Use row-equlibrization
           &512   =512: Use column-equilibrization
             Updated after each call.
             The LU-factorization is stored in the first N*N elements 
             of WORK and the first N elements of IPIVOT
             The QR-factorization is stored in the first N*N+2*N elements
             of WORK and and IPIVOT.
	     The equilib information is stored in the 2*N last elements
	     of LWORK */


/*  INFO    (output) INTEGER */
/*          = 0:  successful exit (A regular) */
/*          = 1:  successful exit (A singular, AX=B compatible) */
/*          = 2:  error exit (A singular, AX=B not compatible). */


/* Author: A. Varga, M. Otter, DLR Oberpfaffenhofen, Feb. 1998. */
/* Modified to reuse factorization: H. Olsson, Dynasim, Feb. 1999. */
/* Modified to use row/column equilibrization: H. Olsson, Dynasim, July 1999 */
/* Modified to not modify A: H. Olsson, Dynasim, Jan. 2000*/

/* ====================================================================== 
*/

/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Executable Statements .. */

    /* Parameter adjustments */
    --iwork;
    --work;
    --ipivot;
    --b;
    a_dim1 = *lda;
    a_offset = a_dim1 + 1;
    a -= a_offset;


    /* Function Body */
	rank_max = 0;
 L15:
	*info = 0;
    n1 = *n + 1;
    n2 = (*n << 1) + 1;
	workTau = *n * *n +1;
	workTau2 = workTau + *n;
	workWork = workTau2 + *n;
	workB= workWork+*n;
    iw = *n * *n + 1;

    iscale = iw + 4* *n;
/*     Compute the LU factorization of A. */

    factjob = *fact &1;
    eqjob = *fact &2;

	for(i__=1;i__<=*n;++i__) {
	    volatile double v=b[i__];
#if 0
#elif 1
			if (v!=v || v>=DBL_MAX || v<=-DBL_MAX) {				
				char str[100];
				sprintfC(str,"Rhs component %d of system %d is %e",i__,*sysnr,v);
				DymosimMessage(str);
#endif
			}
	}
				
	if ((*factor & 3) == 0) {

/*         Save A and then use copy */

/* CC       CALL DCOPY( N*N, A, 1, WORK, 1 ) */
        k = 0;
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
          i__2 = *n;
          for (i__ = 1; i__ <= i__2; ++i__) {
            ++k;
            work[k] = a[i__ + j * a_dim1];
/* L11: */
	    
          }
/* L12: */
        }
      
	/* Determine row and column scaling for equilibrization */
      iequilib=0;
      ierr=0;
      if (eqjob) {dgeequ_(n, n, &a[a_offset], lda, &work[iscale], &work[iscale+*n], 
	      &eqrow, &eqcol, &eqmax, &ierr);}
      if (eqjob && ierr==0) {
	      /* equilib the matrix */
		dlaqge_(n, n, &work[1], n, &work[iscale], &work[iscale+*n], 
			&eqrow, &eqcol, &eqmax, &equed);

		iequilib=(equed=='B')? 512|256 : (equed=='C')? 512 : (equed=='R') ? 256 :0;
      }
      if (factjob == 0 && ierr==0) { /* No need to try LU-factorization if equilibrization failed */
		  ierr=0;


/*         Compute the LU factorization of A. */
	
	   dgetrf_(n, n, &work[1], n, &ipivot[1], &ierr);

/*         Estimate reciprocal condition number and check for singular
ity. */


        if (ierr == 0) {
          anorm = dlange_("1", n, n, &work[1], n, &work[iw], 1);
          dgecon_("1", n, &work[1], n, &anorm, &rcond, &work[iw], &iwork[1]
            , &info2, 1);

/*            real workspace 4*N; integer workspace N */
/* tau-vectors are not yet used and can thus be used as work-space*/

          if (rcond <= *tol) {
            ierr = 1;
	  }
        }
        if (ierr == 0) {
          *factor=1 | iequilib;
        } else {

		/*DymosimMessageInt_("Size of equation system: ",n,24);*/
/*            Singular A, use QR decomposition. */
          factjob =1;
        }
      }
      if (factjob == 1) {
		  /* Restore work to equilibrated A */
		  if (ierr != 0) {
			  /* LU-factorization failed: have to restore*/
			  k = 0;
			  i__1 = *n;
			  for (j = 1; j <= i__1; ++j) {
				  i__2 = *n;
				  for (i__ = 1; i__ <= i__2; ++i__) {
					  ++k;
					  work[k]=a[i__ + j * a_dim1];
					  /* L11: */
					  
				  }
				  /* L12: */
			  }
			  if (iequilib) /* Equilibrate once more: should give same results */
				  dlaqge_(n, n, &work[1], n, &work[iscale], &work[iscale+*n], 
				  &eqrow, &eqcol, &eqmax, &equed);
		  }
		  ierr = 0;

/*        Compute the QR factorization with column pivoting of A: */
/*          A * P = Q * R */
        
        i__1 = *n;
        for (i__ = 1; i__ <= i__1; ++i__) {
          ipivot[i__] = 0;
/* L15: */
	}
	dgeqpf_(n, n, &work[1], n, &ipivot[1], &work[workTau], &work[workTau2], &
		ierr);

/*        workspace 3*N (second tau is not yet used).
	         Details of Householder rotations stored */
/*        in WORK(1:N). */

 

/*        Determine RANK using incremental condition estimation */
/* Uses work(2*n+1:3*n) and work(3*n+1:4*n) as work area */
/* And R-part of QR-factorization */

	ismin = workWork;
	ismax = workWork+*n;
	work[ismin] = 1.;
	work[ismax] = 1.;
	smax = (d__1 = work[1], Dymola_abs(d__1));
	smin = smax;
	if (rank_max>0) {rank=rank_max;goto L25;}
	if ((d__1 = work[1], Dymola_abs(d__1)) == 0.) {
	    rank = 0;
	} else {
	    rank = 1;
	}
	rank_max = rank;

L20:
	if (rank_max < *n) {
	    i__ = rank_max + 1;
	    dlaic1_(&c__2, &rank_max, &work[ismin], &smin, &work[(i__-1) * *n + 1], 
		    &work[(i__-1) * *n + i__], &sminpr, &s1, &c1);
	    dlaic1_(&c__1, &rank_max, &work[ismax], &smax, &work[(i__-1) * *n + 1], 
		    &work[(i__-1) * *n + i__], &smaxpr, &s2, &c2);

		/* Hard lower limit for the residual */
		/* This could be decreased even further */
	    if (smaxpr * DBL_EPSILON <= sminpr) {
			if (smaxpr * *tol <=sminpr && rank_max==rank) {
				++rank;
			}
		i__1 = rank_max;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    work[ismin + i__ - 1] = s1 * work[ismin + i__ - 1];
		    work[ismax + i__ - 1] = s2 * work[ismax + i__ - 1];
/* L30: */
		}
		work[ismin + rank_max] = c1;
		work[ismax + rank_max] = c2;
		smin = sminpr;
		smax = smaxpr;
		++rank_max;
		goto L20;
	    }
	}


/*        Logically partition R = [ R11 R12 ] */
/*                                [  0  R22 ] */
/*        where R11 = R(1:RANK,1:RANK) */

/*        [R11,R12] = [ T11, 0 ] * Y */
	/* After estimating rank*/
L25:
	if (rank < *n) {
	    dtzrqf_(&rank, n, &work[1], n, &work[workTau2], &ierr);
	}
        *factor=2 | iequilib;

        ipivot[(*n)+1]=rank;
/* L80: */
      }
    }
    iequilib = *factor;
	/* Remember b */
    rank=ipivot[(*n)+1];
	rank_def=((*factor&3)==2 && rank<*n);
	if (rank_def) {
		for( i__1 = 1; i__1 <= *n; ++i__1) work[workB+i__1-1]=b[i__1];
	}
    if (iequilib&256) {
	    /* Apply row scaling */
	    for( i__1 = 1; i__1 <= *n; ++i__1) b[i__1] *= work[i__1-1+iscale];
    }
    if ((*factor & 3)==1) {
/*            Solve the system A*X = B, overwriting B with X. */

	    dgetrs_("No transpose", n, &c__1, &work[1], n, &ipivot[1], &b[1], 
		    n, &info2, 12);

    } else if ((*factor & 3)==2) {        
      *info = 1;
      i__1 = *n;


        for (i__ = 1; i__ <= i__1; ++i__) {
          iwork[i__] = 0;
      }

     /*        B := Q' * B */

	dorm2r_("Left", "Transpose", n, &c__1, n, &work[1], n, &work[workTau],
		 &b[1], n, &c1, &info2, 4, 9);
        rank=ipivot[(*n)+1];

 
/*        Check compatibility of the system */

	if (rank >= *n) {
	    *info = 0;
	}
        i__1 = *n;
	for (i__ = rank + 1; i__ <= i__1; ++i__) {
	    if ((d__1 = b[i__], Dymola_abs(d__1)) > *tol) {
		/**info = 2; Removed by Hans Olsson, Dynasim 2000-03-24.*/
	    /* Instead we test the residual for sufficiently small numbers*/
		/* return 0; */
	    }
	    iwork[i__] = 1;
	    b[i__] = 0.;
/* L40: */
	}
	if (rank == 0) {
	    return 0; /* No need to scale if everything is zero */
	}
 

/*        Details of Householder rotations stored in WORK(work) */

/*        B(1:RANK) := inv(T11) * B(1:RANK) */

	dtrsm_("Left", "Upper", "No transpose", "Non-unit", &rank, &c__1, &
		c_b263, &work[1], n, &b[1], n, 4, 5, 12, 8);

/*        B := Y' * B */

	if (rank < *n) {
	    i__1 = rank;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		i__2 = *n - rank + 1;
		dlatzm_("Left", &i__2, &c__1, &work[i__ + (rank + 1-1) * *n], 
			n, &work[workTau2+ i__-1], &b[i__], &b[rank + 1], n, &c1, 4)
			;
/* L50: */
	    }
	}

/*        B := P * B */

	--n2;
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    work[workWork + i__-1] = 1.;
/* L60: */
	}
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    if (work[workWork + i__-1] == 1.) {
		if (ipivot[i__] != i__) {
		    k = i__;
		    t1 = b[k];
		    t2 = b[ipivot[k]];
		    it1 = iwork[k];
		    it2 = iwork[ipivot[k]];
L70:
		    b[ipivot[k]] = t1;
		    iwork[ipivot[k]] = it1;
		    work[workWork + k-1] = 0.;
		    t1 = t2;
		    it1 = it2;
		    k = ipivot[k];
		    t2 = b[ipivot[k]];
		    it2 = iwork[ipivot[k]];
		    if (ipivot[k] != i__) {
			goto L70;
		    }
		    b[i__] = t1;
		    iwork[i__] = it1;
		    work[workWork + k-1] = 0.;
		}
	    }
        }
    }
    if (iequilib&512) {
	    /* Apply column scaling */
	    for( i__1 = 1; i__1 <= *n; ++i__1) b[i__1] *= work[i__1-1+iscale+*n];
    }
	if (rank_def) {
		c1=-1;
		c2=1;
		i__=1;
		dcopy_(n, &work[workB], &i__, &work[workWork], &i__);
		dgemv_("N",n,n,&c1,&a[a_offset],lda,&b[1],&i__,&c2,&work[workWork],&i__,1);
		for(i__=1;i__<=*n;++i__) if (c1=work[workWork+i__-1],Dymola_abs(c1)>*tol) {
			if (rank<rank_max) {
				/* Try to increase rank, restore, and re-test */
				/* This should be done in a more optimal way, but assume it does not happen so often */
				for( i__1 = 1; i__1 <= *n; ++i__1) b[i__1]=work[workB+i__1-1];
				*factor =0; /* Will refactor */
				goto L15;
			} else {
			char str[100];
			sprintfC(str,"Residual component %d of system %d is %e",i__,*sysnr,c1);
			DymosimMessage(str);
			if (!(fEvent && fEvent[2])) *info=2;
			}
		}
	} else if (1) {
		for(i__=1;i__<=*n;++i__) {
		    volatile double v=b[i__];
#if 0
#elif 1
			if (v!=v || v>=DBL_MAX || v<=-DBL_MAX) {				
				char str[200];
				if (varnames) {
					const char*extra="";
					if (strlen(varnames[i__-1])>140) extra="...";
					sprintfC(str,"In system %d %.140s%s = %e",*sysnr,varnames[i__-1],extra,v);
				} else 
					sprintfC(str,"Solution component %d of system %d is %e",i__,*sysnr,v);
				DymosimMessage(str);
				if (!(fEvent && fEvent[2])) *info=2;
			}
#else
			if (v>=1e30 || v<-1e30) {
				char str[100];
				sprintfC(str,"Solution component %d of system %d is %e",i__,*sysnr,v);
				DymosimMessage(str);
			}
#endif
		}
	}
    return 0;

/*     End of DYMLQR2 */

} /* dymlqr2_ */
/* Subroutine */ LIBDS_API int dymmdp_(void)
{


/* Set machine dependent constants in common /DYMCB/: */
/*     EPSMCH    machine precision */
/*               (may be determined by function EPSLON from EISPACK) */
/*     GIANT     largest magnitude */

/* Libraries used: */
/*  EISPACK:  epslon */




    
    return 0;

} /* dymmdp_ */


/* inJacobian_ = 0 during normal simulation, 
   inJacobian_ = 1 (or any non-zero value) during computation of Jacobian for ODE-solver*/


/* dymnl_WithNominal */
/* The function call before computing the Jacobian must be at f(x) and Jacobians at f(x+dx) */

/* -----------------------------------------------------------------------
 */
/*            INFO =  0  Improper input parameters (wrong initialization) 
*/
/*            INFO =  1  Algorithm estimates that the relative error */
/*                       between X and the solution is at most TOL */
/*                       (regular end). */
/*            INFO =  2  Number of calls to FCN has reached or exceeded */
/*                       100*(N+1) for IOPT=1 or 200*(N+1) for IOPT=2. */
/*            INFO =  3  TOL is too small.  No further improvement in the 
*/
/*                       approximate solution X is possible. */
/*            INFO =  4  Iteration is not making good progress. */
/*         ENDIF */
/* DUM   : OUT, DOUBLE( (N**2+13*N)/2 ) */
/*         work array of length LDUM = 20 + (N**2+13*N)/2 */

/* LDUM  : IN, INTEGER */
/*         Provided length of the work array DUM. */
/*         LDUM must not be less than 20 + (N**2+13*N)/2. */

/* IDUM  : IN, INTEGER */
/*         work array of length 20. */

/* Nominal: IN, DOUBLE */
/*         Used for determining suitable scaling of x.*/

/* Method: */
/* DYMNL is a modification of the Powell Hybrid method.  Two of */
/* its main characteristics involve the choice of the correction as */
/* a convex combination of the Newton and scaled gradient directions, */
/* and the updating of the Jacobian by the rank-1 method of Broyden. */
/* The choice of the correction guarantees (under reasonable conditions) 
*/
/* global convergence for starting points far from the solution and a */
/* fast rate of convergence. The Jacobian is calculated at the starting */
/* point by either the user-supplied code or a forward-difference */
/* approximation, but it is not recalculated until the rank-1 method */
/* fails to produce satisfactory progress. */

/* The time required by DYMNL to solve a given problem depends on N, */
/* the behavior of the functions, the accuracy requested, and the */
/* starting point. The number of arithmetic operations needed by DYMNL */
/* is about 11.5*(N**2) to process each evaluation of the functions and */
/* 1.3*(N**3) to process each evaluation of the Jacobian (if IOPT = 1). */

/* This code is the combination of the MINPACK codes (Argonne) HYBRD1 and 
*/
/* HYBRJ1. */

/* References: */
/*   Gramlich, G., Ein nichtlinearer Gleichungssytemloeser mit "reverse- 
*/
/*      communication". DLR TR ER59-91, 1991 */
/*   Powell, M., A Hybrid Method for Nonlinear Equations. */
/*      Numerical Methods for Nonlinear Algebraic Equations, */
/*      P. Rabinowitz, Editor, Gordon and Breach, 1970 */

/* Life cycle: */
/* 1980 MAR  Burton S. Garbow, Kenneth E. Hillstrom, Jorge J. More, */
/*           Argonne National Laboratory. Minpack Project. */
/* 1991 OKT  G. Gramlich, DLR FF-DF: reverse communication interface */
/* 1993 SEP  D. Joos, DLR FF-DR    : iopt=3 introduced */
/* 1997 AUG  M. Otter, DLR FF-DR   : termination criteria improved, */
/*                                   multiple usage at the same time */
/*                                   (internal state is stored in */
/*                                   work arrays). */
/* 2002 FEB H. Olsson, Dynasim. Revised to have nominal value */

DYMOLA_STATIC void dymnl_WithNominal2(integer*infrev, integer iopt, integer n,doublereal*x,
					   doublereal*fvec,doublereal*fjac,doublereal tol,doublereal toldesired,
					   integer*info,doublereal*dum,integer ldum,integer*idum,
					   doublereal const*nominalx,const char*const*varnames,
					   integer*ipivots,int printEvent,int inJacobian) {
    /* Initialized data */

    static const doublereal factor = 100.;
    static const doublereal one = 1.;
    static const doublereal zero = 0.;

    /* Local variables */
    integer j, ml, lr, mu, ibeg, nfev, njev;
    doublereal xtol, xtoldesired;
    integer jopt1, jopt2;
    integer ldfjac, maxfev;

/* -----------------------------------------------------------------------
 */
/*            INFO =  0  Improper input parameters (wrong initialization) 
*/
/*            INFO =  1  Algorithm estimates that the relative error */
/*                       between X and the solution is at most TOL */
/*                       (regular end). */
/*            INFO =  2  Number of calls to FCN has reached or exceeded */
/*                       100*(N+1) for IOPT=1 or 200*(N+1) for IOPT=2. */
/*            INFO =  3  TOL is too small.  No further improvement in the 
*/
/*                       approximate solution X is possible. */
/*            INFO =  4  Iteration is not making good progress. */
/*         ENDIF */
/* DUM   : OUT, DOUBLE( (N**2+15*N)/2 ) */
/*         work array of length LDUM = 20 + (N**2+15*N)/2 */

/* LDUM  : IN, INTEGER */
/*         Provided length of the work array DUM. */
/*         LDUM must not be less than 20 + (N**2+15*N)/2. */

/* IDUM  : IN, INTEGER */
/*         work array of length 20. */

/* Method: */
/* DYMNL is a modification of the Powell Hybrid method.  Two of */
/* its main characteristics involve the choice of the correction as */
/* a convex combination of the Newton and scaled gradient directions, */
/* and the updating of the Jacobian by the rank-1 method of Broyden. */
/* The choice of the correction guarantees (under reasonable conditions) 
*/
/* global convergence for starting points far from the solution and a */
/* fast rate of convergence. The Jacobian is calculated at the starting */
/* point by either the user-supplied code or a forward-difference */
/* approximation, but it is not recalculated until the rank-1 method */
/* fails to produce satisfactory progress. */

/* The time required by DYMNL to solve a given problem depends on N, */
/* the behavior of the functions, the accuracy requested, and the */
/* starting point. The number of arithmetic operations needed by DYMNL */
/* is about 11.5*(N**2) to process each evaluation of the functions and */
/* 1.3*(N**3) to process each evaluation of the Jacobian (if IOPT = 1). */

/* This code is the combination of the MINPACK codes (Argonne) HYBRD1 and 
*/
/* HYBRJ1. */

/* References: */
/*   Gramlich, G., Ein nichtlinearer Gleichungssytemloeser mit "reverse- 
*/
/*      communication". DLR TR ER59-91, 1991 */
/*   Powell, M., A Hybrid Method for Nonlinear Equations. */
/*      Numerical Methods for Nonlinear Algebraic Equations, */
/*      P. Rabinowitz, Editor, Gordon and Breach, 1970 */

/* Life cycle: */
/* 1980 MAR  Burton S. Garbow, Kenneth E. Hillstrom, Jorge J. More, */
/*           Argonne National Laboratory. Minpack Project. */
/* 1991 OKT  G. Gramlich, DLR FF-DF: reverse communication interface */
/* 1993 SEP  D. Joos, DLR FF-DR    : iopt=3 introduced */
/* 1997 AUG  M. Otter, DLR FF-DR   : termination criteria improved, */
/*                                   multiple usage at the same time */
/*                                   (internal state is stored in */
/*                                   work arrays). */

/* Libraries required: */
/* BLAS1:    DNRM2 */

/* Common Blocks: */
/* The subroutines DYMNL* need common /DYMCB/ for the values of */
/*  the machine precision EPSMCH and largest magnitude GIANT. */

/* Example: */

/* The problem is to determine the values of X(1) and X(2), */
/* which solve the nonlinear system of equations */

/*             X(1)*X(2) - X(2)**3 - 1 = 0 */
/*             X(1)**2*X(2) + X(2) - 5 = 0 */

/* ********** */

/*       PROGRAM TEST */
/* C */
/* C     DRIVER FOR DYMNL EXAMPLE. */
/* C */
/*       INTEGER J,N,INFREV,IOPT,INFO,LDUM,IDUM(20) */
/*       PARAMETER (N=2,LDUM=35) */
/*       DOUBLE PRECISION TOL */
/*       DOUBLE PRECISION X(N),FVEC(N),DUM(LDUM),FJAC(N,N) */
/* C */
/*       DOUBLE PRECISION DNRM2, FNORM */
/* C */
/* C     Initialize machine dependent constants */
/*       CALL DYMMDP */
/* C */
/*       IOPT = 1 */
/*       INFREV=-1 */
/*       X(1)=1.D0 */
/*       X(2)=1.D0 */
/* C */
/* C     SET TOL TO ZERO, HENCE THE SQUARE ROOT OF THE MACHINE PRECISION 
*/
/* C     IS ASSUMED INTERNALY. */
/* C */
/*       TOL = 0.D0 */
/* C */
/* 10    CONTINUE */
/*       IF(INFREV .GE. 1 ) THEN */
/*          FVEC(1)=X(1)*X(2)-X(2)**3-1.D0 */
/*          FVEC(2)=X(1)**2*X(2)+X(2)-5.D0 */
/*       END IF */
/*       IF(INFREV .EQ. 2 ) THEN */
/*          FJAC(1,1)=X(2) */
/*          FJAC(1,2)=X(1)-3.*X(2)**2 */
/*          FJAC(2,1)=2.*X(1)*X(2) */
/*          FJAC(2,2)=X(1)**2+1.D0 */
/*       END IF */
/*      CALL DYMNL(INFREV,IOPT,N,X,FVEC,FJAC,TOL,INFO,DUM,LDUM,IDUM,0.0D0)
*/
/*       IF (INFREV .GE. 1) GO TO 10 */
/*       FNORM = DNRM2(N,FVEC,1) */
/*       WRITE (6,1000) FNORM,INFO,(X(J),J=1,N) */
/*  1000 FORMAT (5X,' FINAL L2 NORM OF THE RESIDUALS',E15.7 // */
/*      *        5X,' INFO PARAMETER',16X,I10 // */
/*      *        5X,' FINAL APPROXIMATE SOLUTION' // (5X,3E15.7)) */
/*       STOP */
/*       END */

/* ************ */

/*      FINAL L2 NORM OF THE RESIDUALS  0.1434642E-10 */

/*      INFO PARAMETER                         1 */

/*      FINAL APPROXIMATE SOLUTION */

/*      0.2000000E+01  0.1000000E+01 */

/* ************ */

/* Results obtained with different compilers or machines may be slightly 
*/
/* different. */
/* -----------------------------------------------------------------------
 */




/*     .. executable statements .. */
    /* Parameter adjustments */
    --fvec;
    --x;
    --fjac;
    --dum;
    --idum;

	if (!infrev || !info) return; /* Null-pointer check */
    /* Function Body */

/* ***first executable statement  dymnl */
    ibeg = 20;
    if (*infrev < 1 && !(*infrev<=-20 && *infrev>=-40)) {
		/*        ... initial call */
		if (nominalx==0) {
			for (j = 1; j <= n; ++j) {
				dum[ibeg + j] = one;
			}
		} else {
			for(j = 1; j  <= n; ++j) {
				dum[ibeg + j] = nominalx[j-1];
			}
		}
    }


	if (iopt == 1) {
		jopt1 = 1;
		jopt2 = 0;
	} else if (iopt == 2) {
		jopt1 = 2;
		jopt2 = 0;
	} else if (iopt == 3) {
		jopt1 = 2;
		jopt2 = 1;
	} else {
		*info = 0;
		goto L30;
	}

    maxfev = (n + 1) * 100;
    if (iopt == 2) {
		maxfev *= 2;
    }
	{
		double sqrtn=sqrt(n);
		xtol = tol*sqrtn;
		if (xtol <= zero) {
			xtol = sqrt(DBL_EPSILON);
		}
		xtoldesired = toldesired*sqrtn;
		if (xtoldesired > xtol)
			xtoldesired=xtol;
	}
    ml = n - 1;
    mu = n - 1;
    lr = n * (n + 1) / 2;
    ldfjac = n;
    dymnl1_Fast(infrev, jopt1, n, &x[1], &fvec[1], &fjac[1], ldfjac, xtol, xtoldesired,
		maxfev, ml, mu, nominalx, &dum[ibeg + 1], &dum[ibeg+n*6+1],factor, info, &
	    nfev, &njev, &dum[ibeg + n * 7 + 1], lr, &dum[ibeg + n + 1], &
	    dum[ibeg + n * 2+ 1], &dum[ibeg + n * 3 + 1], &dum[ibeg + 
	    n * 4 + 1], &dum[ibeg + n * 5 + 1], jopt2, &dum[1], &idum[1],varnames,ipivots,printEvent,inJacobian);
    if (*info == 5) {
		*info = 4;
    }
L30:
    return ;
}
DYMOLA_STATIC void dymnl_WithNominal(integer*infrev, integer iopt, integer n,doublereal*x,
					   doublereal*fvec,doublereal*fjac,doublereal tol,doublereal toldesired,
					   integer*info,doublereal*dum,integer ldum,integer*idum,
					   doublereal const*nominalx) {
	dymnl_WithNominal2(infrev, iopt, n,x,
					   fvec,fjac,tol,toldesired,
					   info,dum,ldum,idum,
					   nominalx,0,0,0,0);
}

/* Subroutine */ DYMOLA_STATIC int dymnl_(infrev, iopt, n, x, fvec, fjac, tol, toldesired, info, dum, 
	ldum, idum, epsfcn)
integer *infrev, *iopt;
const integer *n;
doublereal *x, *fvec, *fjac;
const doublereal *tol, *toldesired;
integer *info;
doublereal *dum;
integer*ldum, *idum;
const doublereal *epsfcn;
{
	dymnl_WithNominal2(infrev,*iopt,*n,x,fvec,fjac,*tol,*toldesired,info,dum,*ldum,idum,
		(const doublereal*)(0),(const char*const*)(0),(integer*)(0),0,0);
	return 0;

/* -----------------------------------------------------------------------
 */
/*            INFO =  0  Improper input parameters (wrong initialization) 
*/
/*            INFO =  1  Algorithm estimates that the relative error */
/*                       between X and the solution is at most TOL */
/*                       (regular end). */
/*            INFO =  2  Number of calls to FCN has reached or exceeded */
/*                       100*(N+1) for IOPT=1 or 200*(N+1) for IOPT=2. */
/*            INFO =  3  TOL is too small.  No further improvement in the 
*/
/*                       approximate solution X is possible. */
/*            INFO =  4  Iteration is not making good progress. */
/*         ENDIF */
/* DUM   : OUT, DOUBLE( (N**2+13*N)/2 ) */
/*         work array of length LDUM = 20 + (N**2+13*N)/2 */

/* LDUM  : IN, INTEGER */
/*         Provided length of the work array DUM. */
/*         LDUM must not be less than 20 + (N**2+13*N)/2. */

/* IDUM  : IN, INTEGER */
/*         work array of length 20. */

/* EPSFCN: IN, DOUBLE */
/*         Used in determining a suitable step for the forward-difference 
*/
/*         approximation. This approximation assumes that the relative */
/*         errors in the functions are of the order of epsfcn. If epsfcn 
*/
/*         is less than the machine precision, it is assumed that the */
/*         relative errors in the functions are of the order of the */
/*         machine precision. If iopt=1, then epsfcn can be ignored */
/*         (treat it as a dummy argument). */

/* Method: */
/* DYMNL is a modification of the Powell Hybrid method.  Two of */
/* its main characteristics involve the choice of the correction as */
/* a convex combination of the Newton and scaled gradient directions, */
/* and the updating of the Jacobian by the rank-1 method of Broyden. */
/* The choice of the correction guarantees (under reasonable conditions) 
*/
/* global convergence for starting points far from the solution and a */
/* fast rate of convergence. The Jacobian is calculated at the starting */
/* point by either the user-supplied code or a forward-difference */
/* approximation, but it is not recalculated until the rank-1 method */
/* fails to produce satisfactory progress. */

/* The time required by DYMNL to solve a given problem depends on N, */
/* the behavior of the functions, the accuracy requested, and the */
/* starting point. The number of arithmetic operations needed by DYMNL */
/* is about 11.5*(N**2) to process each evaluation of the functions and */
/* 1.3*(N**3) to process each evaluation of the Jacobian (if IOPT = 1). */

/* This code is the combination of the MINPACK codes (Argonne) HYBRD1 and 
*/
/* HYBRJ1. */

/* References: */
/*   Gramlich, G., Ein nichtlinearer Gleichungssytemloeser mit "reverse- 
*/
/*      communication". DLR TR ER59-91, 1991 */
/*   Powell, M., A Hybrid Method for Nonlinear Equations. */
/*      Numerical Methods for Nonlinear Algebraic Equations, */
/*      P. Rabinowitz, Editor, Gordon and Breach, 1970 */

/* Life cycle: */
/* 1980 MAR  Burton S. Garbow, Kenneth E. Hillstrom, Jorge J. More, */
/*           Argonne National Laboratory. Minpack Project. */
/* 1991 OKT  G. Gramlich, DLR FF-DF: reverse communication interface */
/* 1993 SEP  D. Joos, DLR FF-DR    : iopt=3 introduced */
/* 1997 AUG  M. Otter, DLR FF-DR   : termination criteria improved, */
/*                                   multiple usage at the same time */
/*                                   (internal state is stored in */
/*                                   work arrays). */

/* Libraries required: */
/* BLAS1:    DNRM2 */

/* Common Blocks: */
/* The subroutines DYMNL* need common /DYMCB/ for the values of */
/*  the machine precision EPSMCH and largest magnitude GIANT. */

/* Example: */

/* The problem is to determine the values of X(1) and X(2), */
/* which solve the nonlinear system of equations */

/*             X(1)*X(2) - X(2)**3 - 1 = 0 */
/*             X(1)**2*X(2) + X(2) - 5 = 0 */

/* ********** */

/*       PROGRAM TEST */
/* C */
/* C     DRIVER FOR DYMNL EXAMPLE. */
/* C */
/*       INTEGER J,N,INFREV,IOPT,INFO,LDUM,IDUM(20) */
/*       PARAMETER (N=2,LDUM=35) */
/*       DOUBLE PRECISION TOL */
/*       DOUBLE PRECISION X(N),FVEC(N),DUM(LDUM),FJAC(N,N) */
/* C */
/*       DOUBLE PRECISION DNRM2, FNORM */
/* C */
/* C     Initialize machine dependent constants */
/*       CALL DYMMDP */
/* C */
/*       IOPT = 1 */
/*       INFREV=-1 */
/*       X(1)=1.D0 */
/*       X(2)=1.D0 */
/* C */
/* C     SET TOL TO ZERO, HENCE THE SQUARE ROOT OF THE MACHINE PRECISION 
*/
/* C     IS ASSUMED INTERNALY. */
/* C */
/*       TOL = 0.D0 */
/* C */
/* 10    CONTINUE */
/*       IF(INFREV .GE. 1 ) THEN */
/*          FVEC(1)=X(1)*X(2)-X(2)**3-1.D0 */
/*          FVEC(2)=X(1)**2*X(2)+X(2)-5.D0 */
/*       END IF */
/*       IF(INFREV .EQ. 2 ) THEN */
/*          FJAC(1,1)=X(2) */
/*          FJAC(1,2)=X(1)-3.*X(2)**2 */
/*          FJAC(2,1)=2.*X(1)*X(2) */
/*          FJAC(2,2)=X(1)**2+1.D0 */
/*       END IF */
/*      CALL DYMNL(INFREV,IOPT,N,X,FVEC,FJAC,TOL,INFO,DUM,LDUM,IDUM,0.0D0)
*/
/*       IF (INFREV .GE. 1) GO TO 10 */
/*       FNORM = DNRM2(N,FVEC,1) */
/*       WRITE (6,1000) FNORM,INFO,(X(J),J=1,N) */
/*  1000 FORMAT (5X,' FINAL L2 NORM OF THE RESIDUALS',E15.7 // */
/*      *        5X,' INFO PARAMETER',16X,I10 // */
/*      *        5X,' FINAL APPROXIMATE SOLUTION' // (5X,3E15.7)) */
/*       STOP */
/*       END */

/* ************ */

/*      FINAL L2 NORM OF THE RESIDUALS  0.1434642E-10 */

/*      INFO PARAMETER                         1 */

/*      FINAL APPROXIMATE SOLUTION */

/*      0.2000000E+01  0.1000000E+01 */

/* ************ */

/* Results obtained with different compilers or machines may be slightly 
*/
/* different. */
/* -----------------------------------------------------------------------
 */
} /* dymnl_ */

DYMOLA_STATIC void dymnl1_Fast(integer*infrev, const integer iopt, const integer n, doublereal*x, 
			 doublereal*fvec, doublereal*fjac, const integer ldfjac, const doublereal xtol, 
			 const doublereal xtoldesired, const integer maxfev, const integer ml, 
			 const integer mu, const doublereal*nominalx, doublereal*diag, 
			 doublereal *diag2,
			 const doublereal factor, integer*info, integer*nfev, integer*njev, 
			 doublereal *r__, const integer lr,doublereal* qtf, 
			 doublereal * wa1, doublereal*wa2, doublereal*wa3, doublereal*wa4, 
			 const integer jopt, doublereal*dsave, integer*isave,const char*const*varnames,
			 integer*ipivots,int printEvent,int inJacobian)
{
    /* Initialized data */

    static const doublereal one
			 = 1.;
    static const doublereal p1
			 = .1;
    static const doublereal p5
			 = .5;
    static const doublereal p001
			 = .001;
    static const doublereal p0001
			 = 1e-4;
    static const doublereal zero
			 = 0.;

    /* System generated locals */
    integer fjac_dim1, fjac_offset;
    doublereal d__1, d__2;

    /* Local variables */
    integer i__, j, l, jm1, iwa[1];
    doublereal sum;
    logical sing;
    integer iter, irev;
    doublereal temp;
    doublereal delta;
    logical jeval;
    integer ncsuc;
    doublereal ratio, fnorm, pnorm, fnorm1;
    integer nslow1, nslow2, ncfail;
    doublereal actred, prered;
    integer ifirst;
    doublereal ftolabs, xtolabs;
	logical debugNonLinear;
	logical debugNonLinearDetails;

/*    factor is a positive input variable used in determining the */
/*      initial step bound.  this bound is set to the product of */
/*      factor and the euclidean norm of diag*x if nonzero, or else to */
/*      factor itself.  in most cases factor should lie in the */
/*      interval (.1,100.).  100. is a generally recommended value. */

/*    info is an integer output variable. info is set as follows. */

/*      info = 0  improper input parameters. */

/*      info = 1  relative error between two consecutive iterates is */
/*                at most xtol. */

/*      info = 2  number of calls to fcn has reached or exceeded */
/*                maxfev. */

/*      info = 3  xtol is too small.  no further improvement in the */
/*                approximate solution x is possible. */

/*      info = 4  iteration is not making good progress, as measured */
/*                by the improvement from the last five jacobian */
/*                evaluations. */

/*      info = 5  iteration is not making good progress, as measured */
/*                by the improvement from the last ten iterations. */

/*      sections 4 and 5 contain more details about info. */

/*    nfev is an integer output variable set to the number of */
/*      function evaluations. */

/*    njev is an integer output variable set to the number of */
/*      jacobian evaluations. (if iopt=2, then njev is set zero.) */

/*    r is an output array of length lr which contains the upper */
/*      triangular matrix produced by the qr factorization of the */
/*      final approximate jacobian, stored rowwise. */

/*    lr is a positive integer input variable not less than */
/*      (n*(n+1))/2. */

/*    qtf is an output array of length n which contains the vector */
/*      (q transpose)*fvec. */

/*    wa1, wa2, wa3, and wa4 are work arrays of length n. */

/*    dsave is a double precision work array of length 20. */
/*    isave is an integer work array of length 20. */

/* subprograms called */
/* ------------------ */
/* dymnl-supplied    dymnl2,dymnl3,dymnl4,dymnl5,dymnl6,dymnl7 */
/* blas1-supplied    dnrm2 */
/* fortran-supplied  min0,dmax1,dmin1,dabs */

/* author            guenter m. gramlich, d. joos */
/* date              23.05.91, oberpfaffenhofen */
/* version           1.0       23.09.93 */
/* -----------------------------------------------------------------------
 */
/* Changed substantially 2002-02-11 by Hans Olsson Dynasim */

/* diag2 is used for storing residual*/
/* The code automatically rescales residuals */
/* This is necessary in combination with tearing, where */
/* very badly conditioned systems of equations can occur */

/* If nominalx is non-nil it represents nominal values of*/
/* the variables. If nil we assume that nominal=1 for all variables */
/**/ 



/* C */
/*     .. executable statements .. */

	{
		debugNonLinear=printEvent && !!(printEvent&(1<<8));
		debugNonLinearDetails = debugNonLinear && (printEvent&(1<<7));
	}
    /* Parameter adjustments */
    --wa4;
    --wa3;
    --wa2;
    --wa1;
    --qtf;
    --diag;
    --fvec;
    --x;
    fjac_dim1 = ldfjac;
    fjac_offset = fjac_dim1 + 1;
    fjac -= fjac_offset;
    --r__;
    --dsave;
    --isave;

    /* Function Body */

/* ***first executable statement  dymnl1 */

    xtolabs = xtol;
    ftolabs = xtol;

#ifdef fastnl
#define	actred (dsave[1])
#define delta (dsave[2])
#define fnorm (dsave[3])
#define fnorm1 (dsave[4])
#define pnorm (dsave[5])
#define prered (dsave[6])
#define ratio (dsave[7])
#define sum (dsave[8])
#define temp (dsave[9])
#define xnorm (dsave[10])
	/* dsave 11..13 used by dymnl5 */
	/* dsave 14 used for norm of Jacobian */

	/* Handling of scalar systems: */
	/* 1. Localize an interval enclosing a sign-change (can sort of be generalized to non-scalar) */
	/* 2. Update f(dsave[18])<=0, f(fdsave[19])>=0 */
	/*    Active if isave[20]&1   isave[20]&2   */
	/* 3. If new x will leave dsave[18].. dsave[19] or new x < 0.1*|dsave[18]-save[19]| => bisect. */
#define	i__ (isave[3])
#define iter (isave[4])
#define j (isave[6])
#define jm1 (isave[7])
#define l (isave[8])
#define ncfail (isave[9])
#define ncsuc (isave[10])
#define nslow1 (isave[11])
#define nslow2 (isave[12])
#define irev (isave[13])
#define ifirst (isave[14])
#define jeval (isave[15])
#define sing (isave[16])
#endif
    if (*infrev < 1) {
		/*        ... initial call; initialize local variables */
		actred = 0.;
		delta = 0.;
		fnorm = 0.;
		fnorm1 = 0.;
		pnorm = 0.;
		prered = 0.;
		ratio = 0.;
		sum = 0.;
		temp = 0.;
		*nfev = 0;
		*njev = 0;
		i__ = 0;
		iter = 0;
		l = 0;
		ncfail = 0;
		ncsuc = 0;
		nslow1 = 0;
		nslow2 = 0;
		irev = 0;
		ifirst = 0;
		jeval = FALSE_;
		sing = FALSE_;
		isave[19]=0;
		isave[20]=0;
#ifdef DYMOSIM
		{
			extern int Check5(char*key);
			if (!Check5("")) {
				DymosimMessage("");
				*info=5;
				goto L370real;
			}
		}
#endif
		if (!(*infrev <= -20 && *infrev>=-40)) {
			j = 0;
			{int i;for(i=1;i<=n;++i)
				diag2[i-1]=1;}
			dsave[14]=0;
			iwa[0] = 0;
			jm1 = 0;
		} else {
			j = isave[6];
			jm1 = isave[7];
			iwa[0] = isave[5];
			if (*infrev<-30) {*infrev+=10;if (n!=1) isave[20]=1;}
		}
    } else {
#if 0
		DymosimMessageDouble("F:",fvec[1]);
#endif
		/*        ... get values of local variables from last call */
#ifndef fastnl
		actred = dsave[1];
		delta = dsave[2];
		fnorm = dsave[3];
		fnorm1 = dsave[4];
		pnorm = dsave[5];
		prered = dsave[6];
		ratio = dsave[7];
		sum = dsave[8];
		temp = dsave[9];
		/* dsave 11..13 used by dymnl5 */
		/* dsave 14 used for norm of Jacobian */
		/* dsave 15-17 used for singular systems */
#endif
		*nfev = isave[1];
		*njev = isave[2];
#ifndef fastnl
		i__ = isave[3];
		iter = isave[4];
#endif
		iwa[0] = isave[5];
#ifndef fastnl
		j = isave[6];
		jm1 = isave[7];
		l = isave[8];
		ncfail = isave[9];
		ncsuc = isave[10];
		nslow1 = isave[11];
		nslow2 = isave[12];
		irev = isave[13];
		ifirst = isave[14];
		if (isave[15] == 0) {
			jeval = FALSE_;
		} else {
			jeval = TRUE_;
		}
		if (isave[16] == 0) {
			sing = FALSE_;
		} else {
			sing = TRUE_;
		}
#endif
		/* isave 17..18 used by dymnl5 */
		/* isave 19 used for singular systems and zero start values */
    }
#if defined(DYMOSIM) || defined(GODESS)
	{
		extern int ds_res_int_poll_fast(int,integer n,const doublereal*sol,const char*const*varnames);
		if (ds_res_int_poll_fast(*nfev,n,&x[1],varnames)) {
			*info=5;
			*nfev=-1;
			goto L370a;
		}
	}
#endif
	{
		debugNonLinear=!!(printEvent&(1<<8));
		debugNonLinearDetails = debugNonLinear && (printEvent&(1<<7));
	}
#if DEBUG_NL
	DymosimMessageInt("Infrev",*infrev);
#endif
    switch ((int)*infrev) {
	case 1:  goto L30; /* got initial value, also used if initial had singular Jac */
	case 2:  goto L50; /* got Jacobian */
	case 3:  goto L70; /* got column for Jacobian */
	case 4:  goto L260;/* got new value in inner loop */
	case 5:  goto L5; /* got value for line-search */
	case 6:  goto L6;
	case 11:case 12:case 13:case 14:case 15:case 16: case 17:case 18:case 19:goto L260;
	case 88: goto L88; /* for treating singular Jacobian */

	case 101: *infrev=1; 
		/* Failed to evaluate rhs, set err to max*/
		/* Failed to evaluate original */
		++*nfev;fnorm=DBL_MAX; 
		if (isave[19]==0) {
			if (debugNonLinearDetails) {
				 DymosimMessage("Trying to modify values.");
			}
			for(j=1;j<=n;++j) {
				double dx;
				dx=fabs(x[j]);
				if (nominalx!=0) dx+=nominalx[j-1];
				else dx+=1;
				x[j]+=dx*1e-8;
			}
			isave[19]=-1;
			goto L999;
		} else if (isave[19]==-1) {
			/* Failed to get any residual at all*/
			for(j=1;j<=n;++j)
				fvec[j]=DBL_MAX;
			*info = 5;
			goto L370; 
		}
		goto L40;
	case 102: *infrev=2; 
		goto L70; /* Failed to evalaute analytical Jacobian, try numeric one */
	case 103: *infrev=3; 
		/* Failed during evaluation of numeric Jacobian */ 
			  if (debugNonLinearDetails) {
				  DymosimMessage("Error for numeric Jacobian");
			  }
	case 104: *infrev=4;     
		/* Failed to evaluate rhs, set err to max*/
			  {
				  for (j = 1; j <= n; ++j) {
					  wa2[j] = wa4[j];
					  wa4[j] = fvec[j]*diag2[j-1];
					  fvec[j] = wa2[j];
					  wa2[j] = x[j];
					  x[j] = wa3[j];
					  /* L270: */
				  }
			  }
			  if (debugNonLinearDetails) {
				  DymosimMessage("Handled error by setting error to max");
			  }
		++(*nfev);fnorm1=DBL_MAX; goto L260a;
	case 105:
		/* Failed for debug line search. Give up */
		*infrev=5;
		goto L370real;
	case 106:
		goto L6b;
	case 188: 
		/* Failed to evaluate rhs, set err to max*/
		*infrev=88; ++*nfev;fnorm=DBL_MAX; goto L88a;
	default:
		/* First call continue on */
		break;
    }

    *nfev = 0;
    *njev = 0;

/* check the input parameters for errors. */

/* ...exit */
    *info = 0;
    if (iopt < 1 || iopt > 2 || n <= 0 || xtol < zero || maxfev <= 0 || 
	    ml < 0 || mu < 0 || factor <= zero || ldfjac < n || lr < n *
	     (n + 1) / 2) {
	goto L370;
    }

    if (nominalx !=0) {
		for (j = 1; j <= n; ++j) {
			/*   ** .........exit */
			if (nominalx[j-1] <= zero) {
				goto L370;
			}
		/* L10: */
		}
    }

    /* Note: Added special code 2001-10-17 by Hans:   */
	/* Modified 2003-04-14 */
    /* If we start with infrev between -10 and -40 we */
    /* assume that the residue is already given.      */
    /* If infrev is furthermore between -11 and -40   */
    /* we also assume that the Jacobian is available. */
	/* If infrev is -11-i or -20-i we use i iterations */
/* evaluate the function at the starting point */
/* and calculate its norm. */
	/* If it is -30-i we use i iterations, and forbid new Jacobians, that is handled by setting isave[20] */

/* ...return */

/* initialize iteration counter and monitors. */
	*nfev = 0;
    iter = 1;
    *info = 10;
    ncsuc = 0;
    ncfail = 0;
    nslow1 = 0;
    nslow2 = 0;
    ifirst = 1;

	if (*infrev >=-40 && *infrev<= -10) goto L30;
    *infrev = 1;
    goto L999;

L30:
    ++*nfev;
 
/* beginning of the outer loop. */

/* note that we proceed even with zero residual, we will otherwise get failures */
/* already have residual */
L40:
/* begin block permitting ...exits to 90 */
	if (n==1) {
		if (fvec[1]<0) {
			dsave[18]=x[1];
			isave[20]|=1;
		} else if (fvec[1]>0) {
			dsave[19]=x[1];
			isave[20]|=2;
		}
	}
    jeval = TRUE_;

	if (*infrev <= -20 && *infrev>=-30) {
		/* Scale residual and compute new norm */
		for (i__ = 1; i__ <= n; ++i__) 
			fvec[i__]*=diag2[i__-1];
		fnorm = dnrm2_Fast1(n,&fvec[1]);
		for (i__ = 1; i__ <= n; ++i__) {
			double sumtemp;
			sumtemp=0;
			for(j = 1;j <= n;++j) {
				sumtemp+=fjac[i__*ldfjac+j]*fvec[j];
			}
			qtf[i__]=sumtemp;
		}
		for (j = 1; j <= n; ++j) {
			wa3[j] = diag[j] * x[j];
			/* L110: */
		}
		delta = factor; 
		/* Goto inner loop */
		if (*infrev==-20) {
			*infrev=2;
		} else {
			*infrev=-10-*infrev;
		}
		goto L230;
	}
/* calculate the jacobian matrix. */
#if 0
	goto L60;
#endif
    if (iopt == 2 || isave[19]==-1) {
	goto L60;
    }

/* user supplies jacobian */

/* ...return */
    if (*infrev >=-20 && *infrev< -10) {
		if (*infrev==-10 || *infrev==-11 || *infrev==-20)
			*infrev=2;
		else 
			*infrev=-*infrev-1;
		goto L50;
	}
	*infrev = 2;
    goto L999;

L50:
    ++(*njev);
#ifdef CHECK_JACOBIAN_APPROXIMATION
	{
		goto L65b;
	}
#endif
    goto L80;

L60:
#ifdef CHECK_JACOBIAN_APPROXIMATION
	goto L65b;
#endif
/* code approximates the jacobian */

    if (jopt == 1 && ifirst == 1) {
/* for the first time an approx. of the jacobian is already provided 
*/
	ifirst = 0;
	goto L80;
    }
#ifdef CHECK_JACOBIAN_APPROXIMATION
L65b:
#endif
    irev = 0;
	if (debugNonLinearDetails) {
		DymosimMessage("Compute numeric Jacobian");
	}
L70:
    dymnl5_Fast(&irev, n, &x[1], &fvec[1], &fjac[fjac_offset], ldfjac, ml, mu, 
	    nominalx, &wa1[1], &wa2[1], &dsave[11], &isave[17], nslow2>=6,debugNonLinear,varnames,
		iopt==1);
    if (irev == 1) {
/*   ** ...return */
	*infrev = 3;
	goto L999;

    }
/* Computing MIN */
    *nfev += Dymola_min(ml+mu+1,n);
L80:
	if (n==1) {
		if (!(isave[20]&8)) {
			if (fjac[1+1]>0) 
				isave[20]|=(8|16); /* Positive Jacobian */
			else if (fjac[1+1]<0) 
				isave[20]|=8; /* Negative Jacobian */
		} else if ((isave[20]&3)!=3) {
			/* Known sign of Jacobian, but no interval with sign-change */
			/* The Jacobian has switched sign, but not the residual. That indicates a singularity */
			if (fjac[1+1]!=0 && (fjac[1+1]>0 != !!(isave[20]&16))) {
				goto L359b;
			}
		}

	}
	/* Scale rows of new Jacobian to compensate for ill-conditioned Jacobians*/
	{
		int nrResLess=0;
	{
		int i,j;
#if DEBUG_NL
		if (LocalTime>0.550215) 
		DymosimMessageDoubleMatrix("J",&(fjac[1+1*fjac_dim1]),n,n,fjac_dim1);
#endif
		for(i=1;i<=n;++i) {
			double Jsum=0;
			if (nominalx!=0) {
				for (j=1; j<= n; ++j) {
					Jsum+=Dymola_abs(fjac[i+j*fjac_dim1])*(Dymola_abs(x[j])+nominalx[j-1]);
				}
			} else {
				for (j=1; j<= n; ++j) {
					Jsum+=Dymola_abs(fjac[i+j*fjac_dim1])*(Dymola_abs(x[j])+1);
				}
			}
			if (debugNonLinearDetails) {
				char str[100];
				sprintfC(str,"Residual scale row %d, J_sum=%g",i,Jsum);
				DymosimMessage(str);
			} else if (debugNonLinear&&(Jsum==0)) {
				char str[100];
				sprintfC(str,"Row %d of Jacobian is zero",i);
				DymosimMessage(str);
			}
			if (Jsum==0) {
				/* Count number of Jacobian zero-rows with exactly zero residual */
				/* If these correspond to the loss of singularity the system is consistently overdetermined */
				if (fvec[i]!=0) nrResLess=n+1;
				else nrResLess++;
			}
			Jsum += 1 ; /* Guard against singularity */
			{
				double invScale=1.0 / Jsum;
				diag2[i-1]=invScale;
				for(j=1;j<=n;++j) {
					fjac[i+j*fjac_dim1]*=invScale;
				}
			}
		}
	}
	/* Also scale residual */
	{
		int i;
		for(i=1;i<=n;++i) fvec[i]*=diag2[i-1];
	}
	fnorm = dnrm2_Fast1(n, &fvec[1]);
	if (fnorm!=fnorm || !(fnorm>=0)) fnorm=DBL_MAX;


/* compute the qr factorization of the scaled jacobian. */

	if (ipivots!=0)
		dymnl7_Fast(n, n, &fjac[fjac_offset], ldfjac, ipivots!=0, ipivots, n, &wa1[1], &
	    wa2[1], &wa3[1],varnames, &x[1], printEvent);
	else
    dymnl7_Fast(n, n, &fjac[fjac_offset], ldfjac, c_false, iwa, c__1, &wa1[1], &
	    wa2[1], &wa3[1],varnames, &x[1], printEvent);

	if (*infrev>=0) {
		/* HO 2000-06-27 */
		/* Handle singular matrix. */
		int is_singular;
		int i;
		is_singular=0;
		for(i=1;i<=n;i++) if (wa1[i]==0) is_singular++;
		if (debugNonLinearDetails||(debugNonLinear&&is_singular)) {
			if (is_singular)
				DymosimMessage("Singular Jacobian matrix");
			else {
				double dmi = 1.0;
				double dma = 1.0;
				for(i=1;i<=n;i++) {
					double wabs;
					wabs=fabs(wa1[i]);
					if (i==1 || wabs>dma) dma=wabs;
					if (i==1 || wabs<dmi) dmi=wabs;
				}
				DymosimMessageDouble("Condition estimate of Jacobian matrix",dma/dmi);
			}
		}

		if (is_singular && (fnorm!=0 && nrResLess!=is_singular) && isave[19]!=-1 && !(n==1 && (isave[20]&3)==3)) {
			/* Ignore the case where a scalar system has singular jacobian and we have both upper and lower*/
			/* HO 2000-08-07 only care about singular systems with non-zero residual */
			if (is_singular) {
				/* Singular Jacobian. How to proceed? */
				/* Ignoring would not be meaningful */
				/* Attempt to recompute Jacobian at different point */
				if (isave[19]<=0) {
					/* First time*/
					dsave[17]=fnorm;
					dsave[15]=1e-5;
					isave[19]=0;
				} else {
					x[isave[19]]-=dsave[16];
				}
				isave[19]++;
#if 0
				DymosimMessageInt("Isave:",isave[19]);
#endif
				if (isave[19]>n) {
					isave[19]=1;
					if (dsave[15]>0)
						dsave[15]=-dsave[15];
					else 
						dsave[15]=-dsave[15]*(fabs(dsave[15])<1.5e-4 ? 10 : 2);
					/* Sequence:   1e-5, 1e-4, 1e-3, 2e-3, 4e-3 8e-3, 16e-3, 32e-3, 64e-3, 128e-3 */
					/* Previously: 1e-5, 1e-4, 1e-3, 1e-2, 1e-1 */
					if (fabs(dsave[15])>.2) {
						isave[19]=-1; /* Go back to original once more */
						if (iopt != 2 && debugNonLinear)
							DymosimMessage("Using Numeric Jacobian");
					}
#if 0
					DymosimMessageDouble("Dsave:",dsave[15]);
#endif
				}
				if (isave[19]>=0) {
					dsave[16]=Dymola_max(1.0,fabs(x[isave[19]]))*dsave[15];
					x[isave[19]]+=dsave[16];
				}
				*infrev=1;
				goto L999;
			}
		}
	}
	}
#if 0
	DymosimMessageDouble("X: ", x[1]);
#endif
	if (isave[19]!=0 && isave[19]!=-1 && fnorm>dsave[17]) {
		/* Got non-singular matrix, but unfortunately the error is larger */
		/* We must now fix this because we will otherwise take another step that leads us back into trouble*/
		dsave[16]*=0.5;
		x[isave[19]]-=dsave[16]; /* Start midway */
		*infrev=88;
		if (debugNonLinear) {
			DymosimMessage("Have found non-singular matrix, but larger error");
		}
		goto L999;
	}
	goto L89;
L88:
	++*nfev;
	{
		int i;
		for(i=1;i<=n;++i) fvec[i]*=diag2[i-1];
	}
    fnorm = dnrm2_Fast1(n, &fvec[1]);
	if (fnorm!=fnorm || !(fnorm>=0)) fnorm=DBL_MAX;
L88a:
#if 0
	DymosimMessageDouble("New line search ",fnorm);
	DymosimMessageDouble("Saved norm",dsave[17]);
	DymosimMessageDouble("X: ", x[1]);
	DymosimMessageDouble("F: ", fvec[1]);
#endif
	if (fnorm>=dsave[17]) {
		if (fabs(dsave[16])>=1e-14*Dymola_max(1.0,x[isave[19]])) {
			dsave[16]*=0.5;
			if (fnorm>dsave[17])
				x[isave[19]]-=dsave[16];
			else
				x[isave[19]]+=dsave[16];
			*infrev=88;
#if 0
			DymosimMessageInt("Infrev",*infrev);
#endif
			goto L999;
		}
		/* Done, now we need a new Jacobian*/
		isave[19]=-1; /* Hope that it is not singular again*/
		{
			int i;
			for(i=1;i<=n;++i) 
				fvec[i]/=diag2[i-1];
		}
		goto L40;
	}
L89:
	isave[19]=0;

/* on the first iteration and if mode is 1, scale according */
/* to the norms of the columns of the initial jacobian. */

/* ...exit */
    if (iter != 1) {
	goto L120;
    }
    for (j = 1; j <= n; ++j) {
		if (nominalx) {
			diag[j] = 1.0/(Dymola_abs(x[j])+nominalx[j-1]);
		} else {
			diag[j] = 1.0/(Dymola_abs(x[j])+1);
		}
		/* L90: */
    }
/* on the first iteration, initialize the step bound delta. */

    delta = factor; 
	if (isave[19]!=0)
		delta=Dymola_min(delta,fabs(dsave[15]/10));
	/* Introduced by Hans Olsson, Dynasim 2000-03-24. Simpler than code below */
#if 0
	delta = *factor * xnorm;
/* Otter, July 1997: previously, this was */
/*      IF (delta.EQ.zero) delta = factor */
/*   However, when the starting point x is near zero, this */
/*   will result in a VERY TINY step at the beginning, which */
/*   is of course unreasonable. Therefore, the absolute tolerance */
/*   on x is used for a check, that the first step size does not */
/*   become too small (default of factor=100). */
    if (delta <= xtolabs) {
	delta = *factor;
    }
#endif
/* end Otter */

L120:

/* form (q transpose)*fvec and store in qtf. */

    for (i__ = 1; i__ <= n; ++i__) {
		qtf[i__] = fvec[i__];
		/* L130: */
    }
    for (j = 1; j <= n; ++j) {
		if (fjac[j + j * fjac_dim1] == zero) {
			goto L160;
		}
		sum = zero;
		for (i__ = j; i__ <= n; ++i__) {
			sum += fjac[i__ + j * fjac_dim1] * qtf[i__];
			/* L140: */
		}
		temp = -sum / fjac[j + j * fjac_dim1];
		for (i__ = j; i__ <=n; ++i__) {
			qtf[i__] += fjac[i__ + j * fjac_dim1] * temp;
			/* L150: */
		}
L160:
		/* L170: */
		;
    }

/* copy the triangular factor of the qr factorization into r. */

    sing = FALSE_;
    for (j = 1; j <= n; ++j) {
		l = j;
		for (i__ = 1; i__ <= j-1; ++i__) {
			r__[l] = fjac[i__ + j * fjac_dim1];
			l = l + n - i__;
			/* L180: */
		}
		r__[l] = wa1[j];
		if (wa1[j] == zero) {
			sing = TRUE_;
		}
		/* L200: */
    }

/* accumulate the orthogonal factor in fjac. */

    dymnl6_Fast(n, n, &fjac[fjac_offset], ldfjac, &wa1[1]);

/* rescale if necessary. */

    for (j = 1; j <= n; ++j) {
/* Computing MAX */
		double newDiag=0;
		if (nominalx!=0) {
			newDiag=1.0/(Dymola_abs(x[j])+nominalx[j-1]);
		} else {
			newDiag=1.0/(Dymola_abs(x[j])+1);
		}
		if ((nslow2<6)||newDiag<diag[j] )
			diag[j]=newDiag;
		/* remove the 1|| to only rescale if increasing, makes more sense close to singularities */
		/* Modified to nslow2<6 Hans 2004-08-20. To handle Dynasim #1245 and not influence most cases */
/* L210: */
    }

/* beginning of the inner loop. */

L230:

/* determine the direction p. in wa1. */

    dymnl4_Fast(n, &r__[1], lr, &diag[1], &qtf[1], delta, &wa1[1], &wa2[1], &wa3[
	    1],ipivots,inJacobian);
	if (n==1 && (isave[20]&3)==3) {
		double lower;
		double upper;
		lower=dsave[18];
		upper=dsave[19];
		if (lower>upper) {
			double xxx;
			xxx=lower;
			lower=upper;
			upper=xxx;
		}
		if (x[1]-wa1[1]<lower || x[1]-wa1[1]>upper) {
			/* Technically speaking we should integrate this much more fully, etc */
			/* However, it worked for several examples */
			if (debugNonLinearDetails) {
				DymosimMessage("Using bisection");
			}
			if (upper-lower<1e-13 && (isave[20]&64)==0) {
				/* Too small interval - must be a discontinuity at the edge, Dynasim #7044*/
				if (x[1]-wa1[1]<lower) 
					wa2[1] = lower;
				else 
					wa2[1] = upper;
				isave[20]|=64; /* Only do it once*/
			} else {
				wa2[1] = (upper+lower)/2;
				isave[20]&=~64; /* Allow other code to be activated */
			}
			/* Ensure that we exactly hit limit. Important since test for out-of-bounds use exact limits*/
			wa1[1] = wa2[1] -x[1];
			wa3[1] = diag[1] * wa1[1];
			isave[20]|=4;
			goto LhaveNewP;
		}
	}

/* store the direction p and x + p. calculate the norm of p. */
#if DEBUG_NL
	if (LocalTime>0.550215) 
		DymosimMessageDoubleMatrix("x",&x[1],1,n,1);
	if (LocalTime>0.550215) 
		DymosimMessageDoubleMatrix("-dx",&wa1[1],1,n,1);
#endif
    for (j = 1; j <= n; ++j) {
		wa1[j] = -wa1[j];
		wa2[j] = x[j] + wa1[j];
		wa3[j] = diag[j] * wa1[j];
		/* L240: */
    }
LhaveNewP:
#if 0
	DymosimMessageDouble("p:",-wa1[1]);
	DymosimMessageDouble("x:",wa2[1]);
#endif
    pnorm = dnrm2_Fast1(n, &wa3[1]);

/* on the first iteration, adjust the initial step bound. */

    if (iter == 1) {
	delta = Dymola_min(delta,pnorm);
    }

/* evaluate the function at x + p and calculate its norm. */

/* ...return */
    for (j = 1; j <= n; ++j) {
		wa3[j] = x[j];
		x[j] = wa2[j];
		wa4[j] = fvec[j];
		/* L250: */
    }
	if (*infrev<10 || *infrev>19)
		*infrev = 4;
	else {
		--*infrev;
		if (*infrev==10) {
			*infrev = 0;
			*info=1;
		}
	}
    goto L999;

L260:
	if (n==1) {
		if (fvec[1]<0) {
			dsave[18]=x[1];
			isave[20]|=1;
		} else if (fvec[1]>0) {
			dsave[19]=x[1];
			isave[20]|=2;
		}
	}
#if DEBUG_NL
	if (LocalTime>0.550215) 
		DymosimMessageDoubleMatrix("f",&fvec[1],1,n,1);
#endif

    for (j = 1; j <= n; ++j) {
	wa2[j] = wa4[j];
	wa4[j] = fvec[j]*diag2[j-1];
	fvec[j] = wa2[j];
	wa2[j] = x[j];
	x[j] = wa3[j];
/* L270: */
    }
    ++(*nfev);
    fnorm1 = dnrm2_Fast1(n, &wa4[1]);

	if (debugNonLinearDetails)
			DymosimMessageDouble("Scaled residual",fnorm1);
	if (fnorm1!=fnorm1 || !(fnorm1>=0)) fnorm1=DBL_MAX;
L260a:
	/* test for successful iteration. */
    if (n==1 && (isave[20]&7)==7) {

		/* Bisected. Always accept if fnorm1<DBL_MAX*/
		isave[20]&=~4;
		/* Always accept if old value is not inside new limits */
		if (fnorm1<DBL_MAX && fnorm1>=fnorm) {
			double lower=dsave[18];
			double upper=dsave[19];
			if (lower>upper) {
				double xxx=upper;
				upper=lower;
				lower=xxx;
			}
			if (pnorm>diag[1]*(upper-lower))
				pnorm = diag[1]*(upper-lower); /* Always. */
			if (wa2[1]<lower || wa2[1]>upper) {
				goto L320; /* Always accept it - even if no progress */
			}
		}
		/* Note: It is important to have fnorm1==fnorm handled by this in order to deal with singular Jacobian */
	}
/* compute the scaled actual reduction. */

    actred = -one;
    if (fnorm1 < fnorm) {
/* Computing 2nd power */
	d__1 = fnorm1 / fnorm;
	actred = one - d__1 * d__1;
    }

/* compute the scaled predicted reduction. */
/* Compute: Q'(f+A*wa1)=Q'(f+QRP'*wa1)=Q'f+RP'*wa1, special formula for P'*wa1 */


    l = 1;
    for (i__ = 1; i__ <= n; ++i__) {
	sum = zero;
	for (j = i__; j <= n; ++j) {
		sum += r__[l] * (ipivots? wa1[ipivots[j-1]]: wa1[j]);
	    ++l;
/* L280: */
	}
	wa3[i__] = qtf[i__] + sum;
/* L290: */
    }
    temp = dnrm2_Fast1(n, &wa3[1]);
    prered = zero;
    if (temp < fnorm) {
/* Computing 2nd power */
	d__1 = temp / fnorm;
	prered = one - d__1 * d__1;
    }

/* compute the ratio of the actual to the predicted */
/* reduction. */

    ratio = zero;
    if (prered > zero) {
	ratio = actred / prered;
    }
	if (debugNonLinearDetails) {
		if (prered>=0) {
			DymosimMessageDouble("Old Residual",fnorm);
			DymosimMessageDouble("Predicted relative decrease",(1-sqrt(1-prered)));
			if (fnorm1>=DBL_MAX)
				DymosimMessage("Failed to evaluate");
			else
				DymosimMessageDouble("Actual relative decrease",1-fnorm1/fnorm);
		}
	}
/* update the step bound. */

    if (ratio >= p1) {
	goto L300;
    }
    ncsuc = 0;
    ++ncfail;
	delta = p5 * delta; 
    goto L310;

L300:
    ncfail = 0;
    ++ncsuc;
    if (ratio >= p5 || ncsuc > 1) {
/* Computing MAX */
	d__1 = delta, d__2 = pnorm / p5;
	delta = Dymola_max(d__1,d__2);
    }
    if ((d__1 = ratio - one, Dymola_abs(d__1)) <= p1) {
	delta = pnorm / p5;
    }
L310:


    if (ratio < p0001) {
	goto L330;
    }
L320:
/* successful iteration. update x, fvec, and their norms. */

    for (j = 1; j <= n; ++j) {
		x[j] = wa2[j];
		wa2[j] = diag[j] * x[j];
		fvec[j] = wa4[j];
		/* L320: */
    }
    fnorm = fnorm1;
    ++iter;
    *info = 10;
L330:

/* determine the progress of the iteration. */

    ++nslow1;
    if (actred >= p001) {
	nslow1 = 0;
    }
    if (jeval) {
	++nslow2;
    }
    if (actred >= p1) {
	nslow2 = 0;
    }

/* test for convergence. */

/* Otter, July 1997: original convergence test: */
/*  if (delta .le. xtol*xnorm .or. fnorm .eq. zero) info = 1 */

/*  where xnorm: norm(x) */
/*        fnorm: norm(residue) */
/*        xtol : termination occurs when the relative error between two */
/*               consecutive (scaled) iterates is at most xtol */
/*       delta: bound on the maximum allowed step, i.e., on norm(x_new - x_old).*/

/*  The original convergence test fails, if x=0 is the solution, */
/*  because in such a case delta must be zero. */
/*  Furthermore, the solver may terminate even if the residue */
/*  is large, i.e., when a local minimum is reached. */
/*  In such a case, an error message should be given. */

/*  The assumption has to be made that the equations are well */
/*  scaled (a small residue is desired). Otherwise, a local minimum */
/*  cannot be distinguished from the solution. */

/*  The new convergence criteria requires, that the residue is */
/*  small enough AND that the relative error between two */
/*  consecutive iterates is also small enough (otherwise, the */
/*  solution vector may not reach the required accuracy). */
/*  If "delta" is used in the test, problems occur for linear */
/*  systems because for linear systems the solution is reached */
/*  in one Newton step. However, "delta" is only modestly reduced */
/*  (delta=0.1*delta). As a result, the solver reports */
/*  "not making good progress". Therefore, norm(x_new - x_old) = pnorm */
/*  has to be used in the convergence test. */

#if 0
    if (fnorm <= ftolabs) {
		doublereal pnormmax;
/*        -- Residue is small enough. */
/*        -- Check, whether the required precision in x is */
/*        -- reached. If not, improve the solution further */
	if (xnorm > xtolabs) {
	    pnormmax = *xtol * xnorm;
	} else {
	    pnormmax = xtolabs;
	}

	if (pnorm <= pnormmax || delta <= pnormmax) {
	    *info = 1;
	}
    }
/* end Otter */
#else
	/* Alternative by Hans, 1999-11-12 */
#if 0
	if (debugNonLinearDetails) {
		/* Test to see that it works */
		char str[200];
		sprintfC(str,"%d xtol %.10g %.10g Pnorm %.10g delta %.10g xnorm %.10g fnorm %.10g",
			inJacobian_&&*inJacobian_,
			xtol,xtoldesired,pnorm,delta,pnorm,fnorm);
		dymosimmessage_(str,strlen(str));
	}
#endif
	if ((pnorm<=xtoldesired && fnorm<=xtoldesired)) { 
		
		/* Duplicated below with xtol instead of xtoldesired */
		/* (inJacobian_ && *iopt!=2) || Immediate return if in Jacobian and user-supplied Jacobian (i.e. non-approximate) */
		/* Still not good enough */
		/* dsave[14] computed by dymnl5 is the norm of the Jacobian (actually sum of absolute values). */
		/* If the system of equations is consistent we should have f(x)=J*delta x
		   and thus |f|<=|J|*|delta x| in the solution*/
		/* Only use pnorm because of discussion above */
		/* The scaling of residuals ensure that |J|=1 */
		*info = 1;
	}
#endif

/* .........exit */
    if (*info == 1) {
	goto L370;
    }

/* tests for termination and stringent tolerances. */

    if (*nfev >= maxfev) {
	*info = 2;
    }
/* Computing MAX */
    d__1 = p1 * delta;
    if (p1 * Dymola_max(d__1,pnorm) <= DBL_EPSILON) {
		if (debugNonLinear) {
			DymosimMessageDouble("Cannot reach precision since step bound is",p1 * Dymola_max(d__1,pnorm));
		}
		*info = 3;
    }
    if (nslow2 == 10) {
		if (debugNonLinear) {
			DymosimMessage("Too many Jacobian evaluations due to slow iterations");
		}
	*info = 4;
    }
    if (nslow1 == 45) { /* Increased 1999-12-03 and 2001-09-13 by Hans to avoid premature stops*/
		if (debugNonLinear) {
			DymosimMessage("Too many slow iterations with no progress");
		}
	*info = 5;
    }
    if (*info>=2 && *info<=5 && (pnorm<=xtol && fnorm<=xtol)) {
		/* Not as good as desired but good enough, better using it than failing*/
		/* Introduced by H. Olsson 2000-06-29 */
		*info=1;
	}
/* .........exit */
    if (*info >= 1 && *info <= 5) {
	goto L370;
    }

/* criterion for recalculating jacobian */

/* ...exit */
    if (ncfail ==3) {
		if (n!=1 && isave[20]) {
			*info = 4; /* Cannot compute Jacobian */
			goto L370;
		}
	goto L360;
    }

/* calculate the rank one modification to the jacobian */
/* and update qtf if necessary. */
#if 0
	if (debugNonLinear) {

		/* Test to see that it works */
		char str[200];
		sprintfC(str,"%d ratio %.10g actred %.10g fnorm %.10g fnorm1 %.10g",
			inJacobian_&&*inJacobian_,ratio,actred,fnorm,fnorm1);
		dymosimmessage_(str,strlen(str));
	}
#endif
	if ((actred>-1 || fnorm1<10*fnorm) && fnorm>5e-12) {
#if 0
		if (fnorm<1e-6 && ratio >= p0001) {
			for(j=1;j<=n;++j) {
				double sum=0.0;
				for(i__ = 1; i__ <=n;++i__) sum += fjac[i__ + j * fjac_dim1] * wa4[i__];
				qtf[j]=sum;
			}
			goto L230;
		}
#endif
		/* Only update Jacobian if sensible value. */
		/* Note that fnorm might already have been assigned the value of */
		/* fnorm1 if actred>=0 */

		/* Note that fjac is q as a dense matrix */

		/* Not influenced by permutations: */
		/* In: wa4=f_new */
		/* In: wa3=Q'*f_old+RP'*delta_x [already computed] */
		/* Out: wa2 = (Q'*f_new-(Q'*f_old+RP'*delta_x))/pnorm */
		/* Out: wa1 diagonal scaling of wa1 (=delta_x) */
		for (j = 1; j <= n; ++j) {
			sum = zero;
			for (i__ = 1; i__ <=n; ++i__) {
				sum += fjac[i__ + j * fjac_dim1] * wa4[i__];
			}
			wa2[j] = (sum - wa3[j]) / pnorm;
			wa1[j] = diag[j] * (diag[j] * wa1[j] / pnorm);
			if (ratio >= p0001) {
				qtf[j] = sum;
			}
		}
		
		/* compute the qr factorization of the updated jacobian. */

		/* Note that dymnl3_Fast talks about lower-trapezoidal matrices.*/
		/* */
		/* Thus it treats r as r^T, and solves the transposed problem.*/

		/*Therefore dymnl3_Fast computes a transformation such that:*/
		/*(r^T+wa1*wa2^T)*q is lower triangular, i.e.*/
		/* q^T*(r+wa2*wa1^T) is upper triangular. */

		/*This implies that diagonal scaling of x is a column scaling of r,*/
		/*and the scaling is equivalent to multiplying elements of r by the same scaling.*/
		/**/
		/* A permutation matrix can be seen as replacing r by rp^T in the formula, i.e. */
		/* seeing that without permutation r=q^T*A and thus we operate on */
		/* q^T*(q^T*A+wa2*wa1^T) */ 
		/* With QR with permutations A*P=q*r i.e. q^T*A=r*P^T => */
		/* q^T*(r*P^T+wa2*wa1^T) */
		/* Multiplying from right with P gives: */
		/* q^T*(r+wa2*wa1^T*P)=q^T*(r+wa2*(P^T*wa1)^T */
		dymnl3_Fast(n, n, &r__[1], lr, &wa1[1], &wa2[1], &wa3[1], &sing, ipivots);
		if (n==1 && (isave[20]&3)!=3 && (nslow1>4)) {
			goto L359; /* We have not gotten anywhere. Try a simple line search */
		}
		if (sing) {
			goto L360; /* New Jacobian since it became singular */
		}
		dymnl2_Fast(n, n, &fjac[fjac_offset], ldfjac, &wa2[1], &wa3[1]);
		dymnl2_Fast(c__1, n, &qtf[1], c__1, &wa2[1], &wa3[1]);
		
	}
/* end of the inner loop. */

    jeval = FALSE_;
    goto L230;

	/* Special case for scalar */
L359:
	
	/* restore scaling of residuals */
	{
		int i;
		for(i=1;i<=n;++i) {
			fvec[i]/=diag2[i-1];
		}
	}
L359b:
	if (debugNonLinearDetails) {
			DymosimMessage("Line search to find interval with sign change");
	}
	/* This should be examined in more detail */
	/* But is only used for scalar non-linear equations that have problem with convergence. */
	/* Line search to find better value */
	for (j = 1; j <= n; ++j) {
		wa3[j] = x[j];
		x[j] = x[j]+wa1[j];
		wa4[j] = fvec[j];
		/* L250: */
    }
	/* Perform line search to find problem */
	*infrev = 6;
	goto L999;
	
L6:
	/* Found new value */
	if (n==1) {
		if (fvec[1]==0) {
			*info=1;
			goto L370;
		} else if (fvec[1]<0) {
			dsave[18]=x[1];
			isave[20]|=1;
		} else if (fvec[1]>0) {
			dsave[19]=x[1];
			isave[20]|=2;
		}
	}
	/* We can then use bisection when necessary */
	if ((isave[20]&3)==3) {
		goto L40;
	}

L6b: 
	/* This evaluation was not good enough. Continue with line-search */
	for(j=1;j<=n;++j) {
		wa1[j]=-1.2*wa1[j];
		if (fabs(wa1[j])>1e30*(nominalx ? nominalx[j-1] : 1)) {
			goto L6c;
		}
		x[j]=wa3[j]+wa1[j];
	}
	*infrev = 6;
	goto L999;
L6c:
	/* No progress and too far from start. Restore to original. */
	for(j=1;j<=n;++j) {
		x[j]=wa3[j];
		fvec[j]=wa4[j]*diag2[j-1];
	}
	*info=5;
	if (n==1)
		wa1[1]=1e-3*(fabs(x[1])+(nominalx?nominalx[0]:0));
	goto L370;
L360:

/* end of the outer loop. */

	/* restore scaling of residuals */
	{
		int i;
		for(i=1;i<=n;++i) {
			fvec[i]/=diag2[i-1];
		}
	}
    goto L40;

L370:
	if (!debugNonLinear || *info==1) {
		goto L370real;
	}
L370a:
	/* Decrease scale */
	for(j=1;j<=n;++j) wa1[j]*=0.001;
	DymosimMessage("Line search: DX-norm scaled-residua-norm residual (unscaled)");
	if (n==1 && (isave[20]&3)==3) {
		char str[200];
		double lower;
		double upper;
		DymosimMessage("Had interval with sign-change and still failed");
		lower=dsave[18];
		upper=dsave[19];
		sprintfC(str,"Interval: [%g %g]", lower, upper);
		DymosimMessage(str);
	}
	dymosimmessagedoubles_("Search direction ",&wa1[1],&n,strlen("Search direction"));
	DymosimMessage("To investigate the properties of the function, you can plot the ");
	DymosimMessage("function in the search direction by pasting the following ");
	DymosimMessage("commands in the Dymola command window:");
	DymosimMessage("  Amat={<...>};");
	DymosimMessage("  plotArray(Amat[:,1],Amat[:,2],-1);");
	DymosimMessage("If the graph has discontinuities, local minima above zero,");
	DymosimMessage(" and/or knees this explains the problem.");
	DymosimMessage("");
	DymosimMessage("Amat={");
	if (*nfev!=-1) {
	{
		/* Write initial ones */
		int more=1;
		char str[200];
		char*buf;
		buf=str+sprintfC(str,"{ %.16g, %.16g",0.0, fnorm);
		for(j=1;j<=n;++j) {
			buf=buf+sprintfC(buf,", %.16g",fvec[j]);
			if (j%5==3) {
				DymosimMessage(str);
				buf=str;
			}
		}
		buf=buf+sprintfC(buf,"%s",more?"},":"}};");
		DymosimMessage(str);
	}
	*nfev=0;
	}
L370b:
	for (j = 1; j <= n; ++j) {
		wa3[j] = x[j];
		x[j] = x[j]+wa1[j];
		wa4[j] = fvec[j];
		/* L250: */
    }
	/* Perform line search to find problem */
	*infrev = 5;

	goto L999;
L5:
   for (j = 1; j <= n; ++j) {
	   double d;
	   d=wa4[j];
	wa4[j] = fvec[j]*diag2[j-1];
	wa2[j] = x[j];
	x[j] = wa3[j];
	fvec[j]=d;
/* L270: */
    }
    ++(*nfev);
    fnorm1 = dnrm2_Fast1(n, &wa4[1]);
	if (*nfev==0) fnorm=fnorm1;
	{
		double wnorm;
		int more;
		char str[200];
		char*buf;
		wnorm=dnrm2_Fast1(n, &wa1[1]);
		more=(fnorm1<fnorm*1.1 || *nfev<20 || (n==1 && fnorm1<fnorm*1000)) && wnorm<1e100;
		if (more) {
			for(j=1;j<=n;++j) 
				wa1[j]*=-1.1;
		}
		if (*nfev%2) wnorm=-wnorm;
		buf=str+sprintfC(str,"{ %.16g, %.16g",wnorm, fnorm1);
		for(j=1;j<=n;++j) {
			buf=buf+sprintfC(buf,", %.16g",wa4[j]);
			if (j%5==3) {
				DymosimMessage(str);
				buf=str;
			}
		}
		buf=buf+sprintfC(buf,"%s",more?"},":"}};");
		DymosimMessage(str);
		if (more)
			goto L370b;
	}
	DymosimMessage("plotArray(Amat[:,1],Amat[:,2],-1);");
	DymosimMessage("");
L370real:
	/* end: restore scaling of residuals for error messages */
	{
		int i;
		for(i=1;i<=n;++i) {
			fvec[i]/=diag2[i-1];
		}
	}
/* termination. */

    *infrev = 0;
    goto L999;



/* ... save local variables in work array before return */
L999:
#ifndef fastnl
    dsave[1] = actred;
    dsave[2] = delta;
    dsave[3] = fnorm;
    dsave[4] = fnorm1;
    dsave[5] = pnorm;
    dsave[6] = prered;
    dsave[7] = ratio;
    dsave[8] = sum;
    dsave[9] = temp;
#endif
    isave[1] = *nfev;
    isave[2] = *njev;
#ifndef fastnl
    isave[3] = i__;
    isave[4] = iter;
#endif
    isave[5] = iwa[0];
#ifndef fastnl
    isave[6] = j;
    isave[7] = jm1;
    isave[8] = l;
    isave[9] = ncfail;
    isave[10] = ncsuc;
    isave[11] = nslow1;
    isave[12] = nslow2;
    isave[13] = irev;
    isave[14] = ifirst;
    if (jeval) {
	isave[15] = 1;
    } else {
	isave[15] = 0;
    }
    if (sing) {
	isave[16] = 1;
    } else {
	isave[16] = 0;
    }
#endif
#ifdef fastnl
#undef	actred
#undef delta 
#undef fnorm 
#undef fnorm1
#undef pnorm 
#undef prered
#undef ratio 
#undef sum 
#undef temp
#undef xnorm
	/* dsave 11..13 used by dymnl5 */
	/* dsave 14 used for norm of Jacobian */

#undef	i__ 
#undef iter
#undef j 
#undef jm1
#undef l 
#undef ncfail 
#undef ncsuc 
#undef nslow1 
#undef nslow2 
#undef irev 
#undef ifirst 
#undef jeval 
#undef sing 
#endif
    return;

/* last card of subroutine dymnl1. */

} /* dymnl1_ */


DYMOLA_STATIC void dymnl2_Fast(const integer m, const integer n, doublereal* a, const integer lda, const doublereal* v, const doublereal* w)
{
    /* Initialized data */

    static const doublereal one = 1.;

    /* System generated locals */
    integer a_dim1, a_offset;


    /* Local variables */
    integer i__, j;
    doublereal cos__, sin__;

/* ***begin prologue  dymnl2 */

/* subroutine dymnl2 */

/* given an m by n matrix a, this subroutine computes a*q where */
/* q is the product of 2*(n - 1) transformations */

/* gv(n-1)*...*gv(1)*gw(1)*...*gw(n-1) */

/* and gv(i), gw(i) are givens rotations in the (i,n) plane which */
/* eliminate elements in the i-th and n-th planes, respectively. */
/* q itself is not given, rather the information to recover the */
/* gv, gw rotations is supplied. */

/* the subroutine statement is */

/* subroutine dymnl2(m,n,a,lda,v,w) */

/* where */

/* m is a positive integer input variable set to the number */
/* of rows of a. */

/* n is a positive integer input variable set to the number */
/* of columns of a. */

/* a is an m by n array. on input a must contain the matrix */
/* to be postmultiplied by the orthogonal matrix q */
/* described above. on output a*q has replaced a. */

/* lda is a positive integer input variable not less than m */
/* which specifies the leading dimension of the array a. */

/* v is an input array of length n. v(i) must contain the */
/* information necessary to recover the givens rotation gv(i) */
/* described above. */

/* w is an input array of length n. w(i) must contain the */
/* information necessary to recover the givens rotation gw(i) */
/* described above. */

/* subroutines called */

/* fortran-supplied ... dabs,dsqrt */

/* argonne national laboratory. minpack project. march 1980. */
/* burton s. garbow, kenneth e. hillstrom, jorge j. more */
/* ***end prologue  dymnl2 */
/*     .. executable statements .. */
    /* Parameter adjustments */
    --w;
    --v;
    a_dim1 = lda;
    a_offset = a_dim1 + 1;
    a -= a_offset;

    /* Function Body */

/* apply the first set of givens rotations to a. */

/* ***first executable statement  dymnl2 */

    for (j = n-1; j >= 1; --j) {
		if (fabs(v[j])>one) {
			cos__ = one / v[j];
			/* Computing 2nd power */
			sin__ = sqrt(one - cos__ * cos__);
		} else {
			sin__ = v[j];
			/* Computing 2nd power */
			cos__ = sqrt(one - sin__ * sin__);
		}
		for (i__ = 1; i__ <= m; ++i__) {
			const double aij=a[i__ + j * a_dim1];
			const double ain=a[i__ + n * a_dim1];
			const double temp = cos__ * aij - sin__ * ain;
			a[i__ + n * a_dim1] = sin__ * aij + cos__ * ain;
			a[i__ + j * a_dim1] = temp;
		}
    }

/* apply the second set of givens rotations to a. */

    for (j = 1; j <= n-1; ++j) {
		if (fabs(w[j])>one) {
			cos__ = one / w[j];
			/* Computing 2nd power */
			sin__ = sqrt(one - cos__ * cos__);
		} else {
			sin__ = w[j];
			/* Computing 2nd power */
			cos__ = sqrt(one - sin__ * sin__);
		}
		for (i__ = 1; i__ <= m; ++i__) {
			const double aij=a[i__ + j * a_dim1];
			const double ain=a[i__ + n * a_dim1];
			const double temp = cos__ *aij  + sin__ * ain;
			a[i__ + n * a_dim1] = -sin__ * aij + cos__ * ain;
			a[i__ + j * a_dim1] = temp;
		}
    }
    return ;

/* last card of subroutine dymnl2. */

} /* dymnl2_ */


DYMOLA_STATIC void dymnl3_Fast(const integer m, const integer n, doublereal*s, const integer ls, 
				 const doublereal*u, doublereal*v, doublereal*w, logical*sing,
				 const integer*ipivots)
{
    /* Initialized data */

    static const doublereal one = 1.;
    static const doublereal p5 = .5;
    static const doublereal p25 = .25;
    static const doublereal zero = 0.;

    /* System generated locals */
    doublereal d__1, d__2;

    /* Local variables */
    integer i__, j, l, jj;
    doublereal tan__;
    doublereal cos__, sin__, tau, temp, cotan;

/* ***begin prologue  dymnl3 */

/* subroutine dymnl3 */

/* given an m by n lower trapezoidal matrix s, an m-vector u, */
/* and an n-vector v, the problem is to determine an */
/* orthogonal matrix q such that */

/* (s + u*v^t )*q */

/* is again lower trapezoidal. */

/* this subroutine determines q as the product of 2*(n - 1) */
/* transformations */

/* gv(n-1)*...*gv(1)*gw(1)*...*gw(n-1) */

/* where gv(i), gw(i) are givens rotations in the (i,n) plane */
/* which eliminate elements in the i-th and n-th planes, */
/* respectively. q itself is not accumulated, rather the */
/* information to recover the gv, gw rotations is returned. */

/* the subroutine statement is */

/* subroutine dymnl3(m,n,s,ls,u,v,w,sing) */

/* where */

/* m is a positive integer input variable set to the number */
/* of rows of s. */

/* n is a positive integer input variable set to the number */
/* of columns of s. n must not exceed m. */

/* s is an array of length ls. on input s must contain the lower */
/* trapezoidal matrix s stored by columns. on output s contains */
/* the lower trapezoidal matrix produced as described above. */

/* ls is a positive integer input variable not less than */
/* (n*(2*m-n+1))/2. */

/* u is an input array of length m which must contain the */
/* vector u. */

/* v is an array of length n. on input v must contain the vector */
/* v. on output v(i) contains the information necessary to */
/* recover the givens rotation gv(i) described above. */

/* w is an output array of length m. w(i) contains information */
/* necessary to recover the givens rotation gw(i) described */
/* above. */

/* sing is a logical output variable. sing is set true if any */
/* of the diagonal elements of the output s are zero. otherwise */
/* sing is set false. */

/* subprograms called */

/* fortran-supplied ... dabs,dsqrt */

/* argonne national laboratory. minpack project. march 1980. */
/* burton s. garbow, kenneth e. hillstrom, jorge j. more, */
/* john l. nazareth */
/* ***end prologue  dymnl3 */


/*     .. executable statements .. */
    /* Parameter adjustments */
    --w;
    --u;
    --v;
    --s;

    /* Function Body */

/* ***first executable statement  dymnl3 */

/* initialize the diagonal element pointer. */

    jj = n * (n+1)/2; /* *n * ((*m << 1) - *n + 1) / 2 - (*m - *n); */

/* move the nontrivial part of the last column of s into w. */

    l = jj;
    for (i__ = n; i__ <= m; ++i__) {
		w[i__] = s[l];
		++l;
    }

/* rotate the vector v into a multiple of the n-th unit vector */
/* in such a way that a spike is introduced into w. */

    for (j=n-1; j>=1; --j) {
		jj -= m - j + 1;
		w[j] = zero;
		if (v[j] != zero) {
			
			/*   determine a givens rotation which eliminates the */
			/*   j-th element of v. */
			
			if (((d__1 = v[n], Dymola_abs(d__1)) >= (d__2 = v[j], Dymola_abs(d__2)))) {
				tan__ = v[j] / v[n];
				/* Computing 2nd power */
				d__1 = tan__;
				cos__ = p5 / sqrt(p25 + p25 * (d__1 * d__1));
				sin__ = cos__ * tan__;
				tau = sin__;
			} else {
				cotan = v[n] / v[j];
				/* Computing 2nd power */
				d__1 = cotan;
				sin__ = p5 / sqrt(p25 + p25 * (d__1 * d__1));
				cos__ = sin__ * cotan;
				tau = one;
				if (Dymola_abs(cos__) * DBL_MAX > one) {
					tau = one / cos__;
				}
			}
			
			/*   apply the transformation to v and store the information */
			/*   necessary to recover the givens rotation. */
			
			v[n] = sin__ * v[j] + cos__ * v[n];
			v[j] = tau;
			
			/*   apply the transformation to s and extend the spike in w. */
			
			l = jj;
			for (i__ = j; i__ <= m; ++i__) {
				temp = cos__ * s[l] - sin__ * w[i__];
				w[i__] = sin__ * s[l] + cos__ * w[i__];
				s[l] = temp;
				++l;
			}
		}
    }

/* add the spike from the rank 1 update to w. */

    for (i__ = 1; i__ <= m; ++i__) {
		w[i__] += v[n] * (ipivots ? u[ipivots[i__-1]] : u[i__]);
    }
	
	/* eliminate the spike. */
	
    *sing = FALSE_;
    for (j = 1; j <= n-1; ++j) {
		if (w[j] != zero) {
			
			/*   determine a givens rotation which eliminates the */
			/*   j-th element of the spike. */
			
			if ((d__1 = s[jj], Dymola_abs(d__1)) >= (d__2 = w[j], Dymola_abs(d__2))) {
				tan__ = w[j] / s[jj];
				/* Computing 2nd power */
				d__1 = tan__;
				cos__ = p5 / sqrt(p25 + p25 * (d__1 * d__1));
				sin__ = cos__ * tan__;
				tau = sin__;
			} else {
				cotan = s[jj] / w[j];
				/* Computing 2nd power */
				d__1 = cotan;
				sin__ = p5 / sqrt(p25 + p25 * (d__1 * d__1));
				cos__ = sin__ * cotan;
				tau = one;
				if (Dymola_abs(cos__) * DBL_MAX > one) {
					tau = one / cos__;
				}
			}
			/*   apply the transformation to s and reduce the spike in w. */
			
			l = jj;
			for (i__ = j; i__ <= m; ++i__) {
				temp = cos__ * s[l] + sin__ * w[i__];
				w[i__] = -sin__ * s[l] + cos__ * w[i__];
				s[l] = temp;
				++l;
			}
			
			/*   store the information necessary to recover the */
			/*   givens rotation. */
			
			w[j] = tau;
			
		}
		/*   test for zero diagonal elements in the output s. */
		
		if (s[jj] == zero) {
			*sing = TRUE_;
		}
		jj += m - j + 1;
    }

/* move w back into the last column of the output s. */

    l = jj;
    for (i__ = n; i__ <= m; ++i__) {
		s[l] = w[i__];
		++l;
    }
    if (s[jj] == zero) {
		*sing = TRUE_;
    }
    return ;

/* last card of subroutine dymnl3. */

} /* dymnl3_ */


DYMOLA_STATIC void dymnl4_Fast(const integer n, const doublereal* r__, const integer lr, const doublereal*diag, 
				 const doublereal*qtb, const doublereal delta, doublereal*x, doublereal*wa1, doublereal*wa2,
				 const integer*ipivots,
				 int inJacobian)
{
    /* Initialized data */

    static const doublereal one = 1.;
    static const doublereal zero = 0.;

    /* System generated locals */
    doublereal d__1, d__2, d__3, d__4;


    /* Local variables */
    integer i__, j, k, l, jj, jp1;
    doublereal sum, temp;
    doublereal alpha, bnorm, gnorm, qnorm, sgnorm;
	int divideBy=0;

/* ***begin prologue  dymnl4 */

/* ********* */

/* subroutine dymnl4 */

/* given an m by n matrix a, an n by n nonsingular diagonal */
/* matrix d, an m-vector b, and a positive number delta, the */
/* problem is to determine the convex combination x of the */
/* gauss-newton and scaled gradient directions that minimizes */
/* (a*x - b) in the least squares sense, subject to the */
/* restriction that the euclidean norm of d*x be at most delta. */

/* this subroutine completes the solution of the problem */
/* if it is provided with the necessary information from the */
/* qr factorization of a. that is, if a = q*r, where q has */
/* orthogonal columns and r is an upper triangular matrix, */
/* then dymnl4 expects the full upper triangle of r and */
/* the first n components of (q transpose)*b. */

/* the subroutine statement is */

/* subroutine dymnl4(n,r,lr,diag,qtb,delta,x,wa1,wa2) */

/* where */

/* n is a positive integer input variable set to the order of r. */

/* r is an input array of length lr which must contain the upper */
/* triangular matrix r stored by rows. */

/* lr is a positive integer input variable not less than */
/* (n*(n+1))/2. */

/* diag is an input array of length n which must contain the */
/* diagonal elements of the matrix d. */

/* qtb is an input array of length n which must contain the first */
/* n elements of the vector (q transpose)*b. */

/* delta is a positive input variable which specifies an upper */
/* bound on the euclidean norm of d*x. */

/* x is an output array of length n which contains the desired */
/* convex combination of the gauss-newton direction and the */
/* scaled gradient direction. */

/* wa1 and wa2 are work arrays of length n. */

/* ipivots is an optional pivoting matrix as returned by dymnl7_fast */
/* Additional comments on pivoting*/
/* A=QRP' */
/* min |Ax-b|, |d'*x|<=delta */
/* min |QRP'x-b|, |d'*x|<=delta */
/* y=P'*x x=P*y */
/* min |QRy-b|, |d'*P*y|<=delta */
/* min |Ry-Q'b|, |(P'*d)'*y|<=delta */
/* We compute y in x-vector and then go to returnSection */

/* subprograms called */

/* fortran-supplied ... dabs,dmax1,dmin1,dsqrt */

/* argonne national laboratory. minpack project. march 1980. */
/* burton s. garbow, kenneth e. hillstrom, jorge j. more */
/* ***end prologue  dymnl4 */


/*     .. executable statements .. */
    /* Parameter adjustments */
    --wa2;
    --wa1;
    --x;
    --qtb;
    --diag;
    --r__;

    /* Function Body */

/* ***first executable statement  dymnl4 */

/* first, calculate the gauss-newton direction. */

    jj = n * (n + 1) / 2 + 1;
    for (k = 1; k <= n; ++k) {
		j = n - k + 1;
		jp1 = j + 1;
		jj -= k;
		l = jj + 1;
		sum = zero;

		for (i__ = jp1; i__ <= n; ++i__) {
			sum += r__[l] * x[i__];
			++l;
		}

		temp = r__[jj];
		if (temp == zero) {
			
			l = j;
			for (i__ = 1; i__ <= j; ++i__) {
				/* Computing MAX */
				d__2 = temp, d__3 = (d__1 = r__[l], Dymola_abs(d__1));
				temp = Dymola_max(d__2,d__3);
				l = l + n - i__;
				/* L30: */
			}
			temp = DBL_EPSILON * temp;
			if (temp == zero) {
				temp = DBL_EPSILON;
			}
		}
		if (divideBy>0) {
			/* Ensure consistent down-scaling */
			double t2;
			int i;
			t2=qtb[j];
			for(i=0;i<divideBy;++i) t2*=1e-50;
			x[j]=(t2-sum)/temp;
		} else
			x[j] = (qtb[j] - sum) / temp;
		if (x[j]>1e200 || x[j]<-1e200) {
			++divideBy;
			for(i__=j;i__<=n;++i__) x[j]*=1e-50; /* Scale it down to avoid overflow*/
		}
    }

/* test whether the gauss-newton direction is acceptable. */

    for (j = 1; j <= n; ++j) {
		wa1[j] = zero;
		wa2[j] = (ipivots ? diag[ipivots[j-1]] : diag[j]) * x[j];
		/* L60: */
    }
    qnorm = dnrm2_Fast1(n, &wa2[1]);
    if (qnorm <= delta || inJacobian) {
		goto returnSection;
    }
/* note: early return if acceptable */

/* the gauss-newton direction is not acceptable. */
/* next, calculate the scaled gradient direction. */

    l = 1;
    for (j = 1; j <= n; ++j) {
		temp = qtb[j];
		for (i__ = j; i__ <= n; ++i__) {
			wa1[i__] += r__[l] * temp;
			++l;
		}
		wa1[j] /= (ipivots ? diag[ipivots[j-1]] : diag[j]);
    }

/* calculate the norm of the scaled gradient and test for */
/* the special case in which the scaled gradient is zero. */

    gnorm = dnrm2_Fast1(n, &wa1[1]);
    sgnorm = zero;
    alpha = delta / qnorm;
    if (gnorm == zero) {
	goto L120;
    }

/* calculate the point along the scaled gradient */
/* at which the quadratic is minimized. */

    for (j = 1; j <= n; ++j) {
	wa1[j] = wa1[j] / gnorm / (ipivots ? diag[ipivots[j-1]] : diag[j]);
/* L90: */
    }
    l = 1;
    for (j = 1; j <= n; ++j) {
	sum = zero;
	for (i__ = j; i__ <= n; ++i__) {
	    sum += r__[l] * wa1[i__];
	    ++l;
/* L100: */
	}
	wa2[j] = sum;
/* L110: */
    }
    temp = dnrm2_Fast1(n, &wa2[1]);
    sgnorm = gnorm / temp / temp;

/* test whether the scaled gradient direction is acceptable. */

    alpha = zero;
    if (sgnorm >= delta) {
	goto L120;
    }

/* the scaled gradient direction is not acceptable. */
/* finally, calculate the point along the dogleg */
/* at which the quadratic is minimized. */

    bnorm = dnrm2_Fast1(n, &qtb[1]);
    temp = bnorm / gnorm * (bnorm / qnorm) * (sgnorm / delta);
/* Computing 2nd power */
    d__1 = sgnorm / delta;
/* Computing 2nd power */
    d__2 = temp - delta / qnorm;
/* Computing 2nd power */
    d__3 = delta / qnorm;
/* Computing 2nd power */
    d__4 = sgnorm / delta;
    temp = temp - delta / qnorm * (d__1 * d__1) + sqrt(d__2 * d__2 + (one - 
	    d__3 * d__3) * (one - d__4 * d__4));
/* Computing 2nd power */
    d__1 = sgnorm / delta;
    alpha = delta / qnorm * (one - d__1 * d__1) / temp;
L120:

/* form appropriate convex combination of the gauss-newton */
/* direction and the scaled gradient direction. */

    temp = (one - alpha) * Dymola_min(sgnorm,delta);
    for (j = 1; j <= n; ++j) {
		x[j] = temp * wa1[j] + alpha * x[j];
    }
returnSection:
	/* Pivot result if necessary */

	if (ipivots) {
		/* First copy result to work-area */
		for(j=1;j<=n;++j) wa1[j]=x[j];
		/* Then compute P*x as indicated in dymnl7_fast */
		for(j=1;j<=n;++j) x[ipivots[j-1]]=wa1[j];
	}
    return ;

/* last card of subroutine dymnl4. */

} /* dymnl4_ */


DYMOLA_STATIC void dymnl5_Fast(integer*irev, const integer n, doublereal *x, doublereal *fvec, doublereal *fjac, 
			 const integer ldfjac, const integer ml, const integer mu, const doublereal *nominalx, 
	doublereal *wa1, doublereal *wa2, doublereal* dsave, integer* isave,logical useCentral,
	logical printDetails,const char*const*varnames,int check)
{
    /* Initialized data */

    static const doublereal zero = 0.;

    /* System generated locals */
    integer fjac_dim1, fjac_offset;
    doublereal d__1;

    /* Local variables */
    doublereal h__;
    integer i__, j, k;
    doublereal eps, temp;
    integer msum;


/* ----------------------------------------------------------------------- */
/*     this subroutine computes a forward-difference approximation */
/*     to the n by n jacobian matrix associated with a specified */
/*     problem of n functions in n variables. if the jacobian has */
/*     a banded form, then function evaluations are saved by only */
/*     approximating the nonzero terms. */

/*       irev is a integer input and output variable for driving reverse */
/*         communication */

/*       n is a positive integer input variable set to the number */
/*         of functions and variables. */

/*       x is an input array of length n. */

/*       fvec is an input array of length n which must contain the */
/*         functions evaluated at x. */

/*       fjac is an output n by n array which contains the */
/*         approximation to the jacobian matrix evaluated at x. */

/*       ldfjac is a positive integer input variable not less than n */
/*         which specifies the leading dimension of the array fjac. */

/*       ml is a nonnegative integer input variable which specifies */
/*         the number of subdiagonals within the band of the */
/*         jacobian matrix. if the jacobian is not banded, set */
/*         ml to at least n - 1. */

/*       nominal is an input variable used in determining a suitable */
/*         step length for the forward-difference approximation. */

/*       mu is a nonnegative integer input variable which specifies */
/*         the number of superdiagonals within the band of the */
/*         jacobian matrix. if the jacobian is not banded, set */
/*         mu to at least n - 1. */

/*       wa1,wa2 and wa3 are work arrays of length n. if ml + mu + 1 is a */
/*         least n, then the jacobian is considered dense, and wa2 is */
/*         not referenced. */

/*       dsave,isave are work arrays to save the internal state of the */
/*         subroutine for the next call. */

/* subprograms called */
/* ------------------ */
/* fortran-supplied  dmax1,dabs,dsqrt */

/* author            Guenter M. Gramlich, D. Joos */
/* date              23.05.91, Oberpfaffenhofen */
/* version           1.0       23.09.93 */
/* ----------------------------------------------------------------------- */


/*     .. executable statements .. */
    /* Parameter adjustments */
    --wa2;
    --wa1;
    --fvec;
    --x;
    fjac_dim1 = ldfjac;
    fjac_offset = fjac_dim1 + 1;
    fjac -= fjac_offset;
    --dsave;
    --isave;

    /* Function Body */

/* ***first executable statement  dymnl5 */

    msum = ml + mu + 1;

    if (*irev == 1) {
		/*       ... get last internal state */
		temp = dsave[1];
		h__ = dsave[2];
		eps = dsave[3];
		j = isave[1];
		k = isave[2];
		/*       ... jump to appropriate place */
		if (msum < n) {
			goto L110;
		}
		goto L30;
    } else {
		/*       ... initialize internal state */
		temp = 0.;
		h__ = 0.;
		eps = sqrt(DBL_EPSILON);
		if (useCentral) eps = 1e-6;
		j = 0;
		k = 0;
    }
	
    if (msum < n) {
		goto L70;
    }

/* computation of dense approximate jacobian. */

    j = 1;
    for (i__ = 1; i__ <= n; ++i__) {
	wa1[i__] = fvec[i__];
/* L10: */
    }
L20:
	temp = x[j];
	if (nominalx!=0)
		h__ = eps * Dymola_max(nominalx[j-1],Dymola_abs(temp));
	else
		h__ = eps * Dymola_max(1, Dymola_abs(temp));
	if (h__ == zero) {
		h__ = eps;
	}
	if (k==1)
		x[j] = temp - h__;
	else
		x[j] = temp + h__;
    *irev = 1;
    goto L999;
L30:
    x[j] = temp;
	if (k==0) {
		if (printDetails && check && !useCentral) {
			double rSum;
			rSum=0;
			for (i__ = 1; i__ <= n; ++i__) rSum+=Dymola_abs(fjac[i__ + j * fjac_dim1]);
			for (i__ = 1; i__ <= n; ++i__) {
				double analyticD, numericD;
				analyticD=fjac[i__ + j * fjac_dim1];
				numericD=(fvec[i__] - wa1[i__]) / h__;
			if (Dymola_abs(analyticD-numericD)>1e-3*rSum) {
				char str[600];
				/* First is j = which variable, second is i which equation */
				sprintfC(str,"Jacobian element(%d %.400s,%d)= is %g [analytic] or %g [sum=%g]",j,
					(varnames && varnames[j-1])?varnames[j-1]:"",i__,analyticD,numericD,rSum);
				DymosimMessage(str);
			}
			}
		}
		for (i__ = 1; i__ <= n; ++i__) {
			fjac[i__ + j * fjac_dim1] = (fvec[i__] - wa1[i__]) / h__;
		}
		if (useCentral) {
			k=1;
			goto L20;
		}
	} else {
		double rSum;
		rSum=0;
		for (i__ = 1; i__ <= n; ++i__) rSum+=Dymola_abs(fjac[i__ + j * fjac_dim1]);
		for (i__ = 1; i__ <= n; ++i__) {
			double d;
			double fwd_d;
			d= - (fvec[i__] - wa1[i__]) / h__;
			fwd_d = fjac[i__ + j * fjac_dim1];
			if (printDetails && Dymola_abs(d-fwd_d)>1e-3*(rSum)) {
				char str[600];
				/* First is j = which variable, second is i which equation */
				sprintfC(str,"Jacobian element(%d %.400s,%d)= is %g or %g [sum=%g]",j,
					(varnames && varnames[j-1])?varnames[j-1]:"",i__,fwd_d,d,rSum);
				DymosimMessage(str);
			}	
			fjac[i__ + j * fjac_dim1] =fwd_d*0.5+0.5*d; /* Use central difference */
		}
		k=0;
	}
	if (j == n) {
		goto L50;
	}
    ++j;
    goto L20;
L50:
    for (i__ = 1; i__ <= n; ++i__) {
	fvec[i__] = wa1[i__];
/* L60: */
    }
    goto L160;

L70:

/* computation of banded approximate jacobian. */

    k = 1;
    for (j = 1; j <= n; ++j) {
	wa1[j] = fvec[j];
/* L80: */
    }

L90:
    for (j = k; msum < 0 ? j >= n : j <= n; j += msum) {
	wa2[j] = x[j];
	h__ = eps * (d__1 = wa2[j], Dymola_abs(d__1));
	if (h__ == zero) {
	    h__ = eps;
	}
	x[j] = wa2[j] + h__;
/* L100: */
    }
/* ...return */
    *irev = 1;
    goto L999;

L110:
    for (j = k; msum < 0 ? j >= n : j <= n; j += msum) {
	x[j] = wa2[j];
	h__ = eps * (d__1 = wa2[j], Dymola_abs(d__1));
	if (h__ == zero) {
	    h__ = eps;
	}
	for (i__ = 1; i__ <= n; ++i__) {
	    fjac[i__ + j * fjac_dim1] = zero;
	    if (i__ >= j - mu && i__ <= j + ml) {
		fjac[i__ + j * fjac_dim1] = (fvec[i__] - wa1[i__]) / h__;
	    }
/* L120: */
	}
/* L130: */
    }
    if (k == msum) {
	goto L140;
    }
    ++k;
    goto L90;

L140:
    for (i__ = 1; i__ <= n; ++i__) {
	fvec[i__] = wa1[i__];
/* L150: */
    }
L160:
    *irev = 0;
    goto L999;

/* ... save internal state */
L999:
    dsave[1] = temp;
    dsave[2] = h__;
    dsave[3] = eps;
    isave[1] = j;
    isave[2] = k;
    return;

/* last card of subroutine dymnl5. */

} /* dymnl5_ */


DYMOLA_STATIC void dymnl6_Fast(const integer m, const integer n, doublereal* q, const integer ldq, doublereal*wa)
{
    /* Initialized data */

    static const doublereal one = 1.;
    static const doublereal zero = 0.;

    /* System generated locals */
    integer q_dim1, q_offset;

    /* Local variables */
    integer i__, j, k, l;
    doublereal sum, temp;
    integer minmn;

/* ***begin prologue  dymnl6 */

/* subroutine dymnl6 */

/* this subroutine proceeds from the computed qr factorization of */
/* an m by n matrix a to accumulate the m by m orthogonal matrix */
/* q from its factored form. */

/* the subroutine statement is */

/* subroutine dymnl6(m,n,q,ldq,wa) */

/* where */

/* m is a positive integer input variable set to the number */
/* of rows of a and the order of q. */

/* n is a positive integer input variable set to the number */
/* of columns of a. */

/* q is an m by m array. on input the full lower trapezoid in */
/* the first min(m,n) columns of q contains the factored form. */
/* on output q has been accumulated into a square matrix. */

/* ldq is a positive integer input variable not less than m */
/* which specifies the leading dimension of the array q. */

/* wa is a work array of length m. */

/* subprograms called */

/* fortran-supplied ... min0 */

/* argonne national laboratory. minpack project. march 1980. */
/* burton s. garbow, kenneth e. hillstrom, jorge j. more */
/* ***end prologue  dymnl6 */
/*     .. executable statements .. */
    /* Parameter adjustments */
    --wa;
    q_dim1 = ldq;
    q_offset = q_dim1 + 1;
    q -= q_offset;

    /* Function Body */

/* zero out upper triangle of q in the first min(m,n) columns. */

/* ***first executable statement  dymnl6 */

    minmn = Dymola_min(m,n);
    for (j = 2; j <= minmn; ++j) {
		for (i__ = 1; i__ <= j-1; ++i__) {
			q[i__ + j * q_dim1] = zero;
		}
    }

/* initialize remaining columns to those of the identity matrix. */

    for (j = n+1; j <= m; ++j) {
		for (i__ = 1; i__ <= m; ++i__) {
			q[i__ + j * q_dim1] = zero;
		}
		q[j + j * q_dim1] = one;
    }

/* accumulate q from its factored form. */

    for (l = 1; l <= minmn; ++l) {
		k = minmn - l + 1;
		for (i__ = k; i__ <= m; ++i__) {
			wa[i__] = q[i__ + k * q_dim1];
			q[i__ + k * q_dim1] = zero;
		}
		q[k + k * q_dim1] = one;
		if (wa[k] != zero) {
			double wak_inverse;
			wak_inverse=1/wa[k];
			for (j = k; j <= m; ++j) {
				sum = zero;
				for (i__ = k; i__ <= m; ++i__) {
					sum += q[i__ + j * q_dim1] * wa[i__];
				}
				temp = sum *wak_inverse;
				for (i__ = k; i__ <= m; ++i__) {
					q[i__ + j * q_dim1] -= temp * wa[i__];
				}
			}
		}
    }
    return;

/* last card of subroutine dymnl6. */

} /* dymnl6_ */


DYMOLA_STATIC void dymnl7_Fast(const integer m, const integer n, doublereal* a, const integer lda, 
				 const logical pivot, 
			 integer*ipvt, const integer lipvt, doublereal*sigma, doublereal*acnorm, 
	doublereal*wa,const char*const*varnames, doublereal*x,int printEvent) 
{
    /* Initialized data */

    static const doublereal one = 1.;
    static const doublereal p05 = .05;
    static const doublereal zero = 0.;
	static const doublereal smallPivot =1e-9;

    /* System generated locals */
    integer a_dim1, a_offset;
    doublereal d__1, d__2, d__3;


    /* Local variables */
    integer i__, j, k, jp1;
    doublereal sum;
    integer kmax;
    doublereal temp;
    integer minmn;
    doublereal ajnorm;
	integer wasSmall=0;
	integer firstTime=1;


/* this subroutine uses householder transformations with column */
/* pivoting (optional) to compute a qr factorization of the */
/* m by n matrix a. that is, dymnl7 determines an orthogonal */
/* matrix q, a permutation matrix p, and an upper trapezoidal */
/* matrix r with diagonal elements of nonincreasing magnitude, */
/* such that a*p = q*r. the householder transformation for */
/* column k, k = 1,2,...,min(m,n), is of the form */

/* t */
/* i - (1/u(k))*u*u */

/* where u has zeros in the first k-1 positions. the form of */
/* this transformation and the method of pivoting first */
/* appeared in the corresponding linpack subroutine. */

/* the subroutine statement is */

/* subroutine dymnl7(m,n,a,lda,pivot,ipvt,lipvt,sigma,acnorm,wa) */

/* where */

/* m is a positive integer input variable set to the number */
/* of rows of a. */

/* n is a positive integer input variable set to the number */
/* of columns of a. */

/* a is an m by n array. on input a contains the matrix for */
/* which the qr factorization is to be computed. on output */
/* the strict upper trapezoidal part of a contains the strict */
/* upper trapezoidal part of r, and the lower trapezoidal */
/* part of a contains a factored form of q (the non-trivial */
/* elements of the u vectors described above). */

/* lda is a positive integer input variable not less than m */
/* which specifies the leading dimension of the array a. */

/* pivot is a logical input variable. if pivot is set .true., */
/* then column pivoting is enforced. if pivot is set .false., */
/* then no column pivoting is done. */

/* ipvt is an integer output array of length lipvt. ipvt */
/* defines the permutation matrix p such that a*p = q*r. */
/* column j of p is column ipvt(j) of the identity matrix. */
/* if pivot is .false., ipvt is not referenced. */

/* For a permutation matrix we have: transpose(P)*P=I [GvL] */
/* For this representation it is easy to compute: */

/* transpose(P)*x_i =x[ipvt[i]] or x[ipvt[i-1]] if C-array for ipvt. */
/* x'*P*y=(P'*x)'*y where P'*x is easy to compute */

/* To compute y=P*x we initialize the output and then: */
/* y[ipvt[i]]=x[i] */

/* Note: We use a C-array in other routines to make it easier to handle optional pivoting. */

/* lipvt is a positive integer input variable. if pivot is */
/* .false., then lipvt may be as small as 1. if pivot is */
/* .true., then lipvt must be at least n. */

/* sigma is an output array of length n which contains the */
/* diagonal elements of r. */

/* acnorm is an output array of length n which contains the */
/* norms of the corresponding columns of the input matrix a. */
/* if this information is not needed, then acnorm can coincide */
/* with sigma. */

/* wa is a work array of length n. if pivot is .false., then wa */
/* can coincide with sigma. */

/* subprograms called */

/* fortran-supplied ... dmax1,dsqrt,min0 */

/* minpack. version of december 1978. */
/* burton s. garbow, kenneth e. hillstrom, jorge j. more */
/* ***end prologue  dymnl7 */


/*     .. executable statements .. */
    /* Parameter adjustments */
    --wa;
    --acnorm;
    --sigma;
    a_dim1 = lda;
    a_offset = a_dim1 + 1;
    a -= a_offset;
    --ipvt;

    /* Function Body */

/* ***first executable statement  dymnl7 */

/* compute the initial column norms and initialize several arrays. */

    for (j = 1; j <= n; ++j) {
	acnorm[j] = dnrm2_Fast1(m, &a[j * a_dim1 + 1]);
	sigma[j] = acnorm[j];
	wa[j] = sigma[j];
	if (pivot) {
	    ipvt[j] = j;
	}
	wasSmall=wasSmall || (acnorm[j]<smallPivot);
/* L10: */
    }


/* reduce a to r with householder transformations. */

    minmn = Dymola_min(m,n);
    for (j = 1; j <= minmn; ++j) {
		if (pivot) {
			
			/*   bring the column of largest norm into the pivot position. */
			
			kmax = j;
			for (k = j; k <= n; ++k) {
				/* Commented code: Only pivot if necessary */
				if (sigma[k] > sigma[kmax] /*&& (kmax!=j || sigma[kmax]==0)*/) {
					kmax = k;
				}
				
			}
			if (sigma[kmax]<smallPivot) {
				/* Too small */
				
				if (!!(printEvent&(1<<8))) {
					if (firstTime) {
					  /* This message will appear as heading to the other ones */
					  DymosimMessage("The non-linear system is difficult to solve, this could either be\n"
						"due to a dependency between equations, or due to bad start-values.");
					  firstTime=0;
					}

					if (!wasSmall) {
					 /* Simple analysis of cause. If wasSmall is true at least one column was small before QR */
					 /* Since wasSmall is false none of the columns were small before start of QR and thus */
					 /* this special message, indicating that one has to look for dependency between variables. */
					  DymosimMessage("Note: None of the columns are small on their own - the problem only occurred after eliminating other dependencies between variables.");
					 wasSmall=1;
					}
					{
					  char str[2000];
					  sprintfC(str, "Variable %d(%.400s)=%g has small pivot %g and is hard to solve%s", ipvt[kmax],
						  varnames?varnames[ipvt[kmax]-1]:"?", x[ipvt[kmax]-1], sigma[kmax], (j<minmn)?",":".");
					  DymosimMessage(str);
					 }
				}
			}
			if (kmax != j) {
				/* Swap kmax and j */
				
				for (i__ = 1; i__ <= m; ++i__) {
					temp = a[i__ + j * a_dim1];
					a[i__ + j * a_dim1] = a[i__ + kmax * a_dim1];
					a[i__ + kmax * a_dim1] = temp;
				}
				sigma[kmax] = sigma[j];
				wa[kmax] = wa[j];
				k = ipvt[j];
				ipvt[j] = ipvt[kmax];
				ipvt[kmax] = k;
			}
		}
		
		/*   compute the householder transformation to reduce the */
		/*   j-th column of a to a multiple of the j-th unit vector. */
		
		ajnorm = dnrm2_Fast1(m-j+1, &a[j + j * a_dim1]);
		if (ajnorm != zero) {
			if (a[j + j * a_dim1] < zero) {
				ajnorm = -ajnorm;
			}
			for (i__ = j; i__ <= m; ++i__) {
				a[i__ + j * a_dim1] /= ajnorm;
				/* L50: */
			}
			a[j + j * a_dim1] += one;
			
			/*   apply the transformation to the remaining columns */
			/*   and update the norms. */
			
			jp1 = j + 1;
			for (k = jp1; k <= n; ++k) {
				sum = zero;
				for (i__ = j; i__ <= m; ++i__) {
					sum += a[i__ + j * a_dim1] * a[i__ + k * a_dim1];
				}
				temp = sum / a[j + j * a_dim1];
				for (i__ = j; i__ <= m; ++i__) {
					a[i__ + k * a_dim1] -= temp * a[i__ + j * a_dim1];
				}
				if (!(!pivot || sigma[k] == zero)) {
					
					temp = a[j + k * a_dim1] / sigma[k];
					/* Computing MAX */
					/* Computing 2nd power */
					d__3 = temp;
					d__1 = zero, d__2 = one - d__3 * d__3;
					sigma[k] *= sqrt((Dymola_max(d__1,d__2)));
					/* Computing 2nd power */
					d__1 = sigma[k] / wa[k];
					if (p05 * (d__1 * d__1) <= DBL_EPSILON) {
						
						sigma[k] = dnrm2_Fast1(m - j, &a[jp1 + k * a_dim1]);
						wa[k] = sigma[k];
					}
				}
				;
			}
		}
		sigma[j] = -ajnorm;
		/* L110: */
    }
    return ;

/* last card of subroutine dymnl7. */

} /* dymnl7_ */

DYMOLA_STATIC void dymnlinf1_(const integer*sysnr, const integer*n, const doublereal*sol, const doublereal*res, 
				const integer*nfunc, const integer*njac, const integer*nmax,const char*const*varnames)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Local variables */
    integer i__;
    doublereal norm;


/* Print messages after successful termination of DymNL */

/* SysNr: IN, number of nonlinear equation system */
/* N    : IN, argument N of DymNL */
/* Sol  : IN, argument X of DymNL (solution vector) */
/* Res  : IN, argument FVEC of DymNL (residue) */
/* nFunc: IN, accumulated number of calculations of residue */
/* nJac : IN, accumulated number of calculations of symbolic */
/*        Jacobian */
/* nmax : IN, dimension of nFunc/nJac */

/* ..................................................................... */



/* ..................................................................... */

    /* Parameter adjustments */
    --res;
    --sol;
    --njac;
    --nfunc;

    /* Function Body */
    norm = 0.;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	norm += (d__1 = res[i__], Dymola_abs(d__1));
/* L10: */
    }
    dymosimmessagedouble_("  Solution: Infinity-" "norm of residue = ", &norm, 
	    39);
    if (*sysnr <= *nmax) {
	if (nfunc[*sysnr] > 0) {
	    if (njac[*sysnr] <= 0) {
		DymosimMessageInt_("            Acc. number of residue calculations:", &nfunc[*sysnr], 48);
	    } else {
		DymosimMessageInt_("            Acc. number of residue     calc.:", &nfunc[*sysnr], 51);
	    }
	}

	if (njac[*sysnr] > 0) {
	    DymosimMessageInt_("            Acc. number of symbolic Jacobian calc.:", &njac[*sysnr], 51);
	}
    }

    dymnlinf2_(*n, &sol[1],varnames);

	dymosimmessage_("Residual:", 9);

	dymnlinf2_(*n, &res[1], 0);


    return ;
} /* dymnlinf1_ */




DYMOLA_STATIC void dymnlinf2_(integer n, const doublereal*sol,const char*const*varnames)
{
    /* Local variables */
    integer i__;


/* Print solution vector of DymNL */

/* N    : IN, argument N of DymNL */
/* Sol  : IN, argument X of DymNL (solution vector) */

/* ..................................................................... */



/* ..................................................................... */

    /* Parameter adjustments */
    --sol;

	if (varnames==0) {
		dymosimmessagedoubles_("", &sol[1], &n, 0);
	} else {
		for(i__=1;i__<=n;i__++) {
			char str[200];
			int i;
			for(i=0;varnames[i__-1][i]!='\0';i++) {
				str[i]=varnames[i__-1][i];
				if (i>=140) {
					/* truncate string */
					strcpy(str+i,"...");i+=3;
					break;
				}
			}
			strcpy(str+i," = ");i+=3;
			dymosimmessagedouble_(str,sol+i__,i);
		}
	}
    return ;
} /* dymnlinf2_ */



DYMOLA_STATIC void dymnlerr_(const doublereal*time, const integer*sysnr, const integer*info, const integer*n, 
			   const integer*iopt, const doublereal*tol, const doublereal*sol, const doublereal*res, 
	const integer*nfunc, const integer*njac, const integer*nmax,const char*const*varnames)

{
    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Local variables */
    integer i__, nfcn;
    doublereal norm;


/* Print error message, if nonlinear system of equations could */
/* not be solved */

/* Time : IN, actual time instant */
/* SysNr: IN, number of system of equations */
/* Info : IN, argument INFO of DymNL */
/* N    : IN, argument N of DymNL */
/* iOpt : IN, argument iOpt of DymNL */
/* Tol  : IN, argument TOL of DymNL */
/* Sol  : IN, argument X of DymNL (solution vector) */
/* Res  : IN, argument FVEC of DymNL (residue) */
/* nFunc: IN, accumulated number of calculations of residue */
/* nJac : IN, accumulated number of calculations of symbolic */
/*        Jacobian */
/* nmax : IN, dimension of nFunc/nJac */

/* ..................................................................... */



/* ..................................................................... */


    /* Parameter adjustments */
    --res;
    --sol;
    --njac;
    --nfunc;

    /* Function Body */
	dymosimmessageSev_(1, "ERROR: Failed to solve non-linear system using Newton solver.", 65);
	dymosimmessageSev_(1, "To get more information: Turn on Simulation/Setup/Debug/Nonlinear solver diagnostics/Details",
		100);
    dymosimmessagedoubleSev_(1, "Solution to systems of equations not found at time =", time, 52);

    DymosimMessageIntSev_(1, "   Nonlinear system of equations number = ", sysnr, 
	    42);

    norm = 0.;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	norm += (d__1 = res[i__], Dymola_abs(d__1));
/* L10: */
    }
	if (norm>=DBL_MAX)
		  dymosimmessageSev_(1, "   Failed to evaluate any residual. Try to give better start-values.",68);
	else
		  dymosimmessagedoubleSev_(1, "   Infinity-" "norm of residue = ", &norm, 30);


    if (*info == 0) {
	dymosimmessageSev_(1, "   Improper input parameters (wrong initialization)",
		 51);

    } else if (*info == 2) {
	if (*iopt == 2) {
	    nfcn = (*n + 1) * 200;
	} else {
	    nfcn = (*n + 1) * 100;
	}
	dymosimmessageSev_(1, "   Number of calls to nonlinear solver DymNL has reached or exceeded", 68);
	DymosimMessageIntSev_(1, "   the maximum allowed number of function calls= ", &nfcn, 50);

    } else if (*info == 3) {
	dymosimmessagedoubleSev_(1, "   The tolerance of ", tol, 20);
	dymosimmessageSev_(1, "   is too small. No further improvement", 39);
	dymosimmessageSev_(1, "   in the approximate solution is possible.", 43);

    } else if (*info == 4) {
	dymosimmessageSev_(1, "   Iteration is not making good progress.", 41);

    } else {
	DymosimMessageIntSev_(1, "   Argument INFO = ", info, 19);
	dymosimmessageSev_(1, "   of solver DymNL is wrong.", 28);
    }




    if (*sysnr <= *nmax) {
	if (nfunc[*sysnr] > 0) {
	    if (njac[*sysnr] <= 0) {
		DymosimMessageIntSev_(1, "   Accumulated number of residue calculations:", &nfunc[*sysnr], 46);
	    } else {
		DymosimMessageIntSev_(1, "   Accumulated number of residue       calc.:", &nfunc[*sysnr], 49);
	    }
	}

	if (njac[*sysnr] > 0) {
	    DymosimMessageIntSev_(1, "   Accumulated number of symbolic Jacobian calc.:", &njac[*sysnr], 49);
	}
    }


    dymosimmessageSev_(1, "   Last values of solution vector:", 34);
    dymnlinf2_(*n, &sol[1], varnames);

	if (norm<DBL_MAX) {
		dymosimmessageSev_(1, "   Last values of residual vector:", 34);
		dymnlinf2_(*n, &res[1], 0);
		dymosimmessageSev_(1, " ", 1);
	}
    return;
} /* dymnlerr_ */

/* Subroutine */ DYMOLA_STATIC int dymnon_(qinfrev, qiopt, qnnl, qsol, qres, qjac, qtol, 
	qinfo, qd, qi, printevent, qnlnr, time, qnlfunc, qnljac, qnlmax, 
	ierr___)
integer *qinfrev, *qiopt, *qnnl;
doublereal *qsol, *qres, *qjac, *qtol;
integer *qinfo;
doublereal *qd;
integer *qi, *printevent, *qnlnr;
doublereal *time;
integer *qnlfunc, *qnljac, *qnlmax, *ierr___;
{
    static integer idum[20];







    /* Parameter adjustments */
    --qsol;
    --qres;
    qjac -= 11;
    --qd;
    --qnlfunc;
    --qnljac;

    /* Function Body */
    if (*qinfrev == -1) {
	if (*printevent&(1<<6)) {
	    dymosimmessage_(" ", 1);
	    DymosimMessageInt_("Solving nonlinear system of equations number = ", qnlnr, 47);
	    dymosimmessagedouble_("        iteratively at Time = ", time, 30)
		    ;
	}
    }
    if ((*printevent&(1<<7)) && *qinfrev != 3) {
	dymnlinf2_(*qnnl, &qsol[1], 0);
    }
    if (*qinfrev > 0) {
		++qnlfunc[*qnlnr];
    
		if (*qinfrev == 2) 
			++qnljac[*qnlnr];
	}
    dymnl_(qinfrev, qiopt, qnnl, &qsol[1], &qres[1], &qjac[11], qtol, qtol, qinfo, &
	    qd[1], qi, idum, &c_b362);

    if (*qinfrev <= 0) {
	if (*qinfo != 1) {
		if ((*printevent&(1<<6)) || !(*printevent&(1<<10))) 
	    dymnlerr_(time, qnlnr, qinfo, qnnl, qiopt, qtol, &qsol[1], &qres[
		    1], &qnlfunc[1], &qnljac[1], qnlmax, 0);
	    *ierr___ = 1;
	} else if ((*printevent&(1<<6)) && *qinfo == 1) {
	    dymnlinf1_(qnlnr, qnnl, &qsol[1], &qres[1], &qnlfunc[1], &
		    qnljac[1], qnlmax, 0);
	}
    }
    return 0;
} /* dymnon_ */


/* Subroutine */ DYMOLA_STATIC int dymnon4_(qinfrev, qiopt, qnnl, qsol, qres, qjac, qtol, qtoldesired, 
	qinfo, qd, qi, idum, printevent, qnlnr, time, qnlfunc, qnljac, qnlmax, varnames,
	ierr___)
integer *qinfrev, *qiopt, *qnnl;
doublereal *qsol, *qres, *qjac, *qtol, *qtoldesired;
integer *qinfo;
doublereal *qd;
integer *qi, *printevent, *qnlnr;
doublereal *time;
integer *qnlfunc, *qnljac, *qnlmax, *ierr___;
integer *idum; /* integer array of length 20 */
const char*const*varnames;
{







    /* Parameter adjustments */
    --qnlfunc;
    --qnljac;

    /* Function Body */
    if (*qinfrev == -1) {
	if ((*printevent&(1<<6))) {
	    dymosimmessage_(" ", 1);
	    DymosimMessageInt_("Solving nonlinear system of equations number = ", qnlnr, 47);
	    dymosimmessagedouble_("        iteratively at Time = ", time, 30)
		    ;
	}
    }
    if ((*printevent&(1<<7)) && *qinfrev != 3) {
	dymnlinf2_(*qnnl, qsol,varnames);
    }
    if (*qinfrev > 0) {
		++qnlfunc[*qnlnr];
    
		if (*qinfrev == 2) 
			++qnljac[*qnlnr];
	}
    dymnl_(qinfrev, qiopt, qnnl, qsol, qres, qjac, qtol, qtoldesired, qinfo, 
	    qd, qi, idum, &c_b362);

    if (*qinfrev <= 0) {
	if (*qinfo != 1) {
		if ((*printevent&(1<<6)) || !(*printevent&(1<<10))) 
	    dymnlerr_(time, qnlnr, qinfo, qnnl, qiopt, qtol, qsol, qres, &qnlfunc[1], &qnljac[1], qnlmax, varnames);
	    *ierr___ = 1;
	} else if ((*printevent&(1<<6)) && *qinfo == 1) {
	    dymnlinf1_(qnlnr, qnnl, qsol, qres, &qnlfunc[1], &
		    qnljac[1], qnlmax, varnames);
	}
    }
    return 0;
} /* dymnon4_ */




LIBDS_API void dymnon6_(int*qinfrev, int qiopt, int qnnl, double*qsol,double* qres,double* qjac,
					 double qtol, double qtoldesired,
	int*qinfo, double*qd, int lqd, int*idum, int printevent, int inJacobian, int qnlnr, double time, 
	int*qnlfunc, int*qnljac, int qnlmax, const char*const*varnames,const double*nom,int*ierr___) {
#ifdef DEBUG_NL
	LocalTime=0;
#endif
	if (*qinfrev < 0) {
		if ((printevent&(1<<6))) {
			dymosimmessage_(" ", 1);
			DymosimMessageInt_("Solving nonlinear system of equations number = ", &qnlnr, 47);
			dymosimmessagedouble_("        iteratively at Time = ", &time, 30)
				;
		}
    }
    if ((printevent&(1<<7)) && *qinfrev != 3 && *qinfrev!=5) {
		dymnlinf2_(qnnl, qsol,varnames);
		if (*qinfrev!=-1) {
			dymosimmessage_("Residual:",9);
			dymnlinf2_(qnnl, qres,0);
		}
		dymosimmessage_(" ",1);
    }
    if (*qinfrev > 0) {
		++qnlfunc[qnlnr-1];
    
		if (*qinfrev == 2) 
			++qnljac[qnlnr-1];
	} else if (*qinfrev >=-40 && *qinfrev<=-10) {
		++qnlfunc[qnlnr-1];
		if (*qinfrev <= -11 && *qinfrev>-20)
			++qnljac[qnlnr-1];
	}
    dymnl_WithNominal2(qinfrev, qiopt, qnnl, qsol, qres, qjac, qtol, qtoldesired, qinfo, 
		qd, lqd, idum, nom, varnames,idum+20, printevent, inJacobian);
	
    if (*qinfrev <= 0) {
		if (*qinfo != 1) {
			if ((printevent&(1<<6)) || !(printevent&(1<<10))) 
			dymnlerr_(&time, &qnlnr, qinfo, &qnnl, &qiopt, &qtol, qsol, qres, qnlfunc, 
				qnljac, &qnlmax, varnames);
			*ierr___ = 1;
		} else if ((printevent&(1<<6)) && *qinfo == 1) {
			dymnlinf1_(&qnlnr, &qnnl, qsol, qres, qnlfunc,qnljac, &qnlmax, varnames);
		}
    }
}
LIBDS_API void dymnon5_(int*qinfrev, int qiopt, int qnnl, double*qsol,double* qres,double* qjac,
					 double qtol, double qtoldesired,
	int*qinfo, double*qd, int lqd, int*idum, int printevent, int qnlnr, double time, 
	int*qnlfunc, int*qnljac, int qnlmax, const char*const*varnames,const double*nom,int*ierr___) {
		dymnon6_(qinfrev,qiopt,qnnl,qsol,qres,qjac,qtol,qtoldesired,qinfo,qd,lqd,idum,printevent,0,
			qnlnr,time,qnlfunc,qnljac,qnlmax,varnames,nom,ierr___);
	}






DYMOLA_STATIC int dymnon3_(integer*qinfrev, integer*qiopt, integer*qnnl, double*qsol, double*qres, double*qjac, double*qtol,
			 double*qtoldesired, integer*qinfo, double*qd, integer*lqd, integer*idum, integer*printevent, 
			 integer*qnlnr, double*time, integer*qnlfunc, integer*qnljac, integer*qnlmax,integer*ierr___)
{
	return dymnon4_(qinfrev,qiopt,qnnl,qsol,qres,qjac,qtol,qtoldesired,qinfo,qd,lqd,idum,printevent,
		qnlnr,time,qnlfunc,qnljac,qnlmax,0,ierr___);
}
DYMOLA_STATIC int dymnon2_(integer*qinfrev, integer*qiopt, integer*qnnl, double*qsol,double* qres,double* qjac,double* qtol, 
					integer*qinfo, double*qd, integer*lqd, integer*idum, integer*printevent, integer*qnlnr, 
					double*time, 
					integer*qnlfunc, integer*qnljac, integer*qnlmax, integer*ierr___)
{
	return dymnon3_(qinfrev, qiopt, qnnl, qsol, qres, qjac, qtol, qtol, 
	qinfo, qd, lqd, idum, printevent, qnlnr, time, qnlfunc, qnljac, qnlmax, 
	ierr___); /* Backwards compatibility.*/
}
/* Subroutine */ DYMOLA_STATIC int handleevent_(relation, m, cp, cn, c__, init, printevent, 
	eps, t, anyevent, relation_len)
const char *relation;
integer*m;
doublereal *cp, *cn, *c__;
logical *init;
integer *printevent;
doublereal *eps, *t;
logical *anyevent;
ftnlen relation_len;
{
    /* System generated locals */
    const char* a__1[3];
    integer i__1[3];
    doublereal d__1;

    /* Builtin functions */

    /* Local variables */
    integer mold;
    doublereal delta, epsequ;
    char textline[200];


/* EVENTMAN.FOR */

/* Implements decision table for events. */

/* CopyRight: 		(C) Dynasim AB, 1993-2001 */
/* Author:  		Hilding Elmqvist */
/* Date:    		November 23, 1997 */
/* Version: 		1.9 */
/* Revisions: */
/*   May 14, 1993	Disabling three-valued logic. */
/*   May 23, 1993	Introduced BooleanChanged. */
/*   May 31, 1993        Removed debug printouts. */
/*   July 28, 1993       Changed meaning of PrintEvent */
/*   September 30, 1993  Introduced call to disstr */
/*   November 9, 1993    Introduced hysteresis. */
/*   November 12, 1993   Rewritten - simplified */
/*   March 20, 1994      EpsEqu introduced for equality test (Otter) */
/*   October 24, 1996    Simplified code for output */
/*   November 23, 1997   Fortran code adapted to C code. */
/*                       EPSLON replaced by epsmch (more efficient). */
/*                      Real number (0.9) replaced by double number (0.9D0)*/
/*                       (Otter). */

/* PrintEvent: */
/*   >=5: Log < > == */

/* ---------------------------------------------------------- */
/* 		C < -Delta	-Delta <= C 	C > Delta */
/* 				and C <= Delta */

/* 		M = -1		M = 0		M = 1 */

/* ---------------------------------------------------------- */


    delta = *eps / 1e5;
    *cp = *eps;
    *cn = -(*eps);
    epsequ = DBL_EPSILON * 10.;
    if ((d__1 = *c__ - *cp, Dymola_abs(d__1)) <= epsequ) {
	*cp = *eps * .9;
    } else if ((d__1 = *c__ - *cn, Dymola_abs(d__1)) <= epsequ) {
	*cn = *eps * -.9;
    }
    if (*init) {
	*m = -2;
    }
    mold = *m;
    if (*c__ >= -delta && *c__ <= delta && *m != 0) {
	*m = 0;
	if ((*printevent&(1<<1))) {
/* Writing concatenation */
	    i__1[0] = 11, a__1[0] = "Expression ";
	    i__1[1] = relation_len, a__1[1] = relation;
	    i__1[2] = 8, a__1[2] = " became ";
	    s_cat(textline, a__1, i__1, &c__3, 200);
	    dymosimmessagedouble2_(textline, c__, " ( == 0).", 200, 9);
	}
    } else if (*c__ > delta && *m != 1) {
	*m = 1;
	if ((*printevent&(1<<1))) {
/* Writing concatenation */
	    i__1[0] = 11, a__1[0] = "Expression ";
	    i__1[1] = relation_len, a__1[1] = relation;
	    i__1[2] = 8, a__1[2] = " became ";
	    s_cat(textline, a__1, i__1, &c__3, 200);
	    dymosimmessagedouble2_(textline, c__, " ( > 0).", 200, 8);
	}
    } else if (*c__ < -delta && *m != -1) {
	*m = -1;
	if ((*printevent&(1<<1))) {
/* Writing concatenation */
	    i__1[0] = 11, a__1[0] = "Expression ";
	    i__1[1] = relation_len, a__1[1] = relation;
	    i__1[2] = 8, a__1[2] = " became ";
	    s_cat(textline, a__1, i__1, &c__3, 200);
	    dymosimmessagedouble2_(textline, c__, " ( < 0).", 200, 8);
	}
    }
    if (*m != mold && ! (*init)) {
	*anyevent = TRUE_;
    }
    return 0;
} /* handleevent_ */

DYMOLA_STATIC void handleevent2(const char*rele,const char*sube,logical*ql,double*qp,double*qn,double eps,double qrel,logical init,int printevent,
			 logical*anyevent,logical ltz,logical invres,logical eventiter)
{
	/* New event logic for Modelica, 1999-08-25, Hans Olsson, Dynasim */
	/* In Modelica we can separate x>0,x>=0,x<0,x<=0 and we thus only need */
	/* one crossing-function and two states (e.g. x>0 or not x>0) */
	/* Furthermore the epsilon is only used if necesary, i.e. if we start close to zero. */
	/* qp is the epsilon we use and qn is the maximum allow epsilon (gives the sign). */
	double useqrel;
	logical res;
	useqrel = (eventiter) ? qrel - *qp : qrel;
	res=ltz ? (useqrel<0) : (useqrel>0);
	*qn = (res != ltz) ? - eps : eps;
	if (invres) res = !res;
	if (*ql != res || init) {
		if (!init) *anyevent=1;
		*qp=(Dymola_abs(qrel) < 100*eps || eventiter) ? *qn : *qn*1e-6;
		if ((printevent&(1<<1))) {
			char str[1000];
			sprintfC(str,"Expression %.400s became %s ( %.400s = %g )",rele,res?"true":"false",sube,qrel);
			DymosimMessage(str);
		}
	} else {
		if (!eventiter && Dymola_abs(qrel)<eps) *qp = *qn;
	}
	*ql=res;
}

LIBDS_API void handleevent3(const char*relation,integer*m,doublereal*cp,doublereal*cn,doublereal c__,logical init,
						 integer printevent, doublereal eps, doublereal t,logical*anyevent,
						 logical eventiter,doublereal largeeps,logical convergenceproblem)
/* Subroutine */ 
{
    /* System generated locals */
    const char* a__1[3];
    integer i__1[3];
    /*doublereal d__1;*/

    /* Builtin functions */

    /* Local variables */
    integer mold;
    doublereal delta, epsequ;
    char textline[200];


/* EVENTMAN.FOR */

/* Implements decision table for events. */

/* CopyRight: 		(C) Dynasim AB, 1993-2001 */
/* Author:  		Hilding Elmqvist */
/* Date:    		December 3, 1999 */
/* Version: 		1.10 */
/* Revisions: */
/*   May 14, 1993	Disabling three-valued logic. */
/*   May 23, 1993	Introduced BooleanChanged. */
/*   May 31, 1993        Removed debug printouts. */
/*   July 28, 1993       Changed meaning of PrintEvent */
/*   September 30, 1993  Introduced call to disstr */
/*   November 9, 1993    Introduced hysteresis. */
/*   November 12, 1993   Rewritten - simplified */
/*   March 20, 1994      EpsEqu introduced for equality test (Otter) */
/*   October 24, 1996    Simplified code for output */
/*   November 23, 1997   Fortran code adapted to C code. */
/*                       EPSLON replaced by epsmch (more efficient). */
/*                      Real number (0.9) replaced by double number (0.9D0)*/
/*                       (Otter). */
/*   December 3, 1999    Introduced eventiter (Hans Olsson) rewrote and renamed it*/

/* PrintEvent: */
/*   >=5: Log < > == */

/* ---------------------------------------------------------- */
/* 		C < -Delta	-Delta <= C 	C > Delta */
/* 				and C <= Delta */

/* 		M = -1		M = 0		M = 1 */

/* ---------------------------------------------------------- */


    delta = eps / 1e5;
    if (init && !eventiter) {
		*m = -2;
    }
    mold = *m;
	/* new code: use one-sided crossing functions during event iterations */
	/* This introduces a necessary hysteresis */
	       if (mold == -1 && fabs(c__)<largeeps && eventiter && fabs(*cp)>eps && c__ <      *cp) {
		*cp = *cn =    fabs(*cp);
	} else if (mold ==  1 && fabs(c__)<largeeps && eventiter && fabs(*cp)>eps && c__ >      *cn) {
		*cp = *cn = -  fabs(*cp);
	/* dangerous code: try to increase epsilon during event iterations */
	/* Necessary for some hydraulics problem */
	} else if (eventiter && mold == -1 && c__ < largeeps && c__ > *cp && convergenceproblem) {
		double neweps;
		neweps=Dymola_min(Dymola_max(eps,largeeps),10*fabs(*cp));
		*cp = *cn =  neweps;
		if ((printevent&(1<<1))) {
			a__1[0]="Warning: Event epsilon of ";i__1[0]=strlen(a__1[0]);
			a__1[1]=relation; i__1[1]=strlen(relation);
			a__1[2]=" increased to "; i__1[2]=strlen(a__1[2]);
			s_cat(textline, a__1, i__1, &c__3, 200);
			dymosimmessagedouble_(textline,&neweps,200);
		}
		if (c__ > *cp) {
			*m =  1;
			*cp = *cn = - neweps;
		}
	} else if (eventiter && mold ==  1 && c__ >-largeeps && c__ < *cn && convergenceproblem) {
		double neweps;
		neweps=Dymola_min(Dymola_max(eps,largeeps),10*fabs(*cp));
		*cp = *cn = -neweps;
		if ((printevent&(1<<1))) {
			a__1[0]="Warning: Event epsilon of ";i__1[0]=strlen(a__1[0]);
			a__1[1]=relation; i__1[1]=strlen(relation);
			a__1[2]=" increased to "; i__1[2]=strlen(a__1[2]);
			s_cat(textline, a__1, i__1, &c__3, 200);
			dymosimmessagedouble_(textline,&neweps,200);
		}
		if (c__ < *cn) {
			*m =  -1;
			*cp = *cn =   neweps;
		}
	} else {
		if (eventiter && fabs(*cp)> eps) {
			/* keep previous increase*/
			*cp = fabs(*cp);
			*cn = -(*cp);
		} else {
			/* old code */
			*cp = eps;
			*cn = -(eps);
			epsequ = DBL_EPSILON * 10.;
			if (fabs(c__ - *cp) <= epsequ) {
			*cp = eps * .9;
			} else if (fabs(c__ - *cn) <= epsequ) {
			*cn = eps * -.9;
			}
		}
		if (c__ <= delta && c__ >= - delta)
			*m = 0;
		else if (c__ > delta)
			*m = 1;
		else if (c__ < - delta) 
			*m = -1;
	}
	if ((printevent&(1<<1)) && *m != mold && (!init||(printevent&(1<<2)))) {
		/* Common routines for writing result*/
	    i__1[0] = 11, a__1[0] = "Expression ";
	    i__1[1] = strlen(relation), a__1[1] =relation ;
	    i__1[2] = 8, a__1[2] = " became ";
	    s_cat(textline, a__1, i__1, &c__3, 200);
    if (*m==0) {
/* Writing concatenation */
	    dymosimmessagedouble2_(textline, &c__, " ( == 0).", 200, 9);
    } else if (*m == 1) {
/* Writing concatenation */
	    dymosimmessagedouble2_(textline, &c__, " ( > 0).", 200, 8);
    } else if (*m == -1) {
/* Writing concatenation */
	    dymosimmessagedouble2_(textline, &c__, " ( < 0).", 200, 8);
    }
	}
    if (*m != mold && ! init) {
	*anyevent = TRUE_;
    }
} /* handleevent_ */

LIBDS_API logical handleevent4S(
						  const char*rele,
						  const char*sube,
						  logical*ql,
						  double*qp,
						  double*qn,
						  double eps,
						  double qrel,
						  double*qz,
						  logical init,
						  logical event,
						  logical idemand4,
						  int printevent,
						  logical*anyevent,
						  logical ltz,
						  logical invres,
						  logical*eventiter,
						  double largeeps,
						  logical convergenceproblem,double t) {
	double useqrel;
	logical res;
	logical neg;
#ifdef GODESS
	useqrel = qrel;
	res = ltz ? (useqrel<0) : (useqrel>0);
	neg = res != ltz;
	if (event||init) 
		*qp = *qn = neg ? - eps : eps;
	if (invres) 
		res = !res;
	*ql = res;
	qz[0]=qz[1]=0.1;
	return *ql;
#else
	/* if *ql & 2 events are enabled. Only care about event iteration in that case. */
	useqrel = (*eventiter && (event||init) && (*ql&2)) ? qrel - *qp : qrel;
	res = ltz ? (useqrel<0) : (useqrel>0);
	neg = res != ltz;
#ifdef DYMOSIM
	if (init) {
			extern int Check5(char*key);
			if (!Check5("")) {
				DymosimError("");
				return 0;
			}
		}
#endif
	if (event||init) 
		*qn = neg ? - eps : eps;
	if (invres) 
		res = !res;
	if ((idemand4) && (!(*ql&(2|4)))) {
		if ((*ql&1)==res) {
			if (fabs(qrel)>0.1) {
				/* Did not have events enabled and far from turning point */
				*ql=2 | res;
				*qp = *qn = neg ? -eps : eps; /* Remove any previous increase */
			}
		} else {
			*ql|=4; 
			/* Do not turn it on before the next event*/
			/* The problem is events that come and go */
			/*  */

/*    Sketches explaining the problem.							    */
/*    At some point in time we get a value far from the crossing.   */
/*  Since smooth crossing may appear and reappear it could			*/
/*    be any of the following:										*/
/* 																    */
/*  1															    */
/*               *													*/
/*             /													*/
/*          /														*/
/* ------------------												*/
/*       /															*/
/* 																	*/
/* 2																*/
/*                   *												*/
/*                 /												*/
/*               /													*/
/* 																	*/
/* -------------------												*/
/* 																	*/
/* 																	*/
/* 																	*/
/* 3																*/
/* 																	*/
/* 																	*/
/*                 *												*/
/* 			     /													*/
/* 	 		    /													*/
/* -------------/													*/
/* --------------------------										*/
/* 																	*/
/* 4																*/
/*                            *										*/
/*                           /										*/
/* /--\          /-\        /										*/
/* ----\-------/----\------/										*/
/*       \----/      \-----											*/
/* 																	*/
/* 																	*/
/* 																	*/
/* From the first evaluation at the end of a step					*/
/* we cannot tell cases 1 and 2 apart.								*/
/* If we enable rootfinding in the wrong state						*/
/* it will look as though we have already passed the root           */
/* which will confuse the logic.									*/
/* 																	*/
/* 																	*/
/* The bit '1<<2=4' indicate that we potentially have				*/
/* this problem. A future solution for cases 3 &4 would be          */
/* to turn off this bit after every accepted step (accepted         */
/* means that event handling is complete)							*/

		}
	}
	if ((event && ((*ql&2) || fabs(qrel)>0.1) && res!=(*ql&1)) || init) {
		if (!init) {
			*anyevent=1;
		}
		if (*eventiter) {
			if (Dymola_abs(qrel)<largeeps && convergenceproblem) {
				double neweps;
				neweps=Dymola_min(Dymola_max(eps,largeeps),Dymola_max(eps,10*fabs(*qp)));
				if (neweps>eps && (printevent&(1<<1))) {
					char str[1000];
					sprintfC(str,"Warning: Event epsilon of %.400s increased to %g.\n",rele,neweps);
					DymosimMessage(str);
				}
				*qp = neg ? -neweps : neweps;
			} else {
				double d;
				d=fabs(*qp);
				if (d<eps) d=eps;
				/* keep old increase (if any) */
				*qp = neg ? -d : d;
			}
		} else {
			*qp=(Dymola_abs(qrel) < 100*eps) ? *qn : *qn*1e-6;
		}
		if ((printevent&(1<<1)) && (!init || (printevent&(1<<2)) && !(eventiter && (*ql&1)==res))) {
			char str[1000];
			sprintfC(str,"Expression %.400s became %s ( %.400s = %g )",rele,res?"true":"false",sube,qrel);
			DymosimMessage(str);
		}
		if (*eventiter)
			*ql = res | (*ql&2); /* Keep events if iteration */
		else {
			*ql = res;
		}
		if (fabs(qrel)>0.1) *ql|=2; /* Enable events if far from it */
	} else {
		if ((!*eventiter)&&event) {
			if (fabs(qrel)<=0.1) {
				*ql&=1; /* Disable events */
			}
			/* if (eventiter): Keep increase of qp. */
			if (Dymola_abs(qrel)<eps) 
				*qp = *qn; /* Reasonable value close to zero-crossing */
			else if (Dymola_abs(qrel)>largeeps && Dymola_abs(qrel) > 100*eps) 
				*qp = *qn * 1e-6; /* Decrease it if far from zero-crossing */
			/* otherwise: keep old increase */
		}
		if (!(*ql&2)) {
			/* Not event generating => noEvent */
			if (event)
				*ql=res; /* noEvent logic */
			else
				*ql=res|(*ql&4); /* noEvent but keep disabling */
		}
	} 
	*eventiter=1;

	qz[0]=0.1;
	if (*ql&2) {
		/* Enable event */
		if ((*ql&1)!=(ltz==invres)) {
			qz[0]=-qrel+fabs(*qp);
		} else {
			qz[0]=qrel+fabs(*qp);
		}
	}
	qz[1]=qz[0];
	return (*ql&1);
#endif
}
 
LIBDS_API void handleevent5(const char*rele,const char*sube,logical*ql,double*qp,double*qn,double eps,double qrel,logical init,int printevent,
			 logical*anyevent,logical*anyevent2,logical ltz,logical invres,logical eventiter,doublereal largeeps,logical convergenceproblem)
{
	/* New event logic for Modelica, 1999-08-25, Hans Olsson, Dynasim */
	/* In Modelica we can separate x>0,x>=0,x<0,x<=0 and we thus only need */
	/* one crossing-function and two states (e.g. x>0 or not x>0) */
	/* Furthermore the epsilon is only used if necesary, i.e. if we start close to zero. */
	/* qp is the epsilon we use and qn is the maximum allow normal epsilon (gives the sign). */
	/* 1999-12-03 Allowed even further increase of qp during event iterations, but give warning in these cases. */
	double useqrel;
	logical res;
	logical neg;
	useqrel = eventiter ? qrel - *qp : qrel;
	res=ltz ? (useqrel<0) : (useqrel>0);
	neg=res != ltz;
	*qn = neg ? - eps : eps;
#ifdef DYMOSIM
	if (init) {
			extern int Check5(char*key);
			if (!Check5("")) {
				DymosimError("");
				return ;
			}
		}
#endif
	if (invres) 
		res = !res;
	if (*ql != res || init) {
		if (!init) {
			*anyevent=1;
			if (anyevent2) *anyevent2=1;
		}
		if (eventiter) {
			if (Dymola_abs(qrel)<largeeps && convergenceproblem) {
				double neweps;
				neweps=Dymola_min(Dymola_max(eps,largeeps),Dymola_max(eps,10*fabs(*qp)));
				if (neweps>eps && (printevent&(1<<1))) {
					char str[1000];
					sprintfC(str,"Warning: Event epsilon of %.400s increased to %g.\n",rele,neweps);
					DymosimMessage(str);
				}
				*qp = neg ? -neweps : neweps;
			} else {
				double d;
				d=fabs(*qp);
				if (d<eps) d=eps;
				/* keep old increase (if any) */
				*qp = neg ? -d : d;
			}
		} else {
			*qp=(Dymola_abs(qrel) < 100*eps) ? *qn : *qn*1e-6;
		}
		if ((printevent&(1<<1)) && (!init || (printevent&(1<<2)) && !(eventiter && *ql==res))) {
			char str[1000];
			sprintfC(str,"Expression %.400s became %s ( %.400s = %g )",rele,res?"true":"false",sube,qrel);
			DymosimMessage(str);
		}
	} else if (!eventiter) {
		if (Dymola_abs(qrel)<eps) 
			*qp = *qn; /* Reasonable value close to zero-crossing */
		else if (Dymola_abs(qrel)>largeeps && Dymola_abs(qrel) > 100*eps) 
			*qp = *qn * 1e-6; /* Decrease it if far from zero-crossing */
		/* otherwise: keep old increase */
	} /* if (eventiter): Keep increase of qp. */
	*ql=res;
}
LIBDS_API void handleevent4(const char*rele,const char*sube,logical*ql,double*qp,double*qn,double eps,double qrel,logical init,int printevent,
			 logical*anyevent,logical ltz,logical invres,logical eventiter,doublereal largeeps,logical convergenceproblem)
{
	handleevent5(rele, sube, ql, qp, qn, eps, qrel, init, printevent, anyevent, 0, ltz, invres, eventiter,largeeps, convergenceproblem);
}


LIBDS_API void roundEvent5(const char*expr,double (*func)(double),double qrel,double*val,double*fuzz,double eps,
						logical posi,logical negi,logical reli,
						logical init,logical printevent,logical *anyevent,logical *anyevent2,logical eventiter,double leps,logical slow)
{
	/* For now: ignore eventiter. */
	double res;
	res= (*func)(qrel);
	*fuzz=eps;
	if (res!=*val) {
		if (!init) {
			*anyevent = 1;
			if (anyevent2) *anyevent2=1;
		}
		if ((printevent&(1<<1)) && (!init || (printevent&(1<<2)))) {
			char str[1000];
			sprintfC(str,"Expression %.400s became %g ( unrounded = %.12g )",
				expr,res,qrel);
			DymosimMessage(str);
		}
	}
	*val = res;
}
LIBDS_API void roundEvent4(const char*expr,double (*func)(double),double qrel,double*val,double*fuzz,double eps,
						logical posi,logical negi,logical reli,
						logical init,logical printevent,logical *anyevent,logical eventiter,double leps,logical slow)
{
	roundEvent5(expr,func,qrel,val,fuzz,eps,posi,negi,reli,init,printevent,anyevent,0,eventiter,leps,slow);
}

/* Subroutine */ LIBDS_API int booleanchanged_(varname, var, newvar, printevent, 
	anyevent, varname_len)
const char *varname;
logical *var, *newvar;
integer *printevent;
logical *anyevent;
ftnlen varname_len;
{
    /* System generated locals */
    const char* a__1[4];
    integer i__1[4];

    /* Builtin functions */

    /* Local variables */
    char newval[6], textline[200];


/* Logs changes to boolean state variables. */


    if (*var != *newvar) {
	*anyevent = TRUE_;
	if (*printevent&(1<<1)) {
	    if (*newvar) {
		s_copy(newval, "true.", 6, 5);
	    } else {
		s_copy(newval, "false.", 6, 6);
	    }
/* Writing concatenation */
	    i__1[0] = 9, a__1[0] = "Variable ";
	    i__1[1] = varname_len, a__1[1] = varname;
	    i__1[2] = 12, a__1[2] = " changed to ";
	    i__1[3] = 6, a__1[3] = newval;
	    s_cat(textline, a__1, i__1, &c__4, 200);
	    dymosimmessage_(textline, 200);
	}
    }
    return 0;
} /* booleanchanged_ */


#ifndef GODESS
/* Subroutine */ DYMOLA_STATIC int dgeequ_(const integer*m, const integer*n, doublereal *a, const integer *
	lda, doublereal *r, doublereal *c, doublereal *rowcnd, doublereal *
	colcnd, doublereal *amax, integer *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       March 31, 1993   


    Purpose   
    =======   

    DGEEQU computes row and column scalings intended to equilibrate an   
    M-by-N matrix A and reduce its condition number.  R returns the row   
    scale factors and C the column scale factors, chosen to try to make   
    the largest element in each row and column of the matrix B with   
    elements B(i,j)=R(i)*A(i,j)*C(j) have absolute value 1.   

    R(i) and C(j) are restricted to be between SMLNUM = smallest safe   
    number and BIGNUM = largest safe number.  Use of these scaling   
    factors is not guaranteed to reduce the condition number of A but   
    works well in practice.   

    Arguments   
    =========   

    M       (input) INTEGER   
            The number of rows of the matrix A.  M >= 0.   

    N       (input) INTEGER   
            The number of columns of the matrix A.  N >= 0.   

    A       (input) DOUBLE PRECISION array, dimension (LDA,N)   
            The M-by-N matrix whose equilibration factors are   
            to be computed.   

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= max(1,M).   

    R       (output) DOUBLE PRECISION array, dimension (M)   
            If INFO = 0 or INFO > M, R contains the row scale factors   
            for A.   

    C       (output) DOUBLE PRECISION array, dimension (N)   
            If INFO = 0,  C contains the column scale factors for A.   

    ROWCND  (output) DOUBLE PRECISION   
            If INFO = 0 or INFO > M, ROWCND contains the ratio of the   
            smallest R(i) to the largest R(i).  If ROWCND >= 0.1 and   
            AMAX is neither too large nor too small, it is not worth   
            scaling by R.   

    COLCND  (output) DOUBLE PRECISION   
            If INFO = 0, COLCND contains the ratio of the smallest   
            C(i) to the largest C(i).  If COLCND >= 0.1, it is not   
            worth scaling by C.   

    AMAX    (output) DOUBLE PRECISION   
            Absolute value of largest matrix element.  If AMAX is very   
            close to overflow or very close to underflow, the matrix   
            should be scaled.   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value   
            > 0:  if INFO = i,  and i is   
                  <= M:  the i-th row of A is exactly zero   
                  >  M:  the (i-M)-th column of A is exactly zero   

    ===================================================================== 
  


       Test the input parameters.   

    
   Parameter adjustments   
       Function Body */
    /* System generated locals */
    integer /*a_dim1, a_offset,*/ i__1, i__2;
    doublereal d__1, d__2 /*, d__3*/;
    /* Local variables */
    integer i, j;
    doublereal rcmin, rcmax;
    doublereal bignum, smlnum;
    const integer mp=*m,np=*n,ldap=*lda; 


#define R(I) r[(I)-1]
#define C(I) c[(I)-1]

#define A(I,J) a[(I)-1 + ((J)-1)* ( ldap)]
	
	bignum=1./DBL_MIN;
	smlnum=DBL_MIN;	/* Get machine constants. */

    *info = 0;
    if (mp < 0) {
	*info = -1;
    } else if (np < 0) {
	*info = -2;
    } else if (ldap < Dymola_max(1,mp)) {
	*info = -4;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DGEEQU", &i__1, 6);
	return 0;
    }

/*     Quick return if possible */

    if (mp == 0 || np == 0) {
	*rowcnd = 1.;
	*colcnd = 1.;
	*amax = 0.;
	return 0;
    }


/*     Compute row scale factors. */

    i__1 = mp;
    for (i = 1; i <= mp; ++i) {
	R(i) = 0.;
/* L10: */
    }

/*     Find the maximum element in each row. */

    i__1 = np;
    for (j = 1; j <= np; ++j) {
	i__2 = mp;
	for (i = 1; i <= mp; ++i) {
/* Computing MAX Modified by Hans*/
	    doublereal Aij;
	    Aij=A(i,j);
	    if (Aij<0) Aij=-Aij;
	    if (Aij>R(i)) R(i)=Aij;
/* L20: */
	}
/* L30: */
    }

/*     Find the maximum and minimum scale factors. */

    rcmin = bignum;
    rcmax = 0.;
    i__1 = mp;
    for (i = 1; i <= mp; ++i) {
/* Computing MAX */
	d__1 = rcmax, d__2 = R(i);
	rcmax = Dymola_max(d__1,d__2);
/* Computing MIN */
	d__1 = rcmin, d__2 = R(i);
	rcmin = Dymola_min(d__1,d__2);
/* L40: */
    }
    *amax = rcmax;

    if (rcmin == 0.) {

/*        Find the first zero scale factor and return an error code. */

	i__1 = mp;
	for (i = 1; i <= mp; ++i) {
	    if (R(i) == 0.) {
		*info = i;
		return 0;
	    }
/* L50: */
	}
    } else {

/*        Invert the scale factors. */

	i__1 = mp;
	for (i = 1; i <= mp; ++i) {
/* Computing MIN   */
/* Computing MAX */
	    d__2 = R(i);
	    d__1 = Dymola_max(d__2,smlnum);
	    R(i) = 1. / Dymola_min(d__1,bignum);
/* L60: */
	}

/*        Compute ROWCND = min(R(I)) / max(R(I)) */

	*rowcnd = Dymola_max(rcmin,smlnum) / Dymola_min(rcmax,bignum);
    }

/*     Compute column scale factors */

    i__1 = np;
    for (j = 1; j <= np; ++j) {
	C(j) = 0.;
/* L70: */
    }

/*     Find the maximum element in each column,   */
/*       assuming the row scaling computed above. */

    i__1 = np;
    for (j = 1; j <= np; ++j) {
	i__2 = mp;
	for (i = 1; i <= mp; ++i) {
/* Computing MAX: Modified by Hans */
	    doublereal Aij;
	    Aij=A(i,j);
	    if (Aij<0) Aij=-Aij;
	    Aij*=R(i);
	    if (Aij>C(j)) C(j)=Aij;
/* L80: */
	}
/* L90: */
    }

/*     Find the maximum and minimum scale factors. */

    rcmin = bignum;
    rcmax = 0.;
    i__1 = np;
    for (j = 1; j <= np; ++j) {
/* Computing MIN */
	d__1 = rcmin, d__2 = C(j);
	rcmin = Dymola_min(d__1,d__2);
/* Computing MAX */
	d__1 = rcmax, d__2 = C(j);
	rcmax = Dymola_max(d__1,d__2);
/* L100: */
    }

    if (rcmin == 0.) {

/*        Find the first zero scale factor and return an error code. */

	i__1 = np;
	for (j = 1; j <= np; ++j) {
	    if (C(j) == 0.) {
		*info = mp + j;
		return 0;
	    }
/* L110: */
	}
    } else {

/*        Invert the scale factors. */

	i__1 = np;
	for (j = 1; j <= np; ++j) {
/* Computing MIN   
   Computing MAX */
	    d__2 = C(j);
	    d__1 = Dymola_max(d__2,smlnum);
	    C(j) = 1. / Dymola_min(d__1,bignum);
/* L120: */
	}

/*        Compute COLCND = min(C(J)) / max(C(J)) */

	*colcnd = Dymola_max(rcmin,smlnum) / Dymola_min(rcmax,bignum);
    }

    return 0;

/*     End of DGEEQU */

} /* dgeequ_ */
#undef R
#undef C
#undef A

/* Subroutine */ DYMOLA_STATIC int dlaqge_(const integer*m, const integer*n, doublereal *a, const integer *
	lda, doublereal *r, doublereal *c, doublereal *rowcnd, doublereal *
	colcnd, doublereal *amax, char *equed)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       February 29, 1992   


    Purpose   
    =======   

    DLAQGE equilibrates a general M by N matrix A using the row and   
    scaling factors in the vectors R and C.   

    Arguments   
    =========   

    M       (input) INTEGER   
            The number of rows of the matrix A.  M >= 0.   

    N       (input) INTEGER   
            The number of columns of the matrix A.  N >= 0.   

    A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)   
            On entry, the M by N matrix A.   
            On exit, the equilibrated matrix.  See EQUED for the form of 
  
            the equilibrated matrix.   

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= max(M,1).   

    R       (input) DOUBLE PRECISION array, dimension (M)   
            The row scale factors for A.   

    C       (input) DOUBLE PRECISION array, dimension (N)   
            The column scale factors for A.   

    ROWCND  (input) DOUBLE PRECISION   
            Ratio of the smallest R(i) to the largest R(i).   

    COLCND  (input) DOUBLE PRECISION   
            Ratio of the smallest C(i) to the largest C(i).   

    AMAX    (input) DOUBLE PRECISION   
            Absolute value of largest matrix entry.   

    EQUED   (output) CHARACTER*1   
            Specifies the form of equilibration that was done.   
            = 'N':  No equilibration   
            = 'R':  Row equilibration, i.e., A has been premultiplied by 
  
                    diag(R).   
            = 'C':  Column equilibration, i.e., A has been postmultiplied 
  
                    by diag(C).   
            = 'B':  Both row and column equilibration, i.e., A has been   
                    replaced by diag(R) * A * diag(C).   

    Internal Parameters   
    ===================   

    THRESH is a threshold value used to decide if row or column scaling   
    should be done based on the ratio of the row or column scaling   
    factors.  If ROWCND < THRESH, row scaling is done, and if   
    COLCND < THRESH, column scaling is done.   

    LARGE and SMALL are threshold values used to decide if row scaling   
    should be done based on the absolute size of the largest matrix   
    element.  If AMAX > LARGE or AMAX < SMALL, row scaling is done.   

    ===================================================================== 
  


       Quick return if possible   

    
   Parameter adjustments   
       Function Body */
    /* System generated locals */
    integer /*a_dim1, a_offset,*/ i__1, i__2;
    /* Local variables */
    integer i, j;
    doublereal large;
	doublereal small;
	doublereal cj;
    const integer mp=*m,np=*n,ldap=*lda; 

#define R(I) r[(I)-1]
#define C(I) c[(I)-1]

#define A(I,J) a[(I)-1 + ((J)-1)* ( ldap)]
	
	large=1. / (DBL_MIN / (DBL_EPSILON*FLT_RADIX));
	small=DBL_MIN / (DBL_EPSILON*FLT_RADIX);

    if (mp <= 0 || np <= 0) {
	*(unsigned char *)equed = 'N';
	return 0;
    }

    if (*rowcnd >= .1 && *amax >= small && *amax <= large) {
		
		/*        No row scaling */
		
		if (*colcnd >= .1) {
			
			/*           No column scaling */
			
			*(unsigned char *)equed = 'N';
		} else {
			
			/*           Column scaling */
			
			i__1 = np;
			for (j = 1; j <= np; ++j) {
				cj = C(j);
				i__2 = mp;
				for (i = 1; i <= mp; ++i) {
					A(i,j) *= cj ;
					/* L10: */
				}
				/* L20: */
			}
			*(unsigned char *)equed = 'C';
		}
    } else if (*colcnd >= .1) {
		
		/*        Row scaling, no column scaling */
		
		i__1 = np;
		for (j = 1; j <= np; ++j) {
			i__2 = mp;
			for (i = 1; i <= mp; ++i) {
				A(i,j) *= R(i);
				/* L30: */
			}
			/* L40: */
		}
		*(unsigned char *)equed = 'R';
    } else {
		
		/*        Row and column scaling */
		
		i__1 = np;
		for (j = 1; j <= np; ++j) {
			cj = C(j);
			i__2 = mp;
			for (i = 1; i <= mp; ++i) {
				A(i,j) *= cj * R(i) ;
				/* L50: */
			}
			/* L60: */
		}
		*(unsigned char *)equed = 'B';
    }

    return 0;

/*     End of DLAQGE */

} /* dlaqge_ */
#undef R
#undef C
#undef A

#endif


#endif
   /* DSblock model generated by Dymola from Modelica model sine
 Dymola Version 2016 (64-bit), 2015-04-15 translated this at Thu Oct 06 14:44:22 2016

   */
#ifndef DYN_MULTINSTANCE
#define DYN_MULTINSTANCE 1
#endif

#include <matrixop.h>
/* Declaration of C-structs */
/* Prototypes for functions used in model */
/* Codes used in model */
/* DSblock C-code: */

#define NX_    0
#define NX2_   0
#define NU_    0
#define NY_    1
#define NW_    6
#define NP_    6
#define NPS_   0
#define MAXAuxStr_   0
#define MAXAuxStrLen_   500
#define NHash1_ -173639171
#define NHash2_ -1922225477
#define NHash3_ 2121431217
#define NI_    0
#define NRelF_ 0
#define NRel_  0
#define NTim_  1
#define NSamp_ 0
#define NCons_ 0
#define NA_    1
#define SizePre_ 0
#define SizeEq_ 0
#define SizeDelay_ 0
#define QNLmax_ 0
#define MAXAux 0
#define NrDymolaTimers_ 0
#define NWhen_ 0
#define NCheckIf_ 0
#define NGlobalHelp_ 0
#define NGlobalHelpI_ 0
#ifndef NExternalObject_
#define NExternalObject_ 0
#endif
#include <moutil.c>
PreNonAliasDef(0)
PreNonAliasDef(1)
PreNonAliasDef(2)
PreNonAliasDef(3)
PreNonAliasDef(4)
PreNonAliasDef(5)
#if !defined(DYM2CCUR)
 DYMOLA_STATIC const char*modelName="sine";
#endif
DYMOLA_STATIC const char*usedLibraries[]={0};
DYMOLA_STATIC const char*dllLibraryPath[]={0};
DYMOLA_STATIC const char*default_dymosim_license_filename=
 "c:/users/ridouae/appdata/roaming/dynasim/dymola.lic";
DYMOLA_STATIC const char*GUIDString="{9f1d2555-743d-4e30-a6f9-f678fa064938}";
DYMOLA_STATIC const double cvodeTolerance=1E-005;
#include <dsblock1.c>

/* Define variable names. */

#define Sections_

TranslatedEquations

InitialSection
DYNX(W_,5) = 3.141592653589793;
BoundParameterSection
DYNX(W_,0) = DYNX(DP_,1);
DYNX(W_,1) = DYNX(DP_,2);
DYNX(W_,2) = DYNX(DP_,3);
DYNX(W_,3) = DYNX(DP_,4);
DYNX(W_,4) = DYNX(DP_,5);
InitialSection
InitialSection
InitialStartSection
InitialSection
DefaultSection
InitializeData(0)
InitialSection
InitialSection2
DYNX(W_,3) = DYNX(DP_,4);
DYNX(W_,4) = DYNX(DP_,5);
DYNX(W_,0) = DYNX(DP_,1);
DYNX(W_,1) = DYNX(DP_,2);
DYNX(W_,2) = DYNX(DP_,3);
DYNX(Y_,0) = DYNX(W_,3)+(IF LessTime(DYNX(W_,4), 0) THEN 0 ELSE DYNX(W_,0)*sin(
  6.283185307179586*DYNX(W_,1)*(DYNTime-DYNX(W_,4))+DYNX(W_,2)));
InitialSection2
Init_=false;InitializeData(2);Init_=true;
EndInitialSection

OutputSection
DYNX(Y_,0) = DYNX(W_,3)+(IF LessTime(DYNX(W_,4), 0) THEN 0 ELSE DYNX(W_,0)*sin(
  6.283185307179586*DYNX(W_,1)*(DYNTime-DYNX(W_,4))+DYNX(W_,2)));

DynamicsSection

AcceptedSection1

AcceptedSection2

DefaultSection
InitializeData(1)
EndTranslatedEquations

#include <dsblock6.c>

PreNonAliasNew(0)
StartNonAlias(0)
DeclareVariable("sine.amplitude", "Amplitude of sine wave", 0.0, 0.0,0.0,0.0,0,513)
DeclareVariable("sine.freqHz", "Frequency of sine wave [Hz]", 1, 0.0,0.0,0.0,0,513)
DeclareVariable("sine.phase", "Phase of sine wave [rad|deg]", 0.0, 0.0,0.0,0.0,0,513)
DeclareVariable("sine.offset", "Offset of output signal", 0.0, 0.0,0.0,0.0,0,513)
DeclareVariable("sine.startTime", "Output = offset for time < startTime [s]", \
0.0, 0.0,0.0,0.0,0,513)
DeclareAlias2("sine.y", "Connector of Real output signal", "y", 1, 3, 0, 0)
DeclareVariable("sine.pi", "", 3.141592653589793, 0.0,0.0,0.0,0,2561)
DeclareOutput("y", "", 0, 0.0, 0.0,0.0,0.0,0,512)
DeclareParameter("pi", "", 0, 3.141592653589793, 0.0,0.0,0.0,0,560)
DeclareParameter("amp", "Sinus Amplite", 1, 1, 0.0,0.0,0.0,0,560)
DeclareParameter("sineFreq", "Sinus frequency [Hz]", 2, 0.4, 0.0,0.0,0.0,0,560)
DeclareParameter("phase", "Sinus phase [rad]", 3, 0, 0.0,0.0,0.0,0,560)
DeclareParameter("sineOffset", "Sinus offset", 4, 0, 0.0,0.0,0.0,0,560)
DeclareParameter("start", "Sinus startTime [s]", 5, 0, 0.0,0.0,0.0,0,560)
EndNonAlias(0)

#define DymolaHaveUpdateInitVars 1
#include <dsblock5.c>

DYMOLA_STATIC void UpdateInitVars(double*time, double* X_, double* XD_, double* U_, double* DP_, int IP_[], Dymola_bool LP_[], double* F_, double* Y_, double* W_, double QZ_[], double duser_[], int iuser_[], void*cuser_[],struct DYNInstanceData*did_) {
static Real initStore[1];
}
StartDataBlock
EndDataBlock
