


#if defined(__INTEL_LLVM_COMPILER ) || defined(__INTEL_COMPILER )

#include "f77func.h"
#ifndef ___INTEGER_definition___
typedef int INTEGER; /* for fortran integer */
#define ___INTEGER_definition___
#endif

#ifndef ___logical_definition___
typedef long int logical;
#define ___logical_definition___
#endif

#ifndef ___dcomplex_definition___
typedef struct { double r,i; } dcomplex;
#define ___dcomplex_definition___
#endif

#define MKL_Complex16 dcomplex
//#define INTEGER int
#include <mkl.h>
#include <mkl_blas.h>
#include <mkl_lapack.h>
#include <mkl_scalapack.h>


//prototype for ELPA
void F77_NAME(solve_evp_real,SOLVE_EVP_REAL)(const int*, const int*, double*, const int*, double*, double*, const int*, const int*, const int*, const int*);
void F77_NAME(solve_evp_complex,SOLVE_EVP_COMPLEX)(const int*, const int*, dcomplex*, const int*, dcomplex*, dcomplex*, const int*, const int*, const int*, const int*);
void F77_NAME(elpa_solve_evp_real_2stage_double_impl,ELPA_SOLVE_EVP_COMPLEX_2STAGE_DOUBLE_IMPL)(const int*, int*, double*, const int*, double*, double*, const int*, const int*, const int*, const int*, const int*, const int*);
void F77_NAME(elpa_solve_evp_complex_2stage_double_impl,ELPA_SOLVE_EVP_REAL_2STAGE_DOUBLE_IMPL)(const int*, int*, dcomplex*, const int*, double*, dcomplex*, const int*, const int*, const int*, const int*, const int*, const int*);
void F77_NAME(pdgemm,PDGEMM)(const char*,const char*,const int*,const int*,const int*,double*,double*,int*,int*,int*,double*,int*,int*,int*,double*,double*,int*,int*,int*);
void F77_NAME(pzgemm,PZGEMM)(const char*,const char*,const int*,const int*,const int*,dcomplex*,dcomplex*,int*,int*,int*,dcomplex*,int*,int*,int*,dcomplex*,dcomplex*,int*,int*,int*);

//prototyoe for ScaLAPACK
int Csys2blacs_handle(MPI_Comm comm);
void sl_init_(int* icontext, int* nprow, int* npcolumn);
void blacs_get_(int* icontext, const int* what, int* val);

    //void blacs_gridinit_(int* ConTxt, const char* layout, const int* nprow, const int* npcol);
void Cblacs_gridinit(int* ConTxt, const char* layout, const int nprow, const int npcol);

void blacs_gridinfo_(int* icontext, int* nprow, int* npcolumn, int* myrow, int* mycolumn);
void blacs_gridmap_(int* icontext, int* usermap, int* ldumap, int* nprow, int* npcolumn);


int numroc_(const int* n, const int* nb, const int* iproc,
        const int* isrcproc, const int* nprocs);

int indxg2p_(const int* n, const int* nb, const int* iproc,
                const int* isrcproc, const int* nprocs);

        //note: meaning of element of descinit_
        //parameter( block_cyclic_2d = 1, dlen_ = 9, dtype_ = 1,
        //           ctxt_ = 2, m_ = 3, n_ = 4, mb_ = 5, nb_ = 6,
        //                       rsrc_ = 7, csrc_ = 8, lld_ = 9 )
void descinit_(int* desc, const int* m, const int* n,
        const int* mb, const int* nb, const int* irsrc,
        const int* icsrc, const int* ictxt, const int* lld,
        int* info);

void blacs_exit_(int* cont);


void blacs_gridexit_(int* icontext);
void Cblacs_gridexit(int icontext);
void Cfree_blacs_system_handle(int bhandle2);
void Cblacs_barrier(int icontext, const char* ch);

#else


#include "f77func.h"
#ifdef blaswrap
#define zgemm_ f2c_zgemm
#define zcopy_ f2c_zcopy
#endif

#ifndef ___INTEGER_definition___
typedef int INTEGER; /* for fortran integer */
#define ___INTEGER_definition___ 
#endif

#ifndef ___logical_definition___
typedef long int logical;
#define ___logical_definition___ 
#endif

#ifndef ___dcomplex_definition___
typedef struct { double r,i; } dcomplex;
#define ___dcomplex_definition___ 
#endif

void zhegv_(INTEGER *itype, char *jobz, char *uplo,
            INTEGER *n, dcomplex *a, INTEGER *lda,
            dcomplex *b, INTEGER *ldb,
            double *w, dcomplex *work,
            INTEGER *lwork, double *rwork, INTEGER *info);

int zgesvd_(char *jobu, char *jobvt, INTEGER *m, INTEGER *n, 
	    dcomplex *a, INTEGER *lda,
	    double *s, dcomplex *u, 
	    INTEGER *ldu, dcomplex *vt,
	    INTEGER *ldvt, dcomplex *work, 
	    INTEGER *lwork, double *rwork, INTEGER *info);

int dgesv_(INTEGER *n, INTEGER *nrhs, double *a, INTEGER *lda, INTEGER *ipiv, double *b,
           INTEGER *ldb, INTEGER *info);

void dsyev_(char *JOBZ, char *UPLO,INTEGER *N,double *A,INTEGER *LDA,double *W,
             double *WORK,INTEGER *LWORK, INTEGER *INFO);
void dsbev_(char *JOBZ, char *UPLO, INTEGER *N, INTEGER *KD, double *AB, INTEGER *LDAB,
         double *W, double *Z,
        INTEGER * LDZ,double *WORK, INTEGER *INFO );
void dgesvd_(char *JOBU, char *JOBVT, INTEGER *M, INTEGER *N, double *A, INTEGER *LDA, double *S, double *U,
        INTEGER *LDU, double *VT, INTEGER *LDVT, double *WORK, INTEGER *LWORK, INTEGER *INFO );

void dsygv_(INTEGER *itype, char *jobz, char *uplo, INTEGER *n, 
            double *a, INTEGER *lda, double *b, INTEGER *ldb, 
	    double *w, double *work, INTEGER *lwork, INTEGER *info);

void dsyevx_(char *JOBZ, char *RANGE, char *UPLO, INTEGER *N, double *A, INTEGER *LDA, 
          double *VL, double *VU,
          INTEGER *IL, INTEGER *IU, double *ABSTOL, INTEGER *M, double *W, double *Z, 
           INTEGER *LDZ, double *WORK,
          INTEGER *LWORK, INTEGER *IWORK, INTEGER *IFAIL, INTEGER *INFO);
void dsyevd_(char *JOBZ, char *UPLO, INTEGER *N, double *A, INTEGER *LDA, double *W, double *WORK, 
         INTEGER *LWORK, INTEGER *IWORK,
                       INTEGER  *LIWORK, INTEGER *INFO);
void dsyevr_(char *JOBZ, char *RANGE, char *UPLO, INTEGER *N, double *A, INTEGER *LDA, double *VL, double *VU, 
          INTEGER *IL, INTEGER *IU,
                       double  *ABSTOL, INTEGER *M, double *W, double *Z,INTEGER *LDZ, INTEGER *ISUPPZ, double *WORK, 
        INTEGER *LWORK, INTEGER *IWORK, INTEGER *LIWORK, INTEGER *INFO);

int dstevx_(char *JOBZ, char *RANGE, INTEGER *N, double *D, double *E, double *VL, double *VU, INTEGER *IL, INTEGER *IU, double *ABSTOL,
                        INTEGER *M, double *W, double *Z, INTEGER *LDZ, double *WORK, INTEGER *IWORK, INTEGER *IFAIL, INTEGER *INFO);

void dgetrf_(INTEGER *M, INTEGER *N,double *A, INTEGER *LDA, INTEGER *IPIV, INTEGER *INFO);
void dgetri_( INTEGER *N, double *A, INTEGER *LDA, INTEGER *IPIV, double *WORK, INTEGER *LWORK, INTEGER *INFO);

int dsysv_(char *UPLO, INTEGER *N, INTEGER *NRHS, double *A, INTEGER *LDA, INTEGER *IPIV, double *B, INTEGER *LDB, double *WORK,
                       INTEGER *LWORK, INTEGER *INFO);

int dstedc_(char *COMPZ , INTEGER *N , double *D , double *E , double *Z , 
       INTEGER *LDZ , double *WORK , INTEGER *LWORK , INTEGER *IWORK , INTEGER *LIWORK , 
       INTEGER *INFO);

int dstegr_(char *JOBZ , char *RANGE , INTEGER *N , double *D , double *E , 
       double *VL , double *VU , INTEGER *IL , INTEGER *IU , double *ABSTOL , 
       INTEGER *M , double *W , double *Z , INTEGER *LDZ , INTEGER *ISUPPZ , 
       double *WORK , INTEGER *LWORK , INTEGER *IWORK , INTEGER *LIWORK , INTEGER *INFO
       );

int dsteqr_(char *compz, INTEGER *n, double *d__, 
	    double *e, double *z__, INTEGER *ldz, double *work, INTEGER *info);

int dgeev_(char *jobvl, char *jobvr, INTEGER *n, double *a, INTEGER *lda, 
           double *wr, double *wi, double *vl, 
	   INTEGER *ldvl, double *vr, INTEGER *ldvr, double *work, 
	   INTEGER *lwork, INTEGER *info);

void zheev_(char *JOBZ, char *UPLO, INTEGER *N, dcomplex *A, INTEGER *LDA, double *W, dcomplex *WORK, INTEGER *LWORK, 
       double *RWORK, INTEGER *INFO);
void zheevx_(char *JOBZ, char *RANGE, char *UPLO, INTEGER *N, dcomplex *A, INTEGER *LDA, double *VL, double *VU, INTEGER *IL, INTEGER *IU,
                        double *ABSTOL, INTEGER *M, double *W, dcomplex *Z, INTEGER *LDZ, dcomplex *WORK, INTEGER *LWORK, double *RWORK,
                        INTEGER *IWORK, INTEGER *IFAIL, INTEGER *INFO);


void zheev_(char *JOBZ, char *UPLO, int *N, dcomplex *A, int *LDA, double *W, dcomplex *WORK, int *LWORK, double *RWORK,
                       int *INFO);
void zheevx_(char *JOBZ, char *RANGE, char *UPLO, int *N, dcomplex *A, int *LDA, double *VL, double *VU, int *IL, int *IU,
                        double *ABSTOL, int *M, double *W, dcomplex *Z, int *LDZ, dcomplex *WORK, int *LWORK, double *RWORK,
                        int *IWORK, int *IFAIL, int *INFO);

int dgemm_(char *transa, char *transb, int *m, int *n, int *k, 
           double *alpha, double *a, int *lda,
	   double *b, int *ldb, double *beta, double *c__,
	   int *ldc);

void zgemm_(char* TRANSA, char* TRANSB, int * M, int * N,int *K, dcomplex *alpha, 
         dcomplex *A, int *LDA, dcomplex *B, int*LDB, dcomplex *beta, dcomplex *C, int *LDC);
void zgetrf_(int *m, int *n, dcomplex *a,int *lda,int *ipvt, int *info );
void zgetri_(int *n,dcomplex *a,int *lda, int *ipvt, dcomplex *work, int *lwork, int *info);
int zcopy_(int *n, dcomplex *zx, int *incx, dcomplex *zy, int *incy);

void zgtsv_(INTEGER *n, INTEGER *nrhs, dcomplex *dl, 
	    dcomplex *d__, dcomplex *du, dcomplex *b, INTEGER *ldb,
	    INTEGER *info);

/* added by mari 08.12.2014 */
void zgesv_(INTEGER *n, INTEGER *nrhs, dcomplex *a, INTEGER *lda,
	    INTEGER *ipiv, dcomplex *b, INTEGER *ldb, INTEGER *info);

void dpotrf_(char *uplo, INTEGER *n, double *a, INTEGER *lda, INTEGER *info);
void dpotri_(char *uplo, INTEGER *n, double *a, INTEGER *lda, INTEGER *info);
void dggevx_(char *balanc, char *jobvl, char *jobvr, char *sense,
             INTEGER *n, double *a, INTEGER *lda, double *b, 
	     INTEGER *ldb, double *alphar, double *alphai, double *beta,
             double *vl, INTEGER *ldvl, double *vr, INTEGER *ldvr, 
	     INTEGER *ilo, INTEGER *ihi, double *lscale, double *rscale, 
	     double *abnrm, double *bbnrm, double *rconde, double *rcondv,
             double *work, INTEGER *lwork, INTEGER *iwork, logical *bwork,
             INTEGER *info);

void dsytrd_(char *uplo, INTEGER *n, double *a, INTEGER *lda, double *d__, double *e, 
             double *tau, double *work, INTEGER *lwork, INTEGER *info);


typedef enum {CblasRowMajor=101, CblasColMajor=102} CBLAS_ORDER;
typedef enum {CblasNoTrans=111, CblasTrans=112, CblasConjTrans=113} CBLAS_TRANSPOSE;

/* MKL integer types for LP64 and ILP64 */
#if (!defined(__INTEL_COMPILER)) & defined(_MSC_VER)
    #define MKL_INT64 __int64
    #define MKL_UINT64 unsigned __int64
#else
    #define MKL_INT64 long long int
    #define MKL_UINT64 unsigned long long int
#endif
 
#ifdef MKL_ILP64

/* MKL ILP64 integer types */
#ifndef MKL_INT
    #define MKL_INT MKL_INT64
#endif
#ifndef MKL_UINT
    #define MKL_UINT MKL_UINT64
#endif
#define MKL_LONG MKL_INT64

#else

/* MKL LP64 integer types */
#ifndef MKL_INT
    #define MKL_INT int
#endif
#ifndef MKL_UINT
    #define MKL_UINT unsigned int
#endif
#define MKL_LONG long int

#endif

void cblas_dgemm(const  CBLAS_ORDER Order, const  CBLAS_TRANSPOSE TransA,
                 const  CBLAS_TRANSPOSE TransB, const MKL_INT M, const MKL_INT N,
                 const MKL_INT K, const double alpha, const double *A,
                 const MKL_INT lda, const double *B, const MKL_INT ldb,
                 const double beta, double *C, const MKL_INT ldc);


void DGEMM(const char *transa, const char *transb, const MKL_INT *m, const MKL_INT *n, const MKL_INT *k,
           const double *alpha, const double *a, const MKL_INT *lda, const double *b, const MKL_INT *ldb,
           const double *beta, double *c, const MKL_INT *ldc);

void dgemm(const char *transa, const char *transb, const MKL_INT *m, const MKL_INT *n, const MKL_INT *k,
           const double *alpha, const double *a, const MKL_INT *lda, const double *b, const MKL_INT *ldb,
           const double *beta, double *c, const MKL_INT *ldc);
#endif
