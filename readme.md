# OpenMX patch for Intel OneAPI (icx) compiler

## はじめに

IntelコンパイラがOneAPIとして再編され、メインのC/C++, FortranコンパイラがLLVMベースに変わってしばらく経ちました。従来のicc等が使えなくなり、コンパイルが通らなくなったアプリに関しての悲鳴がちらほら聞こえてきます。

そこで、色々とお世話になった密度汎関数理論計算コードOpenMX[^1]について、OneAPI(icx,ifx)を使ってコンパイルする方法を纏めることにしました。

[^1]: OpenMX本家：[https://www.openmx-square.org/](https://www.openmx-square.org/)

パッチも配布しますので、後半をお読みください。

## 想定環境


- OS: Linux
- OneAPI: version 2025.1 (たぶんそれ以前でも大丈夫です)
- IntelMKL：BLAS, LAPACK, ScaLAPACK用に用います
- IntelMPI：MPI用に使います。

## 修正箇所一覧

### LapackおよびELPA用ヘッダー追加

次のコードをファイルの冒頭に追加。
```
#include "lapack_prototypes.h"
```

対象ファイル

- `Cluster_DFT_Col.c`
- `Cluster_DFT_NonCol.c`
- `Cluster_DFT_Optical_ScaLAPACK.c`
- `DFT.c`
- `Divide_Conquer_LNO.c`
- `EGAC_DFT.c`
- `Eigen_PHH.c`
- `TRAN_DFT.c`
- `TRAN_DFT_NC.c`
- `Mixing_H.c`
- `Mixing_V.c`
- 

### 文字列操作ヘッダーの追加

次のコードをファイルの冒頭に追加。通常使う`string.h`とは異なるので注意
```
#include <strings.h>
```

対象ファイル

- `Cluster_DFT_Col.c`
- `Cluster_DFT_LNO.c`
- `Cluster_DFT_NonCol.c`
- `Cluster_DFT_ON2.c`
- `Cluster_DFT_Optical_ScaLAPACK.c`
- `DFT.c`
- `Generate_Wannier.c`
- `Inputtools.c`
- `MD_pac.c`
- `Mulliken_Charge.c`
- `Occupation_Number_LDA_U.c`

  

### 文字列操作ヘッダーの追加（２）

次のコードをファイルの冒頭に追加。今度は通常使う`string.h`を追加
```
#include <string.h>
```
対象ファイル

- `DFT.c`
- `Force.c`
- `TRAN_Channel_Output.c`

### 関数定義の追加

次のコードをファイルのヘッダーインクルード後あたりに追加
```
void Eigen_lapack_d(double **a, double *ko, int n0, int EVmax);
```

対象ファイル

- `Divide_Conquer_LNO.c`

### 関数定義の追加（２）

次のコードをファイルのヘッダーインクルード後あたりに追加
```
int Lapack_LU_Zinverse(int , dcomplex *);
```

対象ファイル

- `EGAC_DFT.c`

### Lapack関数の宣言をMKLに置き換え

ファイルの中身の冒頭と末尾に、以下のように`#if`～`#endif`行の行を追加。またELPA関連の関数宣言も同時に追加。
```
//(ファイル先頭)
//(追加開始)
#if defined(__INTEL_LLVM_COMPILER ) || defined(__INTEL_COMPILER )
#include "f77_name.h"
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
//(追加終了)

(...元のファイルの中身...)

//(追加開始)
#endif  //ファイル末尾
//(追加終了)
```

対象ファイル

- `lapack_prototypes.h`

### 関数宣言・定義の修正

関数`myselect`の宣言部と定義部の戻り値と引数の型を修正
```
//(ファイル冒頭)
//static logical myselect(double wr, double wi);   //before//
int myselect(const double *, const double *);      //after//

//(ファイル後半)
//static logical myselect(double wr, double wi)   //before//
int myselect(const double *, const double *)      //after//
```

対象ファイル

- `LNO.c`

  

### ompヘッダーのインクルード

次のコードをファイルの冒頭に追加。
```
#include <omp.h>
```

対象ファイル

- `Mixing_V.c`

### unistdヘッダーのインクルード

次のコードをファイルの冒頭に追加。
```
#include <unistd.h>
```

対象ファイル

- `OutData_Binary.c`
- `TRAN_Channel_Output.c`

### INTEGER型のエイリアス化

次のコードをファイルの冒頭に追加。
```
#define INTEGER int
```

対象ファイル

- `TRAN_Calc_CurrentDensity.c`

### fsync()関数のコメントアウト

以下の二行をコメントアウト
```
           //fd = fileno(fp);
           //fsync(fd);
```

対象ファイル

- `TRAN_CDen_Main.c`

※古い関数のため。コードを読む限りはコメントアウトしても特に問題にならないと思われる。

  

### makefile修正

makefileに記載されているコンパイルオプション等は次のようにした
```
#MKLROOT = /opt/intel/mkl  ##環境変数にあればここでは消してしまう
CC = mpiicx -Wl,--allow-multiple-definition -Df77 -D_DEFAULT_SOURCE \
  -std=c11 -restrict -O3 -flto -qmkl -xCORE-AVX512 -qopt-zmm-usage=high \
  -w2 -wd3180 -wd161 -qopt-report=max -qopt-report-file=opt-report.txt \
  -I${MKLROOT}/include -I${MKLROOT}/include/fftw \
  -Wimplicit-function-declaration
FC = mpiifx -assume nounderscore \
  -O3 -qmkl -xCORE-AVX512 -qopt-zmm-usage=high -I${MKLROOT}/include
LIB= -L${MKLROOT}/lib/intel64 -lmkl_scalapack_lp64 -lmkl_intel_lp64 \
  -lmkl_intel_thread -lmkl_core -lmkl_blacs_intelmpi_lp64 -lpthread \
  -lifcore
```
また、MKLROOTに関しては環境変数でセットされていればそちらを使うようにしてmakefileからコメントアウトして良い。

## パッチ

上記の微修正を施したパッチとしてgithubに配置しました。ライセンスは本家を踏襲してGPL3になります。

[https://github.com/atsushi-m-ito/openmx-patch-oneapi/](https://github.com/atsushi-m-ito/openmx-patch-oneapi/)

ここのReleaseからダウンロードできます。

### パッチの適用方法

このパッチの対象は、OpenMX version 3.9.9 になります。事前に本家のソースコード(ver3.9)に、本家のパッチ(3.9.9)を当てておいてください。

その上で以下のコマンドを参考に、今回のパッチ(ファイル名は仮に`v1.tar.gz`)を当ててください。
```
cd openmx3.9/source/
tar xvfz v1.tar.gz --strip-components 1

#以下コンパイル
make -j install
```

  
