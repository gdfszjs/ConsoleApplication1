//////////////////////////////////////////////////////////////////////
// ***************************************************************
//  PBCG          version:  v1.0          date: 09/20/2010
//  -------------------------------------------------------------
//  Name     : PBCG.H
//  Purpose  : 预条件双共轭梯度法求解Ax=b
//  Author   : Chen Wenfei
//  -------------------------------------------------------------
//  Copyright (C) 2010 - All Rights Reserved.
// ***************************************************************
//  调用方法:
//     const int       ITOL = 1, ITMAX = 75;
//     const double    TOL  = 1.0e-9;
//     double          err  = 0.0;
//     int             iter = 0;
//     PBCG pbcgSolver(sa_d, ija_d);
//     pbcgSolver.linbcg(x, b, NP, iter, err, ITOL, TOL, ITMAX);
// ***************************************************************
#pragma once
#ifndef _INC_PBCG_H
#define _INC_PBCG_H

//////////////////////////////////////////////////////////////////////////
//macro define
//------------------------------------------------------------------------
//收敛方式(1,2,3,4)
#ifndef DEF_ITOL_DEFAULT
#define DEF_ITOL_DEFAULT 1
#endif
//------------------------------------------------------------------------
//迭代精度
#ifndef DEF_TOL_DEFAULT
#define DEF_TOL_DEFAULT 1.0e-9
#endif
//------------------------------------------------------------------------
//最大迭代次数
#ifndef DEF_ITMAX_DEFAULT
#define DEF_ITMAX_DEFAULT 75
#endif

//////////////////////////////////////////////////////////////////////////
// MACRO that forbids to use the copy ctor & assign statement
//------------------------------------------------------------------------
#ifndef DISALLOW_COPY_AND_ASSIGN
#define DISALLOW_COPY_AND_ASSIGN(ClassName)   \
    ClassName(const ClassName&);              \
    ClassName& operator=(const ClassName&)
#endif

//////////////////////////////////////////////////////////////////////////
class PBCG
{
public://ctor & dtor
    PBCG(const double *sa, const unsigned long *ija);
    ~PBCG(void);

public://interface functions
    void linbcg(double x[],                     //解向量(out)
        const double b[],                       //方程右式向量(in)
        const unsigned long n,                  //方阵维数
        int &iter,                              //实际迭代次数
        double &err,                            //实际迭代精度
        int itol = DEF_ITOL_DEFAULT,            //收敛方式(1,2,3,4)
        double tol = DEF_TOL_DEFAULT,           //所需迭代精度
        int itmax = DEF_ITMAX_DEFAULT) const;   //所需迭代次数

public:
    void sprstx(double b[], const double x[], const unsigned long n) const;
    void sprsax(double b[], const double x[], const unsigned long n) const;

private:
    void atimes(double b[], const double x[], const unsigned long n, const int itrnsp) const;
    void asolve(double x[], const double b[], const unsigned long n, const int itrnsp) const;
    double snrm(const double sx[], const unsigned long n, const int itol = 1) const;

private://行索引稀疏矩阵
    const double          *_sa;
    const unsigned long   *_ija;

private://forbid to use the copy ctor & assign statement
    DISALLOW_COPY_AND_ASSIGN(PBCG);
};

//////////////////////////////////////////////////////////////////////////
//global functions
//------------------------------------------------------------------------
void sprsin(double *const sa,
            unsigned long *const ija,       //行索引稀疏矩阵(out)
            const unsigned long nmax,       //向量长度一般为2*n*n+1
            const double *const a,          //稀疏系数方阵A(in)
            const unsigned long n,          //方阵的维数n
            const double thresh = 0.0);     //阈值
//------------------------------------------------------------------------
void sprsout(double *const a,               //稀疏系数方阵A(out)
             const unsigned long n,         //方阵维数n
             const double *const sa,
             const unsigned long *const ija);//行索引稀疏矩阵(in)
//------------------------------------------------------------------------
void CreateLapMat(double *const sa,
                  unsigned long *const ija, //行索引稀疏矩阵(out)
                  const int iWidth,         //源图像的宽度
                  const int iHeight);       //源图像的高度

#endif // !_INC_PBCG_H
