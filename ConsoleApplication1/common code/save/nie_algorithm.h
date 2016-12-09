/*******************************

	一般算法的集合
	2008/4/1 by Yongwei Nie

********************************/


#pragma once
#include "nie_type_def.h"

namespace NIE_ALGORITHM
{
class Algorithm
{
public:
	//插值
	static void interp_3d_double_u(MAT_DOU_3D & uf,MAT_C_DOU_3D &uc);
	static void interp_3d_double_d(MAT_DOU_3D & uf,MAT_C_DOU_3D &uc);
	static void interp_3d_p3t_u(MAT_3D_POINT_3D_T &uf,MAT_C_3D_POINT_3D_T &uc);
	static void interp_3d_p3t_d(MAT_3D_POINT_3D_T &uf,MAT_C_3D_POINT_3D_T &uc);
	//copy
	template <class T> 
	static void copymat3d(Matrix3D<T> & lhs,const Matrix3D<T> & rhs);
};

}


void NIE_ALGORITHM::Algorithm::interp_3d_double_u(MAT_DOU_3D &uf, MAT_C_DOU_3D &uc)
{
	int nnc = uc.dim1();
	int mmc = uc.dim2();
	int kkc = uc.dim3();
	
	int nnf = nnc*2-1;
	int mmf = mmc*2-1;
	int kkf = kkc*2-1;

	for(int z = 0;z<nnc;z++)
		for(int y = 0;y<mmc;y++)
			for(int x = 0;x<kkc;x++)
				uf[z*2][y*2][x*2] = uc[z][y][x];

	for(int zf = 0;zf<nnf;zf+=2)
		for(int xf = 0;xf<kkf;xf+=2)
			for(int yf = 1;yf<mmf-1;yf+=2)
				uf[zf][yf][xf] = 0.5*(uf[zf][yf-1][xf] + uf[zf][yf+1][xf]);

	for(int zf = 0;zf<nnf;zf+=2)
		for(int xf = 1;xf<kkf-1;xf+=2)
			for(int yf = 0;yf<mmf;yf++)
				uf[zf][yf][xf] = 0.5*(uf[zf][yf][xf-1] + uf[zf][yf][xf+1]);
	
	for(int zf = 1;zf<nnf-1;zf+=2)
		for(int yf = 0;yf<mmf;yf+=2)
			for(int xf = 0;xf<kkf;xf+=2)
				uf[zf][yf][xf] = 0.5*(uf[zf-1][yf][xf] + uf[zf+1][yf][xf]);

	for(int zf = 1;zf<nnf-1;zf+=2)
		for(int xf = 0;xf<kkf;xf+=2)
			for(int yf = 1;yf<mmf-1;yf+=2)
				uf[zf][yf][xf] = 0.5*(uf[zf][yf-1][xf] + uf[zf][yf+1][xf]);

	for(int zf = 1;zf<nnf-1;zf+=2)
		for(int xf = 1;xf<kkf-1;xf+=2)
			for(int yf = 0;yf<mmf;yf++)
				uf[zf][yf][xf] = 0.5*(uf[zf][yf][xf-1] + uf[zf][yf][xf+1]);
}

void NIE_ALGORITHM::Algorithm::interp_3d_double_d(MAT_DOU_3D &uf, MAT_C_DOU_3D &uc)
{
	int nnc = uc.dim1();
	int mmc = uc.dim2();
	int kkc = uc.dim3();
	
	int nnf = (nnc-1)/2 + 1;
	int mmf = (mmc-1)/2 + 1;
	int kkf = (kkc-1)/2 + 1;

	for(int i = 0;i<nnf;i++)
		for(int j = 0;j<mmf;j++)
			for(int k = 0;k<kkf;k++)
				uf[i][j][k] = uc[i*2][j*2][k*2];
}

void NIE_ALGORITHM::Algorithm::interp_3d_p3t_u(MAT_3D_POINT_3D_T &uf, MAT_C_3D_POINT_3D_T &uc)
{
	int nnf = uf.dim1();
	int mmf = uf.dim2();
	int kkf = uf.dim3();
	int nnc = uc.dim1();
	int mmc = uc.dim2();
	int kkc = uc.dim3();

	MAT_DOU_3D d1(0.0,nnf,mmf,kkf);
	MAT_DOU_3D d2(0.0,nnc,mmc,kkc);

	for(int p = 0;p<3;p++)
	{
		for(int i = 0;i<nnc;i++)
			for(int j = 0;j<mmc;j++)
				for(int k = 0;k<kkc;k++)
					d2[i][j][k] = uc[i][j][k].val[p];
		interp_3d_double_u(d1,d2);
		for(int i = 0;i<nnf;i++)
			for(int j = 0;j<mmf;j++)
				for(int k = 0;k<kkf;k++)
					uf[i][j][k].val[p] = d1[i][j][k];
	}
}

void NIE_ALGORITHM::Algorithm::interp_3d_p3t_d(MAT_3D_POINT_3D_T &uf, MAT_C_3D_POINT_3D_T &uc)
{
	int nnf = uf.dim1();
	int mmf = uf.dim2();
	int kkf = uf.dim3();
	int nnc = uc.dim1();
	int mmc = uc.dim2();
	int kkc = uc.dim3();

	MAT_DOU_3D d1(0.0,nnf,mmf,kkf);
	MAT_DOU_3D d2(0.0,nnc,mmc,kkc);

	for(int p = 0;p<3;p++)
	{
		for(int i = 0;i<nnc;i++)
			for(int j = 0;j<mmc;j++)
				for(int k = 0;k<kkc;k++)
					d2[i][j][k] = uc[i][j][k].val[p];
		interp_3d_double_d(d1,d2);
		for(int i = 0;i<nnf;i++)
			for(int j = 0;j<mmf;j++)
				for(int k = 0;k<kkf;k++)
					uf[i][j][k].val[p] = d1[i][j][k];
	}
}


template <class T>
void NIE_ALGORITHM::Algorithm::copymat3d(Matrix3D<T> &lhs, const Matrix3D<T> &rhs)
{
	int nn = lhs.dim1();
	int mm = lhs.dim2();
	int kk = lhs.dim3();
	for(int i = 0;i<nn;i++)
		for(int j = 0;j<mm;j++)
			for(int k = 0;k<kk;k++)
				lhs[i][j][k] = rhs[i][j][k];
}