/************************************

	class EM3D 用于视频的纹理合成

	2008/4/1 by Yongwei Nie

	完成于2008/4/3

************************************/
#pragma once
#include "nie_algorithm.h"
#include <vector>
#include <cmath>
//#include <iostream>
//using namespace std;
#include "anncell.h"

using namespace NIE_ALGORITHM;

typedef std::vector<MAT_3D_POINT_3D_T *> VP_MAT_3D_POINT_3D_T;   //储存MAT_3D_POINT_3D_T类型指针的vector
typedef std::vector<int> V_INT;							        //储存整形的vector
typedef std::vector<V_INT> V_V_INT;								//储存V_INT的vector
typedef const std::vector<V_INT> C_V_V_INT;

namespace NIE_EM
{
class EM3D
{
public:
	//视屏EM纹理合成，由输入matsrc合成matres（输出）
	static void em(MAT_C_3D_POINT_3D_T & matsrc,MAT_3D_POINT_3D_T & matres,int mr,int hh,C_V_V_INT&vvi);
	
	//创建ANN树
	static ANNCell*ann(MAT_C_3D_POINT_3D_T & mat,int siz,int hh);
	
	//m-step 新的zp
	static void zp(MAT_INT_3D & pi,MAT_DOU_3D & pd,MAT_C_3D_POINT_3D_T & mat,ANNCell * pc,int siz,int hh,
		bool _r = false);
	
	//能量
	static double dist(MAT_C_DOU_3D & lhs,int siz,int hh);

	//分配内存
	static void zpid(MAT_INT_3D *&mati,MAT_DOU_3D *& matd,MAT_C_3D_POINT_3D_T& mat,int siz,int hh);

	//e-step
	static void smooth(MAT_3D_POINT_3D_T & mat,MAT_INT_3D & mati,MAT_DOU_3D & matd,int siz,int hh,ANNCell * pc);

	//反复EM
	static void core(MAT_3D_POINT_3D_T & mat,MAT_INT_3D & mati,MAT_DOU_3D & matd,int siz,int hh,ANNCell * pc);

	//static void display(MAT_C_INT_3D & mat);
};

}

//void NIE_EM::EM3D::display(MAT_C_INT_3D & mat)
//{
//	int nn = mat.dim1();
//	int mm = mat.dim2();
//	int kk = mat.dim3();
//	for(int i = 0;i<nn;i++)
//		for(int j = 0;j<mm;j++)
//			for(int k = 0;k<kk;k++)
//				cout<<mat[i][j][k]<<",";
//	cout<<"********************************"<<endl;
//}


ANNCell * NIE_EM::EM3D::ann(MAT_C_3D_POINT_3D_T & mat,int siz,int hh)
{
	
	int nn,mm,kk,nw,mw,kw;
	nn = mat.dim1();mm = mat.dim2();kk = mat.dim3();
	int w1 = hh/4;
	int w2 = siz/4;
	nw = nn+1-hh;mw=mm+1-siz;kw=kk+1-siz;
	ANNCell * pc = new ANNCell(siz*siz*hh*3,nw*mw*kw);
	for(int iig = 0;iig < nw;iig++)
		for(int jjg = 0;jjg < mw;jjg++)
			for(int kkg = 0;kkg <kw;kkg++)
				for(int iil = iig;iil<iig+hh;iil++)
					for(int jjl = jjg;jjl<jjg+siz;jjl++)
						for(int kkl = kkg;kkl<kkg+siz;kkl++)
							for(int p = 0;p<3;p++)
								pc->readInCoordinates(mat[iil][jjl][kkl].val[p]);
	pc->buildKd_tree();
	return pc;
}

void NIE_EM::EM3D::zp(MAT_INT_3D & pi,MAT_DOU_3D & pd,MAT_C_3D_POINT_3D_T &mat,ANNCell * pc,int siz,int hh,
							  bool _r)
{
	int	nnw = pi.dim1();
	int mmw = pi.dim2();
	int kkw = pi.dim3();
	int w1 = siz/4;
	int w2 = hh/4;
	
	if(_r)
	{
		int total = pc->getpn();
		for(int i = 0;i<nnw;i++)
			for(int j = 0;j<mmw;j++)
				for(int k = 0;k<kkw;k++)
				{
					pi[i][j][k] = rand()%total;
					pd[i][j][k] = 1.0;
				}
	}
	else
	{
		for(int i = 0;i<nnw;i++)
			for(int j = 0;j<mmw;j++)
				for(int k = 0;k<kkw;k++)
				{
					double * _poi = new double[siz*siz*hh*3];
					int nb,mb,kb;
					nb = i*w2;mb=j*w1;kb=k*w1;
					for(int ii = 0;ii<hh;ii++)
						for(int jj = 0;jj<siz;jj++)
							for(int kk = 0;kk<siz;kk++)
								for(int p = 0;p<3;p++)
									_poi[ii*siz*siz*3+jj*siz*3+kk*3+p] = mat[nb+ii][mb+jj][kb+kk].val[p];
					double * di = new double[1];
					di[0] = 0.0;
					pi[i][j][k]= pc->findNearest(_poi,di);
					pd[i][j][k] = di[0];
					delete [] di;
					delete [] _poi;
				}
	}
}

double NIE_EM::EM3D::dist(MAT_C_DOU_3D &lhs,int siz,int hh)
{
	int nn,mm,kk;
	nn = lhs.dim1();
	mm = lhs.dim2();
	kk = lhs.dim3();
	double rs = 0.0;
	for(int i = 0;i<nn;i++)
		for(int j = 0;j<mm;j++)
			for(int k = 0;k<kk;k++)
				rs += lhs[i][j][k];
	return (rs/(nn*mm*kk*siz*siz*hh*3));   //normalize
}

void NIE_EM::EM3D::zpid(MAT_INT_3D *&mati,MAT_DOU_3D *& matd,MAT_C_3D_POINT_3D_T& mat,int siz,int hh)
{
	int w1 = siz/4;
	int w2 = hh /4;
	int nn = mat.dim1();
	int mm = mat.dim2();
	int kk = mat.dim3();
	int nnw = (nn-hh)/w2+1;
	int mmw = (mm-siz)/w1+1;
	int kkw = (kk-siz)/w1+1;
	mati = new MAT_INT_3D(0,nnw,mmw,kkw);
	matd = new MAT_DOU_3D(0.0,nnw,mmw,kkw);
}

void NIE_EM::EM3D::smooth(MAT_3D_POINT_3D_T &mat,MAT_INT_3D & mati,MAT_DOU_3D & matd,int siz,int hh, ANNCell *pc)
{
	ANNpointArray apa = pc->getPointArray();
	int in = mati.dim1();
	int im = mati.dim2();
	int ik = mati.dim3();
	int w1 = siz/4;
	int w2 = hh/4;
	int nn = (in-1)*w2+hh;
	int mm = (im-1)*w1+siz;
	int kk = (ik-1)*w1+siz;
	for(int i = 0;i<nn;i++)
		for(int j = 0;j<mm;j++)
			for(int k = 0;k<kk;k++)
			{
				int ii,jj,kk,_ii,_jj,_kk;
				ii = i / w2;
				jj = j / w1;
				kk = k / w1;
				_ii = i % w2;
				_jj = j % w1;
				_kk = k % w1;
				double total_i = 0.0;
				POINT_3D_T total_d(0.0);
				for(int iiw = ii-3;iiw<=ii;iiw++)
					for(int jjw = jj-3;jjw<=jj;jjw++)
						for(int kkw = kk-3;kkw<=kk;kkw++)
							if(iiw>=0&&jjw>=0&&kkw>=0&&iiw<in&&jjw<im&&kkw<ik)
							{
								for(int p = 0;p<3;p++)
									total_d.val[p] += ((apa[mati[iiw][jjw][kkw]])[((ii-iiw)*w2+_ii)*siz*siz*3 + ((jj-jjw)
									*w1+_jj)*siz*3 + ((kk-kkw)*w1+_kk)*3 + p]) * matd[iiw][jjw][kkw];
								total_i += matd[iiw][jjw][kkw];
							}
				for(int p = 0;p<3;p++)
				{
					total_d.val[p] /= total_i;
					mat[i][j][k].val[p] = total_d.val[p];
				}

			}
}


void NIE_EM::EM3D::core(MAT_3D_POINT_3D_T &mat,MAT_INT_3D & mati,MAT_DOU_3D & matd,int siz,int hh,ANNCell *pc)
{
	double globaldist = 0.0;
	while(1)
	{
		smooth(mat,mati,matd,siz,hh,pc);
		zp(mati,matd,mat,pc,siz,hh);
		double _d = dist(matd,siz,hh);
 		if(fabs(_d-globaldist) < 0.01)
			break;
		else globaldist = _d;
	}
}

void NIE_EM::EM3D::em(MAT_C_3D_POINT_3D_T &matsrc, MAT_3D_POINT_3D_T &matres,int mr,int hh,C_V_V_INT&vvi)
{
	int srcn = matsrc.dim1();
	int srcm = matsrc.dim2();
	int srck = matsrc.dim3();
	int resn = matres.dim1();
	int resm = matres.dim2();
	int resk = matres.dim3();

	VP_MAT_3D_POINT_3D_T v3t(mr);
	VP_MAT_3D_POINT_3D_T v3t_r(mr);
	v3t[mr-1] = new MAT_3D_POINT_3D_T(POINT_3D_T(0.0),srcn,srcm,srck);
	v3t_r[mr-1] = new MAT_3D_POINT_3D_T(POINT_3D_T(0.0),resn,resm,resk);
	Algorithm::copymat3d(*v3t[mr-1],matsrc);
	Algorithm::copymat3d(*v3t_r[mr-1],matres);
	int tmr = mr-1;
	int tn = srcn ,tm = srcm ,tk = srck;
	int rn = resn ,rm = resm ,rk = resk;
	while(--tmr >= 0)
	{	
		tn = (tn-1)/2 + 1;tm = (tm-1)/2 + 1;tk = (tk-1)/2+1;
		rn = (rn-1)/2 + 1;rm = (rm-1)/2 + 1;rk = (rk-1)/2+1;
		v3t[tmr] = new MAT_3D_POINT_3D_T(POINT_3D_T(0.0),tn,tm,tk);
		v3t_r[tmr] = new MAT_3D_POINT_3D_T(POINT_3D_T(0.0),rn,rm,rk);
		Algorithm::interp_3d_p3t_d(*v3t[tmr],*v3t[tmr+1]);
		Algorithm::interp_3d_p3t_d(*v3t_r[tmr],*v3t_r[tmr+1]);
	}
	for(int i = 0;i<mr;i++)
	{
		int sca = (int)vvi[i].size();
		for(int j = 0;j<sca;j++)
		{
			int siz = vvi[i][j];
			ANNCell * pcell = ann(*v3t[i],siz,hh);
			MAT_INT_3D *mati;
			MAT_DOU_3D *matd;
			zpid(mati,matd,*v3t_r[i],siz,hh);
			if(i==0&&j==0) zp(*mati,*matd,*v3t_r[i],pcell,siz,hh,true);
			core(*v3t_r[i],*mati,*matd,siz,hh,pcell);
			delete pcell;delete mati;delete matd;
			pcell = NULL;mati = NULL;matd = NULL;
		}
		Algorithm::interp_3d_p3t_u(*v3t_r[i+1],*v3t_r[i]);
	}

	Algorithm::copymat3d(matres,*v3t_r[mr-1]);
	for(int i = 0;i<mr;i++)
	{
		delete v3t[i];
		delete v3t_r[i];
	}
	v3t.clear();
	v3t_r.clear();
}
