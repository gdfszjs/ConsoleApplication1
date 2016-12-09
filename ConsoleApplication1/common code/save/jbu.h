/************

	Joint Bilateral Filter
	
	class Jbu2D对二维彩色图像进行JBU插值，没有对插值方法进行加速
	2008/3/27 by Yongwei Nie

	class Jbu3D对三维彩色视频进行JBU插值
	2008/3/29 by Yongwei Nie

************/

#ifndef _JBU_H_
#define _JBU_H_
#pragma once 

#include <cmath>
#include "nie_type.h"
#include "datastructure.h"

using namespace NIE_TYPE;
using namespace NIE_DATASTRUCTURE;

typedef nie_point_3d_t POINT_3D_T;

typedef Matrix2D<POINT_3D_T> MAT_2D_POINT_3D_T;
typedef const Matrix2D<POINT_3D_T> MAT_C_2D_POINT_3D_T;

typedef Matrix3D<POINT_3D_T> MAT_3D_POINT_3D_T;
typedef const Matrix3D<POINT_3D_T> MAT_C_3D_POINT_3D_T;

namespace NIE_JBU
{

class Jbu2D
{
public:

/**********************************
hig  是 high resolution image I~
sol  是 low  resolution solution S
res  是 求得的结果
r    是 模板的边长
sigd 是 domain参数
sigr 是 range 参数
***********************************/

	void jbu(MAT_C_2D_POINT_3D_T & hig,MAT_C_2D_POINT_3D_T & sol,MAT_2D_POINT_3D_T & res,const int r,
		double sigd,double sigr);

};

class Jbu3D
{
public:
/*********************************
hig  是 high resolution video I~
sol  是 low  resolution solution S
res  是 求得的结果
r    是 模板的边长
sigd 是 domain参数
sigr 是 range 参数
***********************************/

	void jbu(MAT_C_3D_POINT_3D_T & hig,MAT_C_3D_POINT_3D_T & sol,MAT_3D_POINT_3D_T &res,const int r,
		double sigd,double sigr);
};

}

/*

	#########################################
	#########################################
	#########################################
	#####			实     现		    #####
	#########################################
	#########################################
	#########################################
*/


void NIE_JBU::Jbu2D::jbu(MAT_C_2D_POINT_3D_T &hig, MAT_C_2D_POINT_3D_T &sol, MAT_2D_POINT_3D_T &res, 
						 const int r, double sigd, double sigr)
{
	double sigd2 = sigd * sigd;
	double sigr2 = sigr * sigr;
	int    _r    = r / 2;
	//得到hig的行数
	int hr = hig.getnrows();
	//得到hig的列数
	int hc = hig.getncols();
	//得到sol的行数
	int sr = sol.getnrows();
	//得到sol的列数
	int sc = sol.getncols();
	//res的行数等于hig的行数
	int rr = res.getnrows();
	//res的列数等于hig的列数
	int rc = res.getncols();

	//res和sol之间的长和宽的比例
	double rscale = (double)sr/rr;
	double cscale = (double)sc/rc;
	//hig和sol之间的长和宽的比例
	double _rscale = (double)sr/hr;
	double _cscale = (double)sc/hc;

	int irb,icb;
	
	int ybegin = 0,yend = rr,xbegin = 0,xend = rc;
	while(ybegin*rscale<=_r-1)
		ybegin++;
	while(xbegin*cscale<=_r-1)
		xbegin++;
	while(yend*rscale > sr - _r)
		yend--;
	while(xend*cscale > sc - _r)
		xend--;

	for(int y = ybegin ;y< yend;y++)
		for(int x = xbegin;x< xend;x++)
		{
			POINT_3D_T p3tup(0.0),p3tdown(0.0),p3tres(0.0);
			double dpr = y*rscale,dpc = x*cscale;
			if( dpr == (int)dpr ) irb = -_r;
			else irb = -_r + 1;
			if( dpc == (int)dpc ) icb = -_r;
			else icb = -_r + 1;
			for(int j = (int)dpr + irb;j<= (int)dpr + _r;j++)
				for(int i = (int)dpc + icb;i<= (int)dpc + _r;i++)
				{
					double dpq2 = (dpr-j)*(dpr-j) + (dpc-i)*(dpc-i);
					POINT_3D_T rangediffer2(0.0);
					for(int k = 0;k<3;k++)
						rangediffer2.val[k] = (hig[(int)(dpr/_rscale)][(int)(dpc/_cscale)].val[k] - hig[(int)(j/_rscale)][(int)(i/_cscale)].val[k]) *
						(hig[(int)(dpr/_rscale)][(int)(dpc/_cscale)].val[k] - hig[(int)(j/_rscale)][(int)(i/_cscale)].val[k]);
					for(int k = 0;k<3;k++)
					{
						double tmp = exp(-0.5 * ( dpq2/sigd2 + rangediffer2.val[k]/sigr2));
						p3tup.val[k] += sol[j][i].val[k] * tmp;
						p3tdown.val[k] += tmp;
					}
				}
			
			for(int k = 0;k<3;k++)
				p3tres.val[k] = p3tup.val[k] / p3tdown.val[k];
			res[y][x] = p3tres;
		}
}


void NIE_JBU::Jbu3D::jbu(MAT_C_3D_POINT_3D_T &hig, MAT_C_3D_POINT_3D_T &sol, MAT_3D_POINT_3D_T &res, const int r,
						 double sigd, double sigr)
{
	double sigd2 = sigd*sigd;
	double sigr2 = sigr*sigr;
	int _r = r/2;
	int hnn = hig.dim1();
	int hmm = hig.dim2();
	int hkk = hig.dim3();
	int snn = sol.dim1();
	int smm = sol.dim2();
	int skk = sol.dim3();
	int rnn = res.dim1();
	int rmm = res.dim2();
	int rkk = res.dim3();

	//sol与res之间的长和宽的比例
	double nscale = (double)snn/rnn;
	double mscale = (double)smm/rmm;
	double kscale = (double)skk/rkk;
	//sol与hig之间的长和宽的比例
	double _nscale = (double)snn/hnn;
	double _mscale = (double)smm/hmm;
	double _kscale = (double)skk/hkk;
	
	int nbegin(0),nend(rnn),mbegin(0),mend(rmm),kbegin(0),kend(rkk);
	while(nbegin * nscale <= _r - 1) nbegin++;
	while(nend * nscale > snn - _r)  nend--;
	while(mbegin * mscale <= _r - 1) mbegin++;
	while(mend * mscale > smm - _r)  mend--;
	while(kbegin * kscale <= _r - 1) kbegin++;
	while(kend * kscale > skk - _r)  kend--;

	int inb,imb,ikb;

	for(int n = nbegin;n<nend;n++)
		for(int m = mbegin;m<mend;m++)
			for(int k = kbegin;k<kend;k++)
			{
				POINT_3D_T p3tup(0.0),p3tdown(0.0),p3tres(0.0);
				double dpn = n * nscale,dpm = m * mscale,dpk = k * kscale;
				if(dpn == (int)dpn) inb = -_r;
				else inb = -_r+1;
				if(dpm == (int)dpm) imb = -_r;
				else imb = -_r+1;
				if(dpk == (int)dpk) ikb = -_r;
				else ikb = -_r+1;
				for(int i = (int)dpn + inb;i<(int)dpn + _r;i++)
					for(int j = (int)dpm + imb;j<(int)dpm + _r;j++)
						for(int _k = (int)dpk + ikb;_k<(int)dpk + _r;_k++)
						{
							double dpq2 = (dpn - i)*(dpn-i) + (dpm-j)*(dpm-j) + (dpk-_k)*(dpk-_k);
							POINT_3D_T rangediffer2(0.0);
							for(int kk = 0;kk<3;kk++)
								rangediffer2.val[kk] = (hig[(int)(dpn/_nscale)][(int)(dpm/_mscale)][(int)(dpk/_kscale)].val[kk] - 
								hig[(int)(i/_nscale)][(int)(j/_mscale)][(int)(_k/_kscale)].val[kk]) * (hig[(int)(dpn/_nscale)][(int)(dpm/_mscale)][(int)(dpk/_kscale)].val[kk] - 
								hig[(int)(i/_nscale)][(int)(j/_mscale)][(int)(_k/_kscale)].val[kk]);
							for(int kk = 0;kk<3;kk++)
							{
								double tmp = exp(-0.5*(dpq2/sigd2 + rangediffer2.val[kk] / sigr2));
								p3tup.val[kk] += sol[i][j][_k].val[kk] * tmp;
								p3tdown.val[kk] += tmp;
							}
						}
					for(int kk = 0;kk<3;kk++)
						p3tres.val[kk] = p3tup.val[kk] / p3tdown.val[kk];
					res[n][m][k] = p3tres;
			}
}

#endif /*_JBU_H_*/