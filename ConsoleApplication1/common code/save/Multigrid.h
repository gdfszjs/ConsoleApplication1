/*  
	本程序提供多重网格算法。
	class CSquareMulKer2D 对二维正方形图形（边长必须是2的n次方加1）的泊松方程进行求解。
	2008/3/20 by Yongwei Nie.

	class CRectMulKer2D 对二维长方形图形的泊松方程进行求解，长度和宽度也必须是2的n次方加1
	2008/3/22 by Yongwei Nie
	
	class CMulKer3D 对三维立方体的泊松方程求解，立方体的三个边长必须是2的N次方加1
	2008/3/25 by Yongwei Nie
	
*/

#ifndef _MULTIGRID_H_
#define _MULTIGRID_H_
#include <vector>
#include "datastructure.h"
#include "GaussSeidel.h"

using namespace NIE_DATASTRUCTURE;
using namespace std;

typedef Matrix2D<double> MAT_2D;
typedef const Matrix2D<double> MAT_C_2D;
typedef std::vector<MAT_2D *> VEC_MAT_2D_P;

typedef Matrix3D<double> MAT_3D;
typedef const Matrix3D<double> MAT_C_3D;
typedef std::vector<MAT_3D *> VEC_MAT_3D_P;

namespace NIE_MULTIGRID
{

/**************

对double类型二维数组的多重网格操作
二维数组采用名字空间NIE_DATASTRUCTURE中的Matrix2D进行存储

***************/

class CMulKer2D
{
public:
	//在最粗网格上用迭代法进行计算准确值，这里限定最粗网格拥有9个点
	virtual void slvsml(MAT_2D & u ,MAT_C_2D & rhs) {}
	//将一个vector中数据复制到另外一个中去
	virtual void matcopy(MAT_2D & res ,MAT_C_2D & rhs) {}
	//从一个较粗的网格插值得到较细的网格
	virtual void interp(MAT_2D & uf,MAT_C_2D & uc) {}
	//返回负残差
	virtual void resid(MAT_2D & res, MAT_C_2D & u,MAT_C_2D & rhs) {}
	//半加权限制
	virtual void rstrct(MAT_2D &uc,MAT_C_2D &uf) {}
	//红黑GAUSS-SEIDEL松弛
	virtual void relax(MAT_2D &u, MAT_C_2D &rhs) {}
	//粗细网格插值
	virtual void addint(MAT_2D &uf, MAT_C_2D &uc, MAT_2D &res) {}
	//递归调用多重网格法
	virtual void mg(int j , MAT_2D & u,MAT_C_2D & rhs) {}
	//完全多重网格法求解泊松方程
	virtual void mglin(MAT_2D & u , int ncycle) {}
};

class CSquareMulKer2D : public CMulKer2D//multigrid(2d) kernel 
{		
public:
	//在最粗网格上用迭代法进行计算准确值，这里限定最粗网格拥有9个点
	virtual void slvsml(MAT_2D & u ,MAT_C_2D & rhs);
	//将一个vector中数据复制到另外一个中去
	virtual void matcopy(MAT_2D & res ,MAT_C_2D & rhs);
	//从一个较粗的网格插值得到较细的网格
	virtual void interp(MAT_2D & uf,MAT_C_2D & uc);
	//返回负残差
	virtual void resid(MAT_2D & res, MAT_C_2D & u,MAT_C_2D & rhs);
	//半加权限制
	virtual void rstrct(MAT_2D &uc,MAT_C_2D &uf);
	//红黑GAUSS-SEIDEL松弛
	virtual void relax(MAT_2D &u, MAT_C_2D &rhs);
	//粗细网格插值
	virtual void addint(MAT_2D &uf, MAT_C_2D &uc, MAT_2D &res);
	//递归调用多重网格法
	virtual void mg(int j , MAT_2D & u,MAT_C_2D & rhs);
	//完全多重网格法求解泊松方程
	virtual void mglin(MAT_2D & u , int ncycle);
};

class CRectMulKer2D : public CMulKer2D
{
public:
	//在最粗网格上用迭代法进行计算准确值，这里限定最粗网格拥有9个点
	virtual void slvsml(MAT_2D & u ,MAT_C_2D & rhs);
	//将一个vector中数据复制到另外一个中去
	virtual void matcopy(MAT_2D & res ,MAT_C_2D & rhs);
	//从一个较粗的网格插值得到较细的网格
	virtual void interp(MAT_2D & uf,MAT_C_2D & uc);
	//返回负残差
	virtual void resid(MAT_2D & res, MAT_C_2D & u,MAT_C_2D & rhs);
	//半加权限制
	virtual void rstrct(MAT_2D &uc,MAT_C_2D &uf);
	//红黑GAUSS-SEIDEL松弛
	virtual void relax(MAT_2D &u, MAT_C_2D &rhs);
	//粗细网格插值
	virtual void addint(MAT_2D &uf, MAT_C_2D &uc, MAT_2D &res);
	//递归调用多重网格法
	virtual void mg(int j , MAT_2D & u,MAT_C_2D & rhs);
	//完全多重网格法求解泊松方程
	virtual void mglin(MAT_2D & u , int ncycle);
};

class CMulKer3D
{
public:
	void matcopy(MAT_3D & res,MAT_C_3D & rhs);
	void slvsml(MAT_3D &u,MAT_C_3D &rhs);
	void interp(MAT_3D & uf,MAT_C_3D &uc);
	void resid(MAT_3D &res,MAT_C_3D &u,MAT_C_3D &rhs);
	void rstrct(MAT_3D &uc,MAT_C_3D &uf);
	void relax(MAT_3D &u,MAT_C_3D &rhs);
	void addint(MAT_3D & uf,MAT_C_3D & uc,MAT_3D &res);
	void mg(int j,MAT_3D & u ,MAT_C_3D & rhs);
	void mglin(MAT_3D &u , int ncycle);
};

}



/*
  
  #############################################
  #############################################
  #############################################
  ######                                 ######
  ######   I M P L E M E N T A T I O N   ######
  ######                                 ######
  #############################################
  #############################################
  #############################################
  
*/





/*


  ###########################################

				CSquareMulKer2D

  ###########################################


*/



void NIE_MULTIGRID::CSquareMulKer2D::matcopy(MAT_2D & res ,MAT_C_2D & rhs)
{
	int n = rhs.getnrows();
	for(int i = 0;i<n;i++)
		for(int j = 0;j<n;j++)
			res[i][j] = rhs[i][j];
}

void NIE_MULTIGRID::CSquareMulKer2D::slvsml(MAT_2D & u ,MAT_C_2D & rhs)
{
	double h = 0.5;
	for(int i = 0;i<3;i++)
		for(int j = 0;j<3;j++)
			u[i][j] = 0.0;
	u[1][1] = -h * h * rhs[1][1] / 4.0;
}



void NIE_MULTIGRID::CSquareMulKer2D::interp(MAT_2D & uf,MAT_C_2D & uc)
{
	int ic,iif,jc,jf,nc;

	int nf=uf.getnrows();
	nc=nf/2+1;
	for (jc=0;jc<nc;jc++)
		for (ic=0;ic<nc;ic++) uf[2*ic][2*jc]=uc[ic][jc];
	for (jf=0;jf<nf;jf+=2)
		for (iif=1;iif<nf-1;iif+=2)
			uf[iif][jf]=0.5*(uf[iif+1][jf]+uf[iif-1][jf]);
	for (jf=1;jf<nf-1;jf+=2)
		for (iif=0;iif<nf;iif++)
			uf[iif][jf]=0.5*(uf[iif][jf+1]+uf[iif][jf-1]);

}

void NIE_MULTIGRID::CSquareMulKer2D::resid(MAT_2D &res, MAT_C_2D &u, MAT_C_2D &rhs)
{
	int i,j;
	double h,h2i;

	int n=u.getnrows();
	h=1.0/(n-1);
	h2i=1.0/(h*h);
	for (j=1;j<n-1;j++)
		for (i=1;i<n-1;i++)
			res[i][j] = -h2i*(u[i+1][j]+u[i-1][j]+u[i][j+1]
				+u[i][j-1]-4.0*u[i][j])+rhs[i][j];
	for (i=0;i<n;i++)
		res[i][0]=res[i][n-1]=res[0][i]=res[n-1][i]=0.0;

}

void NIE_MULTIGRID::CSquareMulKer2D::rstrct(MAT_2D &uc,MAT_C_2D &uf)
{
	int ic,iif,jc,jf,ncc;

	int nc=uc.getnrows();
	ncc=2*nc-2;
	for (jf=2,jc=1;jc<nc-1;jc++,jf+=2) {
		for (iif=2,ic=1;ic<nc-1;ic++,iif+=2) {
			uc[ic][jc]=0.5*uf[iif][jf]+0.125*(uf[iif+1][jf]+uf[iif-1][jf]
				+uf[iif][jf+1]+uf[iif][jf-1]);
		}
	}
	for (jc=0,ic=0;ic<nc;ic++,jc+=2) {
		uc[ic][0]=uf[jc][0];
		uc[ic][nc-1]=uf[jc][ncc];
	}
	for (jc=0,ic=0;ic<nc;ic++,jc+=2) {
		uc[0][ic]=uf[0][jc];
		uc[nc-1][ic]=uf[ncc][jc];
	}

}

void NIE_MULTIGRID::CSquareMulKer2D::relax(MAT_2D &u, MAT_C_2D &rhs)
{
	int i,ipass,isw,j,jsw=1;
	double h,h2;

	int n=u.getnrows();
	h=1.0/(n-1);
	h2=h*h;
	for (ipass=0;ipass<2;ipass++,jsw=3-jsw) {
		isw=jsw;
		for (j=1;j<n-1;j++,isw=3-isw)
			for (i=isw;i<n-1;i+=2)
				u[i][j]=0.25*(u[i+1][j]+u[i-1][j]+u[i][j+1]
					+u[i][j-1]-h2*rhs[i][j]);
	}
}

void NIE_MULTIGRID::CSquareMulKer2D::addint(MAT_2D &uf, MAT_C_2D &uc, MAT_2D &res)
{
	int i,j;

	int nf=uf.getnrows();
	interp(res,uc);
	for (j=0;j<nf;j++)
		for (i=0;i<nf;i++)
			uf[i][j] += res[i][j];
}

void NIE_MULTIGRID::CSquareMulKer2D::mg(int j, MAT_2D &u, MAT_C_2D &rhs)
{
	const int NPRE=1,NPOST=1;
	int jpost,jpre,nc,nf;

	nf=u.getnrows();
	nc=(nf+1)/2;
	if (j == 0)
		slvsml(u,rhs);
	else {
		MAT_2D res(nc,nc),v(0.0,nc,nc),temp(nf,nf);
		for (jpre=0;jpre<NPRE;jpre++)
			relax(u,rhs);
		resid(temp,u,rhs);//残差计算 
		rstrct(res,temp);
		mg(j-1,v,res);
		addint(u,v,temp);
		for (jpost=0;jpost<NPOST;jpost++)
			relax(u,rhs);
	}
}

void NIE_MULTIGRID::CSquareMulKer2D::mglin(MAT_2D &u, int ncycle)
{
	int j,jcycle,ng=0,ngrid,nn;
	MAT_2D *uj,*uj1;

	int n=u.getnrows();

	nn=n;
	while (nn >>= 1) ng++;
	VEC_MAT_2D_P  rho(ng);
	nn=n;
	ngrid=ng-1;
	rho[ngrid] = new MAT_2D(nn,nn);
	matcopy(*rho[ngrid],u);
	while (nn > 3) {
		nn=nn/2+1;
		rho[--ngrid]=new MAT_2D(nn,nn);
		rstrct(*rho[ngrid],*rho[ngrid+1]);
	}
	nn=3;
	uj=new MAT_2D(nn,nn);
	slvsml(*uj,*rho[0]);
	for (j=1;j<ng;j++) {
		nn=2*nn-1;
		uj1=uj;
		uj=new MAT_2D(nn,nn);
		interp(*uj,*uj1);
		delete uj1;
		for (jcycle=0;jcycle<ncycle;jcycle++)
			mg(j,*uj,*rho[j]);
	}
	matcopy(u,*uj);
	delete uj;
	for (j=0;j<ng;j++)
		delete rho[j];
}



/*


  ###########################################

				CRectMulKer2D

  ###########################################


*/



void NIE_MULTIGRID::CRectMulKer2D::slvsml(MAT_2D &u, MAT_C_2D &rhs)
{
	int r = u.getnrows();
	int c = u.getncols();
	double h1 = 1.0/(r-1);
	double h2 = 1.0/(c-1);
	for(int i = 0;i<r;i++)
		for(int j = 0;j<c;j++)
			u[i][j] = 0.0;
	int nn;
	if(r == 3)
	{
		nn = c - 2;
		MAT_2D coemat(0.0,nn,nn);
		MAT_2D rmat(nn,1);
		MAT_2D resmat(0.0,nn,1);
		
		for(int i = 0;i<nn;i++)
		{
			coemat[i][i] = -4;
			if(i != 0)	  coemat[i][i-1] = 1;
			if(i != nn - 1 ) coemat[i][i+1] = 1;
		}

		for(int i = 0;i<nn;i++)
			rmat[i][0] = rhs[1][i+1];

		NIE_GAUSSSEIDEL::CGaussSeidel gs(0.000001);
		gs.gs(resmat,coemat,rmat);

		for(int i = 0;i<nn;i++)
			u[1][i+1] = resmat[i][0]*h1*h2;
	}
	else if(c == 3)
	{
		nn = r - 2;
		MAT_2D coemat(0.0,nn,nn);
		MAT_2D rmat(nn,1);
		MAT_2D resmat(0.0,nn,1);

		for(int i = 0;i<nn;i++)
		{
			coemat[i][i] = -4;
			if(i != 0) coemat[i][i-1] = 1;
			if(i != nn-1) coemat[i][i+1] = 1;
		}
		
		for(int i = 0;i<nn;i++)
			rmat[i][0] = rhs[i+1][1];

		NIE_GAUSSSEIDEL::CGaussSeidel gs(0.000001);
		gs.gs(resmat,coemat,rmat);
		for(int i = 0;i<nn;i++)
			u[i+1][1] = resmat[i][0]*h1*h2;
	}
}


void NIE_MULTIGRID::CRectMulKer2D::matcopy(MAT_2D &res,MAT_C_2D &rhs)
{
	int r = res.getnrows();
	int c = res.getncols();
	for(int i = 0;i<r;i++)
		for(int j = 0;j<c;j++)
			res[i][j] = rhs[i][j];
}

void NIE_MULTIGRID::CRectMulKer2D::interp(MAT_2D & uf,MAT_C_2D & uc)
{
	int fr = uf.getnrows();
	int fc = uf.getncols();
	int cr = fr/2 + 1;
	int cc = fc/2 + 1;

	for(int i = 0;i<cr;i++)
		for(int j = 0;j<cc;j++)
			uf[i*2][j*2] = uc[i][j];
	for(int i = 1;i<fr - 1;i+=2)
		for(int j = 0;j<fc;j++)
			uf[i][j] = 0.5 * ( uf[i-1][j] + uf[i+1][j]);
	for(int i = 0;i<fr;i++)
		for(int j = 1;j<fc-1;j+=2)
			uf[i][j] = 0.5 * ( uf[i][j-1] + uf[i][j+1]);
}

void NIE_MULTIGRID::CRectMulKer2D::resid(MAT_2D &res, MAT_C_2D &u, MAT_C_2D &rhs)
{
	int r = u.getnrows();
	int c = u.getncols();
	double h1 = 1.0/(r-1);
	double h2 = 1.0/(c-1);
	double h2i =1.0/(h1 * h2);
	
	for(int i = 1;i<r-1;i++)
		for(int j = 1;j<c-1;j++)
			res[i][j] = -h2i*(u[i-1][j] + u[i+1][j] + u[i][j-1] 
		+ u[i][j+1] - 4.0 * u[i][j]) + rhs[i][j];
	for(int i = 0;i<r;i++)
		res[i][0] = res[i][c-1] = 0.0;
	for(int i = 0;i<c;i++)
		res[0][i] = res[r-1][i] = 0.0;
}

void NIE_MULTIGRID::CRectMulKer2D::rstrct(MAT_2D &uc,MAT_C_2D &uf)
{
	int cr = uc.getnrows();
	int cc = uc.getncols();
	int fr = 2*cr - 2;
	int fc = 2*cc - 2;

	for(int jf = 2 , jc = 1;jc<cc-1;jc++,jf+=2)
		for(int iif =2,ic =1;ic<cr-1;ic++,iif+=2)
			uc[ic][jc] = 0.5 * uf[iif][jf] + 0.125 * ( uf[iif-1][jf] + uf[iif+1][jf] 
		+ uf[iif][jf-1] + uf[iif][jf+1] );
	
	for(int jc = 0,ic = 0;ic < cr ;ic++,jc+=2)
	{
		uc[ic][0] = uf[jc][0];
		uc[ic][cc-1] = uf[jc][fc];
	}
	
	for(int jc= 0, ic = 0;ic < cc;ic++,jc+=2)
	{
		uc[0][ic] = uf[0][jc];
		uc[cr-1][ic] = uf[fr][ic];
	}
}


void NIE_MULTIGRID::CRectMulKer2D::relax(MAT_2D &u, MAT_C_2D &rhs)
{
	int r = u.getnrows();
	int c = u.getncols();
	double h1 = 1.0/(r-1);
	double h2 = 1.0/(c-1);
	double h2i = h1*h2;
	int jsw = 1;
	int isw;
	for(int ipass = 0;ipass < 2;ipass++,jsw=3-jsw){
		isw = jsw;
		for(int j = 1;j<c-1;j++,isw = 3-isw)
			for(int i = isw;i<r-1;i+=2)
				u[i][j] = 0.25 * ( u[i+1][j] + u[i-1][j] + u[i][j+1] + u[i][j-1] - h2i * rhs[i][j] );
	}

}

void NIE_MULTIGRID::CRectMulKer2D::addint(MAT_2D & uf,MAT_C_2D &uc,MAT_2D &res)
{
	int r = uf.getnrows();
	int c = uf.getncols();
	interp(res,uc);
	for(int j  =0;j<c;j++)
		for(int i =0;i<r;i++)
			uf[i][j] += res[i][j];
}

void NIE_MULTIGRID::CRectMulKer2D::mg(int j,MAT_2D & u ,MAT_C_2D &rhs)
{
	const int NPRE = 1 , NPOST = 1;
	int fr = u.getnrows();
	int fc = u.getncols();
	int cr = (fr+1) / 2;
	int cc = (fc+1) / 2;
	if(j==0)
		slvsml(u,rhs);
	else{
		MAT_2D res(cr,cc),v(0.0,cr,cc),temp(fr,fc);
		for(int jpre = 0;jpre< NPRE;jpre++)
			relax(u,rhs);
		resid(temp,u,rhs);
		rstrct(res,temp);
		mg(j-1,v,res);
		addint(u,v,temp);
		for(int jpost = 0;jpost<NPOST;jpost++)
			relax(u,rhs);
	}
}

void NIE_MULTIGRID::CRectMulKer2D::mglin(MAT_2D &u,int ncycle)
{
	int j,jcycle,ng,ngr = 0,ngc = 0,ngrid,nnr,nnc;
	MAT_2D *uj,*uj1;

	int nr=u.getnrows();
	int nc=u.getncols();

	nnr=nr;
	nnc = nc;
	
	int nn = min(nr,nc);

	while (nnr >>= 1) ngr++;
	while (nnc >>= 1) ngc++;

	ng = min(ngr,ngc);

	VEC_MAT_2D_P  rho(ng);
	nnr=nr;
	nnc = nc;
	ngrid=ng-1;
	rho[ngrid] = new MAT_2D(nnr,nnc);
	matcopy(*rho[ngrid],u);
	while (nn > 3) {
		nn=nn/2+1;
		nnr = nnr / 2 + 1;
		nnc = nnc / 2 + 1;
		rho[--ngrid]=new MAT_2D(nnr,nnc);
		rstrct(*rho[ngrid],*rho[ngrid+1]);
	}
	uj=new MAT_2D(nnr,nnc);
	slvsml(*uj,*rho[0]);
	for (j=1;j<ng;j++) {
		nnr = 2 * nnr -1;
		nnc = 2 * nnc -1;
		uj1=uj;
		uj=new MAT_2D(nnr,nnc);
		interp(*uj,*uj1);
		delete uj1;
		for (jcycle=0;jcycle<ncycle;jcycle++)
			mg(j,*uj,*rho[j]);
	}
	matcopy(u,*uj);
	delete uj;
	for (j=0;j<ng;j++)
		delete rho[j];
}




/*


  ###########################################

				CMulKer3D

  ###########################################


*/


void NIE_MULTIGRID::CMulKer3D::matcopy(MAT_3D &res,MAT_C_3D & rhs)
{
	int n = rhs.dim1();
	int m = rhs.dim2();
	int k = rhs.dim3();
	for(int z = 0;z<n;z++)
		for(int y = 0;y<m;y++)
			for(int x = 0;x<k;x++)
				res[z][y][x] = rhs[z][y][x];
}

void NIE_MULTIGRID::CMulKer3D::slvsml(MAT_3D & u,MAT_C_3D & rhs)
{
	int n = u.dim1();
	int m = u.dim2();
	int k = u.dim3();

	double hn = 1.0/(n-1);
	double hm = 1.0/(m-1);
	double hk = 1.0/(k-1);

	for(int z = 0;z<n;z++)
		for(int y = 0;y<m;y++)
			for(int x = 0;x<k;x++)
				u[z][y][x] = 0.0;
	
	int tt = (n-2)*(m-2)*(k-2);

	MAT_2D coemat(0.0,tt,tt);
	MAT_2D rmat(tt,1);
	MAT_2D resmat(0.0,tt,1);
	
	for(int z = 1;z<n-1;z++)
		for(int y = 1;y<m-1;y++)
			for(int x = 1;x<k-1;x++)
			{
				coemat[(z-1) * (m-2)*(k-2) + (y-1) * (k-2) + (x-1)][(z-1) * (m-2)*(k-2) + (y-1) * (k-2) + (x-1)] = -6;
				
				if(y-1 != 0)
					coemat[(z-1) * (m-2)*(k-2) + (y-1) * (k-2) + (x-1)][(z-1) * (m-2)*(k-2) + (y-2) * (k-2) + (x-1)] = 1;

				if(y+1 != m-1)
					coemat[(z-1) * (m-2)*(k-2) + (y-1) * (k-2) + (x-1)][(z-1) * (m-2)*(k-2) + (y) * (k-2) + (x-1)] = 1;
				
				if(x-1 != 0)
					coemat[(z-1) * (m-2)*(k-2) + (y-1) * (k-2) + (x-1)][(z-1) * (m-2)*(k-2) + (y-1) * (k-2) + (x-2)] = 1;
				
				if(x+1 != k-1)
					coemat[(z-1) * (m-2)*(k-2) + (y-1) * (k-2) + (x-1)][(z-1) * (m-2)*(k-2) + (y-1) * (k-2) + (x)] = 1;
				
				if(z-1 != 0)
					coemat[(z-1) * (m-2)*(k-2) + (y-1) * (k-2) + (x-1)][(z-2) * (m-2)*(k-2) + (y-1) * (k-2) + (x-1)] = 1;
				
				if(z+1 != n-1)
					coemat[(z-1) * (m-2)*(k-2) + (y-1) * (k-2) + (x-1)][(z) * (m-2)*(k-2) + (y-1) * (k-2) + (x-1)] = 1;
				
				rmat[(z-1) * (m-2)*(k-2) + (y-1) * (k-2) + (x-1)][0] = rhs[z][y][x];
			}
	NIE_GAUSSSEIDEL::CGaussSeidel gs(0.000001);
	gs.gs(resmat,coemat,rmat);
	
	for(int z = 1;z<n-1;z++)
		for(int y = 1;y<m-1;y++)
			for(int x = 1;x<k-1;x++)
				u[z][y][x] = resmat[(z-1) * (m-2)*(k-2) + (y-1) * (k-2) + (x-1)][0]*hn*hm*hk;
}


void NIE_MULTIGRID::CMulKer3D::interp(MAT_3D &uf,MAT_C_3D & uc)
{
	int nnf = uf.dim1();
	int mmf = uf.dim2();
	int kkf = uf.dim3();

	int nnc = nnf/2 + 1;
	int mmc = mmf/2 + 1;
	int kkc = kkf/2 + 1;

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


void NIE_MULTIGRID::CMulKer3D::resid(MAT_3D & res,MAT_C_3D &u,MAT_C_3D &rhs)
{
	int nn = u.dim1();
	int mm = u.dim2();
	int kk = u.dim3();

	double hn = 1.0/(nn-1);
	double hm = 1.0/(mm-1);
	double hk = 1.0/(kk-1);

	double h2i = 1.0/(hn*hm*hk);

	for(int z = 1;z<nn-1;z++)
		for(int y = 1;y<mm-1;y++)
			for(int x = 1;x<kk-1;x++)
				res[z][y][x] = -h2i*(u[z+1][y][x] + u[z-1][y][x] + u[z][y+1][x] + u[z][y-1][x]
			+ u[z][y][x+1] + u[z][y][x-1] - 6.0 * u[z][y][x]) + rhs[z][y][x];
	
	for(int y = 0;y<mm;y++)
		for(int x = 0;x<kk;x++)
			res[0][y][x] = res[nn-1][y][x] = 0.0;
	for(int z = 0;z<nn;z++)
		for(int y = 0;y<mm;y++)
			res[z][y][0] = res[z][y][kk-1] = 0.0;
	for(int z = 0;z<nn;z++)
		for(int x = 0;x<kk;x++)
			res[z][0][x] = res[z][mm-1][x] = 0.0;
}


void NIE_MULTIGRID::CMulKer3D::rstrct(MAT_3D &uc,MAT_C_3D &uf)
{
	int nnc = uc.dim1();
	int mmc = uc.dim2();
	int kkc = uc.dim3();

	int nnf = 2*nnc - 2;
	int mmf = 2*mmc - 2;
	int kkf = 2*kkc - 2;

	for(int zf = 2,zc=1;zc<nnc-1;zc++,zf+=2)
		for(int yf = 2,yc=1;yc<mmc-1;yc++,yf+=2)
			for(int xf = 2,xc=1;xc<kkc-1;xc++,xf+=2)
				uc[zc][yc][xc] = 0.5 * uf[zf][yf][xf] + 0.083333 * ( uf[zf+1][yf][xf] + uf[zf-1][yf][xf]
			+ uf[zf][yf+1][xf] + uf[zf][yf-1][xf] + uf[zf][yf][xf+1] + uf[zf][yf][xf-1] );

	for(int yf = 0,yc = 0;yc<mmc;yc++,yf+=2)
		for(int xf = 0,xc = 0;xc<kkc;xc++,xf+=2){
			uc[0][yc][xc] = uf[0][yf][xf];
			uc[nnc-1][yc][xc] = uf[nnf][yf][xf];
		}

	for(int zf = 0,zc = 0;zc<nnc;zc++,zf+=2)
		for(int yf = 0,yc = 0;yc<mmc;yc++,yf+=2){
			uc[zc][yc][0] = uf[zf][yf][0];
			uc[zc][yc][kkc-1] = uf[zf][yf][kkf];
		}

	for(int zf = 0,zc = 0;zc<nnc;zc++,zf+=2)
		for(int xf = 0,xc = 0;xc<kkc;xc++,xf+=2){
			uc[zc][0][xc] = uf[zf][0][xf];
			uc[zc][mmc-1][xc] = uf[zf][mmf][xf];
		}
}


void NIE_MULTIGRID::CMulKer3D::relax(MAT_3D &u,MAT_C_3D &rhs)
{
	int nn = u.dim1();
	int mm = u.dim2();
	int kk = u.dim3();

	double hn = 1.0/(nn-1);
	double hm = 1.0/(mm-1);
	double hk = 1.0/(kk-1);
	double h = hn * hm * hk;

	int jsw;
	int isw;
	int ksw = 1;

	for(int ipass = 0;ipass<2;ipass++,ksw = 3-ksw){
		jsw = ksw;
		for(int z = 1;z < nn-1;z++,jsw = 3-jsw){
			isw = jsw;
			for(int x = 1;x<kk-1;x++,isw = 3 - isw)
				for(int y = isw;y<mm-1;y+=2)
					u[z][y][x] = 0.166667*( u[z+1][y][x] + u[z-1][y][x] + u[z][y+1][x] + u[z][y-1][x] + u[z][y][x+1]
				+ u[z][y][x-1] - h * rhs[z][y][x]);
		}
	}

}



void NIE_MULTIGRID::CMulKer3D::addint(MAT_3D &uf, MAT_C_3D &uc, MAT_3D &res)
{
	int nn = uf.dim1();
	int mm = uf.dim2();
	int kk = uf.dim3();

	interp(res,uc);

	for(int z = 0;z<nn;z++)
		for(int y = 0;y<mm;y++)
			for(int x = 0;x<kk;x++)
				uf[z][y][x] += res[z][y][x];
}

void NIE_MULTIGRID::CMulKer3D::mg(int j, MAT_3D &u, MAT_C_3D &rhs)
{
	const int NPRE=1,NPOST=1;
	
	int nf = u.dim1();
	int mf = u.dim2();
	int kf = u.dim3();

	int nc = (nf+1)/2;
	int mc = (mf+1)/2;
	int kc = (kf+1)/2;

	if (j == 0)
		slvsml(u,rhs);
	else {
		MAT_3D res(nc,mc,kc),v(0.0,nc,mc,kc),temp(nf,mf,kf);
		for (int jpre=0;jpre<NPRE;jpre++)
			relax(u,rhs);
		resid(temp,u,rhs);
		rstrct(res,temp);
		mg(j-1,v,res);
		addint(u,v,temp);
		for (int jpost=0;jpost<NPOST;jpost++)
			relax(u,rhs);
	}
}


void NIE_MULTIGRID::CMulKer3D::mglin(MAT_3D &u, int ncycle)
{
	int j,jcycle,ng=0,ngrid;
	MAT_3D *uj,*uj1;

	int nnu=u.dim1();
	int mmu=u.dim2();
	int kku=u.dim3();

	int nnf = nnu;
	int mmf = mmu;
	int kkf = kku;

	int nn = min(nnf,min(mmf,kkf));
	
	while (nn >>= 1) ng++;

	VEC_MAT_3D_P  rho(ng);

	nnf = nnu;
	mmf = mmu;
	kkf = kku;
	
	nn = min(nnf,min(mmf,kkf));
	ngrid=ng-1;

	rho[ngrid] = new MAT_3D(nnf,mmf,kkf);
	matcopy(*rho[ngrid],u);

	while (nn > 3) {
		nn=nn/2+1;
		nnf = nnf/2+1;
		mmf = mmf/2+1;
		kkf = kkf/2+1;
		rho[--ngrid]=new MAT_3D(nnf,mmf,kkf);
		rstrct(*rho[ngrid],*rho[ngrid+1]);
	}

	uj=new MAT_3D(nnf,mmf,kkf);
	slvsml(*uj,*rho[0]);
	for (j=1;j<ng;j++) {
		nnf = 2*nnf-1;
		mmf = 2*mmf-1;
		kkf = 2*kkf-1;
		uj1=uj;
		uj=new MAT_3D(nnf,mmf,kkf);
		interp(*uj,*uj1);
		delete uj1;
		for (jcycle=0;jcycle<ncycle;jcycle++)
			mg(j,*uj,*rho[j]);
	}
	matcopy(u,*uj);
	delete uj;
	for (j=0;j<ng;j++)
		delete rho[j];
}



#endif /* _MULTIGRID_H_ */