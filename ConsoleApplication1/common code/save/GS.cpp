//#include "stdafx.h"
#include "GS.h"
#include <cmath>
GS::GS(void)
{
}
GS::GS(double *&sa_,unsigned long *&ija_,double * &x_, double *&b_,const int n_)
:sa(sa_),ija(ija_),x(x_),b(b_),n(n_)
{
}
GS::~GS(void)
{
}
void GS::Do_GS(int max,double e)
{
	for(int k=0;k<max;k++)
	{
		double temp_e=0.0;
		for(int i=0;i<n;i++)
		{
			double sum=b[i];
			double temp=x[i];
			for(unsigned int j=ija[i];j<ija[i+1];j++)
				sum-=sa[j]*x[ija[j]];
			x[i]=sum/sa[i];
			if(fabs(x[i]-temp)>temp_e)
				temp_e=fabs(x[i]-temp);
		}
		if(temp_e<e)
		{
			//fprintf(stderr,"达到迭代精度 迭代次数为：%d :\n",k);
			return ;
		}
	}
	//fprintf(stderr,"完全迭代 迭代次数为：%d \n",max);
}
