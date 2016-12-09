//#include "stdafx.h"
#include "PCG.h"
#include <cmath>
PCG::PCG(double * &sa_,unsigned long *&ija_,double *&x_,double *&b_,const unsigned long n_,
											const double tol_,const int itmax_)
{
	sa=sa_;
	ija=ija_;
	x=x_;
	b=b_;
	n=n_;
	sb=0;
	tol=tol_;
	itmax=itmax_;
	flag=0;
	nz=0;
}
PCG::PCG(double * &sa_,unsigned long *&ija_,double *&x_,double *&b_,const unsigned long n_,
								const double tol_,const int itmax_,const int flag_,const int nz_)
{
	sa=sa_;
	ija=ija_;
	x=x_;
	b=b_;
	n=n_;
	sb=0;
	tol=tol_;
	itmax=itmax_;
	flag=flag_;
	nz=nz_;
}
PCG::~PCG()
{
	if(sb)
	{
		delete [] sb;
		sb=0;
	}
}
bool PCG::Choldc(double *&sa1,unsigned long *& ija1,const unsigned long n)
{
	unsigned long i,j,k,t,l,m;
	double sum;
//	//////////////////////////////首先进行预处理，使对运算更加稳定
//	for(i=0;i<n;i++)
//	{
//		sum=sa1[i]*sa1[i];
//		for(k=ija1[i];k<ija1[i+1];k++)
//			sum+=sa1[k]*sa1[k];
//		sa1[i]=sqrt(sum);
//	}
//
	for(i=0;i<n;i++)//对每一行分别进行处理
	{
		sum=sa1[i];//从对角线的元素开始
		for(k=ija1[i];ija1[k]<i && k<ija1[i+1];k++)
			sum-=sa1[k]*sa1[k];

		if(sum>0)
			sa1[i]=sqrt(sum);
	//	else 
	//	{
	//		::AfxMessageBox("乔列斯基分解错误！矩阵不符合要求！");
	//		return false;
			
	//	}
		/////////////////对对角线以后的元素处理
	for(t=ija1[i];t<ija1[i+1];t++)
		if(ija1[t]>i)
			break;
		for(j=t;j<ija1[i+1];j++)//第i行，第j个元素
		{
			sum=sa1[j];
			for(k=ija1[i];ija1[k]<i && k<ija1[i+1];k++)
				for( l=ija1[ija1[j]];ija1[l]<i && l<ija1[ija1[j]+1];l++)//列号为ija[j],作为行号
					if(ija1[l]==ija1[k])
					{
						sum-=sa1[k]*sa1[l];//当列值相等时才做累减
					}
			for( m=ija1[ija1[j]];m<ija1[ija1[j]+1];m++)//列号为ija[j],作为行号
				if(ija1[m]==i)
					break;
			sa1[m]=(sa1[i]!=0? sum/sa1[i]:sum);                                                                                                                                                                                                        
		}

	}	
	return true;
}

void PCG::matrix_SSOR(double * &sa,unsigned long * &ija,const double w,const unsigned long n,double * &sa1)
{
    unsigned long i,k;
    for(i=0;i<n;i++)
    {
        for(k=ija[i];k<ija[i+1] && ija[k]<i;k++)
                sa1[k]=+w*sa[k];
        sa1[i]=sa[i]+w*sa[i];
    }
    for(i=0;i<n;i++)
    {
        for(k=ija[i];k<ija[i+1] && ija[k]<i;k++)
               {
                    sa1[k]/=sqrt(sa[ija[k]]);
                    sa1[k]/=sqrt(w*(2-w));
                }    
        sa1[i]/=sqrt(sa[i]);
        sa1[i]/=sqrt(w*(2-w));
    }
}
////函数dsprsax求sa ija 乘以向量组x,结果在b中,n 表示原矩阵为n*n
void PCG::dsprsax(double *sa,unsigned long *ija,double *x,double *b ,unsigned long n)
{
	unsigned long i,k;
		for(i=0;i<n;i++){
			b[i]=sa[i]*x[i];
			for(k=ija[i];k<ija[i+1];k++)
				b[i] += sa[k]*x[ija[k]];
		}
}
void PCG::asolve(double *sa1,unsigned long *ija,double *x,double *b,unsigned long n)
{
	if(flag==0)
	{
		unsigned long i;
		for(i=0;i<n;i++)
			x[i]=(sa[i]!= 0.0?b[i]/sa[i]:b[i]);
	}
if(flag==1 || flag==2)//SSOR 和 Cholsky相同 
{	
	unsigned long i,k;
	double sum;
/*	for(i=0;i<n;i++)
	{
		sum=sa1[i]*sa1[i];
		for(k=ija[i];k<ija[i+1] && ija[k]<i;k++)
			sum+=sa1[k]*sa1[k];
		sa1[i]=sqrt(sum);
	}
*/
	for(i=0;i<n;i++)
	{
		sum=b[i];
		for(k=ija[i];k<ija[i+1] && ija[k]<i;k++)
				sum-=sa1[k]*x[ija[k]];
		x[i]=sum/sa1[i];
	}
	for(i=n-1;i>=0;i--)
	{
		x[i]=x[i]/sa1[i];
			for(k=ija[i];k<ija[i+1] && ija[k]<i;k++)
					x[ija[k]]-=sa1[k]*x[i];
	}
}

}

int PCG::linbcg(void)
{
	if(flag==1 )
	{
		if(nz==0)
			return 0;
			
		if(sb)
		  delete [] sb;
		  
		sb=new double[nz];
		for(unsigned long i=0;i<nz;i++)
			sb[i]=sa[i];
		
		if(!Choldc(sb,ija,n))
			return 0;
	}
	if(flag==2)
	{
		if(nz==0)
			return 0;
			
		if(sb)
		  delete [] sb;
		
		sb=new double[nz];
		const double w=1.3;//超松弛因子取1.3
		matrix_SSOR(sa,ija,w,n,sb);
	}
	
		unsigned long j;
		double ak,akden,bk,err;
		double *p = new double[n];
		double *r = new double[n];
		double *z = new double[n];
		double *w=new double[n];
		double rol;//ρ
		double roll;//ρ'
		double rol0;//ρ0
		//求r=Ax
		dsprsax(sa,ija,x,r,n);
		for(j=0;j<n;j++)     //求r=b-Ax
			r[j] = b[j]-r[j];
		//z0初值 Mz=r
		asolve(sb,ija,z,r,n);
		for(j=0;j<n;j++)//p1=z0
			p[j] = z[j];
		for(rol=0.0,j=0;j<n;j++)//ρ0=rt*z0
			rol+= z[j]*r[j];
		rol0=rol;
		roll=rol;
		iter=1;
		while(iter<=itmax)
		{
			//atimes(n,p,w,0);//求w=Ap
			dsprsax(sa,ija,p,w,n);
			for(akden=0.0,j=0;j<n;j++) //α=ρ/pt*w
				akden += p[j]*w[j];
			ak = roll/akden;

			for(j=0;j<n;j++)//x=x+αp  r=r-αw
			{
				x[j] += ak*p[j];
				r[j] -= ak*w[j];
			}

			//asolve(n,r,z,0);//求解Mz=r
			asolve(sb,ija,z,r,n);
							//求解rol=rt*z
			for(rol=0.0,j=0;j<n;j++)
				rol+= z[j]*r[j];
				bk = rol/roll;//β=ρ/ρ'
			for(j=0;j<n;j++)//p=β*z+p
				p[j] = bk*p[j]+z[j];
		
			roll =rol;
			err=rol/rol0;//ε=ρ/ρ0
			if (err <= tol){
			//	fprintf(stderr,"求解结束，总迭代次数为： %d\n",iter);
				break;
			}
			else
				iter++;
		}
//		if(iter>itmax)
//			;
//				fprintf(stderr,"超过最大迭代次数！迭代完毕！\n");
//		fprintf(stderr,"共轭梯度法求解结束\n");
				//MessageBox("超过最大迭代次数！迭代完毕！",NULL,MB_OK);
		delete [] p;
		delete [] r;
		delete [] z;
		delete [] w;
		
		return iter;
}

