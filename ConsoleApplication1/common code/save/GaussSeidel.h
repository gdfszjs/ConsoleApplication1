//2008/2/22
//Nie yong wei 
//Modified on 2008/3/22 

/***************************
@用高斯-塞德尔方法解线性方程组。注意，只能用来解决拥有小型系数矩阵的线性方程组。
@输入：系数矩阵（用到NIE_DATASTRUCTURE里面的Matrix2D来存储）
	   右端矩阵（同上）
@输出：结果矩阵（同上）
****************************/


#ifndef _GAUSS_SEIDEL_H_
#define _GAUSS_SEIDEL_H_

#include "datastructure.h"
#include <cmath>
using namespace NIE_DATASTRUCTURE;

typedef Matrix2D<double> MAT_2D;
typedef const Matrix2D<double> MAT_C_2D;


namespace NIE_GAUSSSEIDEL{

class CGaussSeidel
{
private:
	double e;
public:
	CGaussSeidel(double _e);
	void gs(MAT_2D & res,MAT_C_2D & coemat,MAT_C_2D & rmat);
};

}

NIE_GAUSSSEIDEL::CGaussSeidel::CGaussSeidel(double _e) : e(_e) { }

void NIE_GAUSSSEIDEL::CGaussSeidel::gs(MAT_2D & res,MAT_C_2D & coemat,MAT_C_2D & rmat)
{
	int n = coemat.getnrows();
	double _e,t;
	while(1)
	{
		_e = 0.0;
		for(int i = 0;i<n;i++)
		{
			t = res[i][0];
			double tmp = 0.0;
			for(int j = 0;j<n;j++)
				if(j != i)
					tmp += coemat[i][j] * res[j][0];
			res[i][0] = (rmat[i][0] - tmp) / coemat[i][i];

			if(fabs(t - res[i][0]) > _e)
				_e = fabs(t-res[i][0]);
		}
		if(_e<e)
			break;
	}
}


 
#endif   /* _GAUSS_SEIDEL_H_ */