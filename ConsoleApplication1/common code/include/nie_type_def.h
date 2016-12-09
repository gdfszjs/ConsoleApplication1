/**********************************

			定义类型
	2008/4/2 by Yongwei Nie 

***********************************/
#pragma once
#include "datastructure.h"
#include "nie_type.h"

using namespace NIE_TYPE;
using namespace NIE_DATASTRUCTURE;

typedef Matrix3D<int>               MAT_INT_3D;
typedef const Matrix3D<int>			MAT_C_INT_3D;

typedef Matrix3D<double>			MAT_DOU_3D;
typedef const Matrix3D<double>		MAT_C_DOU_3D;

typedef nie_point_3d_t				POINT_3D_T;

typedef Matrix3D<POINT_3D_T>		MAT_3D_POINT_3D_T;
typedef const Matrix3D<POINT_3D_T>	MAT_C_3D_POINT_3D_T;