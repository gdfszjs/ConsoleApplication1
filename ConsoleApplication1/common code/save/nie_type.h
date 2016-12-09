/********

	定义常用类型，主要用结构体来进行定义
	2008/3/27
	
	更改于2008/4/1
	1.将不需要的类型去掉了，只留下了nie_point_3d_t

********/
#ifndef _NIE_TYPE_H_
#define _NEI_TYPE_H_

namespace NIE_TYPE
{
struct nie_point_3d_t
{
	double val[4];
	nie_point_3d_t();
	nie_point_3d_t(const double& t);
	nie_point_3d_t(const nie_point_3d_t& rhs);
	nie_point_3d_t& operator=(const nie_point_3d_t& rhs); 
	nie_point_3d_t& operator=(const double & t);
	double& operator[](const int i);
};


}


NIE_TYPE::nie_point_3d_t::nie_point_3d_t(){}

NIE_TYPE::nie_point_3d_t::nie_point_3d_t(const double &t)
{
	val[0] = val[1] = val[2] = val[3] = t;
}

NIE_TYPE::nie_point_3d_t::nie_point_3d_t(const NIE_TYPE::nie_point_3d_t &rhs)
{
	for(int i = 0;i<4;i++)
		val[i] = rhs.val[i];
}

NIE_TYPE::nie_point_3d_t& NIE_TYPE::nie_point_3d_t::operator =(const double &t)
{
	val[0] = val[1] = val[2] = val[3] = t;
	return *this;
}

NIE_TYPE::nie_point_3d_t& NIE_TYPE::nie_point_3d_t::operator =(const NIE_TYPE::nie_point_3d_t &rhs)
{
	for(int i = 0;i<4;i++)
		val[i] = rhs.val[i];
	return *this;
}

double& NIE_TYPE::nie_point_3d_t::operator [](const int i)
{
	return val[i];
}

#endif