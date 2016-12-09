//vec3D.h
#pragma once
#include <cassert>

template <class T> class Vec3D
{
	T x,y,z;
public:
	Vec3D(){ x = y = z = 0;}
	Vec3D(T x,T y,T z) : x(x),y(y),z(z) {}
	Vec3D(const Vec3D<T>& v)
	{
		x = v.x;
		y = v.y;
		z = v.z;
	}

	Vec3D<T>& operator=(const Vec3D<T>& v)
	{
		x = v.x;
		y = v.y;
		z = v.z;
		return *this;
	}
	//////////////////////////////////////////////////////////////////
	//operator+,-,
	friend Vec3D<T> operator+(const Vec3D<T>& v1,const Vec3D<T>& v2)
	{
		Vec3D<T> v;
		v.x = v1.x + v2.x;
		v.y = v1.y + v2.y;
		v.z = v1.z + v2.z;
		return v;
	}

	friend Vec3D<T> operator-(const Vec3D<T>& v1,const Vec3D<T>& v2)
	{
		Vec3D<T> v;
		v.x = v1.x - v2.x;
		v.y = v1.y - v2.y;
		v.z = v1.z - v2.z;
		return v;
	}
	////////////////////////////////////////////////////////////////////
	friend Vec3D<T> operator +(const Vec3D<T>& v ,T t)
	{
		Vec3D<T> temp;
		temp.x = v.x + t;
		temp.y = v.y + t;
		temp.z = v.z + t;
		return temp;
	}

	friend Vec3D<T> operator +(T t,const Vec3D<T>& v)
	{
		Vec3D<T> temp;
		temp.x = v.x + t;
		temp.y = v.y + t;
		temp.z = v.z + t;
		return temp;
	}

	friend Vec3D<T> operator -(const Vec3D<T>& v ,T t)
	{
		Vec3D<T> temp;
		temp.x = v.x - t;
		temp.y = v.y - t;
		temp.z = v.z - t;
		return temp;
	}

	friend Vec3D<T> operator -(T t,const Vec3D<T>& v)
	{
		Vec3D<T> temp;
		temp.x = t - v.x;
		temp.y = t - v.y;
		temp.z = t - v.z;
		return temp;
	}

	///////////////////////////////////////////////////////////////////
	Vec3D<T> operator +=(const Vec3D<T>& v)
	{
		x += v.x;
		y += v.y;
		z += v.z;
		return *this;
	}
	
	Vec3D<T> operator -=(const Vec3D<T>& v)
	{
		x-=v.x;
		y-=v.y;
		z-=v.z;
		return *this;
	}

	Vec3D<T> operator +=(T t)
	{
		x += t;
		y += t;
		z += t;
		return *this;
	}
	
	Vec3D<T> operator -=(T t)
	{
		x-=t;
		y-=t;
		z-=t;
		return *this;
	}


	//* , / 
	friend T operator*(const Vec3D<T>& v1,const Vec3D<T>& v2)
	{
		return (v1.x*v2.x + v1.y * v2.y + v1.z * v2.z);
	}

	friend Vec3D<T> operator*(const Vec3D<T>& v1, T t)
	{
		Vec3D<T> v;
		v.x = v1.x * t;
		v.y = v1.y * t;
		v.z = v1.z * t;
		return v;
	}

	friend Vec3D<T> operator*(T t,const Vec3D<T>& v1)
	{
		Vec3D<T> v;
		v.x = v1.x * t;
		v.y = v1.y * t;
		v.z = v1.z * t;
		return v;
	
	}

	Vec3D<T>& operator*=(T t)
	{
		x *= t;
		y *= t;
		z *= t;
		return *this;
	}

	Vec3D<T>& operator /=(T t)
	{
		assert(t!=0);
		x /= t;
		y /= t;
		z /= t;
		return *this;
	}

	friend Vec3D<T> operator /(const Vec3D<T> & v1,T t)
	{
		Vec3D<T> v;
		v.x = v1.x/t;
		v.y = v1.y/t;
		v.z = v1.z/t;
		return v;
	}

	//operator []
	T operator[](int elem)
	{
		assert(elem>=0 && elem <=2);
		if(elem == 0)
			return x;
		else if(elem == 1)
			return y;
		else
			return z;
	}

};//class Vec3D
