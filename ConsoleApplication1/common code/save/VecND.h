//vecND.h
//2009.4.20

#pragma once
#include <cassert>

template <class T, size_t N> class VecND
{
private:
	T val[N];
public:
	VecND(){ for(int i = 0;i<N;i++) val[N] = 0; }
	VecND(T * v)
	{
		for(int i = 0;i<N;i++)
			val[i] = v[i];
	}
	VecND(const VecND<T,N>& v)
	{
		for(int i = 0;i<N;i++)
			val[i] = v.val[i];
	}

	VecND<T,N>& operator=(const VecND<T,N>& v)
	{
		for(int i = 0;i<N;i++)
			val[i] = v.val[i];
		return *this;
	}
	//////////////////////////////////////////////////////////////////
	//operator+,-,
	friend VecND<T,N> operator+(const VecND<T,N>& v1,const VecND<T,N>& v2)
	{
		VecND<T,N> v;
		for(int i = 0;i<N;i++)
			v.val[i] = v1.val[i] + v2.val[i];
		return v;
	}

	friend VecND<T,N> operator-(const VecND<T,N>& v1,const VecND<T>& v2)
	{
		VecND<T,N> v;
		for(int i = 0;i<N;i++)
			v.val[i] = v1.val[i] - v2.val[i];
		return v;
	}
	////////////////////////////////////////////////////////////////////
	friend VecND<T,N> operator +(const VecND<T,N>& v ,T t)
	{
		VecND<T,N> temp;
		for(int i = 0 ;i<N;i++)
			temp.val[i] = v.val[i] - t;
		return temp;
	}

	friend VecND<T,N> operator +(T t,const VecND<T,N>& v)
	{
		VecND<T,N> temp;
		for(int i = 0;i<N;i++)
			temp.val[i] = t + v.val[i];
		return temp;
	}

	friend VecND<T,N> operator -(const VecND<T,N>& v ,T t)
	{
		VecND<T,N> temp;
		for(int i = 0;i<N;i++)
			temp.val[i] = v.val[i] - t;
		return temp;
	}

	friend VecND<T,N> operator -(T t,const VecND<T,N>& v)
	{
		VecND<T,N> temp;
		for(int i = 0;i<N;i++)
			temp.val[i] = t - v.val[i];
		return temp;
	}

	///////////////////////////////////////////////////////////////////
	VecND<T,N> operator +=(const VecND<T,N>& v)
	{
		for(int i = 0;i<N;i++)
			val[i] += v.val[i];
		return *this;
	}

	VecND<T,N> operator -=(const VecND<T,N>& v)
	{
		for(int i = 0;i<N;i++)
			val[i] -= v.val[i];
		return *this;
	}

	VecND<T,N> operator +=(T t)
	{
		for(int i = 0;i<N;i++)
			val[i] += t;
		return *this;
	}

	VecND<T,N> operator -=(T t)
	{
		for(int i = 0;i<N;i++)
			val[i] -= t;
		return *this;
	}


	//* , / 
	friend T operator*(const VecND<T,N>& v1,const VecND<T,N>& v2)
	{
		T r = 0;
		for(int i = 0;i<N;i++)
			r += v1.val[i] * v2.val[i];
		return r;
	}

	friend VecND<T,N> operator*(const VecND<T,N>& v1, T t)
	{
		VecND<T,N> v;
		for(int i = 0;i<N;i++)
			v.val[i] = v1.val[i] * t;
		return v;
	}

	friend VecND<T,N> operator*(T t,const VecND<T,N>& v1)
	{
		VecND<T,N> v;
		for(int i = 0;i<N;i++)
			v.val[i] = v1.val[i] * t;
		return v;

	}

	VecND<T,N>& operator*=(T t)
	{
		for(int i = 0;i<N;i++)
			val[i] *= t;
		return *this;
	}

	VecND<T,N>& operator /=(T t)
	{
		assert(t!=0);
		for(int i = 0;i<N;i++)
			val[i] /= t;
		return *this;
	}

	friend VecND<T,N> operator /(const VecND<T,N> & v1,T t)
	{
		VecND<T,N> v;
		for(int i = 0;i<N;i++)
			v.val[i] = v1.val[i] / t;
		return v;
	}

	//operator []
	T operator[](int elem)
	{
		assert(elem>=0 && elem <=2);
		return val[elem];
	}

};//class VecND
