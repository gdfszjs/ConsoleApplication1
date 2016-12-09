/*****
	提供对数据结构的c++定义
	2008/3/20

	修改于2008/3/25
	增加了Matrix3D的定义
*****/
#ifndef _DATA_STRUCTURE_H_
#define _DATA_STRUCTURE_H_

namespace NIE_DATASTRUCTURE
{
/****
二维矩阵的定义
****/
template <class T>
class Matrix2D
{
private:
	int nrows;
	int ncols;
	T **v;
public:
	Matrix2D();
	Matrix2D(int r,int c);
	Matrix2D(const T & a,int r,int c);
	Matrix2D(const T * a,int r,int c);
	Matrix2D(const Matrix2D & rhs);
	Matrix2D & operator=(const Matrix2D & rhs);
	Matrix2D & operator=(const T & a);
	inline T * operator[](const int i);
	inline const T * operator[](const int i) const;
	inline int getnrows() const;
	inline int getncols() const;
	~Matrix2D();
};

/*****
三维矩阵的定义
*****/
template <class T>
class Matrix3D
{
private:
	int nn;
	int mm;
	int kk;
	T *** v;
public:
	Matrix3D();
	Matrix3D(int n,int m,int k);
	Matrix3D(const T & a ,int n,int m,int k);
	inline T** operator[](const int i);
	inline const T* const * operator[](const int i)const;
	inline int dim1()const;
	inline int dim2()const;
	inline int dim3()const;
	~Matrix3D();
};


}

template <class T>
NIE_DATASTRUCTURE::Matrix2D<T>::Matrix2D() : nrows(0),ncols(0),v(0){}

template <class T>
NIE_DATASTRUCTURE::Matrix2D<T>::Matrix2D(int r,int c) : nrows(r),ncols(c),v(new T*[r])
{
	v[0] = new T[r*c];
	for(int i = 1;i<r;i++)
		v[i] = v[i-1] + c;
}

template <class T>
NIE_DATASTRUCTURE::Matrix2D<T>::Matrix2D(const T & a ,int r,int c) : nrows(r),ncols(c),v(new T*[r])
{
	v[0] = new T[r*c];
	for(int i = 1;i<r;i++)
		v[i] = v[i-1] + c;
	for(int i = 0;i<r;i++)
		for(int j =0;j<c;j++)
			v[i][j] = a;
}

template <class T>
NIE_DATASTRUCTURE::Matrix2D<T>::Matrix2D(const T *a, int r, int c) : nrows(r),ncols(c),v(new T*[r])
{
	v[0] = new T[r*c];
	for(int i = 1;i<r;i++)
		v[i] = v[i-1] + c;
	for(int i = 0;i<r;i++)
		for(int j = 0;j<c;j++)
			v[i][j] = *a++;
}

template <class T>
NIE_DATASTRUCTURE::Matrix2D<T>::Matrix2D(const NIE_DATASTRUCTURE::Matrix2D<T> &rhs):nrows(rhs.nrows),
ncols(rhs.ncols),v(new T*[rhs.nrows])
{
	v[0] = new T[nrows * ncols];
	for(int i = 1;i<nrows;i++)
		v[i] = v[i-1] + ncols;
	for(int i = 0;i<nrows;i++)
		for(int j = 0;j<ncols;j++)
			v[i][j] = rhs[i][j];
}


template <class T>
NIE_DATASTRUCTURE::Matrix2D<T> &
NIE_DATASTRUCTURE::Matrix2D<T>::operator =(const NIE_DATASTRUCTURE::Matrix2D<T> &rhs)
{
	if(this != &rhs)
	{
		if(nrows != rhs.nrows || ncols != rhs.ncols)
		{
			if(v!=0)
			{
				delete [] ( v[0] );
				delete [] v;
			}
			nrows = rhs.nrows;
			ncols = rhs.ncols;
			
			v = new T*[nrows];
			v[0] = new T[ nrows * ncols ];
			for(int i = 1;i<nrows;i++)
				v[i] = v[i-1] + ncols;
		}

		for(int i = 0 ;i< nrows;i++)
			for(int j = 0;j<ncols;j++)
				v[i][j] = rhs[i][j];
	}
	return *this;
}


template <class T>
NIE_DATASTRUCTURE::Matrix2D<T> &
NIE_DATASTRUCTURE::Matrix2D<T>::operator =(const T &a)
{
	for(int i = 0;i<nrows;i++)
		for(int j = 0;j<ncols;j++)
			v[i][j] = a;
	return *this;
}


template <class T>
inline T * NIE_DATASTRUCTURE::Matrix2D<T>::operator [](const int i)
{
	return v[i];
}

template <class T>
inline const T * NIE_DATASTRUCTURE::Matrix2D<T>::operator [](const int i) const
{
	return v[i];
}

template <class T>
inline int NIE_DATASTRUCTURE::Matrix2D<T>::getnrows() const
{
	return nrows;
}

template <class T>
inline int NIE_DATASTRUCTURE::Matrix2D<T>::getncols() const
{
	return ncols;
}

template <class T>
NIE_DATASTRUCTURE::Matrix2D<T>::~Matrix2D()
{
	if(v!=0)
	{
		delete [] (v[0]);
		delete [] v;
	}
}

template <class T>
NIE_DATASTRUCTURE::Matrix3D<T>::Matrix3D() : nn(0),mm(0),kk(0){}

template <class T>
NIE_DATASTRUCTURE::Matrix3D<T>::Matrix3D(int n, int m, int k) : nn(n),mm(m),kk(k),v(new T**[n])
{
	v[0] = new T*[n*m];
	v[0][0] = new T[n*m*k];
	for(int j = 1;j<m;j++)
		v[0][j] = v[0][j-1] + k;
	for(int i = 1;i<n;i++){
		v[i] = v[i-1] + m;
		v[i][0] = v[i-1][0] + m*k;
		for(int j = 1;j<m;j++)
			v[i][j] = v[i][j-1] + k;
	}
}


template <class T>
NIE_DATASTRUCTURE::Matrix3D<T>::Matrix3D(const T &a, int n, int m, int k) : nn(n),mm(m),kk(k),v(new T**[n])
{
	v[0] = new T*[n*m];
	v[0][0] = new T[n*m*k];
	for(int j = 1;j<m;j++)
		v[0][j] = v[0][j-1] + k;
	for(int i = 1;i<n;i++){
		v[i] = v[i-1] + m;
		v[i][0] = v[i-1][0] + m*k;
		for(int j = 1;j<m;j++)
			v[i][j] = v[i][j-1] + k;
	}

	for(int z = 0;z<n;z++)
		for(int y = 0;y<m;y++)
			for(int x = 0;x<k;x++)
				v[z][y][x] = a;
}

template <class T>
inline T** NIE_DATASTRUCTURE::Matrix3D<T>::operator [](const int i)
{
	return v[i];
}

template <class T>
inline const T * const * NIE_DATASTRUCTURE::Matrix3D<T>::operator [](const int i) const
{
	return v[i];
}

template <class T>
inline int NIE_DATASTRUCTURE::Matrix3D<T>::dim1() const
{
	return nn;
}

template <class T>
inline int NIE_DATASTRUCTURE::Matrix3D<T>::dim2() const
{
	return mm;
}

template <class T>
inline int NIE_DATASTRUCTURE::Matrix3D<T>::dim3() const 
{
	return kk;
}

template <class T>
NIE_DATASTRUCTURE::Matrix3D<T>::~Matrix3D()
{
	if(v!=0)
	{
		delete [] ( v[0][0] );
		delete [] ( v[0] );
		delete [] ( v );
	}
}

#endif /* _DATA_STRUCTURE_H_ */