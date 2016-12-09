#pragma once

//typedef real double;
//typedef int int; 
class GS
{
public:
	GS(void);
	GS(double *&sa_,int *&ija_,double * &x_, double *&b_,const int n_); 
public:
	double *x;
	double *b;
	double *sa;
	int *ija;
	int n;
public:
	~GS(void);
public:
	void Do_GS(int max,double e);
};
