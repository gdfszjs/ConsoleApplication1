#pragma once

class GS
{
public:
	GS(void);
	GS(double *&sa_,unsigned long *&ija_,double * &x_, double *&b_,const int n_); 
public:
	double *x;
	double *b;
	double *sa;
	unsigned long *ija;
	int n;
public:
	~GS(void);
public:
	void Do_GS(int max,double e);
};
