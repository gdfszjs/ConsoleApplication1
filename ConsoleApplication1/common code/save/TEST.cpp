#include "stdafx.h"
#include "TEST.h"
#include <iostream>
#include <typeinfo>
using namespace std;
using namespace TestSuite;

void Test::do_test(bool cond, const std::string &lbl, const char *fname, long lineno)
{
	if(cond)
		succeed_();
	else
		do_fail(lbl,fname,lineno);
}

void Test::do_fail(const std::string &lbl, const char *fname, long lineno)
{
	++nFail;
	if(osptr)
	{
		*osptr<<typeid(*this).name()
			<<"failure: ("<<lbl<<") , "<<fname
			<<" ( line"<<lineno<<")"<<endl;
	}
}

long Test::report()const 
{
	if(osptr)
	{
		*osptr<<"Test \""<<typeid(*this).name()
			<<"\"\n   Passed  "<<nPass
			<<"\tFailed  "<<nFail
			<<endl;
	}
	return nFail;
}

