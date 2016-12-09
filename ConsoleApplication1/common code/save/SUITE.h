//TestSuite : SUITE.H
#pragma once

#include <vector>
#include <stdexcept>
#include "TEST.h"

using std::logic_error;
using std::vector;

#pragma warning(disable: 4290)

namespace TestSuite
{
	class TestSuiteError : public logic_error
	{
	public:
		TestSuiteError(const string & msg = ""):logic_error(msg){}
	};//class TestSuiteError

	class Suite
	{
		string name;
		ostream * osptr;
		vector<Test*> tests;
		void reset();

		Suite(const Suite&);
		Suite& operator=(const Suite&);

	public:
		Suite(const string & name, ostream * osptr = &cout):name(name)
		{
			this->osptr = osptr;
		}

		string getName() const { return name; }
		long getNumPassed() const;
		long getNumFailed() const;

		const ostream* getStream() const { return osptr; }
		void setStream(ostream * osptr) { this->osptr = osptr; }

		void addTest(Test*) throw(TestSuiteError);
		void addSuite(const Suite&);

		void run();
		long report() const;
		void free();
	};//class Suite
}//namespace TestSuite


