#pragma once

#include "ANN.h"

class ANNCell{
public:
	explicit ANNCell(const int& cn,const int& pn);
	~ANNCell();
public:
	void readInCoordinates(const double & cor);
	void addPointToArray();
	void buildKd_tree();
	int findNearest(const ANNpoint &p,double *& distance);
	ANNpointArray& getPointArray(void);
	int getpn()
	{
		return pointNum;
	}
	int getcn()
	{
		return coordinateNum;
	}
private:
	bool isReady;
	int coordinateCount;
	int pointCount;
	int coordinateNum;
	int pointNum;
	ANNpoint commonPoint;
	ANNpointArray arrayOfPoint;
	ANNkd_tree * kdTree;
};
