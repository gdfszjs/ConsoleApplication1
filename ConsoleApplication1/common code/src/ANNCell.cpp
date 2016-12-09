//#include "stdafx.h"
#include "ANNCell.h"

ANNCell::ANNCell(const int& cn,const int& pn)
							:	isReady(false)
							,	coordinateNum(cn)
							,	pointNum(pn)
							,	commonPoint(NULL)
							,	arrayOfPoint(NULL)
							,	kdTree(NULL)
{
	commonPoint = annAllocPt(coordinateNum);
	arrayOfPoint = annAllocPts(pointNum,coordinateNum);
	pointCount = 0;
	coordinateCount = 0;
}

ANNCell::~ANNCell(){
	annDeallocPt(commonPoint);
	annDeallocPts(arrayOfPoint);
	if(kdTree){
		delete kdTree;
	}
	annClose();
}

void ANNCell::readInCoordinates(const double & cor){

	arrayOfPoint[pointCount][coordinateCount] = (ANNcoord)cor;
	
	if(++coordinateCount == coordinateNum){
		coordinateCount = 0;
		if(++pointCount == pointNum)
		{
			isReady = true;
		}
	}
}

void ANNCell::addPointToArray()
{
	
}

void ANNCell::buildKd_tree(){
	if(isReady){
		kdTree = new ANNkd_tree(arrayOfPoint,pointNum,coordinateNum);
	}
	else{
		std::cerr<<"BuildKd_tree() - error : Is not ready!"<<std::endl;
	}
}

int ANNCell::findNearest(const ANNpoint & p,double *& distance){
	if(isReady){
		ANNidxArray index = new ANNidx[1];
		kdTree->annkSearch(p,1,index,distance);
		return index[0];
	}
	else{
		std::cerr<<"findNearest() - error : Is not ready!"<<std::endl;
		return NULL;
	}
}

ANNpointArray& ANNCell::getPointArray(){
	return arrayOfPoint;
}

