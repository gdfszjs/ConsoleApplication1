#pragma once

#include "cv.h"
#include <vector>
using std::vector;



/************************************************************************
class GMM (Gauss Mixture Model)
************************************************************************/
#define GMM_DIM 3
#define GMM_LAMDA 1
#define LAMDA_B 100

class GMM
{
public:

	struct M_Color{ float r, g, b; };

	GMM();
	~GMM();

	void	Clear();

	int		SetNumClusters(int num);
	int		GetNumClusters() const { return m_nNumClusters; }

	void	AddInstance(const float* pInstance);
	int		GetNumInstances() const { return (int)m_vInstances.size(); }
	
	M_Color GetInstance(int i) { assert(i < GetNumInstances()); return m_vInstances[i];}
	bool	Cluster();

	float	Rp(float color[]);

private:
	//////////////////////////////////////////////////////////////////////////
	// source data
	
	vector<M_Color>	m_vInstances;	// instances
	//////////////////////////////////////////////////////////////////////////
	// cluster results
	int		m_nNumClusters;
	int*	m_pIDs;				// IDs of each instance
	vector<double> vecWeight;       //The weight of every cluster
	vector<CvMat*> vecCenter;       //The center of every cluster
	vector<CvMat*> vecInvCov;          //The covariant matrix of every cluster.
	vector<double> vecDet;          //The det of every covariant matrix.
	double * m_pInvCovParam;
	double * m_pPowHalfDet;         //pow(det,0.5);

};