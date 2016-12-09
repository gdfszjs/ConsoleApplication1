#include "stdafx.h"
#include "GMM.h"
#include <cmath>
#include <iostream>
using namespace std;

//#ifdef _DEBUG
//#undef THIS_FILE
//static char THIS_FILE[]=__FILE__;
//#define new DEBUG_NEW
//#endif

/************************************************************************
class GMM:	Gauss Mixture Model                                                                     
************************************************************************/

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

#define DEFAULT_CLUSTER_NUM 5
#define DEFAULT_DONUM 0
GMM::GMM()
{
	m_pIDs = NULL;
	vecWeight.clear();
	vecCenter.clear();
	vecInvCov.clear();
	vecDet.clear();
	m_pInvCovParam = NULL;
	m_pPowHalfDet = NULL;
	SetNumClusters(DEFAULT_CLUSTER_NUM);
}

GMM::~GMM()
{
	m_vInstances.clear();
	Clear();
}



int GMM::SetNumClusters(int num)
{
	Clear();

	
	m_nNumClusters = num;

	return num;
}

void GMM::Clear()
{
	if(m_pIDs){
		delete []m_pIDs;
		m_pIDs = NULL;
	}
	
	vecWeight.clear();
	vecDet.clear();

	int num = (int)vecCenter.size();
	//assert(num == m_nNumClusters);
	
	for(int i = 0;i<num;i++)
	{
		if(vecCenter[i])
			cvReleaseMat(&vecCenter[i]);
		if(vecInvCov[i])
			cvReleaseMat(&vecInvCov[i]);
	}
	vecCenter.clear();
	vecInvCov.clear();
	if(m_pInvCovParam)
		delete [] m_pInvCovParam;
	if(m_pPowHalfDet)
		delete [] m_pPowHalfDet;
}

void GMM::AddInstance(const float* pInstance)
{
	M_Color color;
	memcpy(&color, pInstance, sizeof(M_Color));
	m_vInstances.push_back(color);
}

bool GMM::Cluster()
{
	if(GetNumInstances()<m_nNumClusters) return false;

	Clear();

	int i;
	int iCluster;

	//////////////////////////////////////////////////////////////////////////
	// groups pixels using k-means algorithm implemented in OpenCV 


	// make matrix used in cvKMeans2 function
	CvMat* pSample = cvCreateMat(GetNumInstances(), GMM_DIM, CV_32FC1);                            //(1)pSample temporary
	
	for(i=0; i<GetNumInstances(); i++)
	{
		cvmSet(pSample, i, 0, m_vInstances[i].r);
		cvmSet(pSample, i, 1, m_vInstances[i].g);
		cvmSet(pSample, i, 2, m_vInstances[i].b);
	}

	// allocate space for m_pIDs
	m_pIDs = new int[GetNumInstances()];


	CvMat* pIDMat = cvCreateMat(GetNumInstances(), 1, CV_32SC1);                                   //(2)pIDMat

	cvKMeans2(pSample, m_nNumClusters, pIDMat, cvTermCriteria( CV_TERMCRIT_EPS+CV_TERMCRIT_ITER, 1000, 1e-10 ));

	for(i=0; i<GetNumInstances(); i++){
		m_pIDs[i] = pIDMat->data.i[i];
	}

	cvReleaseMat(&pIDMat);                                                                         //(2) pIDMat is released


	//////////////////////////////////////////////////////////////////////////
	// calculating attributes of GMM

	// statistics: count of instances for each cluster
	int *pNumInstOfCluster = new int[m_nNumClusters];                                              //(3) pNumInstOfCluster temporary
	memset(pNumInstOfCluster, 0, sizeof(int)*m_nNumClusters);
	for(i=0;i<GetNumInstances();i++)
		pNumInstOfCluster[m_pIDs[i]]++;

	
 	//
	//
	//Initializing ...
	//
	//

	vector<double> vecWeightT;
	vector<CvMat *> vecCenterT;
	vector<CvMat *> vecInvCovT;
	vector<double> vecDetT;

	vecWeightT.clear();
	vecCenterT.clear();
	vecInvCovT.clear();
	vecDetT.clear();
	
	for(int i = 0;i<m_nNumClusters;i++)
	{
		vecWeightT.push_back(0);
		vecCenterT.push_back(NULL);
		vecInvCovT.push_back(NULL);
		vecDetT.push_back(NULL);
	}


	m_pInvCovParam = new double[m_nNumClusters*6];
	m_pPowHalfDet = new double[m_nNumClusters];

	
	// set weights according to ratio of counts
	for(int i = 0;i<m_nNumClusters;i++)
		vecWeight.push_back((double)pNumInstOfCluster[i] / (double)GetNumInstances());

	

	// prepare matrices used in OpenCV function cvCalcCovarMatrix
	CvMat** ppVects = new CvMat*[m_nNumClusters];
	for(iCluster=0; iCluster<m_nNumClusters; iCluster++){
		ppVects[iCluster] = (CvMat*)cvAlloc(sizeof(CvMat)*pNumInstOfCluster[iCluster]);
	}
	CvMat*** pppVects = new CvMat**[m_nNumClusters];
	for(iCluster=0; iCluster<m_nNumClusters; iCluster++){
		pppVects[iCluster] = (CvMat**)cvAlloc(sizeof(CvMat*)*pNumInstOfCluster[iCluster]);
	}

	int	*pVectIndex = new int[m_nNumClusters];	// index of current vector for each cluster
	memset(pVectIndex, 0, sizeof(int)*m_nNumClusters);

	CvMat *pCurrentVectMatrix;
	int iIndex;
	for(i=0;i<GetNumInstances();i++)
	{
		iCluster = m_pIDs[i];
		iIndex = pVectIndex[iCluster]++;
		pCurrentVectMatrix = pppVects[iCluster][iIndex] = &ppVects[iCluster][iIndex];
		*pCurrentVectMatrix = cvMat(1, GMM_DIM, CV_32FC1, pSample->data.ptr+i*pSample->step);
	}
	
	// for each cluster, calculate the mean vector and covariant matrices
	for(iCluster=0; iCluster<m_nNumClusters; iCluster++)
	{
		CvMat* pInvCovMat = cvCreateMat(GMM_DIM, GMM_DIM, CV_32FC1);
		CvMat* pCenterMat = cvCreateMat(1, GMM_DIM, CV_32FC1);

		cvCalcCovarMatrix((const CvArr**)pppVects[iCluster], pNumInstOfCluster[iCluster],
			pInvCovMat, pCenterMat, CV_COVAR_NORMAL + CV_COVAR_SCALE);		
		
		/*for(int ii = 0;ii<GMM_DIM;ii++)
		{
			for(int jj = 0;jj <GMM_DIM; jj++)
			{
				cout<<cvmGet(pInvCovMat,ii,jj)<<",";
			}
			cout<<endl;
		}*/

		
		float det = (float)cvmInvert(pInvCovMat, pInvCovMat);
		
		vecInvCov.push_back(pInvCovMat);
		vecDet.push_back(det);

		vecCenter.push_back(pCenterMat);		
	}

	int DONUM = 0;

	//EM algorithm, to optimize GMM's parameters.

	CvMat * pProbability = cvCreateMat(m_nNumClusters,GetNumInstances(),CV_32FC1);
	CvMat * pPosteriorProbability = cvCreateMat(m_nNumClusters,GetNumInstances(),CV_32FC1);
	CvMat * matSample = cvCreateMat(1,3,CV_32FC1);
	CvMat * matDiff = cvCreateMat(1,3,CV_32FC1);
	CvMat * matTemp13 = cvCreateMat(1,3,CV_32FC1);
	CvMat * matTemp11 = cvCreateMat(1,1,CV_32FC1);
	CvMat * matTemp13_2 = cvCreateMat(1,3,CV_32FC1);
	CvMat * matTemp11_2 = cvCreateMat(1,1,CV_32FC1);
	CvMat * matTemp33 = cvCreateMat(3,3,CV_32FC1);
	CvMat * matTemp33_2 = cvCreateMat(3,3,CV_32FC1);
	CvMat * matTemp33_3 = cvCreateMat(3,3,CV_32FC1);
	CvMat * matTemp33_4 = cvCreateMat(3,3,CV_32FC1);

	while(DONUM++ < DEFAULT_DONUM)
	{
		//Compute the probability of every instance of every cluster
		for(int ii = 0;ii<m_nNumClusters;ii++)
		{	
			for(int jj = 0;jj<GetNumInstances();jj++)
			{
				/*((float*)(matSample->data.ptr))[0] = ((float*)(pSample->data.ptr + pSample->step*jj))[0];
				((float*)(matSample->data.ptr))[1] = ((float*)(pSample->data.ptr + pSample->step*jj))[1];
				((float*)(matSample->data.ptr))[2] = ((float*)(pSample->data.ptr + pSample->step*jj))[2];*/


				cvmSet(matSample,0,0,cvmGet(pSample,jj,0));
				cvmSet(matSample,0,1,cvmGet(pSample,jj,1));
				cvmSet(matSample,0,2,cvmGet(pSample,jj,2));

				cvSub(matSample,vecCenter[ii],matDiff);

				cvMatMul(matDiff,vecInvCov[ii],matTemp13);
				cvGEMM(matTemp13,matDiff,1,NULL,0,matTemp11,CV_GEMM_B_T);

				double dd = -0.5 * cvmGet(matTemp11,0,0); 

				double val =( 1.0/(pow(2 * 3.14159265,GMM_DIM/2.0) * pow(vecDet[ii],0.5))) * exp(-0.5 * cvmGet(matTemp11,0,0)); 
				cvmSet(pProbability,ii,jj,val);
		
				
			}
		}
		
		//Compute posterior probability of every instance
		
		for(int ii = 0;ii<m_nNumClusters;ii++)
		{
			for(int jj = 0;jj<GetNumInstances();jj++)
			{
				double ppUp = vecWeight[ii] * cvmGet(pProbability,ii,jj);
				double ppDown = 0.0;
				for(int iii = 0;iii<m_nNumClusters;iii++)
					ppDown += vecWeight[iii] * cvmGet(pProbability,iii,jj);
				cvmSet(pPosteriorProbability,ii,jj,ppUp/ppDown);
			}
		}
		

		//Compute new weight of every cluster
	
		for(int ii = 0;ii<m_nNumClusters;ii++)
		{
			for(int jj = 0;jj<GetNumInstances();jj++)
				vecWeightT[ii] += cvmGet(pPosteriorProbability,ii,jj);
			vecWeightT[ii] /= (double)GetNumInstances();
		}

		
		//Compute new center:
		for(int ii = 0;ii<m_nNumClusters;ii++)
		{
			double down = 0.0;
			CvMat * matRes = cvCreateMat(1,3,CV_32FC1);

			cvmSet(matTemp13,0,0,0.0);
			cvmSet(matTemp13,0,1,0.0);
			cvmSet(matTemp13,0,2,0.0);
			
			for(int jj = 0;jj<GetNumInstances();jj++)
			{
				cvmSet(matSample,0,0,cvmGet(pSample,jj,0));
				cvmSet(matSample,0,1,cvmGet(pSample,jj,1));
				cvmSet(matSample,0,2,cvmGet(pSample,jj,2));
				
				cvmSet(matTemp11,0,0,cvmGet(pPosteriorProbability,ii,jj));

				cvGEMM(matTemp11,matSample,1,matTemp13,1,matTemp13);

				down += cvmGet(pPosteriorProbability,ii,jj);				
			}

			cvmSet(matTemp11_2,0,0,1.0/down);

			cvMatMul(matTemp11_2 ,matTemp13,matRes);

			vecCenterT[ii] = matRes;

		}

		
		//Compute new covariant matrix:

		
		for(int ii = 0;ii<m_nNumClusters;ii++)
		{
			CvMat * matRes = cvCreateMat(3,3,CV_32FC1);
			double down = 0.0;

			for(int iii = 0;iii<3;iii++)
				for(int jjj = 0;jjj<3;jjj++)
					cvmSet(matTemp33,iii,jjj,0.0);

			for(int jj = 0;jj<GetNumInstances();jj++)
			{
				cvmSet(matSample,0,0,cvmGet(pSample,jj,0));
				cvmSet(matSample,0,1,cvmGet(pSample,jj,1));
				cvmSet(matSample,0,2,cvmGet(pSample,jj,2));

				cvSub(matSample,vecCenter[ii],matDiff);

				for(int iii = 0;iii<3;iii++)
					for(int jjj = 0;jjj<3;jjj++)
					{
						if(iii != jjj)
							cvmSet(matTemp33_3,iii,jjj,0.0);
						else
							cvmSet(matTemp33_3,iii,jjj,cvmGet(pPosteriorProbability,ii,jj));
					}


				cvGEMM(matDiff,matDiff,1,NULL,0,matTemp33_2,CV_GEMM_A_T);

				cvGEMM( matTemp33_3,matTemp33_2,1,matTemp33,1,matTemp33);

				down += cvmGet(pPosteriorProbability,ii,jj);
			}

			for(int iii = 0;iii<3;iii++)
				for(int jjj = 0;jjj<3;jjj++)
				{
					if(iii != jjj)
						cvmSet(matTemp33_4,iii,jjj,0.0);
					else
						cvmSet(matTemp33_4,iii,jjj,1.0/down);
				}

			cvMatMul( matTemp33_4,matTemp33,matRes);
			double det = cvInvert(matRes,matRes);
			vecInvCovT[ii] = matRes;
			vecDetT[ii] = det;
		}

		//Delete old
		for(int ii = 0;ii<m_nNumClusters;ii++)
		{
			cvReleaseMat(&vecCenter[ii]);
			cvReleaseMat(&vecInvCov[ii]);
		}

		//Exchange the parameters and reset new:
		for(int ii = 0;ii<m_nNumClusters;ii++)
		{
			vecWeight[ii] = vecWeightT[ii];
			vecDet[ii] = vecDetT[ii];
			vecCenter[ii] = vecCenterT[ii];
			vecInvCov[ii] = vecInvCovT[ii];
			vecWeightT[ii] = 0.0;
			vecDetT[ii] = 0.0;
			vecCenterT[ii] = NULL;
			vecInvCovT[ii] = NULL;
		}

 	}
	
	cvReleaseMat(&matTemp33_4);
	cvReleaseMat(&matTemp33_3);
	cvReleaseMat(&matTemp33_2);
	cvReleaseMat(&matTemp33);

	cvReleaseMat(&matTemp11_2);
	cvReleaseMat(&matTemp13_2);
	cvReleaseMat(&matTemp11);
	cvReleaseMat(&matTemp13);

	cvReleaseMat(&matDiff);
	cvReleaseMat(&matSample);
	cvReleaseMat(&pPosteriorProbability);
	cvReleaseMat(&pProbability);
	
	for(iCluster = 0;iCluster < m_nNumClusters;iCluster++)
	{
		double * pInvCovParam = &m_pInvCovParam[iCluster*6];
		pInvCovParam[0] = 0.5 * cvmGet(vecInvCov[iCluster], 0, 0);
		pInvCovParam[1] = 0.5 * ( cvmGet(vecInvCov[iCluster], 0, 1) + cvmGet(vecInvCov[iCluster], 1, 0) );
		pInvCovParam[2] = 0.5 * ( cvmGet(vecInvCov[iCluster], 0, 2) + cvmGet(vecInvCov[iCluster], 2, 0) );
		pInvCovParam[3] = 0.5 * cvmGet(vecInvCov[iCluster], 1, 1);
		pInvCovParam[4] = 0.5 * ( cvmGet(vecInvCov[iCluster], 1, 2) + cvmGet(vecInvCov[iCluster], 2, 1) );
		pInvCovParam[5] = 0.5 * cvmGet(vecInvCov[iCluster], 2, 2);

		m_pPowHalfDet[iCluster] = pow(vecDet[iCluster],0.5);
		//m_pPowHalfDet[iCluster] = 0.5 * log(vecDet[iCluster]);

	}

 	delete []pNumInstOfCluster;
	for(iCluster=0; iCluster<m_nNumClusters; iCluster++){
		cvFree((void**)&pppVects[iCluster]);
	}
	delete []pppVects;

	for(iCluster=0; iCluster<m_nNumClusters; iCluster++){
		cvFree((void**)&ppVects[iCluster]);
	}
	delete []ppVects;

	cvReleaseMat(&pSample);

	return true;

}

float GMM::Rp(float color[])
{
	float pDiff[3];
	//double res = 0.0;
	double minCost = HUGE;
	for(int i = 0;i<m_nNumClusters;i++)
	{
		for(int ii = 0;ii<3;ii++)
			pDiff[ii] = color[ii] - cvmGet(vecCenter[i],0,ii);

		double *pInvCov = &m_pInvCovParam[i*6];

		double xTVx = pInvCov[0] * pDiff[0] * pDiff[0] + // 0.5 * a11 x^2
			pInvCov[1] * pDiff[0] * pDiff[1] +				// 0.5 * (a12+a21)xy
			pInvCov[2] * pDiff[0] * pDiff[2] +				// 0.5 * (a13+a31)xz
			pInvCov[3] * pDiff[1] * pDiff[1] +				// 0.5 * a22 y^2
			pInvCov[4] * pDiff[1] * pDiff[2] +				// 0.5 * (a23+a32)yz
			pInvCov[5] * pDiff[2] * pDiff[2];				// 0.5 *a33 z^2
	
		/*double val =1.0/ m_pPowHalfDet[i] * exp(-xTVx); 
		res += val * vecWeight[i];*/
		double cost = xTVx + m_pPowHalfDet[i] - log(vecWeight[i]);
		minCost = cost < minCost? cost : minCost;

	}
	return minCost;
	//return -log(res);
	//return res;
}