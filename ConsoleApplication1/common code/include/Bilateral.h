#pragma once

#include <cv.h>
#include <cxcore.h>
#include <highgui.h>
#include <cmath>
#include <iostream>
using namespace std;

namespace Bilateral
{
	static bool 
	ColourImageBilateral8U(const IplImage * srcImage,
				IplImage *& baseImage,
				int coreWidth,				
				double spaceSigma,
				double rangeSigma
				);
	
	static bool 
		ColourImageBilateral32F(const IplImage * srcImage,
		IplImage *& baseImage,
		int coreWidth,				
		float spaceSigma,
		float rangeSigma
		);

	static bool 
		GrayImageBilateral32F(const IplImage * srcImage,
		IplImage *& baseImage,
		int coreWidth,
		float spaceSigma,
		float rangeSigma
		);

	//8u
	static bool
		GrayImageBilateral(const IplImage * srcImage,
		IplImage *& baseImage,
		int coreWidth,				
		double spaceSigma,
		double rangeSigma
		);

	/*
	###########################################################
	###########################################################
	###########################################################
	##############                               ##############
	##############  I M P L E M E N T A T I O N  ##############
	##############                               ##############
	###########################################################
	###########################################################
	###########################################################
	*/


	static bool 
		ColourImageBilateral8U(const IplImage * srcImage,
		IplImage *& baseImage,
		int coreWidth,				
		double spaceSigma,
		double rangeSigma
		)
	{
		if(!srcImage)
			return false;

		int width = srcImage->width;
		int height = srcImage->height;
		int channels = srcImage->nChannels;
		int depth = srcImage->depth;

		/*if(channels != 3)
			return false;*/
		if(depth != 8)
			return false;


		IplImage * iplCIELab	=	cvCreateImage(cvSize(width,height),depth,channels);		
		cvCvtColor(srcImage,iplCIELab,CV_BGR2Lab);
		baseImage				=	cvCreateImage(cvSize(width,height),depth,channels);		

		double * spaceDist = new double[coreWidth * coreWidth];								

		uchar uOut[3] = { 0,0,0 };

		for(int i = 0;i<coreWidth;i++)
			for(int j = 0;j<coreWidth;j++)
				spaceDist[i * coreWidth + j] = sqrt((double)((i - coreWidth/2)*(i-coreWidth/2) + (j - coreWidth/2)*(j - coreWidth/2)));

		for(int i = 0;i<height;i++)
			for(int j = 0;j<width;j++)
			{
				uchar * pCIELabCore		=	&((uchar*)(iplCIELab->imageData + iplCIELab->widthStep * i))[j * channels];
				uchar * pSrcImageCore	=	&((uchar*)(srcImage->imageData + srcImage->widthStep * i))[j * channels];

				CvScalar cs;
				cs.val[0] = cs.val[1] = cs.val[2] = cs.val[3] = 0.0;

				for(int ii = -coreWidth/2;ii<=coreWidth/2;ii++)
					for(int jj = -coreWidth/2;jj<=coreWidth/2;jj++)
					{
						uchar * pCIELab;
						uchar * pSrcImage;
						double rangeDist;
						double dIntegral;

						if(i + ii >= 0 && i + ii < height && j + jj >= 0 && j + jj < width)
						{
							pCIELab		=	&((uchar*)(iplCIELab->imageData + iplCIELab->widthStep * (i+ii)))[(j+jj) * channels];
							pSrcImage	=	&((uchar*)(srcImage->imageData + srcImage->widthStep * (i+ii)))[(j+jj) * channels];

							rangeDist = sqrt((double)((pCIELab[0] - pCIELabCore[0])*(pCIELab[0] - pCIELabCore[0])
								+ (pCIELab[1] - pCIELabCore[1]) * (pCIELab[1] - pCIELabCore[1])
								+ (pCIELab[2] - pCIELabCore[2]) * (pCIELab[2] - pCIELabCore[2])));

							dIntegral = exp(-0.5 * ( spaceDist[(ii + coreWidth /2) * coreWidth +  jj + coreWidth/2] / spaceSigma )
								* ( spaceDist[(ii + coreWidth /2) * coreWidth +  jj + coreWidth/2] / spaceSigma ))
								* exp(-0.5 * ( rangeDist / rangeSigma ) * ( rangeDist / rangeSigma ));
						}
						else 
						{
							pCIELab = uOut;
							pSrcImage = uOut;
							rangeDist = 0;
							dIntegral = 0;
						}

						cs.val[3] += dIntegral;

						for(int k = 0;k<channels;k++)
							cs.val[k] += pSrcImage[k] * dIntegral;
					}

					for(int k = 0;k<channels;k++)
						cs.val[k] /= cs.val[3];
					cvSet2D(baseImage,i,j,cs);
			}

			delete[] spaceDist;
			return true;
	}

	static bool 
		ColourImageBilateral32F(const IplImage * srcImage,
		IplImage *& baseImage,
		int coreWidth,				
		float spaceSigma,
		float rangeSigma
		)
	{
		if(!srcImage)
			return false;

		int width = srcImage->width;
		int height = srcImage->height;
		int channels = srcImage->nChannels;
		int depth = srcImage->depth;

		/*if(channels != 3)
			return false;*/
		if(depth != IPL_DEPTH_32F)
			return false;


		IplImage * iplCIELab	=	cvCreateImage(cvSize(width,height),depth,channels);		
		cvCvtColor(srcImage,iplCIELab,CV_BGR2Lab);
		baseImage				=	cvCreateImage(cvSize(width,height),depth,channels);		

		float * spaceDist = new float[coreWidth * coreWidth];								

		float fOut[3] = { 0.0f,0.0f,0.0f };

		for(int i = 0;i<coreWidth;i++)
			for(int j = 0;j<coreWidth;j++)
				spaceDist[i * coreWidth + j] = sqrt((float)((i - coreWidth/2)*(i-coreWidth/2) + (j - coreWidth/2)*(j - coreWidth/2)));

		for(int i = 0;i<height;i++)
			for(int j = 0;j<width;j++)
			{
				float * pCIELabCore		=	&((float*)(iplCIELab->imageData + iplCIELab->widthStep * i))[j * channels];
				float * pSrcImageCore	=	&((float*)(srcImage->imageData + srcImage->widthStep * i))[j * channels];

				CvScalar cs;
				cs.val[0] = cs.val[1] = cs.val[2] = cs.val[3] = 0.0;

				for(int ii = -coreWidth/2;ii<=coreWidth/2;ii++)
					for(int jj = -coreWidth/2;jj<=coreWidth/2;jj++)
					{
						float * pCIELab;
						float * pSrcImage;
						float rangeDist;
						float dIntegral;

						if(i + ii >= 0 && i + ii < height && j + jj >= 0 && j + jj < width)
						{
							pCIELab		=	&((float*)(iplCIELab->imageData + iplCIELab->widthStep * (i+ii)))[(j+jj) * channels];
							pSrcImage	=	&((float*)(srcImage->imageData + srcImage->widthStep * (i+ii)))[(j+jj) * channels];

							rangeDist = sqrt((float)((pCIELab[0] - pCIELabCore[0])*(pCIELab[0] - pCIELabCore[0])
								+ (pCIELab[1] - pCIELabCore[1]) * (pCIELab[1] - pCIELabCore[1])
								+ (pCIELab[2] - pCIELabCore[2]) * (pCIELab[2] - pCIELabCore[2])));

							dIntegral = exp(-0.5f * ( spaceDist[(ii + coreWidth /2) * coreWidth +  jj + coreWidth/2] / spaceSigma )
								* ( spaceDist[(ii + coreWidth /2) * coreWidth +  jj + coreWidth/2] / spaceSigma ))
								* exp(-0.5f * ( rangeDist / rangeSigma ) * ( rangeDist / rangeSigma ));
						}
						else 
						{
							pCIELab = fOut;
							pSrcImage = fOut;
							rangeDist = 0;
							dIntegral = 0;
						}

						cs.val[3] += dIntegral;

						for(int k = 0;k<channels;k++)
							cs.val[k] += pSrcImage[k] * dIntegral;
					}

					for(int k = 0;k<channels;k++)
						cs.val[k] /= cs.val[3];
					cvSet2D(baseImage,i,j,cs);
			}

			delete[] spaceDist;
			return true;
	}

	static bool GrayImageBilateral32F(const IplImage * srcImage, IplImage *& baseImage, int coreWidth, float spaceSigma, float rangeSigma )
	{
		if(!srcImage)
			return false;
		int w = srcImage->width;
		int h = srcImage->height;
		int nc = srcImage->nChannels;
		int d = srcImage->depth;
		if(nc != 1 || d != IPL_DEPTH_32F)
			return false;

		baseImage = cvCreateImage(cvSize(w,h),d,nc);

		int r = coreWidth /2;
		float sS =-0.5f/(spaceSigma * spaceSigma);
		float sR =-0.5f/(rangeSigma * rangeSigma);

		for(int i = r;i<h-r;i++)
		{
			for(int j = r;j<w-r;j++)
			{
				float * pC = &((float*)(srcImage->imageData + srcImage->widthStep * i))[j];
				float up = 0.0;
				float down = 0.0;

				for(int ii = i - r;ii<=i+r;ii++)
				{
					for(int jj = j - r;jj<=j+r;jj++)
					{
						float * p = &((float*)(srcImage->imageData + srcImage->widthStep * ii))[jj];
						float temp = exp(sR * (p[0] - pC[0]) * (p[0] - pC[0]) + sS * ((ii-i)* (ii-i)+(jj-j)*(jj-j)));
						down += temp;
						up += temp * p[0];
					}
				}
				((float*)(baseImage->imageData + baseImage->widthStep*i))[j] =(float)(up / down);
			}
		}
		return true;




		/*if(!srcImage)
			return false;
		int w = srcImage->width;
		int h = srcImage->height;
		int nc = srcImage->nChannels;
		int d = srcImage->depth;
		if(nc != 1 || d != IPL_DEPTH_8U)
			return false;

		baseImage = cvCreateImage(cvSize(w,h),d,nc);

		int r = coreWidth /2;
		float sS =-0.5f/(spaceSigma * spaceSigma);
		float sR =-0.5f/(rangeSigma * rangeSigma);

		for(int i = r;i<h-r;i++)
		{
			for(int j = r;j<w-r;j++)
			{
				uchar * pC = &((uchar*)(srcImage->imageData + srcImage->widthStep * i))[j];
				float up = 0.0;
				float down = 0.0;

				for(int ii = i - r;ii<=i+r;ii++)
				{
					for(int jj = j - r;jj<=j+r;jj++)
					{
						uchar * p = &((uchar*)(srcImage->imageData + srcImage->widthStep * ii))[jj];
						float temp = exp(sR * (p[0] - pC[0]) * (p[0] - pC[0]) + sS * ((ii-i)* (ii-i)+(jj-j)*(jj-j)));
						down += temp;
						up += temp * p[0];
					}
				}
				((uchar*)(baseImage->imageData + baseImage->widthStep*i))[j] =(uchar)(up / down);
			}
		}
		return true;*/
	}

	static bool 
		GrayImageBilateral(const IplImage * srcImage,
		IplImage *& baseImage,
		int coreWidth,				
		double spaceSigma,
		double rangeSigma
		)
	{
		if(!srcImage)
			return false;

		int width = srcImage->width;
		int height = srcImage->height;
		int channels = srcImage->nChannels;
		int depth = srcImage->depth;

		if(channels != 1)
			return false;
		if(depth != 8)
			return false;

		baseImage				=	cvCreateImage(cvSize(width,height),depth,channels);		

		double * spaceDist = new double[coreWidth * coreWidth];								

		uchar uOut[1] = { 0 };

		for(int i = 0;i<coreWidth;i++)
			for(int j = 0;j<coreWidth;j++)
				spaceDist[i * coreWidth + j] = sqrt((double)((i - coreWidth/2)*(i-coreWidth/2) + (j - coreWidth/2)*(j - coreWidth/2)));

		for(int i = 0;i<height;i++)
			for(int j = 0;j<width;j++)
			{
				uchar * pSrcImageCore	=	&((uchar*)(srcImage->imageData + srcImage->widthStep * i))[j * channels];

				CvScalar cs;
				cs.val[0] = cs.val[1] = cs.val[2] = cs.val[3] = 0.0;

				for(int ii = -coreWidth/2;ii<=coreWidth/2;ii++)
					for(int jj = -coreWidth/2;jj<=coreWidth/2;jj++)
					{
						uchar * pSrcImage;
						double rangeDist;
						double dIntegral;

						if(i + ii >= 0 && i + ii < height && j + jj >= 0 && j + jj < width)
						{
							pSrcImage	=	&((uchar*)(srcImage->imageData + srcImage->widthStep * (i+ii)))[(j+jj) * channels];

							rangeDist = fabs((double)(pSrcImage[0] - pSrcImageCore[0]));

							dIntegral = exp(-0.5 * ( spaceDist[(ii + coreWidth /2) * coreWidth +  jj + coreWidth/2] / spaceSigma )
								* ( spaceDist[(ii + coreWidth /2) * coreWidth +  jj + coreWidth/2] / spaceSigma ))
								* exp(-0.5 * ( rangeDist / rangeSigma ) * ( rangeDist / rangeSigma ));
						}
						else 
						{
							pSrcImage = uOut;
							rangeDist = 0;
							dIntegral = 0;
						}

						cs.val[3] += dIntegral;

						for(int k = 0;k<channels;k++)
							cs.val[k] += pSrcImage[k] * dIntegral;
					}

					for(int k = 0;k<channels;k++)
						cs.val[k] /= cs.val[3];
					cvSet2D(baseImage,i,j,cs);
			}

			delete[] spaceDist;
			return true;
	}

}
