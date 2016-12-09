#pragma once

#include "common.h"
#include <vector>
#include <cmath>

#include "Bilateral.h"
#include "NMessage.h"

//using namespace NMessage;
//using namespace Bilateral;

//#define __DEBUG__

#ifdef __DEBUG__

#include <iostream>
using std::cout;
using std::endl;

#endif

namespace ImageDecompose
{
	// Decompose image into
	// base + detail1 + detail2 + ...
	class ImageDecompose
	{
	private:
		IplImage * src;
		IplImage * base;
		IplImage * log_base;
		vector<IplImage*> detail;
		vector<IplImage*> vecBase;
		vector<IplImage*> log_detail;
		size_t level;
		size_t coreWidth;
	    double contrast;

	public:
		ImageDecompose();
		ImageDecompose(IplImage * s,size_t l,size_t cw,double co);
		~ImageDecompose();
				
		void Decompose();
		void Save();
		IplImage* Enhance(double * w,int n);
		IplImage * GetBase(int n)
		{
			ASSERT(n < (int)vecBase.size());
			return vecBase[n];
		}
	protected:
		template<typename T> inline typename T nMax(const T a, const T b);
		template<typename T> inline typename T nMin(const T a, const T b);
		template<typename T> inline typename T log10(const T x);

		void imageLog10(const IplImage * s, IplImage *& d);
		void imageExp10(const IplImage * s, IplImage *& d,int depth);

		void clear();

		int intensityRange(const IplImage * img);
		int rgbRange(const IplImage* img,int rgb,int & maxV);

		float * spaceSigma(float s);
		float * rangeSigma(float r);

		void genDetail(const IplImage * s,const IplImage * b , IplImage *& d);
		void genDetailDis(const IplImage * s , IplImage *& d);
	};

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
	
	ImageDecompose::ImageDecompose()
	{
		src = NULL;
		base = NULL;
		log_base = NULL;
		detail.clear();
		vecBase.clear();
		level = 0;
		coreWidth = 0;
		contrast = 0.0;
	}

	ImageDecompose::ImageDecompose(IplImage * s,size_t l,size_t cw,double co) : src(s),level(l),coreWidth(cw),contrast(co)
	{
		base = NULL;
		log_base = NULL;
		detail.clear();
		vecBase.clear();
		log_detail.clear();
	}

	ImageDecompose::~ImageDecompose()
	{
		/*if(src)
			cvReleaseImage(&src);
		src = NULL;*/
		clear();
	}

	void ImageDecompose::clear()
	{
		if(base)
			cvReleaseImage(&base);
		
		base = NULL;

		if(log_base)
			cvReleaseImage(&log_base);
		log_base = NULL;

		for(int i = 0;i<(int)detail.size();i++)
		{
			if(detail[i])
				cvReleaseImage(&detail[i]);
			if(log_detail[i])
				cvReleaseImage(&log_detail[i]);
			if(vecBase[i])
				cvReleaseImage(&vecBase[i]);
		}
		
		detail.clear();
		vecBase.clear();
		log_detail.clear();
	}

	template<typename T>
	inline T
		ImageDecompose::nMax(const T a, const T b) { return a>b?a:b; }

	template<typename T>
	inline T
		ImageDecompose::nMin(const T a, const T b) { return a<b?a:b; }

	template<typename T>
	inline T
		ImageDecompose::log10(const T x)
	{
		static const double inv_log_base=1.0/log(10.0);

		return (T)(log((double)x)*inv_log_base);
	}

	void ImageDecompose::imageLog10(const IplImage * s,IplImage *& d)
	{
		if(!d) cvReleaseImage(&d);

		//Float point type after log operation.
		d = cvCreateImage(cvSize(s->width,s->height),IPL_DEPTH_32F,s->nChannels);

		for(int i = 0;i<d->height;i++)
		{
			for(int j = 0;j<d->width;j++)
			{
				CvScalar cs = cvGet2D(s,i,j);
#ifdef __DEBUG__
				if(i == 0 && j == 0)
					cout<<"Before Log: "<<cs.val[0]<<endl;
#endif
				for(int k = 0;k<s->nChannels;k++)
					if(cs.val[k]==0.0)
						cs.val[k] = -10000000.0;
					else
						cs.val[k]=log10(cs.val[k]);
				
				cvSet2D(d,i,j,cs);
#ifdef __DEBUG__
				if(i == 0 && j == 0)
				{
					CvScalar c = cvGet2D(d,i,j);
					cout<<"After Log: "<<c.val[0]<<endl;
				}
#endif
			}
		}
	}

	

	void ImageDecompose::imageExp10(const IplImage * s, IplImage *& d ,int depth)
	{
		if(!d) cvReleaseImage(&d);

		//Float point type after log operation.
		d = cvCreateImage(cvSize(s->width,s->height),depth,s->nChannels);

		for(int i = 0;i<d->height;i++)
		{
			for(int j = 0;j<d->width;j++)
			{
				CvScalar cs = cvGet2D(s,i,j);
#ifdef __DEBUG__
				if(i == 0 && j == 0)
					cout<<"Before Exp: "<<cs.val[0]<<endl;
#endif
				for(int k = 0;k<s->nChannels;k++)
						cs.val[k]=pow(10.0,cs.val[k]);

				cvSet2D(d,i,j,cs);
#ifdef __DEBUG__
				if(i == 0 && j == 0)
				{
					CvScalar c = cvGet2D(d,i,j);
					cout<<"After Exp: "<<c.val[0]<<endl;
				}
#endif
			}
		}
	}

	int ImageDecompose::intensityRange(const IplImage * img)
	{
		if(img->nChannels!=3)
			return -1;
		int maxValue = -10000;
		int minValue = 10000;
		for(int i = 0;i<img->height;i++)
			for(int j = 0;j<img->width;j++)
			{
				CvScalar cs = cvGet2D(img,i,j);
				int intensity =static_cast<int>( (cs.val[0] * 1.0 + cs.val[1] * 40.0 + cs.val[2] * 20.0) / 61.0 ) ;
				if(intensity > maxValue)
					maxValue = intensity;
				if(intensity < minValue)
					minValue = intensity;
			}
			return ( maxValue - minValue );
	}

	int ImageDecompose::rgbRange(const IplImage * img ,int rgb,int & maxV)
	{
		if(img->nChannels!=3)
			return -1;
		int maxValue = -10000;
		int minValue = 10000;
		for(int i = 0;i<img->height;i++)
			for(int j = 0;j<img->width;j++)
			{
				CvScalar cs = cvGet2D(img,i,j);
				if(cs.val[rgb] > maxValue)
					maxValue = (int)cs.val[rgb];
				if(cs.val[rgb] < minValue)
					minValue = (int)cs.val[rgb];
			}
		maxV = maxValue;
		return (maxValue - minValue);
	}

	float * ImageDecompose::spaceSigma(float s)
	{
		float * ss = new float[level];
		ss[0] = s;
		if(level >=2 )
		{
			ss[1] = 1.73205080756888f * ss[0];
			for(size_t l = 2;l < level;l++)
				ss[l] = pow(2.0f,(float)(l - 1)) * ss[l-1];
		}
		return ss;
	}
	
	float * ImageDecompose::rangeSigma(float r)
	{
		float * sr = new float[level];
		sr[0] = r;
		for(size_t l = 1;l < level;l++)
			sr[l] = sr[0] / pow(2.0f,(float)l);
		return sr;
	}

	void ImageDecompose::genDetail(const IplImage * s,const IplImage * b , IplImage *& d)
	{
		if(d)
			cvReleaseImage(&d);
		
		d = cvCreateImage(cvSize(s->width,s->height),s->depth,s->nChannels);
		
		for(int i = 0;i<s->height;i++)
			for(int j = 0;j<s->width;j++)
			{
				CvScalar cs = cvGet2D(s,i,j);
				CvScalar cb = cvGet2D(b,i,j);
				CvScalar cd;
				for(int k = 0;k<s->nChannels;k++)
					 cd.val[k] = cs.val[k] - cb.val[k];
				cvSet2D(d,i,j,cd);
			}
	}
	void ImageDecompose::genDetailDis(const IplImage * s , IplImage *& d)
	{
		d = cvCreateImage(cvSize(s->width,s->height),IPL_DEPTH_8U,s->nChannels);

		double * scaleFactor = new double[s->nChannels];

		for(int k = 0;k<s->nChannels;k++)
		{
			int maxV;
			int range = rgbRange(src,k,maxV);
			const double gamma = contrast /range;
			scaleFactor[k] = 1.0 / (maxV * gamma);
		}

		for(int i = 0;i<s->height;i++)
			for(int j = 0;j<s->width;j++)
			{
				CvScalar cs = cvGet2D(s,i,j);
				for(int k = 0;k<s->nChannels;k++)
					cs.val[k] = nMax(nMin(255.0 * pow(scaleFactor[k] * cs.val[k] , 1.0 / 2.2 ),255.0),0.0);
				cvSet2D(d,i,j,cs);
			}
		delete [] scaleFactor;
	}


	void ImageDecompose::Decompose()
	{
		if(!src)
			throw NMessage::NMessage("ImageDecompose::Decompose() - error : src is NULL");
		if(level == 0)
			throw NMessage::NMessage("ImageDecompose::Decompose() - error : level is 0");
		clear();

		float * spaceS = spaceSigma((float)coreWidth);
		float * rangeS = rangeSigma(intensityRange(src) / 40.0f);
		
		IplImage * log_image = NULL;
		imageLog10(src,log_image);
		
		IplImage * filtered_log_image = NULL;
		IplImage * detail_image = NULL;
		IplImage * exp_detail_image = NULL; 
		IplImage * display_detail_image = NULL;

		for(size_t i = 0;i<level;i++)
		{
			Bilateral::ColourImageBilateral32F(log_image,filtered_log_image,(int)coreWidth,spaceS[i],rangeS[i]);
			/*if(i == 0)
			{
				imageExp10(filtered_log_image,base,IPL_DEPTH_8U);
				CvvImage cvv;
				cvv.CopyOf(base);
				cvv.Save("e:\\base.bmp");
				cvv.Destroy();
			}*/

			imageExp10(filtered_log_image,base,IPL_DEPTH_8U);
			vecBase.push_back(base);

			if(i == level-1)
			{
				imageExp10(filtered_log_image,base,IPL_DEPTH_8U);
				log_base = cvCreateImage(cvSize(filtered_log_image->width,filtered_log_image->height),filtered_log_image->depth,filtered_log_image->nChannels);
				cvCopy(filtered_log_image,log_base);
			}
			genDetail(log_image,filtered_log_image,detail_image);
			log_detail.push_back(detail_image);
			imageExp10(detail_image,exp_detail_image,IPL_DEPTH_32F);
			genDetailDis(exp_detail_image,display_detail_image);
			detail.push_back(display_detail_image);

			cvReleaseImage(&log_image);
			//cvReleaseImage(&detail_image);
			cvReleaseImage(&exp_detail_image);
			detail_image = NULL;
			
			log_image = filtered_log_image;			
		}
		
		if(log_image)
			cvReleaseImage(&log_image);
		if(spaceS)
			delete [] spaceS;
		if(rangeS)
			delete [] rangeS;

		log_image = NULL;
		spaceS = NULL;
		rangeS = NULL;
	}

	void ImageDecompose::Save()
	{
		if(!base)
			throw NMessage::NMessage("ImageDecompose::Save() - error : base is NULL");

		CvvImage c;
		c.CopyOf(base);
		c.Save("E:\\id_base.bmp");
		c.Destroy();

		for(int i = 0;i<(int)detail.size();i++)
		{
			if(!detail[i])
				throw NMessage::NMessage("ImageDecompose::Save() - error : detail[i] is NULL");
			if(!vecBase[i])
				throw NMessage::NMessage("ImageDecompose::Save() - error : vecBase[i] is NULL");
			CvvImage * cvv = new CvvImage();
			CvvImage * cvv2 = new CvvImage();
			cvv->CopyOf(detail[i]);
			cvv2->CopyOf(vecBase[i]);
			char fileName[20];
			char fileName2[20];
			sprintf(fileName,"E:\\id_detail%d.bmp",i);
			sprintf(fileName2,"E:\\id_base%d.bmp",i);
			cvv->Save(fileName);
			cvv2->Save(fileName2);
			cvv->Destroy();
			cvv2->Destroy();
			delete cvv;
			delete cvv2;
		}
	}

	IplImage * ImageDecompose::Enhance(double * w,int n)
	{
		if(!log_base || n != (int)level)
			return NULL;
		int height = log_base->height;
		int width = log_base->width;
		IplImage * ipl = cvCreateImage(cvSize(width,height),IPL_DEPTH_8U,log_base->nChannels);
		for(int i = 0;i<height;i++)
			for(int j = 0;j<width;j++)
			{
				CvScalar cs = cvGet2D(log_base,i,j);
				for(int l = 0;l<(int)level;l++)				
				{
					CvScalar temp = cvGet2D(log_detail[l],i,j);
					for(int k = 0;k<base->nChannels;k++)
						cs.val[k] += temp.val[k] * w[l];
				}
				for(int k = 0;k<base->nChannels;k++)
					cs.val[k] = pow(10.0,cs.val[k]);
				cvSet2D(ipl,i,j,cs);
			}
		return ipl;
	}
}