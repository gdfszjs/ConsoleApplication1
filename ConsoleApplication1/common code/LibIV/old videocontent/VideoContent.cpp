//#include "stdafx.h"
#include "VideoContent.h"
//#include "LibIV.h"

VideoContent::VideoContent(void) :	pContent(NULL),
									m_nWidth(0),
									m_nHeight(0),
									m_nPages(0),
									m_nChannels(0),
									m_dFps(0.0)								
{

}

VideoContent::VideoContent(int width,
						   int height,
						   int page,
						   int channels,
						   double fps
						   )
{
	pContent = NULL;
	create(width,height,page,channels,fps);
}

VideoContent::VideoContent(const VideoContent & vc)
{
	assert(vc.pContent);
	m_nWidth = vc.m_nWidth;m_nHeight = vc.m_nHeight;
	m_nPages = vc.m_nPages;m_nChannels = vc.m_nChannels;
	m_dFps = vc.m_dFps;
	
	pContent = mem_alloc(m_nPages,m_nHeight,m_nWidth);
	
	mem_copy(pContent,vc.pContent,sizeof(CvScalar)*m_nPages*m_nHeight*m_nWidth);
	
}

void VideoContent::create(int width,int height,int page,int channels,double fps)
{
	m_nWidth = width;
	m_nHeight = height;
	m_nPages = page;
	m_nChannels = channels;
	m_dFps = fps;
	
	if(pContent)
		mem_free(&pContent);
	pContent = mem_alloc(m_nPages,m_nHeight,m_nWidth);
	memset(pContent[0][0],0,sizeof(CvScalar) * m_nPages * m_nHeight * m_nWidth);
}

VideoContent & VideoContent::operator=(const VideoContent & vc)
{
	cleanUp();
	m_nWidth = vc.m_nWidth;m_nHeight = vc.m_nHeight;
	m_nPages = vc.m_nPages;m_nChannels = vc.m_nChannels;
	m_dFps = vc.m_dFps;

	pContent = mem_alloc(m_nPages,m_nHeight,m_nWidth);

	mem_copy(pContent,vc.pContent,sizeof(CvScalar)*m_nPages*m_nHeight*m_nWidth);

	return (*this);
}

CvScalar *** VideoContent::mem_alloc(int p,int h,int w)
{
	CvScalar *** a = new CvScalar**[p];
	a[0] = new CvScalar*[p*h];
	a[0][0] = new CvScalar[p*h*w];

	for(int i = 1;i<p;i++)
	{
		a[i] = a[i-1] + h;
		a[i][0] = a[i-1][0] + w * h;
	}
	
	for(int i = 0;i<p;i++)
	{
		for(int j = 1;j<h;j++)
		{
			a[i][j] = a[i][j-1] + w;
		}
	}
	return a;
}

void VideoContent::mem_free(CvScalar **** p)
{
	if(!p || !*p)
		return;
	delete [] (*p)[0][0];
	delete [] (*p)[0];
	delete [] (*p);
	(*p) = NULL;
}

void VideoContent::mem_copy(CvScalar *** dst, CvScalar *** src, unsigned int max_count)
{
	assert(dst&&src);
	memcpy(dst[0][0],src[0][0],max_count);
}


VideoContent::~VideoContent(void)
{
	cleanUp();
}

void VideoContent::cleanUp()
{
	mem_free(&pContent);
	m_nWidth = 0;
	m_nHeight = 0;
	m_nChannels = 0;
	m_nPages = 0;
	m_dFps = 0.0;
	
}



bool VideoContent::readVideo(const char * fileName)
{
	cleanUp();
	
	CvCapture * pCapture = cvCreateFileCapture(fileName);
	if(!pCapture)
		return false;
	
	m_nWidth	=	(int)cvGetCaptureProperty(pCapture,CV_CAP_PROP_FRAME_WIDTH);
	m_nHeight	=	(int)cvGetCaptureProperty(pCapture,CV_CAP_PROP_FRAME_HEIGHT);
	m_nPages	=	(int)cvGetCaptureProperty(pCapture,CV_CAP_PROP_FRAME_COUNT);
	m_dFps		=   cvGetCaptureProperty(pCapture,CV_CAP_PROP_FPS);
	
	pContent = mem_alloc(m_nPages,m_nHeight,m_nWidth);

	IplImage * ipl_tmp = NULL;
	
	for(int i = 0;i<m_nPages;i++)
	{
		cvSetCaptureProperty(pCapture,CV_CAP_PROP_POS_FRAMES,i);
		ipl_tmp = cvQueryFrame(pCapture);
		
		for(int j = 0;j<m_nHeight;j++)
			for(int k = 0;k<m_nWidth;k++)
				pContent[i][j][k] = cvGet2D(ipl_tmp,m_nHeight-j-1,k);
	}

	m_nChannels = ipl_tmp->nChannels;

	cvReleaseCapture(&pCapture);
	return true;
}

bool VideoContent::saveVideo(const char * fileName)
{
	if(!pContent)
		return false;
	CvVideoWriter * pWriter = cvCreateVideoWriter(fileName,0,m_dFps,cvSize(m_nWidth,m_nHeight),1);
	if(!pWriter)
		return false;
	IplImage * ipl_tmp = cvCreateImage(cvSize(m_nWidth,m_nHeight),IPL_DEPTH_8U,m_nChannels);
	//ipl_tmp->origin = 1;

	for(int i = 0;i<m_nPages;i++)
	{
		for(int j = 0;j<m_nHeight;j++)
			for(int k = 0;k<m_nWidth;k++)
				cvSet2D(ipl_tmp,j,k,pContent[i][j][k]);
		
		cvWriteFrame(pWriter,ipl_tmp);
	}
	cvReleaseImage(&ipl_tmp);
	cvReleaseVideoWriter(&pWriter);
	return true;
}

bool VideoContent::saveVideo(const char * fileName,int n)
{
	char _fileName[1024];
	char _cNum[256];
	sprintf(_cNum,"%d",n);
	
	int len = 0;
	int len2;
	while(fileName[len] != '\0')
	{
		len++;
	}
	if(len <= 4)
		return false;
	len2 = len;
	len--;
	while(fileName[len] != '.')
	{
		len--;
	}
	int lenn = 0;
	while(_cNum[lenn] != '\0')
	{
		lenn++;
	}
	
	for(int i = 0;i<len;i++)
		_fileName[i] = fileName[i];
	for(int i = 0;i<lenn;i++)
		_fileName[i+len] = _cNum[i];
	for(int i = len;i<len2;i++)
		_fileName[lenn + i] = fileName[i];

	_fileName[len2 + lenn] = '\0';

	return saveVideo(_fileName);
}

bool VideoContent::saveVideo(const char * fileName,int b,int e)
{
	if(!pContent)
		return false;
	CvVideoWriter * pWriter = cvCreateVideoWriter(fileName,CV_FOURCC('X','V','I','D'),m_dFps,cvSize(m_nWidth,m_nHeight),1);
	if(!pWriter)
		return false;
	IplImage * ipl_tmp = cvCreateImage(cvSize(m_nWidth,m_nHeight),IPL_DEPTH_8U,m_nChannels);
	ipl_tmp->origin = 1;

	for(int i = b;i<=e;i++)
	{
		for(int j = 0;j<m_nHeight;j++)
			for(int k = 0;k<m_nWidth;k++)
				cvSet2D(ipl_tmp,j,k,pContent[i][j][k]);
		cvWriteFrame(pWriter,ipl_tmp);
	}
	cvReleaseImage(&ipl_tmp);
	cvReleaseVideoWriter(&pWriter);
	return true;
}

IplImage * VideoContent::getFrame(int p)
{
	if(!pContent)
		return NULL;
	if(p<0 || p >= m_nPages)
		return NULL;

	IplImage * presult = cvCreateImage(cvSize(m_nWidth,m_nHeight),IPL_DEPTH_8U,m_nChannels);
	
	for(int i = 0;i<m_nHeight;i++)
	{
		for(int j = 0;j<m_nWidth;j++)
			cvSet2D(presult,i,j,pContent[p][i][j]);
	}
	return presult;
}