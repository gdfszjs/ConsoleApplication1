/*
	===========================
	VideoContent	09/10/30
	===========================
	A new version of VideoContent.
	
*/

#pragma once
#include "common.h"

class VideoContent
{
public:
	// ------------------------------------------------------------------------------
	//Initialize a VideoContent. No content.
	VideoContent(void);														
	VideoContent(int width,int height,int page,int channels,double fps);		

	// ------------------------------------------------------------------------------
	//Copy constructor.
	VideoContent(const VideoContent & vc);
	
	// ------------------------------------------------------------------------------
	//Deconstructor
	~VideoContent(void);
public:
	// ------------------------------------------------------------------------------
	//Create the video
	void create(int width,int height,int page,int channels,double fps);
	//Copy operator
	VideoContent & operator=(const VideoContent & vc);

	//Read and save the video
	bool readVideo(const char * fileName);
	bool saveVideo(const char * fileName);
	bool saveVideo(const char * fileName,int n);
	// Save the frames that from b to e
	bool saveVideo(const char * fileName,int b,int e);

	// Get the pth frame
	IplImage* getFrame(int p);
	
	
	// ------------------------------------------------------------------------------
	//Get the video's information
	void getParameters(int& width, int &height,int &page,int& channel ,double& fps) const
	{ width = m_nWidth;height=m_nHeight;page = m_nPages;channel=m_nChannels; fps = m_dFps; }

	int getPage(void)   const { return m_nPages;    }
	int getWidth(void)  const { return m_nWidth;    }
	int getHeight(void) const { return m_nHeight;   }
	double getFps(void) const { return m_dFps;      }
	int getChannels()   const { return m_nChannels; }
	
	// ------------------------------------------------------------------------------
	//The video, a 3-dimension array.
	CvScalar *** pContent;

private:
	void cleanUp(void);

	// ------------------------------------------------------------------------------
	//Memory operations
	CvScalar *** mem_alloc(int p,int h,int w);
	void mem_free(CvScalar **** p);
	void mem_copy(CvScalar *** dst, CvScalar *** src, unsigned int max_count);

private:
	int m_nWidth, m_nHeight, m_nPages, m_nChannels;
	double m_dFps;	
};
