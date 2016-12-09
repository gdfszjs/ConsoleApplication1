#pragma once

#include <fstream>
using std::ofstream;
using std::ifstream;
using std::ios;

class UnifiedMask
{
public:
	
	enum MASK_CONTAINER_TYPE { MCT_MAT_2D_INT = 0, MCT_INT = 1 };

	UnifiedMask();
	~UnifiedMask();

	
	void * Load(MASK_CONTAINER_TYPE type);							//Load mask from files.

	//int * mask
	void Save(const int * pmask,int w,int h);			//Save mask to files.

	//MAT_2D_INT
	void Save(const MAT_2D_INT * pmask,int w,int h);

protected:

	MAT_2D_INT * loadMat2D(char * filename);
	int		   * loadInt(char * filename);
};


UnifiedMask::UnifiedMask()
{
	
}

UnifiedMask::~UnifiedMask()
{
	
}

void UnifiedMask::Save(const int * pmask,int w,int h)
{
	CFileDialog fileSave(FALSE,_T("*.txt"),_T("text"),OFN_HIDEREADONLY | OFN_OVERWRITEPROMPT | OFN_FILEMUSTEXIST,
		_T("Text Files (*.txt)|*.txt"));

	if(fileSave.DoModal() == IDOK)
	{
		CString filePath = fileSave.GetPathName();
		try
		{										
			if(filePath.Right(3).CompareNoCase(_T("txt")))
				throw new CFileException;

			wchar_t * wText = filePath.GetBuffer(0);
			DWORD dwNum = WideCharToMultiByte(CP_OEMCP,NULL,wText,-1,NULL,0,NULL,FALSE);
			char * psText = new char[dwNum];
			WideCharToMultiByte(CP_OEMCP,NULL,wText,-1,psText,dwNum,NULL,FALSE);

			ofstream mout(psText,ios::out);
			if(!mout)
				return;

			mout<<(int)MCT_INT<<' '<<w<<' '<<h<<' ';

			for(int i = 0;i<w*h;i++)
				mout<<pmask[i]<<' ';


			filePath.ReleaseBuffer();
			delete [] psText;


		}
		catch(CFileException * fe)
		{
			AfxMessageBox(_T("无法写入非txt文件"));
			delete fe;
		}
	}
}

void UnifiedMask::Save(const MAT_2D_INT * pmask,int w,int h)
{
	CFileDialog fileSave(FALSE,_T("*.txt"),_T("text"),OFN_HIDEREADONLY | OFN_OVERWRITEPROMPT | OFN_FILEMUSTEXIST,
		_T("Text Files (*.txt)|*.txt"));

	if(fileSave.DoModal() == IDOK)
	{
		CString filePath = fileSave.GetPathName();
		try
		{										
			if(filePath.Right(3).CompareNoCase(_T("txt")))
				throw new CFileException;

			wchar_t * wText = filePath.GetBuffer(0);
			DWORD dwNum = WideCharToMultiByte(CP_OEMCP,NULL,wText,-1,NULL,0,NULL,FALSE);
			char * psText = new char[dwNum];
			WideCharToMultiByte(CP_OEMCP,NULL,wText,-1,psText,dwNum,NULL,FALSE);

			ofstream mout(psText,ios::out);
			if(!mout)
				return;
			
			mout<<(int)MCT_MAT_2D_INT<<' '<<w<<' '<<h<<' ';
			
			for(int i = 0;i<h;i++)
				for(int j = 0;j<w;j++)
					mout<<(*pmask)[i][j]<<' ';


			filePath.ReleaseBuffer();
			delete [] psText;


		}
		catch(CFileException * fe)
		{
			AfxMessageBox(_T("无法写入非txt文件"));
			delete fe;
		}
	}
}

void * UnifiedMask::Load(MASK_CONTAINER_TYPE type)
{
	CFileDialog fileLoad(TRUE,_T("*.txt"),_T("text"),OFN_HIDEREADONLY | OFN_OVERWRITEPROMPT | OFN_FILEMUSTEXIST,
		_T("Text Files (*.txt)|*.txt"));

	if(fileLoad.DoModal() == IDOK)
	{
		CString filePath = fileLoad.GetPathName();
		try
		{										
			if(filePath.Right(3).CompareNoCase(_T("txt")))
				throw new CFileException;
			
			wchar_t * wText = filePath.GetBuffer(0);
			DWORD dwNum = WideCharToMultiByte(CP_OEMCP,NULL,wText,-1,NULL,0,NULL,FALSE);
			char * psText = new char[dwNum];
			WideCharToMultiByte(CP_OEMCP,NULL,wText,-1,psText,dwNum,NULL,FALSE);
	
			ifstream in(psText);
			if(!in)
				return NULL;
			
			if(type == MCT_MAT_2D_INT)
			{
				filePath.ReleaseBuffer();
				delete [] psText;
				return loadMat2D(psText);
			}
			else if(type == MCT_INT)
			{
				filePath.ReleaseBuffer();
				delete [] psText;
				return loadInt(psText);
			}
			else
			{
				filePath.ReleaseBuffer();
				delete [] psText;
				return NULL;
			}
			
		}
		catch(CFileException * fe)
		{
			AfxMessageBox(_T("无法打开非txt文件"));
			return NULL;
			delete fe;
		}
	}
	else
		return NULL;
}

MAT_2D_INT * UnifiedMask::loadMat2D(char * filename)
{
	int type;
	int width;
	int height;
	int a;
	int index = 0;
	ifstream in(filename);
	if(!in)
		return NULL;

	in>>type>>width>>height;
	
	MAT_2D_INT * pmat = new MAT_2D_INT(0,height,width);
	

	while(in>>a)
	{
		(*pmat)[index / width][index % width ] = a;
		index++;
	}

	return pmat;
}

int * UnifiedMask::loadInt(char * filename)
{
	int type;
	int width;
	int height;
	int a;
	int index = 0;
	
	ifstream in(filename);
	if(!in)
		return NULL;
	in>>type>>width>>height;

	int * pint = new int[height * width];
	memset(pint,0,width * height * sizeof(int));
	
	while(in>>a)
		pint[index++] = a;
	return pint;
}
