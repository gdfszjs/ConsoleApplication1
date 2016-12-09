#include "StdAfx.h"
#include "WorkerThread.h"

CWorkerThread::CWorkerThread(WorkerThreadProc _workerproc,int pri,UINT siz,DWORD fla,LPSECURITY_ATTRIBUTES attrs) :
pthread(NULL),workerproc(_workerproc),nPriority(pri),nStackSize(siz),dwCreateFlags(fla),lpSecurityAttrs(attrs){}

CWorkerThread::~CWorkerThread(void)
{
	if(pthread)
	{
		delete pthread;
		pthread = NULL;
	}
}

void CWorkerThread::BeginWorkerThread(LPVOID pParam)
{
	pthread = AfxBeginThread(workerproc,pParam,nPriority,nStackSize,dwCreateFlags,lpSecurityAttrs);
	pthread->m_bAutoDelete = FALSE;
}

long CWorkerThread::ResumeWorkerThread()
{
	if(pthread)
		return (long)pthread->ResumeThread();
	else return -1;
}

long CWorkerThread::SuspendWorkerThread()
{
	if(pthread)
		return (long)pthread->SuspendThread(); 
	else return -1;
}


long CWorkerThread::GetExitCode(void)
{
	if(!pthread)
		return -1;
	else
	{
		DWORD dwExitCode;
		::GetExitCodeThread(pthread->m_hThread,&dwExitCode);
		return (long)dwExitCode;
	}
}