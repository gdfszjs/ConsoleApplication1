/*****************************
	
	工作者线程
	2008/3/30 by Yongwei Nie
	
	更改于2008/3/31

	1. 将CWinThread::m_bAutoDelete设为false
	2. ResumeWorkerThread 和 SuspendWorkerThread 加上返回值 (long)
	3. GetExitCode()获取线程退出的原因，若线程并未结束，返回STILL_ACTIVE
	
*****************************/
#pragma once

class CWorkerThread
{
public:
	typedef UINT (*WorkerThreadProc)(LPVOID pParam);
private:
	CWinThread * pthread;
	WorkerThreadProc workerproc;
	int nPriority;
	UINT nStackSize;
	DWORD dwCreateFlags;
	LPSECURITY_ATTRIBUTES lpSecurityAttrs;
	
public:
	explicit CWorkerThread(WorkerThreadProc _workerproc,int pri = THREAD_PRIORITY_NORMAL,UINT siz = 0,
		DWORD fla = CREATE_SUSPENDED,LPSECURITY_ATTRIBUTES attrs = NULL);
	void BeginWorkerThread(LPVOID pParam);
    long ResumeWorkerThread();
	long SuspendWorkerThread();
	long GetExitCode();
	~CWorkerThread(void);
};
