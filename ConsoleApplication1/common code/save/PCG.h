
#pragma once
//初雨 2008.12 Chu Yu.
class PCG //前置条件共轭梯度 对Ax=b求解 其中A为输入条件为对称、正定的矩阵 
{
private:
	double *sa,*sb,*x,*b;
	unsigned long *ija;								//What's this?
	double tol;//迭代精度
	int itmax;//迭代次数
	int iter;//实际迭代次数
	unsigned long n;//矩阵边长 
	int flag;									//What's this? Precondition factor
	unsigned long nz;//稀疏矩阵中非零元的数目		
public:
	~PCG(void);
	//构造函数参数列表如下：
	//sa,ija 输入的矩阵A，以稀疏形式存储        //How? sa is used to store values and ija is a auxiliary array.
	//x,最终求解的结果
	//b方程右边的矩阵b
	//n A的行数
	//sb ijb 预处理因子
	// tol 迭代精度
	// itmax 最大迭代次数
	//iter 实际迭代次数
	PCG(double * &sa_,unsigned long *&ija_,double *&x_,double *&b_,const unsigned long n_,
											const double tol_,const int itmax_);
	PCG(double * &sa_,unsigned long *&ija_,double *&x_,double *&b_,const unsigned long n_,
								const double tol_,const int itmax_,const int flag_,const int nz_);//当flag为0时，sb，ijb 一定要赋值！
	
	//预处理因子
	//flag=0,为默认的雅可比处理因子，即把主对角元作为预处理因子
	//flag=2，为超松弛预处理因子,sa1 为分解的结果 
	void matrix_SSOR( double * &sa, unsigned long * &ija,const double w,const unsigned long n,double * &sa1);	
	//flag=1 为乔列斯基不完全分解，分解的结果保存在sa1中
	bool Choldc(double *&sa1, unsigned long *& ija1,const unsigned long n);								

	//进行PCG 迭代
	int linbcg();

	//linbcg 中用到的函数 计算b=Ax
	//sa ija 为A ，n为其行数
	void asolve(double *sa,unsigned long *ija,double *x,double *b, unsigned long n);

	///linbcg 用到的函数 ，求解Mx=b 输出为x 其中M为预处理因子，为不完全Cholsky分解的结果，为下三角部分
	void dsprsax(double *sa,unsigned long *ija,double *x,double *b ,unsigned long n);
};
