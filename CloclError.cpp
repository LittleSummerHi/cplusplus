// CloclError.cpp : 定义 DLL 的初始化例程。
//

#include "stdafx.h"
#include "CloclError.h"
#include "stdafx.h"
#include "RTK.h"
#include "struct.h"
#include "Lambda.h"
#include "string.h"
#define MAXDIM     20

struct AttiDetGlobal* pglobal;
TORGOBSDATA *m_ptMMData;
TORGOBSDATA *m_ptS1MData;
TORGOBSDATA *m_ptMMDataNew;
TORGOBSDATA *m_ptS1MDataNew;

/**********************************************************************************************


#ifdef _DEBUG
#define new DEBUG_NEW
#endif

//
//TODO:  如果此 DLL 相对于 MFC DLL 是动态链接的，
//		则从此 DLL 导出的任何调入
//		MFC 的函数必须将 AFX_MANAGE_STATE 宏添加到
//		该函数的最前面。
//
//		例如:
//
//		extern "C" BOOL PASCAL EXPORT ExportedFunction()
//		{
//			AFX_MANAGE_STATE(AfxGetStaticModuleState());
//			// 此处为普通函数体
//		}
//
//		此宏先于任何 MFC 调用
//		出现在每个函数中十分重要。  这意味着
//		它必须作为函数中的第一个语句
//		出现，甚至先于所有对象变量声明，
//		这是因为它们的构造函数可能生成 MFC
//		DLL 调用。
//
//		有关其他详细信息，
//		请参阅 MFC 技术说明 33 和 58。
//

// CCloclErrorApp

BEGIN_MESSAGE_MAP(CCloclErrorApp, CWinApp)
END_MESSAGE_MAP()


// CCloclErrorApp 构造

CCloclErrorApp::CCloclErrorApp()
{
// TODO:  在此处添加构造代码，
// 将所有重要的初始化放置在 InitInstance 中
}


// 唯一的一个 CCloclErrorApp 对象

CCloclErrorApp theApp;


// CCloclErrorApp 初始化

BOOL CCloclErrorApp::InitInstance()
{
CWinApp::InitInstance();

return TRUE;
}









*********************************************************************************************************/

double SDP[MAXNUMSV];//载波相位单差
double R[MAXNUMSV*MAXNUMSV];//twBlock=Rzt，为上三角矩阵
double twBlock[MAXNUMSV*MAXNUMSV];//twBlock=U*phatP
double ZV[180];
double threshold;
/*  wavelen*SDP=c*RecedTime-wavelen*NdTime+v
*****正交变换消除接收机钟差RecedTtime
*****SDP为载波相位小数部分差值
*****NdTime为接收机钟差造成的整周模糊度,zt为其浮点解
*/

void GetMatrixSDP(double *SDP, TORGOBSDATA *m_ptMMDataNew, TORGOBSDATA *m_ptS1MDataNew, double *ZV)
{
	int i = 0;
	int zv = 0;
	for (i = 0; i<m_ptMMDataNew->u1SVBlocks; i++)
	{
		zv = m_ptMMDataNew->tSv[i].u1PID;
		SDP[i] = m_ptS1MDataNew->tSv[i].dPhase - m_ptMMDataNew->tSv[i].dPhase - ZV[zv];//载波相位单差小数部分
	}    //	ZV[zv] = (int)(m_ptS1MDataNew->tSv[i].dPhase - m_ptMMDataNew->tSv[i].dPhase)
}

void CalculatePEk2(double *P, int m)//household变换中消除钟差项时用到的P~矩阵
{
	double a = -1. / (m - sqrt((double)m));
int i = 0, j = 0;
for (i = 0; i<m - 1; i++)
{
	for (j = 0; j<m; j++)
	{
		//~P阵
		P[i*m + j] = a;
		if (j == 0)
			P[i*m + j] = 1 / sqrt((double)m);
		else if (i == (j - 1))
			P[i*m + j] = 1 + a;
	}
}
}
void GetFMatrix(double *F, int m)//P~=FD
{
	int i = 0, j = 0;
	double a = -1. / (m - sqrt((double)m));
	for (i = 0; i<m - 1; i++)
	{
		for (j = 0; j<m - 1; j++)
		{
			F[i*(m - 1) + j] = a;
			if (i == j)
			{
				F[i*(m - 1) + j] = 1 + a;
			}
		}
	}
}

void CalculatePhatP(double *PhatP, double *P, double *SDP, int m)
{
	int i = 0, j = 0;
	for (i = 0; i<m - 1; i++)
	{
		for (j = 0; j<m; j++)
		{
			PhatP[i] += P[i*m + j] * SDP[j];//PhatP=P~y

		}
	}
}

void CalculateQat(double *Qat, double *R, int m)//浮点解的协方差矩阵
{
	double delta = 1 / (2 * pglobal->InitPara.fPhasePrecision1*pglobal->InitPara.fPhasePrecision1);

	int i = 0, j = 0, k = 0; unsigned char bInvert;

	double SZ_[MAXNUMSV*MAXNUMSV];//Z变换矩阵

	memset(SZ_, 0, sizeof(double)*MAXNUMSV*MAXNUMSV);
	memcpy(SZ_, R, sizeof(double)*m*m);

	for (i = 0; i<m; i++)
	{
		for (j = i + 1; j<m; j++)
		{
			//SZ_ =SZ^ 转置
			double temp = SZ_[i*m + j];
			SZ_[i*m + j] = SZ_[j*m + i];
			SZ_[j*m + i] = temp;
		}
	}
	for (i = 0; i<m; i++)
	{
		for (j = 0; j<m; j++)
		{
			for (k = 0; k<m; k++)
			{
				//Q = SZ_*SZ;
				Qat[i*m + j] += delta*SZ_[i*m + k] * R[k*m + j];//浮点解的协方差矩阵
			}
		}
	}
	bInvert = InvertMatrix(Qat, 1e-5, m);
}

void householderQR(double *X, double *Q, int m, int n)
{
	int i = 0, j = 0, k = 0, l = 0;

	double beta = 0.0;
	double v[MAXDIM];
	double x[MAXDIM];
	double B[MAXDIM*MAXDIM];
	double C[MAXDIM*MAXDIM];
	double Hk[MAXDIM*MAXDIM];

	//计算Q阵。
	memset(Q, 0, sizeof(double)*m*m);
	for (i = 0; i<m; i++)
		Q[i*m + i] = 1.0; //Q = I^m.

	//////加一个约束 使行小于列的情况也可以用
	int t;
	if (m>n || m == n)
		t = n;
	else
		t = m;

	for (j = 0; j<t; j++)   //for j=1:n.
	{
		//异常处理。
		if (m - j == 1) continue;

		memset(v, 0, sizeof(double)*MAXDIM);
		memset(x, 0, sizeof(double)*MAXDIM);
		for (i = j; i<m; i++)
		{
			x[i - j] = X[i*n + j];  //A[j:m,j].
		}
		beta = house(x, v, m - j);//用到House变换函数在下面

		//update A[j:m,j:n] = (I^(m-j/*+1*/) - beta*v*vt)A(j:m,j:n).
		//(1) :计算 (I^(m-j/*+1*/) - beta*v*vt).
		memset(B, 0, sizeof(double)*MAXDIM*MAXDIM);
		for (i = 0; i<m - j; i++)
		{
			for (k = 0; k<m - j; k++)
			{
				double a = 0;
				if (i == k)
					a = 1.0000000;
				B[i*(m - j) + k] = a - beta * v[i] * v[k];
			}
		}
		////////
		memset(Hk, 0, sizeof(double)*MAXDIM*MAXDIM);
		for (i = 0; i<j; i++)
		{
			Hk[i*m + i] = 1.0;
		}
		for (i = 0; i<m - j; i++)
		{
			for (k = 0; k<m - j; k++)
			{
				Hk[(i + j)*m + k + j] = B[i*(m - j) + k];
			}
		}
		//(2) 计算C = B *A .
		memset(C, 0, sizeof(double)*MAXDIM*MAXDIM);
		for (i = 0; i<m; i++)
		{
			for (l = 0; l<n; l++)
			{
				for (k = 0; k<m; k++)
				{
					C[i*n + l] += Hk[i*m + k] * X[k*n + l];
				}
			}
		}
		memcpy(X, C, sizeof(double)*m*n);

		//（5）计算Q阵。


		memset(C, 0, sizeof(double)*MAXDIM*MAXDIM);
		for (i = 0; i<m; i++)
		{
			for (k = 0; k<m; k++)
			{
				for (l = 0; l<m; l++)
				{
					C[i*m + k] += Hk[i*m + l] * Q[l*m + k];
				}
			}
		}

		memcpy(Q, C, sizeof(double)*m*m);
	}
}

//----------------------------------------------------------------------//
//	功能描述:	(Householder向量）给定x	& R^n ,本函数计算满足 v(1) =1   //
//				的v &R^n 和 beta & R 使得 P = I^n - beta*v*vt 是正交阵 //
//				且Px = norm_2(x)*e1.									//
//				参见《矩阵计算》第三版 算法 5.1.1.						//
//	返 回 值:	beta													//
//	参    数:	x ，v ： n维向量。										//
//	修订历史:														  	//
//----------------------------------------------------------------------//
double house(const double *x, double *v, int n)
{
	//[1] 对向量X进行加权运算，避免上溢。
	double X[MAXDIM]; double norm = 0; int i = 0;
	double beta = 0; double a = 0;

	memcpy(X, x, sizeof(double)*n);

	for (i = 0; i<n; i++)
	{
		norm += X[i] * X[i];
	}
	norm = sqrt(norm);
	for (i = 0; i<n; i++)
	{
		X[i] = X[i] / norm;
	}

	v[0] = 1.0;  //n ==1 	
	for (i = 1; i<n; i++)
	{
		a += X[i] * X[i];
		v[i] = X[i];
	}

	if (a > 1e-5)
	{
		double c = sqrt(X[0] * X[0] + a);
		if (X[0] <0)
			v[0] = X[0] - c;
		else
			//v[0] = X[0] +c;  
			v[0] = -a / (X[0] + c); //parlett(1971年）

		beta = 2 * v[0] * v[0] / (a + v[0] * v[0]);
		c = v[0];
		for (i = 0; i<n; i++)
		{
			v[i] = v[i] / c;
		}
	}

	return beta;
}

void GetTUVRMatrix(const double *A, const double *T, double *U,  double *R, int m)
{
	int i = 0, j = 0, k = 0;
	//U阵 [m-1,m-1]
	for (i = 0; i<(m - 1)*(m - 1); i++)
	{
		U[i] = T[i];
	}
	//R 阵[m-1,m-1]。
	for (i = 0; i<m - 1; i++)
	{
		for (j = 0; j<m - 1; j++)
		{
			R[i*(m - 1) + j] = A[i*(m - 1) + j];
			if (i>j)
				R[i*(m - 1) + j] = 0;
		}
	}
}
void CalculatetwBlock(double *twBlock,double *U, double *PhatP, int m)
{
	int i = 0, j = 0;
	for (i = 0; i<m - 1; i++)
	{
		for (j = 0; j<m - 1; j++)
		{
			twBlock[i] += U[i*(m - 1) + j] *PhatP[j];
		}
	}
}

//----------------------------------------------------------------------//
//	功能描述:	递归求解x变量。											//
//	返 回 值:	None.													//
//	参    数:	x^n --一维变量	Rx = b									//
//				其中R 为上三角矩阵。									//
//	修订历史:														  	//
//----------------------------------------------------------------------//
void RecursiveX(double *x, const double *b, const double *R, int n)
{
	int i = 0, j = 0;
	double a = 0.0;
	for (i = n - 1; i >= 0; i--)
	{
		a = b[i];
		for (j = i + 1; j <= n - 1; j++)
		{
			a -= x[j] * R[i*n + j];
		}

		x[i] = a / R[i*n + i];
	}
}

////////用QR分解由z反求N
void CalculateDtime(double *SDP, double *z, double *Dtime, double *P, double *F, int m)

{
	int i = 0, j = 0, k = 0;
	double Fz[MAXNUMSV];
	double RpInvert[MAXNUMSV];
	double N[MAXNUMSV];
	double Rp[MAXNUMSV];
	double T[MAXNUMSV];
	double TFz[MAXNUMSV];
	double N1[MAXNUMSV];
	double Up[MAXNUMSV];
	for (i = 0; i<m - 1; i++)
	{
		Fz[i] = 0;
		for (j = 0; j<m - 1; j++)
		{
			Fz[i] += F[i*(m - 1) + j] * z[j];
		}
	}

	householderQR(P, T, m - 1, m);//P进行QR分解****************************考虑维数
	GetTUVRMatrix(P, T, Up, Rp, m);//P分解得到U,R****************************考虑维数

	for (i = 0; i<m; i++)
	{
		Rp[i] = 0;
		for (j = 0; j<m; j++)
		{
			Rp[i] += T[i*m + j] * P[j];
		}
	}

	memcpy(RpInvert, Rp, sizeof(double)* 9);
	InvertMatrix(RpInvert, 1e-4, m);
	for (i = 0; i<m - 1; i++)
	{
		TFz[i] = 0;
		for (j = 0; j<m- 1; j++)
		{
			TFz[i] += T[i * (m - 1) + j] * Fz[j];
		}
	}

	for (i = 0; i<m; i++)
	{
		N[i] = 0;
		for (j = 0; j<m; j++)
		{
			N[i] += RpInvert[i*m + j] * TFz[j];
		}
	}
	N1[i] = fabs(N[i]);
	Dtime[i] = SDP[1] + N1[1];

}








struct ParaDesc defaults[] =
{
	"RTK", "20140715", "RTK定位算法模块"

};
extern "C"__declspec(dllexport)
int GetDefault(struct ParaDesc **pDesc)
{
	*pDesc = defaults;
	return sizeof(defaults) / sizeof(struct ParaDesc);
}

extern "C"__declspec(dllexport)
bool Init(struct AttiDetGlobal* pGlobal)
{
	pglobal = pGlobal;
	printf("子模块 RTK 正在启动。\n");
	return true;
}
extern "C"__declspec(dllexport)
bool UnInit()
{
	return true;
}

extern "C"__declspec(dllexport)
bool Process(void)
{

	int m1;//当前历元的卫星数
	double DD_Phase[MAXNUMSV];//双差相位
	double DD_Pseudo[MAXNUMSV];//双差伪距
	double Nt[MAXNUMSV];//钟差整周绝对值
	double Ntfloat[MAXNUMSV];//当前钟差整周减去初始整周
	double n = sqrt(2.0);
	double threshold_1[MAXNUMSV];
	memset(threshold_1, 0, sizeof(double)*MAXNUMSV);
	double threshold_0;
	double threshold_2[MAXNUMSV];
	memset(threshold_2, 0, sizeof(double)*MAXNUMSV);
	threshold_0 = 10 * (n / pglobal->InitPara.SAT[0].dWAVELEN - 0.01);
	m1 = m_ptMMData->u1SVBlocks;

	/*for (i = 0; i < m1; i++)
	{
	m_ptMMDataNew->tSv[i].dPhase_float = fmod(m_ptMMDataNew->tSv[i].dPhase, 1.0);//fmod函数：计算第一个指定参数除以第二个指定参数的余数。包含在math.h
	m_ptS1MDataNew->tSv[i].dPhase_float = fmod(m_ptS1MDataNew->tSv[i].dPhase, 1.0);//得到载波相位小数部分
	}
	*/

	int o = 0;
	for (int i = 1; i < m1; i++)

	{
		/*DD_Phase[o] = (m_ptMMDataNew->tSv[i].dPhase - m_ptS1MDataNew->tSv[i].dPhase) - (m_ptMMDataNew->tSv[0].dPhase - m_ptS1MDataNew->tSv[0].dPhase);*/
		DD_Phase[o] = (m_ptMMDataNew->tSv[i].dPhase_float - m_ptS1MDataNew->tSv[i].dPhase_float) - (m_ptMMDataNew->tSv[0].dPhase_float - m_ptS1MDataNew->tSv[0].dPhase_float);//双差载波相位
		DD_Pseudo[o] = (m_ptMMDataNew->tSv[i].dPseudo - m_ptS1MDataNew->tSv[i].dPseudo) - (m_ptMMDataNew->tSv[0].dPseudo - m_ptS1MDataNew->tSv[0].dPseudo);//双差伪距（当前伪距差减去初始伪距差）
		Ntfloat[o] = DD_Pseudo[o] / pglobal->InitPara.SAT[0].dWAVELEN - DD_Phase[o];//双差整周
		/*Ntfloat[o] = DD_Phase[o]-DD_Pseudo[o] / pglobal->InitPara.SAT[0].dWAVELEN;*/
		Nt[o] = fabs(Ntfloat[o]);//fabs:求浮点数的绝对值
		threshold_1[o] = (int)(Nt[o] + threshold_0);
		threshold_2[o] = (int)(Nt[o] - threshold_0);//整周的上下界
		o = o + 1;
	}

	///////////////////////[7] 记录新升起的卫星的载波整数部分  减去初始整数
	double SDP[MAXNUMSV];
	int tempVs[MAXNUMSV];
	int VsCurrent[MAXNUMSV];
	memset(tempVs, 0, sizeof(int)*MAXNUMSV);
	memset(VsCurrent, 0, sizeof(int)*MAXNUMSV);
	int i;
	for (i = 0; i < m_ptMMData->u1SVBlocks; i++)
	{
		tempVs[i] = m_ptMMData->tSv[i].u1PID;//基站伪随机噪声码编号
		VsCurrent[i] = m_ptMMData->tSv[i].u1PID;
	}
	memset(SDP, 0, sizeof(double)*MAXNUMSV);
	{
		for (i = 0; i < m1; i++)
		{
			int zv = 0;
			zv = VsCurrent[i];
			ZV[zv] = (int)(m_ptS1MDataNew->tSv[i].dPhase - m_ptMMDataNew->tSv[i].dPhase);//新的载波相位差值取整
		}
	}
	GetMatrixSDP(SDP, m_ptMMDataNew, m_ptS1MDataNew, ZV);//输出SDP为载波相位差小数部分


	////////////////////////////////////////////////////////////////////////////////
	///////////////////////[]得到P F阵
	double P[MAXNUMSV*MAXNUMSV];//household变换中消除钟差项时用到的正交P~矩阵，维数（m-1)m
	memset(P, 0, sizeof(double)*MAXNUMSV*MAXNUMSV);
	CalculatePEk2(P, m1);


	double F[MAXNUMSV*MAXNUMSV];// P~= FD，F维数（m-1)(m-1)
	memset(F, 0, sizeof(double)*MAXNUMSV*MAXNUMSV);
	GetFMatrix(F, m1);//P~=FD




	double PhatP[MAXNUMSV];//PhatP=P~y
	memset(PhatP, 0, sizeof(double)*MAXNUMSV);
	CalculatePhatP(PhatP, P, SDP, m1);

	double Q[MAXNUMSV*MAXNUMSV];//Q[m-1,m-1]//Q应该是转置过的,Q*PhatP=Rzt
	double U[MAXNUMSV* MAXNUMSV];  //U[m-1,m-1,].///QT=(U,V)T
	double R[MAXNUMSV* MAXNUMSV];

	householderQR(F, Q, m1 - 1, m1 - 1);//phatp=Fz,对F进行QR分解
	GetTUVRMatrix(F, Q, U, R, m1);//F分解得到U,R
	CalculatetwBlock(twBlock, U, PhatP, m1);//twBlock=U*phatp

	///////计算浮点解和方差
	double zt[MAXNUMSV];//浮点解
	double Qat[MAXNUMSV*MAXNUMSV];//浮点解的协方差矩阵
	double ztLambda[5 * MAXNUMSV];
	memset(ztLambda, 0, sizeof(double)* 5 * MAXNUMSV);//
	memset(zt, 0, sizeof(double)*MAXNUMSV);
	memset(Qat, 0, sizeof(double)*MAXNUMSV*MAXNUMSV);
	RecursiveX(zt, twBlock, R, m1 - 1); //递归求解x变量。
	CalculateQat(Qat, R, m1 - 1);
	memcpy(ztLambda, zt, sizeof(double)*MAXNUMSV);//ztlambda为浮点解
	//////////////lambda算法计算整周模糊度
	pglobal->pResult->blambda = Lambda(ztLambda, Qat, m1 - 1);//ztlambda-模糊度浮点解，Qat-浮点解方差协方差阵。 S_line-整周维数	

	double Dtime[MAXNUMSV];//接收机钟差
	CalculateDtime(SDP, ztLambda, Dtime, P, F, m1 - 1);
}