// CloclError.cpp : ���� DLL �ĳ�ʼ�����̡�
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
//TODO:  ����� DLL ����� MFC DLL �Ƕ�̬���ӵģ�
//		��Ӵ� DLL �������κε���
//		MFC �ĺ������뽫 AFX_MANAGE_STATE ����ӵ�
//		�ú�������ǰ�档
//
//		����:
//
//		extern "C" BOOL PASCAL EXPORT ExportedFunction()
//		{
//			AFX_MANAGE_STATE(AfxGetStaticModuleState());
//			// �˴�Ϊ��ͨ������
//		}
//
//		�˺������κ� MFC ����
//		������ÿ��������ʮ����Ҫ��  ����ζ��
//		��������Ϊ�����еĵ�һ�����
//		���֣������������ж������������
//		������Ϊ���ǵĹ��캯���������� MFC
//		DLL ���á�
//
//		�й�������ϸ��Ϣ��
//		����� MFC ����˵�� 33 �� 58��
//

// CCloclErrorApp

BEGIN_MESSAGE_MAP(CCloclErrorApp, CWinApp)
END_MESSAGE_MAP()


// CCloclErrorApp ����

CCloclErrorApp::CCloclErrorApp()
{
// TODO:  �ڴ˴���ӹ�����룬
// ��������Ҫ�ĳ�ʼ�������� InitInstance ��
}


// Ψһ��һ�� CCloclErrorApp ����

CCloclErrorApp theApp;


// CCloclErrorApp ��ʼ��

BOOL CCloclErrorApp::InitInstance()
{
CWinApp::InitInstance();

return TRUE;
}









*********************************************************************************************************/

double SDP[MAXNUMSV];//�ز���λ����
double R[MAXNUMSV*MAXNUMSV];//twBlock=Rzt��Ϊ�����Ǿ���
double twBlock[MAXNUMSV*MAXNUMSV];//twBlock=U*phatP
double ZV[180];
double threshold;
/*  wavelen*SDP=c*RecedTime-wavelen*NdTime+v
*****�����任�������ջ��Ӳ�RecedTtime
*****SDPΪ�ز���λС�����ֲ�ֵ
*****NdTimeΪ���ջ��Ӳ���ɵ�����ģ����,ztΪ�両���
*/

void GetMatrixSDP(double *SDP, TORGOBSDATA *m_ptMMDataNew, TORGOBSDATA *m_ptS1MDataNew, double *ZV)
{
	int i = 0;
	int zv = 0;
	for (i = 0; i<m_ptMMDataNew->u1SVBlocks; i++)
	{
		zv = m_ptMMDataNew->tSv[i].u1PID;
		SDP[i] = m_ptS1MDataNew->tSv[i].dPhase - m_ptMMDataNew->tSv[i].dPhase - ZV[zv];//�ز���λ����С������
	}    //	ZV[zv] = (int)(m_ptS1MDataNew->tSv[i].dPhase - m_ptMMDataNew->tSv[i].dPhase)
}

void CalculatePEk2(double *P, int m)//household�任�������Ӳ���ʱ�õ���P~����
{
	double a = -1. / (m - sqrt((double)m));
int i = 0, j = 0;
for (i = 0; i<m - 1; i++)
{
	for (j = 0; j<m; j++)
	{
		//~P��
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

void CalculateQat(double *Qat, double *R, int m)//������Э�������
{
	double delta = 1 / (2 * pglobal->InitPara.fPhasePrecision1*pglobal->InitPara.fPhasePrecision1);

	int i = 0, j = 0, k = 0; unsigned char bInvert;

	double SZ_[MAXNUMSV*MAXNUMSV];//Z�任����

	memset(SZ_, 0, sizeof(double)*MAXNUMSV*MAXNUMSV);
	memcpy(SZ_, R, sizeof(double)*m*m);

	for (i = 0; i<m; i++)
	{
		for (j = i + 1; j<m; j++)
		{
			//SZ_ =SZ^ ת��
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
				Qat[i*m + j] += delta*SZ_[i*m + k] * R[k*m + j];//������Э�������
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

	//����Q��
	memset(Q, 0, sizeof(double)*m*m);
	for (i = 0; i<m; i++)
		Q[i*m + i] = 1.0; //Q = I^m.

	//////��һ��Լ�� ʹ��С���е����Ҳ������
	int t;
	if (m>n || m == n)
		t = n;
	else
		t = m;

	for (j = 0; j<t; j++)   //for j=1:n.
	{
		//�쳣����
		if (m - j == 1) continue;

		memset(v, 0, sizeof(double)*MAXDIM);
		memset(x, 0, sizeof(double)*MAXDIM);
		for (i = j; i<m; i++)
		{
			x[i - j] = X[i*n + j];  //A[j:m,j].
		}
		beta = house(x, v, m - j);//�õ�House�任����������

		//update A[j:m,j:n] = (I^(m-j/*+1*/) - beta*v*vt)A(j:m,j:n).
		//(1) :���� (I^(m-j/*+1*/) - beta*v*vt).
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
		//(2) ����C = B *A .
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

		//��5������Q��


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
//	��������:	(Householder����������x	& R^n ,�������������� v(1) =1   //
//				��v &R^n �� beta & R ʹ�� P = I^n - beta*v*vt �������� //
//				��Px = norm_2(x)*e1.									//
//				�μ���������㡷������ �㷨 5.1.1.						//
//	�� �� ֵ:	beta													//
//	��    ��:	x ��v �� nά������										//
//	�޶���ʷ:														  	//
//----------------------------------------------------------------------//
double house(const double *x, double *v, int n)
{
	//[1] ������X���м�Ȩ���㣬�������硣
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
			v[0] = -a / (X[0] + c); //parlett(1971�꣩

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
	//U�� [m-1,m-1]
	for (i = 0; i<(m - 1)*(m - 1); i++)
	{
		U[i] = T[i];
	}
	//R ��[m-1,m-1]��
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
//	��������:	�ݹ����x������											//
//	�� �� ֵ:	None.													//
//	��    ��:	x^n --һά����	Rx = b									//
//				����R Ϊ�����Ǿ���									//
//	�޶���ʷ:														  	//
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

////////��QR�ֽ���z����N
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

	householderQR(P, T, m - 1, m);//P����QR�ֽ�****************************����ά��
	GetTUVRMatrix(P, T, Up, Rp, m);//P�ֽ�õ�U,R****************************����ά��

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
	"RTK", "20140715", "RTK��λ�㷨ģ��"

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
	printf("��ģ�� RTK ����������\n");
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

	int m1;//��ǰ��Ԫ��������
	double DD_Phase[MAXNUMSV];//˫����λ
	double DD_Pseudo[MAXNUMSV];//˫��α��
	double Nt[MAXNUMSV];//�Ӳ����ܾ���ֵ
	double Ntfloat[MAXNUMSV];//��ǰ�Ӳ����ܼ�ȥ��ʼ����
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
	m_ptMMDataNew->tSv[i].dPhase_float = fmod(m_ptMMDataNew->tSv[i].dPhase, 1.0);//fmod�����������һ��ָ���������Եڶ���ָ��������������������math.h
	m_ptS1MDataNew->tSv[i].dPhase_float = fmod(m_ptS1MDataNew->tSv[i].dPhase, 1.0);//�õ��ز���λС������
	}
	*/

	int o = 0;
	for (int i = 1; i < m1; i++)

	{
		/*DD_Phase[o] = (m_ptMMDataNew->tSv[i].dPhase - m_ptS1MDataNew->tSv[i].dPhase) - (m_ptMMDataNew->tSv[0].dPhase - m_ptS1MDataNew->tSv[0].dPhase);*/
		DD_Phase[o] = (m_ptMMDataNew->tSv[i].dPhase_float - m_ptS1MDataNew->tSv[i].dPhase_float) - (m_ptMMDataNew->tSv[0].dPhase_float - m_ptS1MDataNew->tSv[0].dPhase_float);//˫���ز���λ
		DD_Pseudo[o] = (m_ptMMDataNew->tSv[i].dPseudo - m_ptS1MDataNew->tSv[i].dPseudo) - (m_ptMMDataNew->tSv[0].dPseudo - m_ptS1MDataNew->tSv[0].dPseudo);//˫��α�ࣨ��ǰα����ȥ��ʼα��
		Ntfloat[o] = DD_Pseudo[o] / pglobal->InitPara.SAT[0].dWAVELEN - DD_Phase[o];//˫������
		/*Ntfloat[o] = DD_Phase[o]-DD_Pseudo[o] / pglobal->InitPara.SAT[0].dWAVELEN;*/
		Nt[o] = fabs(Ntfloat[o]);//fabs:�󸡵����ľ���ֵ
		threshold_1[o] = (int)(Nt[o] + threshold_0);
		threshold_2[o] = (int)(Nt[o] - threshold_0);//���ܵ����½�
		o = o + 1;
	}

	///////////////////////[7] ��¼����������ǵ��ز���������  ��ȥ��ʼ����
	double SDP[MAXNUMSV];
	int tempVs[MAXNUMSV];
	int VsCurrent[MAXNUMSV];
	memset(tempVs, 0, sizeof(int)*MAXNUMSV);
	memset(VsCurrent, 0, sizeof(int)*MAXNUMSV);
	int i;
	for (i = 0; i < m_ptMMData->u1SVBlocks; i++)
	{
		tempVs[i] = m_ptMMData->tSv[i].u1PID;//��վα�����������
		VsCurrent[i] = m_ptMMData->tSv[i].u1PID;
	}
	memset(SDP, 0, sizeof(double)*MAXNUMSV);
	{
		for (i = 0; i < m1; i++)
		{
			int zv = 0;
			zv = VsCurrent[i];
			ZV[zv] = (int)(m_ptS1MDataNew->tSv[i].dPhase - m_ptMMDataNew->tSv[i].dPhase);//�µ��ز���λ��ֵȡ��
		}
	}
	GetMatrixSDP(SDP, m_ptMMDataNew, m_ptS1MDataNew, ZV);//���SDPΪ�ز���λ��С������


	////////////////////////////////////////////////////////////////////////////////
	///////////////////////[]�õ�P F��
	double P[MAXNUMSV*MAXNUMSV];//household�任�������Ӳ���ʱ�õ�������P~����ά����m-1)m
	memset(P, 0, sizeof(double)*MAXNUMSV*MAXNUMSV);
	CalculatePEk2(P, m1);


	double F[MAXNUMSV*MAXNUMSV];// P~= FD��Fά����m-1)(m-1)
	memset(F, 0, sizeof(double)*MAXNUMSV*MAXNUMSV);
	GetFMatrix(F, m1);//P~=FD




	double PhatP[MAXNUMSV];//PhatP=P~y
	memset(PhatP, 0, sizeof(double)*MAXNUMSV);
	CalculatePhatP(PhatP, P, SDP, m1);

	double Q[MAXNUMSV*MAXNUMSV];//Q[m-1,m-1]//QӦ����ת�ù���,Q*PhatP=Rzt
	double U[MAXNUMSV* MAXNUMSV];  //U[m-1,m-1,].///QT=(U,V)T
	double R[MAXNUMSV* MAXNUMSV];

	householderQR(F, Q, m1 - 1, m1 - 1);//phatp=Fz,��F����QR�ֽ�
	GetTUVRMatrix(F, Q, U, R, m1);//F�ֽ�õ�U,R
	CalculatetwBlock(twBlock, U, PhatP, m1);//twBlock=U*phatp

	///////���㸡���ͷ���
	double zt[MAXNUMSV];//�����
	double Qat[MAXNUMSV*MAXNUMSV];//������Э�������
	double ztLambda[5 * MAXNUMSV];
	memset(ztLambda, 0, sizeof(double)* 5 * MAXNUMSV);//
	memset(zt, 0, sizeof(double)*MAXNUMSV);
	memset(Qat, 0, sizeof(double)*MAXNUMSV*MAXNUMSV);
	RecursiveX(zt, twBlock, R, m1 - 1); //�ݹ����x������
	CalculateQat(Qat, R, m1 - 1);
	memcpy(ztLambda, zt, sizeof(double)*MAXNUMSV);//ztlambdaΪ�����
	//////////////lambda�㷨��������ģ����
	pglobal->pResult->blambda = Lambda(ztLambda, Qat, m1 - 1);//ztlambda-ģ���ȸ���⣬Qat-����ⷽ��Э������ S_line-����ά��	

	double Dtime[MAXNUMSV];//���ջ��Ӳ�
	CalculateDtime(SDP, ztLambda, Dtime, P, F, m1 - 1);
}