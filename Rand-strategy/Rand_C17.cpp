// DEStd.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"

#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <fstream>
#include <random>
//#include "fc.h"

using namespace std;

//constexpr int D = 100;
//const int NP = 1615;
//constexpr int D = 50;
//const int NP = 1018;
constexpr int D = 30;
const int NP = 724;

const int Ainit = NP * 15 / 10;
const int H = 5;
const int maxFEs = 10000 * D;
const int Gen = 100 * D;
const int runtime = 51;
const double M_PI = 3.14159265358979323846;

double X[NP + Ainit][D];
double f[NP + Ainit];

double V[NP][D];
double U[NP][D];
double solution[D];
double his_mem[H][2];
double ObjValSol;
double GlobalMin;
double GlobalMins[runtime];
double lb;
double ub;
double Fi[NP];
double Cri[NP];
double Mcr[NP];
double MF[NP];
double sucess_F[NP];
double sucess_Cr[NP];
double delta_f[NP];
double Ri[NP];
double sumOfdelta_f;
double pb = 0.085;
double Dis_max, Dis_min;
double mean_Cr, mean_F;

int best;
int worst;
int bestp[NP];
int bestA[Ainit];
int worstA[Ainit];
int dy_NP;
int ANum;
int pbSize;
int sucess, his_k = 1, archive = 0, his_f = 0;
int FEs, FEsmin, FEsmean, modify_count, SucNum;
double mean, Std, SucRate;
double Accept = 0;   //	1.0e-8;	//  

random_device rd;
mt19937 gen(rd());

typedef double (*FunctionCallback)(double sol[D]);

char FILENAME[50] = ("C17_RAND_30.txt");
ofstream fft(FILENAME);
ofstream diversity;

#include "CEC2017.h"


FunctionCallback function = &f1;

void initRandom(int index)
{
	int j;
	double r;

	r = ((double)rand() / (double)(RAND_MAX));				//第一个随机数有可能一样

	for (j = 0; j < D; j++) {
		r = ((double)rand() / (double)(RAND_MAX));
		X[index][j] = r * (ub - lb) + lb;
		solution[j] = X[index][j];
	}
	f[index] = function(solution);   FEs++;
}

void initial()
{
	int i;

	for (i = 0; i < NP; i++) {
		initRandom(i);
	}
	for (int idx = 0; idx < H; idx++)
	{
		his_mem[idx][0] = 0.3;
		his_mem[idx][1] = 0.8;
	}

	//his_mem[H - 1][0] = 0.9;
	//his_mem[H - 1][1] = 0.9;
}

void MemorizeBestSolution()
{
	int i;

	GlobalMin = f[0];
	for (i = 1; i < dy_NP; i++) {
		if (f[i] < GlobalMin)
		{
			GlobalMin = f[i];
			best = i;
		}
	}
}

double gaussrand(double E, double V)
{
	static double V1, V2, S;
	static int phase = 0;
	double X;

	if (phase == 0) {
		do {
			double U1 = (double)rand() / RAND_MAX;
			double U2 = (double)rand() / RAND_MAX;
			V1 = 2 * U1 - 1;
			V2 = 2 * U2 - 1;
			S = V1 * V1 + V2 * V2;
		} while (S >= 1 || S == 0);
		X = V1 * sqrt(-2 * log(S) / S);
	}
	else {
		X = V2 * sqrt(-2 * log(S) / S);
	}
	phase = 1 - phase;
	X = X * V + E;
	return X;
}

double cauchy_dis(double mu, double gamma)
{
	return mu + gamma * tan(M_PI * ((double)rand() / RAND_MAX - 0.5));
}

const double pmax = 0.25;
const double pmin = pmax / 2.0;


#include <algorithm>   // for std::sort, std::partial_sort

void pbestsort()
{
	// 计算 pb 和 pbSize 保持不变
	pbSize = max(2, (int)(pb * dy_NP));

	// 1. 对 bestp 数组（长度 dy_NP）按 f[索引] 升序排序
	//    原逻辑：bestp 初始化为 [0, 1, 2, ..., dy_NP-1]
	for (int i = 0; i < dy_NP; ++i)
		bestp[i] = i;

	std::sort(bestp, bestp + dy_NP,
		[&](int a, int b) { return f[a] < f[b]; });

	// 2. 对 bestA 数组（长度 archive）部分排序，只要求前 60% 元素有序（最小的 60%）
	//    原代码外层循环只到 0.6*archive，实现了部分排序的效果，
	//    但内层循环却遍历整个数组，导致后面 40% 元素可能被交换到前面但并未完全有序。
	//    这里使用 std::partial_sort 更清晰地表达意图。
	//for (int i = 0; i < archive; ++i)
	//	bestA[i] = NP + i;

	//int topK = static_cast<int>(0.4 * archive);
	//if (topK > 0 && topK <= archive) {
	//	std::partial_sort(bestA, bestA + topK, bestA + archive,
	//		[&](int a, int b) { return f[a] < f[b]; });
	//}
	//else if (topK > archive) {
	//	// 理论上不会发生，安全起见：全排序
	//	std::sort(bestA, bestA + archive,
	//		[&](int a, int b) { return f[a] < f[b]; });
	//}
}

void DEMutation()
{
	int i, j, r2_range, indexH, np_range;
	int neighbour1, neighbour2, neighbour3, rand3, temp;
	double r, Fw, Dis_r, Dis_i, fit_norm, fpbest, frand, Fa;

	for (i = 0; i < dy_NP; i++) {

		//r = (double)rand() / (double)RAND_MAX;
		//indexH = (int)(r * H) % H;
		MF[i] = his_mem[his_f][0];
		Mcr[i] = his_mem[his_f][1];
		//MF[i] = his_mem[indexH][0];
		//Mcr[i] = his_mem[indexH][1];

		Fa = cauchy_dis(MF[i], 0.1);
		while (Fa <= 0) {
			Fa = cauchy_dis(MF[i], 0.1);
		}

		if (FEs < 0.25 * maxFEs && Fa < 0.5)
			Fa = 0.5;


		//if (FEs > 0.6 * maxFEs && Fa < 0.6)
		//	Fa = 0.6;


		if (Fa >= 1) Fa = 1;

		Fi[i] = Fa;


		r = (double)rand() / (double)RAND_MAX;
		temp = (int)(r * pbSize) % pbSize;
		neighbour1 = bestp[temp];
		while (neighbour1 == i) {
			r = (double)rand() / (double)RAND_MAX;
			temp = (int)(r * pbSize) % pbSize;
			neighbour1 = bestp[temp];
		}

		if (FEs < 0.5 * maxFEs)
		{
			r = (double)rand() / (double)RAND_MAX;
			np_range = max(3, (int)(0.7 * dy_NP));
			neighbour2 = bestp[(int)(r * np_range) % np_range];
			while (neighbour2 == i || neighbour2 == neighbour1) {
				r = (double)rand() / (double)RAND_MAX;
				neighbour2 = bestp[(int)(r * np_range) % np_range];
			}
		}
		else
		{

			r = (double)rand() / (double)RAND_MAX;
			np_range = max(3, (int)((0.5 - 0.2 * (double)(FEs / maxFEs)) * dy_NP));
			neighbour2 = bestp[(int)(r * np_range) % np_range];
			while (neighbour2 == i || neighbour2 == neighbour1) {
				r = (double)rand() / (double)RAND_MAX;
				neighbour2 = bestp[(int)(r * np_range) % np_range];
			}
		}

		//r = (double)rand() / (double)RAND_MAX;
		//neighbour2 = (int)(r * dy_NP) % dy_NP;
		//while (neighbour2 == i || neighbour2 == neighbour1){
		//	r = (double)rand() / (double)RAND_MAX;
		//	neighbour2 = (int)(r * dy_NP) % dy_NP;
		//}


		//r2_range = dy_NP + 0.4 * archive;
		//r = (double)rand() / (double)RAND_MAX;
		//neighbour3 = (int)(r * r2_range) % r2_range;
		//while (neighbour3 == i || neighbour3 == neighbour1 || neighbour3 == neighbour2){
		//	r = (double)rand() / (double)RAND_MAX;
		//	neighbour3 = (int)(r * r2_range) % r2_range;
		//} 
		//if (neighbour3 >= dy_NP)
		//{
		//	neighbour3 = bestA[neighbour3 - dy_NP];
		//}

		r2_range = dy_NP + archive;
		r = (double)rand() / (double)RAND_MAX;
		neighbour3 = (int)(r * r2_range) % r2_range;
		while (neighbour3 == i || neighbour3 == neighbour1 || neighbour3 == neighbour2) {
			r = (double)rand() / (double)RAND_MAX;
			neighbour3 = (int)(r * r2_range) % r2_range;
		}
		if (neighbour3 >= dy_NP)
		{
			neighbour3 += NP - dy_NP;
		}

		r = (double)rand() / RAND_MAX;
		rand3 = (int)(r * dy_NP) % dy_NP;
		while (rand3 == neighbour3 || rand3 == neighbour2 || rand3 == i || rand3 == neighbour1) {
			r = ((double)rand() / (double)RAND_MAX);
			rand3 = (int)(r * dy_NP) % dy_NP;
		}

		//if (FEs < 0.25 * maxFEs) {
		//	Fw = 0.8 * Fa;
		//}
		//else {
		//	if (FEs < 0.4 * maxFEs) {
		//		Fw = 0.9 * Fa;
		//	}
		//	else {
		//		Fw = 1.5 * Fa;
		//	}
		//}

		//Fw = 0.3 + 0.5 * (double)FEs / maxFEs;


		for (j = 0; j < D; j++) {


			if (FEs < 0.6 * maxFEs)
			{
				V[i][j] = X[rand3][j] + Fa * (X[neighbour1][j] - X[rand3][j]) + Fa * (X[neighbour2][j] - X[neighbour3][j]);
			}
			else
			{
				if (f[neighbour2] < f[neighbour3])
					V[i][j] = X[rand3][j] + Fa * (X[neighbour1][j] - X[rand3][j]) + Fa * (X[neighbour2][j] - X[neighbour3][j]);
				else
					V[i][j] = X[rand3][j] + Fa * (X[neighbour1][j] - X[rand3][j]) + Fa * (X[neighbour3][j] - X[neighbour2][j]);
			}





			while (V[i][j] < lb)
				V[i][j] = (lb + X[i][j]) / 2.0;
			while (V[i][j] > ub)
				V[i][j] = (ub + X[i][j]) / 2.0;

		}
	}
}

double Diversity() {
	int i, j;
	double avg[D], var[D], Dvar, sum;

	// 防止除零（ub==lb 时最大值方差为0）
	double range = ub - lb;
	double Var_max = (range * range) / 4.0;
	if (Var_max == 0.0) return 0.0;   // 所有维无多样性

	// 均值
	for (j = 0; j < D; j++) {
		sum = 0.0;
		for (i = 0; i < dy_NP; i++) sum += X[i][j];
		avg[j] = sum / dy_NP;
	}

	// 方差
	for (j = 0; j < D; j++) {
		sum = 0.0;
		for (i = 0; i < dy_NP; i++) {
			double diff = X[i][j] - avg[j];
			sum += diff * diff;          // 比 pow 更快
		}
		var[j] = sum / dy_NP;
	}

	// 归一化
	sum = 0.0;
	for (j = 0; j < D; j++) sum += var[j] / Var_max;
	Dvar = sum / D;
	return Dvar;
}

void DECrossover()
{
	int i, j;
	int k;
	double r;


	for (i = 0; i < dy_NP; i++) {

		r = ((double)rand() / (double)(RAND_MAX));
		k = (int)(r * D) % D;

		if (Mcr[i] == -1)
			Cri[i] = 0.0;
		else {
			Cri[i] = gaussrand(Mcr[i], 0.1);

			if (FEs < 0.25 * maxFEs)
				Cri[i] = max(Cri[i], 0.5);
			//else if (FEs < 0.7 * maxFEs)
			//	Cri[i] = max(Cri[i], 0.6);


		}

		if (Cri[i] <= 0)
			Cri[i] = 0;
		if (Cri[i] >= 1)
			Cri[i] = 1;

		for (j = 0; j < D; j++) {
			r = ((double)rand() / (double)(RAND_MAX));
			if (r <= Cri[i] || j == k)
			{

				double c = 0.3 * (double)FEs / maxFEs;
				U[i][j] = (1 - c) * V[i][j] + c * X[i][j];
				//U[i][j] = V[i][j];

			}
			else {
				U[i][j] = X[i][j];
			}
		}
	}
}

void DESelection()
{
	int i, j, index, temp, is_sort, pool_size, sel_idx, index_rep;
	double wk, r, k1 = 0.0, k2 = 0.0, t1 = 0.0, t2 = 0.0;
	sumOfdelta_f = 0; sucess = 0;

	is_sort = false;
	for (i = 0; i < dy_NP && FEs < maxFEs; i++) {

		ObjValSol = function(U[i]);		FEs++;
		if (ObjValSol < f[i]) {

			if (archive < ANum)
			{
				index = archive + NP;
				for (int idx = 0; idx < D; idx++) {
					X[index][idx] = X[i][idx];
				}
				f[index] = f[i];
				archive++;
			}
			else
			{

				if (!is_sort)
				{
					for (int m = 0; m < ANum; m++)
					{
						worstA[m] = NP + m;
					}

					for (int m = 0; m < ANum * 0.3; m++)
					{
						for (int n = m + 1; n < ANum; n++)
						{
							if (f[worstA[m]] < f[worstA[n]])
							{
								int temp = worstA[m];
								worstA[m] = worstA[n];
								worstA[n] = temp;

							}
						}
					}
					is_sort = true;
				}

				pool_size = max(1, (int)(ANum * 0.3));
				r = ((double)rand() / (double)RAND_MAX);
				sel_idx = (int)(r * pool_size) % pool_size;

				index_rep = worstA[sel_idx];

				for (int m = 0; m < D; m++) {
					X[index_rep][m] = X[i][m];
				}
				f[index_rep] = f[i];
			}

			sucess_F[sucess] = Fi[i];	sucess_Cr[sucess] = Cri[i];	delta_f[sucess] = f[i] - ObjValSol;
			sumOfdelta_f += delta_f[sucess];
			sucess++;

		}

		if (ObjValSol <= f[i])
		{
			for (j = 0; j < D; j++) {
				X[i][j] = U[i][j];
			}
			f[i] = ObjValSol;	SucNum++;
		}

		if (ObjValSol <= Accept) {
			FEsmin = FEs;
		}

	}

	if (sucess != 0)
	{
		for (int idx = 0; idx < sucess; idx++)
		{

			wk = (double)delta_f[idx] / sumOfdelta_f;
			t1 += (double)wk * sucess_Cr[idx];
			t2 += (double)wk * sucess_Cr[idx] * sucess_Cr[idx];
			k1 += (double)wk * sucess_F[idx] * sucess_F[idx];
			k2 += (double)wk * sucess_F[idx];

		}
		if (his_k != H)
		{

			//if (FEs < 0.3 * maxFEs)
			//{
			//	if (his_mem[his_k][1] == -1 || t1 == 0)
			//	{
			//		his_mem[his_k][1] = -1;
			//	}
			//	else
			//	{
			//		mean_Cr = t2 / t1;
			//		his_mem[his_k][1] = mean_Cr;
			//	}
			//	mean_F = k1 / k2;
			//	his_mem[his_k][0] = mean_F;
			//}
			//else
			//{
			//	if (his_mem[his_k][1] == -1 || t1 == 0)
			//	{
			//		his_mem[his_k][1] = -1;
			//	}
			//	else
			//	{
			//		mean_Cr = t2 / t1;
			//		his_mem[his_k][1] = (mean_Cr + his_mem[his_k][1]) / 2.0;
			//	}
			//	mean_F = k1 / k2;
			//	his_mem[his_k][0] = (mean_F + his_mem[his_f][0])/ 2.0;
			//}

			if (his_mem[his_k][1] == -1 || t1 == 0)
			{
				his_mem[his_k][1] = -1;
			}
			else
			{
				mean_Cr = t2 / t1;
				his_mem[his_k][1] = (mean_Cr + his_mem[his_f][1]) / 2.0;
			}
			mean_F = k1 / k2;
			his_mem[his_k][0] = (mean_F + his_mem[his_f][0]) / 2.0;

		}

		his_k = (his_k + 1) % (H);
		his_f = (his_f + 1) % (H);

	}
	//else {
	//	r = (double)rand() / RAND_MAX;
	//	temp = (int)(r * H) % H;
	//	his_mem[his_k][1] = (his_mem[temp][1] + his_mem[his_k][1]) / 2.0;
	//	his_mem[his_k][0] = (his_mem[temp][0] + his_mem[his_k][0]) / 2.0;

	//}


}

int getWorstIndex(int max_idx)
{
	int w = 0;
	for (int idx = 1; idx < max_idx; idx++)
	{
		if (f[idx] > f[w])
			w = idx;
	}
	return w;

}

void elementEx(int org, int des)
{

	if (org == des)
		return;

	for (int idx = 0; idx < D; idx++)
	{
		double temp = X[org][idx];
		X[org][idx] = X[des][idx];
		X[des][idx] = temp;
	}

	double tmpf = f[org];
	f[org] = f[des];
	f[des] = tmpf;

}
void MemorizeWorstArchive()
{
	int m;
	worst = NP;

	for (m = NP; m < NP + archive; m++) {
		if (f[m] > f[worst]) {
			worst = m;
		}
	}
}

void DEUpdatesize()
{

	int new_NP, NPmin = 5, ANum1, i, MoveNum;
	// 计算新的种群大小
	//new_NP = round((double)(NPmin - NP) / maxFEs * FEs + NP);
	new_NP = int(NP + double((NPmin - NP) * (double)FEs / maxFEs) + 0.5);

	if (new_NP < NPmin) new_NP = NPmin; // 安全检查

	if (new_NP < dy_NP)
	{
		// 循环移除多余个体
		// 每次在当前有效范围内找到最差的，交换到当前有效范围的末尾
		int current_size = dy_NP;
		while (current_size > new_NP)
		{
			int worst_idx = getWorstIndex(current_size);
			elementEx(worst_idx, current_size - 1);
			current_size--;
		}

		dy_NP = new_NP;
	}

	ANum1 = dy_NP * 15 / 10; // 更新存档大小上限
	if (archive >= ANum1)
	{
		MoveNum = archive - ANum1;
		for (i = 0; i < MoveNum; i++)
		{
			MemorizeWorstArchive();
			elementEx(worst, NP + archive - 1);
			archive--;
		}
	}
	ANum = ANum1;

	pb = 0.085 + 0.085 * double((double)FEs / maxFEs);

}

void DEstd(int Fun)
{
	int iter, run, index, i;
	struct Result {
		double FEs;
		double fitG;
		double MuCr;
		double MuF;
		double Diversity;
	};
	Result	Red[100];
	fstream ffit;
	mean = 0;    Std = 0;		index = 1;	SucRate = 0;
	FEsmean = 0; modify_count = 0;
	srand(time(NULL));
	function(solution);

	for (run = 0; run < runtime; run++) {
		FEs = 0;	FEsmin = 0;		SucNum = 0;
		archive = 0;    his_k = 1;   his_f = 0;
		initial();	dy_NP = NP;		ANum = Ainit;
		for (iter = 0; FEs < maxFEs; iter++) {

			pbestsort();
			DEMutation();
			DECrossover();
			DESelection();
			DEUpdatesize();
			MemorizeBestSolution();

			//if (FEs > 0.9 * maxFEs)
			//	cout << 'B';

			if (run == 5 && FEs / Gen == index) {
				Red[index - 1].FEs = (double)FEs / maxFEs;
				Red[index - 1].fitG = GlobalMin;
				Red[index - 1].MuCr = mean_Cr;
				Red[index - 1].MuF = mean_F;
				Red[index - 1].Diversity = Diversity();
				index++;
			}
			if (FEsmin != 0)
				break;
		}
		MemorizeBestSolution();
		if (run == 5 && FEs < maxFEs) {
			Red[index - 1].FEs = FEs;
			Red[index - 1].fitG = GlobalMin;
			Red[index - 1].MuCr = mean_Cr;
			Red[index - 1].MuF = mean_F;
			Red[index - 1].Diversity = Diversity();
			index++;
		}

		cout.width(2);
		cout << run + 1 << ".run:  " << setiosflags(ios::scientific) << setprecision(2) << "  GlobalMin =  " << GlobalMin << "  Bees: " << best + 1
			<< "	FEsmin = " << FEsmin << "	FEs = " << FEs << "	SucRate =  " << double(SucNum) / double(FEs) << endl;
		GlobalMins[run] = GlobalMin;
		mean = mean + GlobalMin;
		SucRate = SucRate + double(SucNum) / double(FEs);
		if (GlobalMin <= Accept) {				//  *** result *** 
			FEsmean = FEsmean + FEsmin;
			modify_count++;
		}
	}

	//  ***    Evolutionary process data collection   EPDC ***
	char Filename[100] = "./Data_Collection/DEstd_Test_Fun_";
	char buf[10];

	itoa(Fun, buf, 10);
	strcat(Filename, buf);
	strcat(Filename, ".txt");

	ffit.open(Filename, ios::out);

	ffit << "           ***      DEstd Algorithm Test Result      ***" << endl << endl;
	ffit << "Function " << Fun << endl;
	ffit << setiosflags(ios::scientific) << setprecision(2);
	ffit << "FEs	" << "GlobalMin		" << "MuCr	" << "MuF	    " << "Diversity" << endl;

	for (i = 0; i < index - 1; i++) {
		ffit << Red[i].FEs << "	     " << Red[i].fitG << "       " << Red[i].MuCr << "        " << Red[i].MuF << "          " << Red[i].Diversity << endl;
	}
	ffit << endl;
	ffit.close();
	cout << endl;
	//  ***    Evolutionary process data collection  END    ***



	if (modify_count != 0)
		FEsmean = FEsmean / modify_count;
	mean = mean / runtime;
	SucRate = SucRate / runtime;
	for (run = 0; run < runtime; run++) {
		Std = Std + pow(GlobalMins[run] - mean, 2);
	}
	Std = sqrt(Std / runtime);
	cout << "Means of " << runtime << " runs:  " << setiosflags(ios::scientific) << setprecision(2) << mean << endl;
	cout << "Std of " << runtime << " runs:  " << setiosflags(ios::scientific) << setprecision(2) << Std << endl;
	cout << "Means of FEs = " << FEsmean << endl;
	cout << "The results for the number of accept is  " << modify_count << endl << endl;

	fft << scientific << setprecision(2) << mean << endl;

}


int main()
{
	FunctionCallback functionArr[30] = { &f1, &f2, &f3, &f4, &f5, &f6, &f7, &f8, &f9, &f10, &f11, &f12, &f13, &f14, &f15,
										&f16, &f17, &f18, &f19, &f20, &f21, &f22, &f23, &f24, &f25, &f26, &f27, &f28, &f29, &f30 };
	//FunctionCallback functionArr[12] = { &f1, &f2, &f3, &f4, &f5, &f6, &f7, &f8, &f9, &f10,&f11,&f12 };
	fstream ffit;
	ffit.open("DEstd_D30CEC17.txt", ios::out);
	ffit << "           ***      DEstd Algorithm Test Result      ***" << endl << endl;

	ffit << "Parameter setting: " << "NP = " << NP << " D = " << D << " MaxFEs = " << maxFEs << "   runtime = " << runtime << endl << endl;
	diversity << "FEs     " << "Diversity";

	int i;
	for (i = 0; i < 30; i++)
	{

		function = functionArr[i];
		gen.seed(rd());
		DEstd(i + 1);
		fileflag = 0;
		ffit << "Function :" << i + 1 << " Means of " << runtime << " runs:  " << setiosflags(ios::scientific) << setprecision(2) << mean;
		ffit << "   Std of " << runtime << " runs:  " << Std << setiosflags(ios::scientific) << setprecision(4) << "    SucRate = " << SucRate << endl;
		ffit << setiosflags(ios::scientific) << setprecision(2) << "Mean( Std ):    " << mean << "(" << Std << ")";
		ffit << "   Means of FEs = " << FEsmean << "	The number of accept(0.00) is  " << modify_count << endl << endl;

	}
	ffit.close();
	return 0;

}



//Dis_r = 0.0;
//for (int j = 0; j < D; j++) {
//	Dis_r += (X[rand3][j] - X[neighbour1][j]) * (X[rand3][j] - X[neighbour1][j]);
//}
//Dis_r = sqrt(Dis_r); // 欧式距离
//
//Dis_i = 0.0;
//for (int j = 0; j < D; j++) {
//	Dis_i += (X[i][j] - X[neighbour1][j]) * (X[i][j] - X[neighbour1][j]);
//}
//Dis_i = sqrt(Dis_i); // 欧式距离
//
//frand = f[rand3];
//fpbest = f[neighbour1];
//fit_norm = (frand - fpbest) / (f[i] - fpbest + 1e-10); // 防止分母0

//ofstream record("Dis.txt");
//record << "Dis_r       "<<"Dis_i      "<<"Devide_Dis        " << "F_r      " << "F_i      " << "D_f       " << "fit_norm		" << "Cr         " << "Fw         " << "Fa			" << endl;



//if(FEs > 0.3*maxFEs)
//	record << Dis_r << "		"<<Dis_i<<"       "<<((double)Dis_i/(Dis_i + Dis_r)) << "       " << f[rand3] << "        " << f[i] << "        " << (f[i] - f[rand3]) << "         " << fit_norm << "			" << Mcr[i] << "          " << Fw << "            " << Fi[i] << endl;
//else
//	record << Dis_r << "		" << Dis_i << "       " << ((double)Dis_i / (Dis_i + Dis_r)) << "       " << f[rand3] << "        " << f[i] << "        " << (f[i] - f[rand3]) << "         " << fit_norm << "			" << Mcr[i] << "          " << Fw << "            " << Fi[i] << endl;