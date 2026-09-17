#include "stdafx.h"
#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <fstream>
#include <random>
using namespace std;

//constexpr int D = 100;
//const int NP = 1615;
//constexpr int D = 50;
//const int NP = 1018;
constexpr int D = 30;
const int NP = 724;
const int Ainit = NP * 7 / 10;
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
double lb, ub;
double Fi[NP];
double Cri[NP];
double Mcr[NP];
double MF[NP];
double sucess_F[NP];
double sucess_Cr[NP];
double delta_f[NP];
double sumOfdelta_f;
double pb = 0.085;
double mean_Cr, mean_F;
int best, worst;
int bestp[NP];
int worstA[Ainit];
int dy_NP, ANum, pbSize;
int sucess, his_k = 1, archive = 0, his_f = 0;
int FEs, FEsmin, FEsmean, modify_count, SucNum;
double mean, Std, SucRate;
double Accept = 0;

random_device rd;
mt19937 gen(rd());

typedef double (*FunctionCallback)(double sol[D]);
char FILENAME[50] = ("C17_RAND_50.txt");
ofstream fft(FILENAME);

#include "CEC2017.h"

FunctionCallback function = &f1;

template<typename T>
T Max(T a, T b) { return a > b ? a : b; }

void initRandom(int index)
{
	int j;
	double r;
	for (j = 0; j < D; j++) {
		r = ((double)rand() / (double)(RAND_MAX));
		X[index][j] = r * (ub - lb) + lb;
		solution[j] = X[index][j];
	}
	f[index] = function(solution); FEs++;
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

void pbestsort()
{
	int i, j, temp;
	pbSize = Max(2, (int)(pb * dy_NP));
	for (i = 0; i < dy_NP; i++)
		bestp[i] = i;
	for (i = 0; i < dy_NP - 1; i++) {
		for (j = i + 1; j < dy_NP; j++) {
			if (f[bestp[i]] > f[bestp[j]]) {
				temp = bestp[i];
				bestp[i] = bestp[j];
				bestp[j] = temp;
			}
		}
	}
}

void DEMutation()
{
	int i, j, r2_range, np_range;
	int neighbour1, neighbour2, neighbour3, rand3, temp;
	double r, Fa;

	for (i = 0; i < dy_NP; i++) {
		MF[i] = his_mem[his_f][0];
		Mcr[i] = his_mem[his_f][1];

		Fa = cauchy_dis(MF[i], 0.1);
		while (Fa <= 0) {
			Fa = cauchy_dis(MF[i], 0.1);
		}
		if (FEs < 0.25 * maxFEs && Fa < 0.5)
			Fa = 0.5;
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

		if (FEs < 0.3 * maxFEs)
		{
			r = (double)rand() / (double)RAND_MAX;
			np_range = Max(3, (int)(0.7 * dy_NP));
			neighbour2 = bestp[(int)(r * np_range) % np_range];
			while (neighbour2 == i || neighbour2 == neighbour1) {
				r = (double)rand() / (double)RAND_MAX;
				neighbour2 = bestp[(int)(r * np_range) % np_range];
			}
		}
		else
		{
			r = (double)rand() / (double)RAND_MAX;
			np_range = Max(3, (int)((0.5 - 0.2 * (double)(FEs / maxFEs)) * dy_NP));
			neighbour2 = bestp[(int)(r * np_range) % np_range];
			while (neighbour2 == i || neighbour2 == neighbour1) {
				r = (double)rand() / (double)RAND_MAX;
				neighbour2 = bestp[(int)(r * np_range) % np_range];
			}
		}

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
		while (rand3 == i || rand3 == neighbour1 || rand3 == neighbour2) {
			r = ((double)rand() / (double)RAND_MAX);
			rand3 = (int)(r * dy_NP) % dy_NP;
		}

		for (j = 0; j < D; j++) {
			if (FEs < 0.3 * maxFEs)
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

void DECrossover()
{
	int i, j, k;
	double r, c;
	for (i = 0; i < dy_NP; i++) {
		r = ((double)rand() / (double)(RAND_MAX));
		k = (int)(r * D) % D;

		if (Mcr[i] == -1)
			Cri[i] = 0.0;
		else {
			Cri[i] = gaussrand(Mcr[i], 0.1);
			if (FEs < 0.25 * maxFEs)
				Cri[i] = Max(Cri[i], 0.5);
		}

		if (Cri[i] <= 0) Cri[i] = 0;
		if (Cri[i] >= 1) Cri[i] = 1;

		c = 0.3 * (double)FEs / maxFEs;
		for (j = 0; j < D; j++) {
			r = ((double)rand() / (double)(RAND_MAX));
			if (r <= Cri[i] || j == k)
			{
				U[i][j] = (1 - c) * V[i][j] + c * X[i][j];
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
		ObjValSol = function(U[i]); FEs++;
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

				pool_size = Max(1, (int)(ANum * 0.3));
				r = ((double)rand() / (double)RAND_MAX);
				sel_idx = (int)(r * pool_size) % pool_size;
				index_rep = worstA[sel_idx];

				for (int m = 0; m < D; m++) {
					X[index_rep][m] = X[i][m];
				}
				f[index_rep] = f[i];
			}

			sucess_F[sucess] = Fi[i]; sucess_Cr[sucess] = Cri[i]; delta_f[sucess] = f[i] - ObjValSol;
			sumOfdelta_f += delta_f[sucess];
			sucess++;
		}

		if (ObjValSol <= f[i])
		{
			for (j = 0; j < D; j++) {
				X[i][j] = U[i][j];
			}
			f[i] = ObjValSol; SucNum++;
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
		his_k = (his_k + 1) % H;
		his_f = (his_f + 1) % H;
	}
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
	new_NP = int(NP + double((NPmin - NP) * (double)FEs / maxFEs) + 0.5);
	if (new_NP < NPmin) new_NP = NPmin;

	if (new_NP < dy_NP)
	{
		int current_size = dy_NP;
		while (current_size > new_NP)
		{
			int worst_idx = getWorstIndex(current_size);
			elementEx(worst_idx, current_size - 1);
			current_size--;
		}
		dy_NP = new_NP;
	}

	ANum1 = dy_NP * 7 / 10;
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
	Result Red[100];
	fstream ffit;
	mean = 0; Std = 0; index = 1; SucRate = 0;
	FEsmean = 0; modify_count = 0;
	srand(time(NULL));
	function(solution);

	for (run = 0; run < runtime; run++) {
		FEs = 0; FEsmin = 0; SucNum = 0;
		archive = 0; his_k = 1; his_f = 0;
		initial(); dy_NP = NP; ANum = Ainit;
		for (iter = 0; FEs < maxFEs; iter++) {
			pbestsort();
			DEMutation();
			DECrossover();
			DESelection();
			DEUpdatesize();
			MemorizeBestSolution();

			if (run == 5 && FEs / Gen == index) {
				Red[index - 1].FEs = (double)FEs / maxFEs;
				Red[index - 1].fitG = GlobalMin;
				Red[index - 1].MuCr = mean_Cr;
				Red[index - 1].MuF = mean_F;
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
			index++;
		}

		cout.width(2);
		cout << run + 1 << ".run:  " << setiosflags(ios::scientific) << setprecision(2) << "  GlobalMin =  " << GlobalMin << "  Bees: " << best + 1
			<< "	FEsmin = " << FEsmin << "	FEs = " << FEs << "	SucRate =  " << double(SucNum) / double(FEs) << endl;
		GlobalMins[run] = GlobalMin;
		mean = mean + GlobalMin;
		SucRate = SucRate + double(SucNum) / double(FEs);
		if (GlobalMin <= Accept) {
			FEsmean = FEsmean + FEsmin;
			modify_count++;
		}
	}

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
	fstream ffit;
	ffit.open("DEstd_D30CEC17.txt", ios::out);
	ffit << "           ***      DEstd Algorithm Test Result      ***" << endl << endl;
	ffit << "Parameter setting: " << "NP = " << NP << " D = " << D << " MaxFEs = " << maxFEs << "   runtime = " << runtime << endl << endl;

	int i;
	for (i = 0; i < 30; i++) {
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