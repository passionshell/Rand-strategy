#pragma once

#include"string.h"

double Os[D], y[D], z[D];
double Mr[D][D];
int SS[D], fileflag;
const double INF = 1.0e99;


void rotatefunc(double* x, double* xrot)
{
	int i, j;
	for (i = 0; i < D; i++) {
		xrot[i] = 0;
		for (j = 0; j < D; j++)
			xrot[i] = xrot[i] + x[j] * Mr[i][j];
	}
}

void read_cec2022_file(char* funcnum, int num)
{
	char FName[45], tempc[5];
	int i, j;
	ifstream ffit;

	strcpy(FName, "../CEC2022_input_data/shift_data_");
	strcat(FName, funcnum);
	strcat(FName, ".txt");
	cout << FName << endl;

	ffit.open(FName);
	if (!ffit)
		cout << "Cannot open input file." << endl;

	for (i = 0; i < D; i++)
		ffit >> Os[i];

	ffit.close();

	strcpy(FName, "../CEC2022_input_data/M_");
	strcat(FName, funcnum);
	strcat(FName, "_D");
	itoa(D, tempc, 10); //cout<< tempc << endl;
	strcat(FName, tempc);
	strcat(FName, ".txt");
	cout << FName << endl;

	ffit.open(FName);
	if (!ffit)
		cout << "Cannot open input file." << endl;

	for (i = 0; i < D; i++) {
		for (j = 0; j < D; j++)
			ffit >> Mr[i][j];
	}
	ffit.close();

	if (num == 3) {
		strcpy(FName, "../CEC2022_input_data/shuffle_data_");
		strcat(FName, funcnum);
		strcat(FName, "_D");
		itoa(D, tempc, 10);
		strcat(FName, tempc);
		strcat(FName, ".txt");
		cout << FName << endl;

		ffit.open(FName);
		if (!ffit)
			cout << "Cannot open input file." << endl;

		for (i = 0; i < D; i++)
			ffit >> SS[i];

		ffit.close();
	}
	fileflag = 1;
}


double f1(double sol[D])								/* zakharov */
{
	int i, num;
	double sh_rate, top = 0.0,sum1= 0.0,sum2 = 0.0;
	char func_num[4] = "1";

	lb = -100, ub = 100;
	num = 2;
	if (fileflag == 0) {
		read_cec2022_file(func_num, num);
	}

	sh_rate = 1;
	for (i = 0; i < D; i++) {
		y[i] = (sol[i] - Os[i]) * sh_rate;
	}

	rotatefunc(y, z);
	for (i = 0; i < D; i++) {
		sum1 = sum1 + z[i] * z[i];
		sum2 = sum2 + 0.5 *(i+1)* z[i];
	}
	top = sum1 + pow(sum2,2) + pow(sum2 ,4);

	return top;
}

double f2(double sol[D])								/* Rosenbrock's */			
{
	int i, num;
	double sh_rate, tmp1,tmp2, top;
	char func_num[4] = "2";

	lb = -100, ub = 100;
	num = 2;
	if (fileflag == 0) 
	{
		read_cec2022_file(func_num, num);
	}

	sh_rate = 2.048 / 100.0;
	for (i = 0; i < D; i++) {
		y[i] = (sol[i] - Os[i]) * sh_rate;
	}

	rotatefunc(y, z);
	top = 0.0, tmp1 = 0.0, tmp2 = 0.0;
	z[0] +=  1.0;    
	for (i = 0; i < D - 1; i++) {
		z[i + 1] +=  1.0;
		tmp1 = z[i]* z[i] - z[i + 1];
		tmp2 = z[i] - 1.0;                                       
		top += (100 * tmp1 * tmp1) + (tmp2 * tmp2);
	}

	return top;
}

double f3(double sol[D])								/* Schwefel's 1.2  */
{
	int i, num;
	double sh_rate, tmp, tmp2,vari, top;
	char func_num[4] = "3";

	lb = -100, ub = 100;
	num = 2;
	if (fileflag == 0)
		read_cec2022_file(func_num, num);

	sh_rate = 1.0;
	for (i = 0; i < D ; i++)
	{
		y[i] = (sol[i] - Os[i]) * sh_rate;
	}

	rotatefunc(y, z);
	top = 0.0,tmp = 0.0,tmp2 = 0.0;
	for (i = 0; i < D - 1; i++)
	{
		z[i] = pow(y[i] * y[i] + y[i + 1] * y[i + 1], 0.5);
		tmp = sin(50.0 * pow(z[i], 0.2));
		top += pow(z[i], 0.5) + pow(z[i], 0.5) * tmp * tmp;
	}
	top = top * top / (D - 1) / (D - 1);

	return top;
}


double f4(double sol[D])								/* Noncontinuous Rastrigin's  */
{
	int i, num;
	double sh_rate, tmp1, tmp2, top;
	char func_num[4] = "4";

	lb = -100, ub = 100;
	num = 2;
	if (fileflag == 0)
		read_cec2022_file(func_num, num);


	for (i = 0; i < D; i++)
	{
		y[i] = sol[i];
		if (fabs(y[i] - Os[i]) > 0.5)
			y[i] = Os[i] + floor(2.0 * (y[i] - Os[i]) + 0.5) / 2.0;
	}

	sh_rate = 5.12 / 100.0;
	for (i = 0; i < D; i++)
		y[i] = (y[i] - Os[i]) * sh_rate;
	 
	rotatefunc(y, z);
	
	tmp1 = 0.0, tmp2 = 0.0, top = 0.0;
	for (i = 0; i < D; i++) 
	{
		tmp1 = z[i] * z[i];
		tmp2 = 10.0 * cos(2.0 * M_PI * z[i]);
		top += (tmp1 - tmp2 + 10.0);
	}

	return top;
}

double f5(double sol[D])									/* Levy */		
{
	int i, num;
	double sh_rate, w[D], term1, term3, sum, sum1, wi, newv, top;
	char func_num[4] = "5";

	lb = -100, ub = 100;
	num = 2;
	if (fileflag == 0)
		read_cec2022_file(func_num, num);

	sh_rate = 1.0;
	for (i = 0; i < D; i++)
		y[i] = (sol[i] - Os[i]) * sh_rate;

	rotatefunc(y, z);

	sum1 = 0.0;
	for (i = 0; i < D; i++)
	{
		w[i] = 1.0 + (z[i] - 0.0)/4.0;                                     //论文是1.0 + (z[i] + 1.0)/4.0
	}

	term1 = pow((sin(M_PI * w[0])), 2);
	term3 = pow((w[D - 1] - 1), 2) * (1 + pow((sin(2 * M_PI * w[D - 1])), 2));

	sum = 0.0;

	for (i = 0; i < D - 1; i++)
	{
		wi = w[i];
		newv = pow((wi - 1), 2) * (1 + 10 * pow((sin(M_PI * wi + 1)), 2));
		sum = sum + newv;
	}

	top = term1 + sum + term3;

	return top;
}

double f6(double sol[D])
{
	int i, tmp, num, cf_num = 3;
	char func_num[4] = "6";
	double alpha,sum_z,r2,sh_rate, top;
	int G[3], G_nx[3];
	double Gp[3] = { 0.4, 0.4, 0.2 };
	lb = -100;	ub = 100;

	tmp = 0;
	for (i = 0; i < cf_num - 1; i++) {
		G_nx[i] = ceil(Gp[i] * D);
		tmp += G_nx[i];
	}
	G_nx[cf_num - 1] = D - tmp;

	G[0] = 0;
	for (i = 1; i < cf_num; i++) {
		G[i] = G[i - 1] + G_nx[i - 1];
	}

	num = 3;
	if (fileflag == 0)
		read_cec2022_file(func_num, num);

	for (i = 0; i < D; i++)
		z[i] = sol[i] - Os[i];

	rotatefunc(z, y);

	for (i = 0; i < D; i++)
		z[i] = y[SS[i] - 1];

	top = 0.0;
	sh_rate = 1.0;											 /* Bent_Cigar */
	z[0] = z[0] * sh_rate;
	top = z[0] * z[0];
	for (i = 1; i < G[1]; i++) {
		z[i] = z[i] * sh_rate;
		top = top + pow(10.0, 6.0) * z[i] * z[i];
	}


	sh_rate = 5.0 / 100.0;									/* HGBat, provdided by Hans - Georg Beyer(HGB) */
	alpha = 1.0 / 4.0;
	r2 = 0.0;
	sum_z = 0.0;
	for (i = G[1]; i < G[2]; i++)
	{
		z[i] = z[i] * sh_rate;
		z[i] = z[i] - 1.0;//shift to orgin
		r2 += z[i] * z[i];
		sum_z += z[i];
	}
	top += pow(fabs(pow(r2, 2.0) - pow(sum_z, 2.0)), 2 * alpha) + (0.5 * r2 + sum_z) / G_nx[1] + 0.5;


	sh_rate = 5.12 / 100.0;	
	for (i = G[2]; i < D; i++) {
		z[i] = z[i] * sh_rate;
		top += (z[i] * z[i] - 10.0 * cos(2.0 * M_PI * z[i]) + 10.0);
	}


	return top;
}

double f7(double sol[D])
{
	int i, tmp, num, cf_num = 6;
	char func_num[4] = "7";
	double alpha, sh_rate,r2,sum_z ,top;
	double temp, tmp1, tmp2, tmp3;
	double sum1, sum2;
	int G[6], G_nx[6];
	double Gp[6] = { 0.1,0.2,0.2,0.2,0.1,0.2 };
	lb = -100;	ub = 100;

	tmp = 0;
	for (i = 0; i < cf_num - 1; i++) {
		G_nx[i] = ceil(Gp[i] * D);
		tmp += G_nx[i];
	}
	G_nx[cf_num - 1] = D - tmp;

	G[0] = 0;
	for (i = 1; i < cf_num; i++) {
		G[i] = G[i - 1] + G_nx[i - 1];
	}

	num = 3;
	if (fileflag == 0)
		read_cec2022_file(func_num, num);

	for (i = 0; i < D; i++)
		z[i] = sol[i] - Os[i];

	rotatefunc(z, y);

	for (i = 0; i < D; i++)
		z[i] = y[SS[i] - 1];

	top = 0.0;										//hgbat function
	sh_rate = 5.0 / 100.0;	alpha = 1.0 / 4.0;	r2 = 0.0;
	sum_z = 0.0;	tmp1 = 0.0;

	for (i = 0; i < G[1]; i++)
	{
		tmp1 = z[i] * sh_rate;
		tmp1 = tmp1 - 1.0;//shift to orgin
		r2 += tmp1 * tmp1;
		sum_z += tmp1;
	}
	top += pow(fabs(pow(r2, 2.0) - pow(sum_z, 2.0)), 2 * alpha) + (0.5 * r2 + sum_z) / G_nx[0] + 0.5;


	sh_rate = 5.0 / 100.0;						   /* Katsuura  */
	sum_z = 1.0;	tmp3 = pow(1.0 * G_nx[1], 1.2);
	double t = 0.0;

	for (i = G[1]; i < G[2]; i++)
	{
		t = z[i] * sh_rate;
		temp = 0.0;
		for (int j = 1; j <= 32; j++)
		{
			tmp1 = pow(2.0, j);
			tmp2 = tmp1 * t;
			temp += fabs(tmp2 - floor(tmp2 + 0.5)) / tmp1;
		}
		sum_z *= pow(1.0 + ((i-G[1]) + 1) * temp, 10.0 / tmp3);
	}
	tmp1 = 10.0 / G_nx[1] / G_nx[1];
	top += sum_z * tmp1 - tmp1;


	sh_rate = 1.0;                                     /* Ackley's  */
	sum1 = 0.0, sum2 = 0.0;

	for (i = G[2]; i < G[3]; i++)
	{
		z[i] = z[i] * sh_rate;
		sum1 += z[i] * z[i];
		sum2 += cos(2.0 * M_PI * z[i]);
	}
	sum1 = -0.2 * sqrt(sum1 / G_nx[2]);
	sum2 /= G_nx[2];
	top += exp(1.0) - 20.0 * exp(sum1) - exp(sum2) + 20.0;



	sh_rate = 5.12 / 100.0;                             //rastrigin's function
	sum1 = 0.0, sum2 = 0.0;


	for (i = G[3]; i < G[4]; i++)
	{
		t = z[i] * sh_rate;
		top = top + ((t * t) - 10.0 * cos(2.0 * M_PI * t) + 10.0);
	}

	
	sh_rate = 1000.0 / 100.0;                        //modified achwefel's function
	sum_z = 0.0;
	tmp1 = 0.0;
	for (i = G[4]; i < G[5]; i++)
	{
		t = z[i] * sh_rate;
		t += 4.209687462275036e+002;
		if (t > 500)
		{
			sum_z -= (500.0 - fmod(t, 500)) * sin(pow(500.0 - fmod(t, 500), 0.5));
			tmp1 = (t - 500.0) / 100;
			sum_z += tmp1 * tmp1 / G_nx[4];
		}
		else if (t < -500)
		{
			sum_z -= (-500.0 + fmod(fabs(t), 500)) * sin(pow(500.0 - fmod(fabs(t), 500), 0.5));
			tmp1 = (t + 500.0) / 100;
			sum_z += tmp1 * tmp1 / G_nx[4];
		}
		else
			sum_z -= t * sin(pow(fabs(t), 0.5));
	}
	sum_z += 4.189828872724338e+002 * G_nx[4];
	top = top + sum_z;




	sh_rate = 1.0;                        //schaffer's F7 function
	sum_z = 0.0;

	tmp1 = 0.0;
	for (i = G[5]; i < D - 1; i++)
	{
		z[i] = pow(z[(i - G[5])] * z[(i - G[5])] + z[(i - G[5]) + 1] * z[(i - G[5]) + 1], 0.5);
		tmp1 = sin(50.0 * pow(z[i], 0.2));
		sum_z += pow(z[i], 0.5) + pow(z[i], 0.5) * tmp1 * tmp1;
	}
	sum_z = sum_z * sum_z / (G_nx[5] - 1) / (G_nx[5] - 1);
	top = top + sum_z;


	return top;
}

double f8(double sol[D])
{
	int i, tmp, num, cf_num = 5;
	char func_num[4] = "8";
	double alpha, sh_rate, top;
	double  tmp1, tmp2, tmp3;
	double temp;
	double sum1, sum2, sum_z, r2;
	int G[5], G_nx[5];
	double Gp[5] = { 0.3,0.2,0.2,0.1,0.2 };
	lb = -100;	ub = 100;

	tmp = 0;
	for (i = 0; i < cf_num - 1; i++) {
		G_nx[i] = ceil(Gp[i] * D);
		tmp += G_nx[i];
	}
	G_nx[cf_num - 1] = D - tmp;

	G[0] = 0;
	for (i = 1; i < cf_num; i++) {
		G[i] = G[i - 1] + G_nx[i - 1];
	}

	num = 3;
	if (fileflag == 0)
		read_cec2022_file(func_num, num);

	for (i = 0; i < D; i++)
		z[i] = sol[i] - Os[i];

	rotatefunc(z, y);

	for (i = 0; i < D; i++)
		z[i] = y[SS[i] - 1];

	top = 0.0;
	sh_rate = 5.0 / 100.0;						   // katsuura function
	tmp1 = 0.0, tmp2 = 0.0, temp = 0.0, sum_z = 1.0;
	tmp3 = pow(1.0 * G_nx[0], 1.2);

	for (i = G[0]; i < G[1]; i++) {

		z[i] = z[i] * sh_rate;
		temp = 0.0;
		for (int j = 1; j <= 32; j++)
		{
			tmp1 = pow(2.0, j);
			tmp2 = tmp1 * z[i];
			temp += fabs(tmp2 - floor(tmp2 + 0.5)) / tmp1;
		}
		sum_z *= pow(1.0 + (i + 1) * temp, 10.0 / tmp3);
	}
	tmp1 = 10.0 / G_nx[0] / G_nx[0];
	top += tmp1 * sum_z - tmp1;

	sh_rate = 5.0 / 100.0;						   // happyCat function
	alpha = 1.0 / 8.0;
	r2 = 0.0;
	sum_z = 0.0;
	for (i = G[1]; i < G[2]; i++) {
		z[i] = z[i] * sh_rate;
		z[i] = z[i] - 1.0;

		r2 += z[i] * z[i];
		sum_z += z[i];
	}
	top += pow(fabs(r2 - G_nx[1]), 2 * alpha) + (0.5 * r2 + sum_z) / G_nx[1] + 0.5;


	sh_rate = 5.0 / 100.0;						   // expanded girewank's plus rosenbrock's function
	tmp1 = 0.0, tmp2 = 0.0, temp = 0.0, sum_z = 0.0;
	z[G[2]] *= sh_rate;
	z[G[2]] += 1.0;
	for (i = G[2]; i < G[3] - 1; i++)
	{
		z[i + 1] *= sh_rate;
		z[i + 1] += 1.0;
		tmp1 = z[i] * z[i] - z[i + 1];
		tmp2 = z[i] - 1.0;
		temp = 100.0 * tmp1 * tmp1 + tmp2 * tmp2;
		sum_z += (temp * temp) / 4000.0 - cos(temp) + 1.0;
	}
	tmp1 = z[G[2] + G_nx[2] - 1] * z[G[2]+G_nx[2] - 1] - z[G[2]];
	tmp2 = z[G[2]+G_nx[2] - 1] - 1.0;
	temp = 100.0 * tmp1 * tmp1 + tmp2 * tmp2;
	top += sum_z + (temp * temp) / 4000.0 - cos(temp) + 1.0;


	sh_rate = 1000.0 / 100.0;						//shwefel_func
	sum_z = 0.0;
	tmp1 = 0.0;
	for (i = G[3]; i < G[4]; i++) {					
		z[i] = z[i] * sh_rate;
		z[i] += 4.209687462275036e+002;

		if (z[i] > 500) {
			sum_z -= (500.0 - fmod(z[i], 500)) * sin(pow(500.0 - fmod(z[i], 500), 0.5));
			tmp1 = (z[i] - 500.0) / 100;
			sum_z += tmp1 * tmp1 / G_nx[3];
		}
		else if (z[i] < -500) {
			sum_z -= (-500.0 + fmod(fabs(z[i]), 500)) * sin(pow(500.0 - fmod(fabs(z[i]), 500), 0.5));
			tmp1 = (z[i] + 500.0) / 100;
			sum_z += tmp1 * tmp1 / G_nx[3];
		}
		else
			sum_z -= z[i] * sin(pow(fabs(z[i]), 0.5));
	}
	sum_z += 4.189828872724338e+002 * G_nx[3];
	top = top + sum_z;


	//sh_rate = 1.0;
	sum1 = sum2 = 0.0;
	for (i = G[4]; i < D; i++) {						// ackley_func (f5)
		sum1 += z[i] * z[i];
		sum2 += cos(2.0 * M_PI * z[i]);
	}
	sum1 = -0.2 * sqrt(sum1 / G_nx[4]);
	sum2 /= G_nx[4];
	top += exp(1.0) - 20.0 * exp(sum1) - exp(sum2) + 20.0;

	return top;
}

double cf_Os[6][D], cf_Mr[6][D][D], cf_SS[6][D];	// cf need shift rotate and shuffle data

void read_cec2022_cf(char* funcnum, int num, int cf_num)		// for Composition Function
{
	char FName[45], tempc[5];
	int i, j, k;
	ifstream ffit;

	strcpy(FName, "../CEC2022_input_data/shift_data_");
	strcat(FName, funcnum);
	strcat(FName, ".txt");
	cout << FName << endl;

	ffit.open(FName);
	if (!ffit)
		cout << "Cannot open input file." << endl;

	for (k = 0; k < cf_num; k++)
	{
		for (i = 0; i < D; i++)
			ffit >> cf_Os[k][i];

		ffit.ignore(numeric_limits<streamsize>::max(), '\n');

	}


	ffit.close();

	strcpy(FName, "../CEC2022_input_data/M_");
	strcat(FName, funcnum);
	strcat(FName, "_D");
	itoa(D, tempc, 10); //cout<< tempc << endl;
	strcat(FName, tempc);
	strcat(FName, ".txt");
	cout << FName << endl;

	ffit.open(FName);
	if (!ffit)
		cout << "Cannot open input file." << endl;

	for (k = 0; k < cf_num; k++)
		for (i = 0; i < D; i++)
			for (j = 0; j < D; j++)
				ffit >> cf_Mr[k][i][j];
	ffit.close();

	if (num == 3) {
		strcpy(FName, "../CEC2022_input_data/shuffle_data_");
		strcat(FName, funcnum);
		strcat(FName, "_D");
		itoa(D, tempc, 10);
		strcat(FName, tempc);
		strcat(FName, ".txt");
		cout << FName << endl;

		ffit.open(FName);
		if (!ffit)
			cout << "Cannot open input file." << endl;

		for (k = 0; k < cf_num; k++)
			for (i = 0; i < D; i++)
				ffit >> cf_SS[k][i];
		ffit.close();
	}
	fileflag = 1;
}

void movedata_cf(int num, int cf_no, int r_flag)		// for Composition Function
{
	int i, j;

	for (i = 0; i < D; i++)
		Os[i] = cf_Os[cf_no][i];

	if (r_flag == 1)
		for (i = 0; i < D; i++)
			for (j = 0; j < D; j++)
				Mr[i][j] = cf_Mr[cf_no][i][j];
	else
		for (i = 0; i < D; i++)
			for (j = 0; j < D; j++)
				if (i == j)
					Mr[i][j] = 1;
				else
					Mr[i][j] = 0;

	if (num == 3)
		for (i = 0; i < D; i++)
			SS[i] = cf_SS[cf_no][i];
}

double cf_cal(double* x, double* delta, double* bias, double* fit, int cf_num)
{
	int i, j;
	double w[6];
	double top, w_max = 0, w_sum = 0;

	for (i = 0; i < cf_num; i++) {
		fit[i] += bias[i];
		w[i] = 0;
		for (j = 0; j < D; j++) {
			w[i] += pow(x[j] - cf_Os[i][j], 2.0);
		}
		if (w[i] != 0)
			w[i] = pow(1.0 / w[i], 0.5) * exp(-w[i] / 2.0 / D / pow(delta[i], 2.0));
		else
			w[i] = INF;
		if (w[i] > w_max)
			w_max = w[i];
	}

	for (i = 0; i < cf_num; i++) {
		w_sum = w_sum + w[i];
	}
	if (w_max == 0) {
		for (i = 0; i < cf_num; i++)
			w[i] = 1;
		w_sum = cf_num;
	}
	top = 0.0;
	for (i = 0; i < cf_num; i++) {
		top = top + w[i] / w_sum * fit[i];
	}

	return top;
}

/* 标记af的函数是组合函数中需要的，但是前面没有定义过的函数, 充当一个辅助函数, auxiliary function
组合函数还是 double f9(double sol[D]), 为了区分开，我就使用了af开头命名 */

double af1(double sol[D])						//elliptic function
{
	int i, num;
	double sh_rate, top;

	sh_rate = 1.0;
	for (i = 0; i < D; i++)
		y[i] = (sol[i] - Os[i]) * sh_rate;

	rotatefunc(y, z);

	top = 0.0;
	for (i = 0; i < D; i++)
	{
		top += pow(10.0, 6.0 * i / (D - 1)) * z[i] * z[i];;
	}

	return top;
}

double af2(double sol[D])					//bent_cigart function
{
	int i, num;
	double sh_rate, top,sum;

	sh_rate = 1.0;
	for (i = 0; i < D; i++)
		y[i] = (sol[i] - Os[i]) * sh_rate;

	rotatefunc(y, z);

	sum = 0.0;
	top = z[0] * z[0];
	for (i = 1; i < D; i++)
		top += pow(10.0, 6.0) * z[i] * z[i];

	return top;

}

double af3(double sol[D])							//discus function
{
	int i, num;
	double sh_rate, top, sum;

	sh_rate = 1.0;
	for (i = 0; i < D; i++)
		y[i] = (sol[i] - Os[i]) * sh_rate;

	rotatefunc(y, z);

	top = pow(10.0,6.0)*z[0]*z[0];
	for (i = 1; i < D; i++)
	{
		top += z[i] * z[i];
	}
	return top;
}

double af4(double sol[D],int r_flag)                 // schwefel's function
{
	int i, num;
	double sh_rate, tmp,top;

	sh_rate = 1000.0 / 100.0;
	

	if (r_flag == 1)
	{
		for (i = 0; i < D; i++)
			y[i] = (sol[i] - Os[i]) * sh_rate;

		rotatefunc(y, z);
	}
	else{
		for (i = 0; i < D; i++)
			z[i] = (sol[i] - Os[i]) * sh_rate;
	}

	top = 0.0;
	for (i = 0; i < D; i++) {

		z[i] += 4.209687462275036e+002;
		if (z[i] > 500) {
			top -= (500.0 - fmod(z[i], 500)) * sin(pow(500.0 - fmod(z[i], 500), 0.5));
			tmp = (z[i] - 500.0) / 100;
			top += tmp * tmp / D;
		}
		else if (z[i] < -500) {
			top -= (-500.0 + fmod(fabs(z[i]), 500)) * sin(pow(500.0 - fmod(fabs(z[i]), 500), 0.5));
			tmp = (z[i] + 500.0) / 100;
			top += tmp * tmp / D;
		}
		else
			top -= z[i] * sin(pow(fabs(z[i]), 0.5));
	}
	top += 4.189828872724338e+002 * D;

	return top;


}

double af5(double sol[D])							//hgbat function
{
	int i, num;
	double sh_rate, top;
	double alpha, r2, sum_z;

	sh_rate = 5.0 / 100.0;
	for (i = 0; i < D; i++)							 
		y[i] = (sol[i] - Os[i]) * sh_rate;			

	rotatefunc(y, z);					

	alpha = 1.0 / 4.0;
	r2 = 0.0; 	sum_z = 0.0;
	top = 0.0;
	for (i = 0; i < D; i++) {
		z[i] = z[i] - 1.0;				
		r2 += z[i] * z[i];
		sum_z += z[i];
	}
	top = pow(fabs(pow(r2, 2.0) - pow(sum_z, 2.0)), 2 * alpha) + (0.5 * r2 + sum_z) / D + 0.5;

	return top;
}

double af6(double sol[D])                         //escaffer6_func
{
	int i, num;
	double sh_rate, top,temp1,temp2;

	sh_rate = 1.0;                       
	for (i = 0; i < D; i++)
		y[i] = (sol[i] - Os[i]) * sh_rate;

	rotatefunc(y, z);

	temp1 = 0.0, temp2 = 0.0,top = 0.0;
	for (i = 0; i < D - 1; i++)
	{
		temp1 = sin(sqrt(z[i] * z[i] + z[i + 1] * z[i + 1]));
		temp1 = temp1 * temp1;
		temp2 = 1.0 + 0.001 * (z[i] * z[i] + z[i + 1] * z[i + 1]);
		top += 0.5 + (temp1 - 0.5) / (temp2 * temp2);
	}
	temp1 = sin(sqrt(z[D - 1] * z[D - 1] + z[0] * z[0]));
	temp1 = temp1 * temp1;
	temp2 = 1.0 + 0.001 * (z[D - 1] * z[D - 1] + z[0] * z[0]);
	top += 0.5 + (temp1 - 0.5) / (temp2 * temp2);

	return top;
}

double af7(double sol[D])                           //griewank function
{
	int i, num;
	double sh_rate, top, tmp,sum;
	
	sh_rate = 600.0 / 100.0;
	for (i = 0; i < D; i++)
	{
		y[i] = (sol[i] - Os[i]) * sh_rate;
	}

	rotatefunc(y, z);

	top = 0.0, tmp = 0.0, sum = 1.0;
	for (i = 0; i < D; i++)
	{
		tmp += z[i] * z[i];
		sum *= cos(z[i] / sqrt(1.0 + i));
	}
	top = 1.0 + tmp / 4000.0 - sum;

	return top;
}

double f9(double sol[D])                         //composition function 1
{
	int i, num, cf_num = 5;
	char func_num[4] = "9";
	double fit[5], top;
	double delta[5] = { 10, 20, 30, 40, 50 };
	double bias[5] = { 0, 200, 300, 100, 400 };

	num = 2;
	if (fileflag == 0)
		read_cec2022_cf(func_num, num, cf_num);

	i = 0;
	movedata_cf(num, i, 1);
	fit[i] = f2(sol);										//rosenbrock_func
	fit[i] = 10000 * fit[i] / 1e+4;

	i = 1;
	movedata_cf(num, i, 1);
	fit[i] = af1(sol);										//ellips_func
	fit[i] = 10000 * fit[i] / 1e+10;

	i = 2;
	movedata_cf(num, i, 1);
	fit[i] = af2(sol);										//bent_cigar_func
	fit[i] = 10000 * fit[i] / 1e+30;

	i = 3;
	movedata_cf(num, i, 1);
	fit[i] = af3(sol);										//discus_func
	fit[i] = 10000 * fit[i] / 1e+10;

	i = 4;
	movedata_cf(num, i, 0);
	fit[i] = af1(sol);										//ellips_func
	fit[i] = 10000 * fit[i] / 1e+10;

	top = cf_cal(sol, delta, bias, fit, cf_num);

	return top;

}

double f10(double sol[D])                   //composition function 2
{
	int i, num, cf_num = 3;
	char func_num[4] = "10";
	double fit[3], top;
	double delta[3] = { 20,10,10 };
	double bias[3] = { 0,200,100 };

	num = 2;
	if (fileflag == 0)
		read_cec2022_cf(func_num, num, cf_num);

	i = 0;
	movedata_cf(num, i, 0); 
	fit[i] = af4(sol,0);						//schwefel_func r_falg = 0 , 不旋转

	i = 1;
	movedata_cf(num, i, 1);
	fit[i] = f4(sol);							//rastrigin_func

	i = 2;
	movedata_cf(num, i, 1);
	fit[i] = af5(sol);							//hgbat_func

	top = cf_cal(sol, delta, bias, fit, cf_num);

	return top;

}

double f11(double sol[D])
{
	int i,num,cf_num = 5;
	char func_num[4] = "11";
	double fit[5], top;
	double delta[5] = { 20,20,30,30,20 };
	double bias[5] = { 0,200,300,400,200 };

	num = 2;
	if (fileflag == 0)
		read_cec2022_cf(func_num, num, cf_num);

	i = 0;
	movedata_cf(num, i, 1);
	fit[i] = af6(sol);											//escaffer6_func
	fit[i] = 10000 * fit[i] / 2e+7;

	i = 1;
	movedata_cf(num, i, 1);
	fit[i] = af4(sol, 1);										//schwefel_funcr

	i = 2;
	movedata_cf(num, i, 1);
	fit[i] = af7(sol);											//griewank_func
	fit[i] = 1000 * fit[i] / 100;

	i = 3;
	movedata_cf(num, i, 1);
	fit[i] = f2(sol);											//rosenbrock_func

	i = 4;
	movedata_cf(num, i, 1);
	fit[i] = f4(sol);											//rastrigin_func
	fit[i] = 10000 * fit[i] / 1e+3;

	top = cf_cal(sol, delta, bias, fit, cf_num);

	return top;
}

double f12(double sol[D])
{
	int i,num, cf_num = 6;
	char func_num[4] = "12";
	double fit[6], top;
	double delta[6] = { 10,20,30,40,50,60 };
	double bias[6] = { 0,300,500,100,400,200 };

	num = 2;
	if (fileflag == 0)
		read_cec2022_cf(func_num, num, cf_num);

	i = 0;
	movedata_cf(num, i, 1);
	fit[i] = af5(sol);												//hgbat_func
	fit[i] = 10000 * fit[i] / 1000;

	i = 1;
	movedata_cf(num, i, 1);
	fit[i] = f4(sol);												//rastrigin_func	
	fit[i] = 10000 * fit[i] / 1e+3;

	i = 2;
	movedata_cf(num, i, 1);
	fit[i] = af4(sol,1);										     //schwefel_func
	fit[i] = 10000 * fit[i] / 4e+3;

	i = 3;
	movedata_cf(num, i, 1);											
	fit[i] = af2(sol);												//bent_cigar_func
	fit[i] = 10000 * fit[i] / 1e+30;

	i = 4;
	movedata_cf(num, i, 1);
	fit[i] = af1(sol);												//ellips_func
	fit[i] = 10000 * fit[i] / 1e+10;

	i = 5;
	movedata_cf(num, i, 1);
	fit[i] = af6(sol);												//escaffer6_func
	fit[i] = 10000 * fit[i] / 2e+7;									

	top = cf_cal(sol, delta, bias, fit, cf_num);

	return top;

}






