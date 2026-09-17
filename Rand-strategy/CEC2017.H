//  CEC17 Test Function Suite 
//

#include "string.h"

double Os[D], y[D], z[D];
double Mr[D][D];
int SS[D], fileflag = 0;

void rotatefunc (double *x, double *xrot)
{
	int i,j;
    for(i = 0; i < D; i++){
        xrot[i] = 0;
		for (j = 0; j < D; j++)
			xrot[i] = xrot[i] + x[j] * Mr[i][j];
    }
}

void read_cec2017_file( char *funcnum, int num )
{ 
	char FName[70], tempc[5];
	int i, j;
	ifstream ffit;
	
	strcpy( FName, "../CEC2017_input_data/shift_data_" );
	strcat( FName, funcnum );
	strcat( FName, ".txt" );
	cout << FName << endl;

	ffit.open( FName );
	if( !ffit )
		cout << "Cannot open input file." << endl;

	for(i = 0; i < D; i++)	
		ffit >> Os[i];

	ffit.close( );

	strcpy( FName, "../CEC2017_input_data/M_" );
	strcat( FName, funcnum );
	strcat( FName, "_D" );
	itoa( D, tempc, 10 ); //cout<< tempc << endl;
	strcat( FName, tempc );
	strcat( FName, ".txt" );
	cout << FName << endl;

	ffit.open( FName );
	if( !ffit )
		cout << "Cannot open input file." << endl;

	for(i = 0; i < D; i++){
		for (j = 0; j < D; j++)	
			ffit >> Mr[i][j];
	}
	ffit.close( );

	if( num == 3 ){
		strcpy( FName, "../CEC2017_input_data/shuffle_data_" );
		strcat( FName, funcnum );
		strcat( FName, "_D" );
		itoa( D, tempc, 10 );
		strcat( FName, tempc );
		strcat( FName, ".txt" );
		cout << FName << endl;

		ffit.open( FName );
		if( !ffit )
			cout << "Cannot open input file." << endl;

		for(i = 0; i < D; i++)	
			ffit >> SS[i];

		ffit.close( );
	}
	fileflag = 1;
}

double f1(double sol[D])			// bent_cigar_func
{
    int i, num;
    double sh_rate, top = 0;
	char func_num[4] = "1";

	lb = -100;	ub = 100;
	num = 2;
	if(fileflag == 0)
		read_cec2017_file( func_num, num );

	sh_rate = 1.0;
	for(i = 0; i < D; i++)							//shift 
		y[i] = (sol[i] - Os[i]) * sh_rate;			//shrink to the original search range

	rotatefunc (y, z);				//rotate

	top = z[0]*z[0];
	for(i = 1; i < D; i++)
		top += pow(10.0, 6.0) * z[i] * z[i];

	return top;
}

double f2(double sol[D])			/* sum of different power */
{
	int i, num;
	double sh_rate, top = 0;
	char func_num[4] = "2";


	lb = -100;	ub = 100;
	num = 2;
	if (fileflag == 0)
		read_cec2017_file(func_num, num);

	sh_rate = 1.0;
	for (i = 0; i < D; i++)							//shift 
		y[i] = (sol[i] - Os[i]) * sh_rate;			//shrink to the original search range

	rotatefunc(y, z);				//rotate

	for (i = 0; i < D; i++) {
		double xi = z[i];
		double newv = pow((abs(xi)), (i + 1));
		top = top + newv;
	}

	return top;
}

double f3(double sol[D])			 /* zakharov */
{
    int i, num;
    double sh_rate, top = 0;
	char func_num[4] = "3";

	lb = -100;	ub = 100;
	num = 2;
	if(fileflag == 0)
		read_cec2017_file( func_num, num );

	sh_rate = 1.0;
	for(i = 0; i < D; i++)							//shift 
		y[i] = (sol[i] - Os[i]) * sh_rate;			//shrink to the original search range

	rotatefunc (y, z);				//rotate

	double sum1 = 0.0;
	double sum2 = 0.0;
	for (i = 0; i < D; i++){
		double xi = z[i];
		sum1 = sum1 + pow(xi, 2);
		sum2 = sum2 + 0.5 * (i+1) * xi;
	}
	
	top = sum1 + pow(sum2,2) + pow(sum2,4);

	return top;
}

double f4(double sol[D])			// rosenbrock_func
{
    int i, num;
    double sh_rate, top = 0;
	char func_num[4] = "4";

	lb = -100;	ub = 100;
	num = 2;
	if(fileflag == 0)
		read_cec2017_file( func_num, num );

	sh_rate = 2.048 / 100.0;
	for(i = 0; i < D; i++)							//shift 
		y[i] = (sol[i] - Os[i]) * sh_rate;			//shrink to the original search range

	rotatefunc (y, z);			//rotate

	z[0] = z[0] + 1.0;
	for(i = 0; i < D-1; i++){
		z[i+1] = z[i+1] + 1.0;
		top += 100 * (z[i] * z[i] - z[i+1]) * (z[i] * z[i] - z[i+1]) + (z[i] - 1.0) * (z[i] - 1.0);
	}

	return top ;
}

double f5(double sol[D])				// rastrigin_func
{
    int i, num;
    double sh_rate, top;
	char func_num[4] = "5";

	lb = -100;	ub = 100;
	num = 2;
	if(fileflag == 0)
		read_cec2017_file( func_num, num );

	sh_rate = 5.12 / 100.0;
	for(i = 0; i < D; i++)							//shift 
		y[i] = (sol[i] - Os[i]) * sh_rate;			//shrink to the original search range

	rotatefunc (y, z);			//rotate

	top = 0.0;
	for (i = 0; i < D; i++)
		top += (z[i] * z[i] - 10.0 * cos( 2.0 * M_PI * z[i]) + 10.0 );

	return top ;
}

double f6(double sol[D])				 /* schaffer_F7_func  Schwefel's 1.2  */
{
    int i, num;
    double sh_rate, top;
	char func_num[4] = "6";
	double tmp;
	lb = -100;	ub = 100;

	num = 2;
	if(fileflag == 0)
		read_cec2017_file( func_num, num );

	sh_rate = 1.0;
	for(i = 0; i < D; i++)							//shift 
		y[i] = (sol[i] - Os[i]) * sh_rate;			//shrink to the original search range

	rotatefunc (y, z);								//rotate

 	top = 0.0;
	for( i= 0; i< D-1; i++)	
	{
		z[i] = pow( y[i] * y[i] + y[i+1] * y[i+1], 0.5);
		tmp = sin( 50.0 * pow(z[i], 0.2));
		top += pow(z[i], 0.5) + pow(z[i], 0.5) * tmp * tmp ;
	}
	top = top * top / (D-1) / (D-1);
	
	return top ;
}

/*  *********   */
double f7(double sol[D])				 /* Lunacek Bi_rastrigin Function */
{
    int i, num;
    double sh_rate, top;
	char func_num[4] = "7";

	lb = -100;	ub = 100;
	num = 2;
	if(fileflag == 0)
		read_cec2017_file( func_num, num );

	top = 0.0;
	double mu0 = 2.5, d = 1.0, s, mu1, tmp, tmp1, tmp2;
	double * tmpx = new double[D];
	s = 1.0 - 1.0 / (2.0 * pow(D + 20.0, 0.5) - 8.2 );
	mu1 = -pow( (mu0 * mu0 - d) / s, 0.5 );

	sh_rate = 10.0 / 100.0;
	for (i = 0; i < D; i++)							//shift 
		y[i] = (sol[i] - Os[i]) * sh_rate;			//shrink to the original search range

	for(i = 0; i < D; i++){
		tmpx[i] = 2 * y[i];
        if (Os[i] < 0.0)
            tmpx[i] *= -1.;
    }
	for (i = 0; i < D; i++){
		z[i] = tmpx[i];
		tmpx[i] += mu0;
	}

    tmp1 =0.0;  tmp2 = 0.0;
	for (i = 0; i < D; i++){
		tmp = tmpx[i] - mu0;
		tmp1 += tmp * tmp;
		tmp = tmpx[i] - mu1;
		tmp2 += tmp * tmp;
	}
	tmp2 *= s;
	tmp2 += d * D;
	tmp = 0.0;

	rotatefunc (y, z);			//rotate
	for ( i = 0;  i < D; i++){
		tmp += cos(2.0 * M_PI * y[i]);
	}	
	if(tmp1 < tmp2)
		top = tmp1;
	else
		top = tmp2;

	top += 10.0 * (D - tmp);
	delete []tmpx;

	return top ;
}

double f8(double sol[D])				 /* step_rastrigin_func  */
{
    int i, num;
    double sh_rate, top;
	char func_num[4] = "8";
	lb = -100;	ub = 100;

	num = 2;
	if(fileflag == 0)
		read_cec2017_file( func_num, num );

	for (i = 0; i < D; i++)
	{
		if(fabs(y[i] - Os[i]) > 0.5)
			y[i] = Os[i] + floor( 2 * (y[i] - Os[i]) + 0.5) / 2;
	}


	sh_rate = 5.12 / 100.0;
	for(i = 0; i < D; i++)							//shift 
		y[i] = (sol[i] - Os[i]) * sh_rate;			//shrink to the original search range

	rotatefunc (y, z);								//rotate

 	top = 0.0;
	for (i = 0; i < D; i++){
		top += ( z[i] * z[i] - 10.0 * cos(2.0 * M_PI * z[i]) + 10.0);
	}
	
	return top;
}

double f9(double sol[D])				 /* Levy function */
{
    int i, num;
    double sh_rate, top;
	char func_num[4] = "9";

	lb = -100;	ub = 100;
	num = 2;
	if(fileflag == 0)
		read_cec2017_file( func_num, num );

	sh_rate = 1.0;
	for(i = 0; i < D; i++)							//shift 
		y[i] = (sol[i] - Os[i]) * sh_rate;			//shrink to the original search range

	rotatefunc (y, z);			//rotate

	top = 0.0;
	double * w = new double[D];

	double sum1 = 0.0;
	for (i = 0; i < D; i++){
	   w[i] = 1.0 + ( z[i] - 1.0) / 4.0;
	}
	
	double term1 = pow( (sin( M_PI * w[0]) ),2);
	double term3 = pow( (w[D-1] - 1), 2) * (1 + pow( (sin(2 * M_PI * w[D-1]) ), 2) );
	
	double sum = 0.0;

	for (i = 0; i < D - 1; i++){
		double wi = w[i];
        double newv = pow( (wi-1), 2) * (1 + 10 * pow( (sin(M_PI * wi + 1) ), 2) );
		sum = sum + newv;
	}
	
	top = term1 + sum + term3;
	delete []w;					 // ADD THIS LINE to free memory! Thanks for Dr. Janez

	return top ;
}

double f10(double sol[D])				// schwefel_func
{
    int i, num;
    double sh_rate, top, tmp;
	char func_num[4] = "10";

	lb = -100;	ub = 100;
	num = 2;
	if( fileflag == 0)
		read_cec2017_file( func_num, num );

	sh_rate = 1000.0 / 100.0;
	for(i = 0; i < D; i++)							//shift 
		y[i] = (sol[i] - Os[i]) * sh_rate;			//shrink to the original search range

	rotatefunc (y, z);			//rotate

	top = 0;
	for (i = 0; i < D; i++){
	
		z[i] += 4.209687462275036e+002;
		if(z[i] > 500){
			top -= (500.0 - fmod(z[i], 500)) * sin(pow(500.0 - fmod(z[i], 500), 0.5));
			tmp = (z[i] - 500.0) / 100;
			top += tmp * tmp / D;
		} else if (z[i] < -500){
			top -= (-500.0 + fmod(fabs(z[i]), 500)) * sin(pow(500.0 - fmod(fabs(z[i]), 500), 0.5));
			tmp = (z[i] + 500.0) / 100;
			top += tmp * tmp / D;
		}else
			top -= z[i] * sin(pow(fabs(z[i]), 0.5));
	}
	top += 4.189828872724338e+002 * D;

	return top ;
}

/*		*** Hybrid Function ***		*/ 

double f11(double sol[D])				// hf01
{
	int i, tmp, num, cf_num = 3;
	char func_num[4] = "11";
    double sh_rate, top;
	int G[3], G_nx[3];
	double Gp[3] = {0.2, 0.4, 0.4};
	lb = -100;	ub = 100;

	tmp = 0;
	for(i = 0; i < cf_num - 1 ; i++){
		G_nx[i] = ceil(Gp[i] * D);
		tmp += G_nx[i];
	}
	G_nx[cf_num-1] = D - tmp;

	G[0]=0;
	for(i = 1; i < cf_num; i++){
		G[i] = G[i-1] + G_nx[i-1];
	}

	num = 3;
	if(fileflag == 0)
		read_cec2017_file( func_num, num );

	for(i = 0; i < D; i++)				//shift 
		z[i] = sol[i] - Os[i];

	rotatefunc (z, y);					//rotate
 
	for(i = 0; i < D; i++)				//shuffle
		z[i] = y[SS[i] - 1];

 	top = 0.0;

	double sum1 = 0.0;					/* zakharov */
	double sum2 = 0.0;
	sh_rate = 1.0;
	for(i = G[0]; i < G[1]; i++){	
		double xi = z[i];
		sum1 = sum1 + pow(xi, 2);
		sum2 = sum2 + 0.5 * (i + 1) * xi;
	}
	top = sum1 + pow(sum2, 2) + pow(sum2, 4);

	double tmp1, tmp2;					/* Rosenbrock's */
	sh_rate = 2.048 / 100;;
	for(i = G[1]; i < G[2]; i++){
		z[i] = z[i] * sh_rate;
	}

	z[G[1]] += 1.0;						//shift to orgin
	for(i = G[1]; i < G[2]; i++){
		z[i+1] += 1.0;					//shift to orgin
		tmp1 = z[i] * z[i] - z[i+1];
		tmp2 = z[i] - 1.0;
		top += 100.0 * tmp1 * tmp1 + tmp2 * tmp2;
	}

	sh_rate = 5.12 / 100.0;				// rastrigin_func (f8, 9)
	for(i = G[2]; i < D; i++){
		z[i] =  z[i] * sh_rate;
		top += (z[i] * z[i] - 10.0 * cos( 2.0 * M_PI * z[i]) + 10.0 );	
	}
	
	return top ;
}

double f12(double sol[D])				// hf02
{
	int i, tmp, num, cf_num = 3;
	char func_num[4] = "12";
    double sh_rate, top;
	int G[3], G_nx[3];
	double Gp[3] = {0.3, 0.3, 0.4};
	lb = -100;	ub = 100;

	tmp = 0;
	for(i = 0; i < cf_num - 1 ; i++){
		G_nx[i] = ceil(Gp[i] * D);
		tmp += G_nx[i];
	}
	G_nx[cf_num-1] = D - tmp;

	G[0] = 0;
	for(i = 1; i < cf_num; i++){
		G[i] = G[i-1] + G_nx[i-1];
	}

	num = 3;
	if(fileflag == 0)
		read_cec2017_file( func_num, num );

	for(i = 0; i < D; i++)				//shift 
		z[i] = sol[i] - Os[i];

	rotatefunc (z, y);					//rotate
 
	for(i = 0; i < D; i++)				//shuffle
		z[i] = y[SS[i] - 1];

	top = 0;
//	sh_rate = 1.0;						/* Ellipsoidal */
	for(i = G[0]; i < G[1]; i++){
       top += pow(10.0, 6.0 * i / (G_nx[0] - 1)) * z[i] * z[i];	   
	}

	sh_rate = 1000.0 / 100.0;			// schwefel_func 
	for(i = G[1]; i < G[2]; i++){
		z[i] = z[i] * sh_rate;
		z[i] += 4.209687462275036e+002;

		if(z[i] > 500){
			top -= (500.0 - fmod(z[i], 500)) * sin(pow(500.0 - fmod(z[i], 500), 0.5));
			tmp = (z[i] - 500.0) / 100;
			top += tmp * tmp / G_nx[1];
		} else if (z[i] < -500){
			top -= (-500.0 + fmod(fabs(z[i]), 500)) * sin(pow(500.0 - fmod(fabs(z[i]), 500), 0.5));
			tmp = (z[i] + 500.0) / 100;
			top += tmp * tmp / G_nx[1];
		}else
			top -= z[i] * sin(pow(fabs(z[i]), 0.5));
	}
	top += 4.189828872724338e+002 * G_nx[1];

//	sh_rate = 1.0;						// bent_cigar_func
	top += z[G[2]] * z[G[2]];
	for(i = G[2] + 1; i < D; i++){
		top += pow(10.0, 6.0) * z[i] * z[i];
	}
	
	return top ;
}

double f13(double sol[D])				// hf03
{
	int i, tmp, num, cf_num = 3;
	char func_num[4] = "13";
    double sh_rate, top;
	int G[3], G_nx[3];
	double Gp[3] = {0.3, 0.3, 0.4};
	lb = -100;	ub = 100;

	tmp = 0;
	for(i = 0; i < cf_num - 1 ; i++){
		G_nx[i] = ceil(Gp[i] * D);
		tmp += G_nx[i];
	}
	G_nx[cf_num-1] = D - tmp;

	G[0] = 0;
	for(i = 1; i < cf_num; i++){
		G[i] = G[i-1] + G_nx[i-1];
	}

	num = 3;
	if(fileflag == 0)
		read_cec2017_file( func_num, num );

	for(i = 0; i < D; i++)				//shift 
		z[i] = sol[i] - Os[i];

	rotatefunc (z, y);					//rotate
 
	for(i = 0; i < D; i++)				//shuffle
		z[i] = y[SS[i] - 1];

//	sh_rate = 1.0;							// bent_cigar_func (f2)
	top = z[G[0]] * z[G[0]];
	for(i = G[0] + 1; i < G[1]; i++){
		top = pow(10.0, 6.0) * z[i] * z[i];
	}

	double tmp1, tmp2;						/* Rosenbrock's */
	sh_rate = 2.048 / 100;
	for(i = G[1]; i < G[2]; i++){
		z[i] = z[i] * sh_rate;
	}

	z[G[1]] += 1.0;						//shift to orgin
	for(i = G[1]; i < G[2]; i++){
		z[i+1] += 1.0;					//shift to orgin
		tmp1 = z[i] * z[i] - z[i+1];
		tmp2 = z[i] - 1.0;
		top += 100.0 * tmp1 * tmp1 + tmp2 * tmp2;
	}

											 /* Lunacek Bi_rastrigin Function */
	double mu0 = 2.5, d = 1.0, s, mu1;
	double * tmpx = new double[D];
	s = 1.0 - 1.0 / (2.0 * pow(D + 20.0, 0.5) - 8.2 );
	mu1 = -pow( (mu0 * mu0 - d) / s, 0.5 );

	sh_rate = 10.0 / 100.0;
	for(i = G[2]; i < D; i++){
		y[i] = z[i] * sh_rate;
	}

	for(i = G[2]; i < D; i++){	
		tmpx[i] = 2 * y[i];
        if (Os[i] < 0.0)
            tmpx[i] *= -1.0;
    }

	for(i = G[2]; i < D; i++){	
		z[i] = tmpx[i];
		tmpx[i] += mu0;
	}

    tmp1 =0.0;  tmp2 = 0.0;
	for(i = G[2]; i < D; i++){	
		tmp = tmpx[i] - mu0;
		tmp1 += tmp * tmp;
		tmp = tmpx[i] - mu1;
		tmp2 += tmp * tmp;
	}
	tmp2 *= s;
	tmp2 += d * G_nx[2];
	tmp = 0.0;

	for(i = G[2]; i < D; i++){	
		tmp += cos(2.0 * M_PI * z[i]);
	}	
	if(tmp1 < tmp2)
		top = tmp1;
	else
		top = tmp2;
		
	top += 10.0 * (G_nx[2] - tmp);
	delete []tmpx;
	return top ;
}

double f14(double sol[D])				// hf04
{
	int i, tmp, num, cf_num = 4;
	char func_num[4] = "14";
    double sh_rate, top;
	int G[4], G_nx[4];
	double Gp[4] = {0.2, 0.2, 0.2, 0.4};
	lb = -100;	ub = 100;

	tmp = 0;
	for(i = 0; i < cf_num - 1 ; i++){
		G_nx[i] = ceil(Gp[i] * D);
		tmp += G_nx[i];
	}
	G_nx[cf_num-1] = D - tmp;

	G[0]=0;
	for(i = 1; i < cf_num; i++){
		G[i] = G[i-1] + G_nx[i-1];
	}

	num = 3;
	if(fileflag == 0)
		read_cec2017_file( func_num, num );

	for(i = 0; i < D; i++)				//shift 
		z[i] = sol[i] - Os[i];

	rotatefunc (z, y);					//rotate
 
	for(i = 0; i < D; i++)				//shuffle
		z[i] = y[SS[i] - 1];

	top = 0;
//	sh_rate = 1.0;							/* Ellipsoidal */
	for(i = G[0]; i < G[1]; i++){
       top += pow(10.0, 6.0 * i / (G_nx[0] - 1)) * z[i] * z[i];	   
	}

    double sum1, sum2;						/* Ackley's  */
    sum1 = 0.0;    sum2 = 0.0;
	for(i = G[1]; i < G[2]; i++){
		sum1 += z[i] * z[i];
		sum2 += cos(2.0 * M_PI * z[i]);
	}
	sum1 = -0.2 * sqrt( sum1 / G_nx[1] );
	sum2 /= G_nx[1];

	top += exp(1.0) - 20.0*exp(sum1) - exp(sum2) + 20.0; // EEEEEEE


	double top1 = 0;						// schaffer_F7_func	
//	sh_rate = 1.0;
 	for(i = G[2]; i < G[3] - 1; i++){
/*
		y[i] = pow(z[i] * z[i] + z[i+1] * z[i+1], 0.5);
		tmp = sin( 50.0 * pow(y[i], 0.2) );
		top1 += pow(z[i], 0.5) + pow(y[i], 0.5) * tmp * tmp ;
*/
		z[i] = pow(y[i] * y[i] + y[i+1] * y[i+1], 0.5);
		tmp = sin( 50.0 * pow(z[i], 0.2) );
		top1 += pow(z[i], 0.5) + pow(z[i], 0.5) * tmp * tmp ;

	}
	top1 = top1 * top1 / (G_nx[2] - 1) /(G_nx[2] - 1);
	top += top1;

	sh_rate = 5.12 / 100.0;					// rastrigin_func (f5)
	for(i = G[3]; i < D; i++){
		z[i] =  z[i] * sh_rate;
		top += (z[i] * z[i] - 10.0 * cos( 2.0 * M_PI * z[i]) + 10.0 );	
	}

	return top ;
}

double f15(double sol[D])				// hf05
{
	int i, tmp, num, cf_num = 4;
	char func_num[4] = "15";
    double sh_rate, top;
	int G[4], G_nx[4];
	double Gp[4] = {0.2, 0.2, 0.3, 0.3};
	lb = -100;	ub = 100;

	tmp = 0;
	for(i = 0; i < cf_num - 1 ; i++){
		G_nx[i] = ceil(Gp[i] * D);
		tmp += G_nx[i];
	}
	G_nx[cf_num-1] = D - tmp;

	G[0]=0;
	for(i = 1; i < cf_num; i++){
		G[i] = G[i-1] + G_nx[i-1];
	}

	num = 3;
	if(fileflag == 0)
		read_cec2017_file( func_num, num );

	for(i = 0; i < D; i++)				//shift 
		z[i] = sol[i] - Os[i];

	rotatefunc (z, y);					//rotate
 
	for(i = 0; i < D; i++)				//shuffle
		z[i] = y[SS[i] - 1];

//	sh_rate = 1.0;							// bent_cigar_func (f2)
	top = z[G[0]] * z[G[0]];
	for(i = G[0] + 1; i < G[1]; i++){
		top = pow(10.0, 6.0) * z[i] * z[i];
	}

	double alpha, r2, sum_z;				// hgbat_func 
	alpha = 1.0 / 4.0;

	sh_rate = 5.0 / 100.0;
	for(i = G[1]; i < G[2]; i++){
		z[i] = z[i] * sh_rate;
	}

	r2 = 0.0;	sum_z=0.0;
	for(i = G[1]; i < G[2]; i++){	
		z[i] = z[i] - 1.0;		//shift to orgin
        r2 += z[i] * z[i];
		sum_z += z[i];
    }
    top += pow( fabs( pow(r2, 2.0) - pow(sum_z, 2.0)), 2 * alpha ) + (0.5 * r2 + sum_z )/ G_nx[1] + 0.5;

	sh_rate = 5.12 / 100.0;					// rastrigin_func (f5)
	for(i = G[2]; i < G[3]; i++){
		z[i] =  z[i] * sh_rate;
		top += (z[i] * z[i] - 10.0 * cos( 2.0 * M_PI * z[i]) + 10.0 );	
	}

	double tmp1, tmp2;						/* Rosenbrock's */
	sh_rate = 2.048 / 100;
	for(i = G[3]; i < D; i++){
		z[i] = z[i] * sh_rate;
	}

	z[G[3]] += 1.0;						//shift to orgin
	for(i = G[3]; i < D - 1; i++){
		z[i+1] += 1.0;					//shift to orgin
		tmp1 = z[i] * z[i] - z[i+1];
		tmp2 = z[i] - 1.0;
		top += 100.0 * tmp1 * tmp1 + tmp2 * tmp2;
	}
	
	return top ;
}

double f16(double sol[D])				// hf06
{
	int i, tmp, num, cf_num = 4;
	char func_num[4] = "16";
    double sh_rate, top;
	int G[4], G_nx[4];
	double Gp[4] = {0.2, 0.2, 0.3, 0.3};
	lb = -100;	ub = 100;

	tmp = 0;
	for(i = 0; i < cf_num - 1 ; i++){
		G_nx[i] = ceil(Gp[i] * D);
		tmp += G_nx[i];
	}
	G_nx[cf_num-1] = D - tmp;

	G[0] = 0;
	for(i = 1; i < cf_num; i++){
		G[i] = G[i-1] + G_nx[i-1];
	}

	num = 3;
	if(fileflag == 0)
		read_cec2017_file( func_num, num );

	for(i = 0; i < D; i++)				//shift 
		z[i] = sol[i] - Os[i];

	rotatefunc (z, y);					//rotate
 
	for(i = 0; i < D; i++)				//shuffle
		z[i] = y[SS[i] - 1];

	top = 0;
	double temp1, temp2;					// escaffer6_func

//	sh_rate = 1.0;
	for(i = G[0]; i < G[1] - 1; i++){
        temp1 = sin(sqrt(z[i] * z[i] + z[i+1] * z[i+1]));
		temp1 = temp1 * temp1;
        temp2 = 1.0 + 0.001 * (z[i] * z[i] + z[i+1] * z[i+1]);
        top += 0.5 + (temp1 - 0.5) / (temp2 * temp2);
    }
    temp1 = sin(sqrt(z[G[1]-1] * z[G[1]-1] + z[G[0]] * z[G[0]]));
	temp1 = temp1 * temp1;
    temp2 = 1.0 + 0.001 * (z[G[1]-1] * z[G[1]-1] + z[G[0]] * z[G[0]]);
    top += 0.5 + (temp1 - 0.5) / (temp2 * temp2);

	double alpha, r2, sum_z;				// hgbat_func 
	alpha = 1.0 / 4.0;

	sh_rate = 5.0 / 100.0;
	for(i = G[1]; i < G[2]; i++){
		z[i] = z[i] * sh_rate;
	}

	r2 = 0.0;	sum_z=0.0;
	for(i = G[1]; i < G[2]; i++){	
		z[i] = z[i] - 1.0;		//shift to orgin
        r2 += z[i] * z[i];
		sum_z += z[i];
    }
    top += pow( fabs( pow(r2, 2.0) - pow(sum_z, 2.0)), 2 * alpha ) + (0.5 * r2 + sum_z )/ G_nx[1] + 0.5;
	
	double tmp1, tmp2;						/* Rosenbrock's */
	sh_rate = 2.048 / 100;
	for(i = G[2]; i < G[3]; i++){
		z[i] = z[i] * sh_rate;
	}

	z[G[2]] += 1.0;						//shift to orgin
	for(i = G[2]; i < G[3]; i++){
		z[i+1] += 1.0;					//shift to orgin
		tmp1 = z[i] * z[i] - z[i+1];
		tmp2 = z[i] - 1.0;
		top += 100.0 * tmp1 * tmp1 + tmp2 * tmp2;
	}
	
	sh_rate = 1000.0 / 100.0;				// schwefel_func 
	for(i = G[3]; i < D; i++){
		z[i] = z[i] * sh_rate;
		z[i] += 4.209687462275036e+002;

		if(z[i] > 500){
			top -= (500.0 - fmod(z[i], 500)) * sin(pow(500.0 - fmod(z[i], 500), 0.5));
			tmp = (z[i] - 500.0) / 100;
			top += tmp * tmp / G_nx[3];
		} else if (z[i] < -500){
			top -= (-500.0 + fmod(fabs(z[i]), 500)) * sin(pow(500.0 - fmod(fabs(z[i]), 500), 0.5));
			tmp = (z[i] + 500.0) / 100;
			top += tmp * tmp / G_nx[3];
		}else
			top -= z[i] * sin(pow(fabs(z[i]), 0.5));
	}
	top += 4.189828872724338e+002 * G_nx[3];
	
	return top ;
}

double f17(double sol[D])				// hf07
{
	int i, j, tmp, num, cf_num = 5;
	char func_num[4] = "17";
    double sh_rate, top;
	int G[5], G_nx[5];
	double Gp[5] = {0.1, 0.2, 0.2, 0.2, 0.3};
	lb = -100;	ub = 100;

	tmp = 0;
	for(i = 0; i < cf_num - 1 ; i++){
		G_nx[i] = ceil(Gp[i] * D);
		tmp += G_nx[i];
	}
	G_nx[cf_num-1] = D - tmp;

	G[0]=0;
	for(i = 1; i < cf_num; i++){
		G[i] = G[i-1] + G_nx[i-1];
	}

	num = 3;
	if(fileflag == 0)
		read_cec2017_file( func_num, num );

	for(i = 0; i < D; i++)				//shift 
		z[i] = sol[i] - Os[i];

	rotatefunc (z, y);					//rotate
 
	for(i = 0; i < D; i++)				//shuffle
		z[i] = y[SS[i] - 1];

	double temp, tmp1, tmp2, tmp3;			/* Katsuura  */
	top = 1.0;
	tmp3 = pow(1.0 * G_nx[0], 1.2);

	sh_rate = 5.0/100.0;
	for(i = G[0]; i < G[1]; i++){
		temp = 0.0;	tmp1 = 1;
		for(j = 1; j <= 32; j++){
			tmp1 = tmp1 * 2.0;
			tmp2 = tmp1 * z[i];
			temp += fabs( tmp2 - floor(tmp2 + 0.5) ) / tmp1;
		}
		top *= pow(1.0 + (i + 1) * temp,10.0 / tmp3);
    }
	tmp1 = 10.0 / G_nx[0]/ G_nx[0];
    top = top * tmp1 - tmp1;

    double sum1, sum2;						/* Ackley's  */
    sum1 = 0.0;    sum2 = 0.0;
	for(i = G[1]; i < G[2]; i++){
		sum1 += z[i] * z[i];
		sum2 += cos(2.0 * M_PI * z[i]);
	}
	sum1 = -0.2 * sqrt( sum1 / G_nx[1] );
	sum2 /= G_nx[1];

	top += exp(1.0) - 20.0*exp(sum1) - exp(sum2) + 20.0; // EEEEEEE

	sh_rate = 5.0/100.0;						// grie_rosen_func
	for(i = G[2]; i < G[3]; i++){
		z[i] = z[i] * sh_rate;
	}

	z[G[2]] = z[G[2]] + 1.0;		//shift to orgin
	for(i = G[2]; i < G[3] - 1; i++){
		z[i+1] += 1.0;				//shift to orgin
		tmp1 = z[i] * z[i] - z[i+1];
		tmp2 = z[i] - 1.0;
        temp = 100.0 * tmp1 * tmp1 + tmp2 * tmp2;
        top += (temp * temp) / 4000.0 - cos(temp) + 1.0; 
    }
	tmp1 = z[G_nx[2]-1] * z[G_nx[2]-1] - z[G[2]];
	tmp2 = z[G_nx[2]-1] - 1.0;
    temp = 100.0 * tmp1 * tmp1 + tmp2 * tmp2;
    top += (temp * temp) / 4000.0 - cos(temp) + 1.0 ;

	sh_rate = 1000.0 / 100.0;				// schwefel_func 
	for(i = G[3]; i < G[4]; i++){
		z[i] = z[i] * sh_rate;
		z[i] += 4.209687462275036e+002;

		if(z[i] > 500){
			top -= (500.0 - fmod(z[i], 500)) * sin(pow(500.0 - fmod(z[i], 500), 0.5));
			tmp = (z[i] - 500.0) / 100;
			top += tmp * tmp / G_nx[3];
		} else if (z[i] < -500){
			top -= (-500.0 + fmod(fabs(z[i]), 500)) * sin(pow(500.0 - fmod(fabs(z[i]), 500), 0.5));
			tmp = (z[i] + 500.0) / 100;
			top += tmp * tmp / G_nx[3];
		}else
			top -= z[i] * sin(pow(fabs(z[i]), 0.5));
	}
	top += 4.189828872724338e+002 * G_nx[3];

	sh_rate = 5.12 / 100.0;					// rastrigin_func (f5)
	for(i = G[4]; i < D; i++){
		z[i] =  z[i] * sh_rate;
		top += (z[i] * z[i] - 10.0 * cos( 2.0 * M_PI * z[i]) + 10.0 );	
	}

	return top ;
}

double f18(double sol[D])				// hf08
{
	int i, tmp, num, cf_num = 5;
	char func_num[4] = "18";
    double sh_rate, top;
	int G[5], G_nx[5];
	double Gp[5] = {0.2, 0.2, 0.2, 0.2, 0.2};
	lb = -100;	ub = 100;

	tmp = 0;
	for(i = 0; i < cf_num - 1 ; i++){
		G_nx[i] = ceil(Gp[i] * D);
		tmp += G_nx[i];
	}
	G_nx[cf_num-1] = D - tmp;

	G[0] = 0;
	for(i = 1; i < cf_num; i++){
		G[i] = G[i-1] + G_nx[i-1];
	}

	num = 3;
	if(fileflag == 0)
		read_cec2017_file( func_num, num );

	for(i = 0; i < D; i++)				//shift 
		z[i] = sol[i] - Os[i];

	rotatefunc (z, y);					//rotate
 
	for(i = 0; i < D; i++)				//shuffle
		z[i] = y[SS[i] - 1];

	top = 0;
	
	//	sh_rate = 1.0;						/* Ellipsoidal */
	for(i = G[0]; i < G[1]; i++){
       top += pow(10.0, 6.0 * i / (G_nx[0] - 1)) * z[i] * z[i];	   
	}

    double sum1, sum2;						/* Ackley's  */
    sum1 = 0.0;    sum2 = 0.0;
	for(i = G[1]; i < G[2]; i++){
		sum1 += z[i] * z[i];
		sum2 += cos(2.0 * M_PI * z[i]);
	}
	sum1 = -0.2 * sqrt( sum1 / G_nx[1] );
	sum2 /= G_nx[1];

	top += exp(1.0) - 20.0 * exp(sum1) - exp(sum2) + 20.0; // EEEEEEE

	sh_rate = 5.0/100.0;						// grie_rosen_func
	for(i = G[2]; i < G[3]; i++){
		z[i] = z[i] * sh_rate;
	}

	sh_rate = 5.12 / 100.0;					// rastrigin_func (f5)
	for(i = G[2]; i < G[3]; i++){
		z[i] =  z[i] * sh_rate;
		top += (z[i] * z[i] - 10.0 * cos( 2.0 * M_PI * z[i]) + 10.0 );	
	}

	double alpha, r2, sum_z;				// hgbat_func 
	alpha = 1.0 / 4.0;

	sh_rate = 5.0 / 100.0;
	for(i = G[3]; i < G[4]; i++){
		z[i] = z[i] * sh_rate;
	}

	r2 = 0.0;	sum_z=0.0;
	for(i = G[3]; i < G[4]; i++){	
		z[i] = z[i] - 1.0;		//shift to orgin
        r2 += z[i] * z[i];
		sum_z += z[i];
    }
    top += pow( fabs( pow(r2, 2.0) - pow(sum_z, 2.0)), 2 * alpha ) + (0.5 * r2 + sum_z )/ G_nx[1] + 0.5;

//	sh_rate = 1.0;							/* discus_func  */
	top += pow(10.0, 6.0) * z[G[4]] * z[G[4]];
	for(i = G[4]; i < D; i++){
		top += z[i] * z[i];
	}	

	return top ;
}

double f19(double sol[D])				// hf09
{
	int i, tmp, num, cf_num = 5;
	char func_num[4] = "19";
    double sh_rate, top, ak, bk;
	int G[5], G_nx[5];
	double Gp[5] = {0.2, 0.2, 0.2, 0.2, 0.2};
	lb = -100;	ub = 100;

	tmp = 0;
	for(i = 0; i < cf_num - 1 ; i++){
		G_nx[i] = ceil(Gp[i] * D);
		tmp += G_nx[i];
	}
	G_nx[cf_num-1] = D - tmp;

	G[0] = 0;
	for(i = 1; i < cf_num; i++){
		G[i] = G[i-1] + G_nx[i-1];
	}

	num = 3;
	if(fileflag == 0)
		read_cec2017_file( func_num, num );

	for(i = 0; i < D; i++)				//shift 
		z[i] = sol[i] - Os[i];

	rotatefunc (z, y);					//rotate
 
	for(i = 0; i < D; i++)				//shuffle
		z[i] = y[SS[i] - 1];

//	sh_rate = 1.0;							// bent_cigar_func (f2)
	top = z[G[0]] * z[G[0]];
	for(i = G[0] + 1; i < G[1]; i++){
		top = pow(10.0, 6.0) * z[i] * z[i];
	}

	sh_rate = 5.12 / 100.0;					// rastrigin_func (f5)
	for(i = G[1]; i < G[2]; i++){
		z[i] =  z[i] * sh_rate;
		top += (z[i] * z[i] - 10.0 * cos( 2.0 * M_PI * z[i]) + 10.0 );	
	}
	
	double temp, tmp1, tmp2;
	sh_rate = 5.0/100.0;					// grie_rosen_func
	for(i = G[2]; i < G[3]; i++){
		z[i] = z[i] * sh_rate;
	}

	z[G[2]] = z[G[2]] + 1.0;		//shift to orgin
	for(i = G[2]; i < G[3] - 1; i++){
		z[i+1] += 1.0;				//shift to orgin
		tmp1 = z[i] * z[i] - z[i+1];
		tmp2 = z[i] - 1.0;
        temp = 100.0 * tmp1 * tmp1 + tmp2 * tmp2;
        top += (temp * temp) / 4000.0 - cos(temp) + 1.0; 
    }
	tmp1 = z[G_nx[2]-1] * z[G_nx[2]-1] - z[G[2]];
	tmp2 = z[G_nx[2]-1] - 1.0;
    temp = 100.0 * tmp1 * tmp1 + tmp2 * tmp2;
    top += (temp * temp) / 4000.0 - cos(temp) + 1.0 ;
	
    int j, k_max;							// weierstrass_func 
    double sum, sum2, a, b;
    
	a = 0.5;	b = 3.0;
    k_max = 20;
 
	sh_rate = 0.5 / 100.0;
	for(i = G[3]; i < G[4]; i++){
		z[i] = z[i] * sh_rate;
	}

	for(i = G[3]; i < G[4]; i++){
		sum = 0.0;		sum2 = 0.0;
		ak = 1;	bk = 1;
		for (j = 0; j <= k_max; j++) {
			sum += ak * cos( 2.0 * M_PI * bk * (z[i] + 0.5) );
			ak = ak * a;	bk = bk * b;
		}
		top += sum;
	}
	sum2 = 0.0;
	ak = 1;	bk = 1;
	for (j = 0; j <= k_max; j++) {
		sum2 += ak * cos(2.0 * M_PI * bk * 0.5);
		ak = ak * a;	bk = bk * b;
	}
	top -= G_nx[3] * sum2;

	double temp1, temp2;					// escaffer6_func

//	sh_rate = 1.0;
	for(i = G[4]; i < D - 1; i++){
        temp1 = sin(sqrt(z[i] * z[i] + z[i+1] * z[i+1]));
		temp1 = temp1 * temp1;
        temp2 = 1.0 + 0.001 * (z[i] * z[i] + z[i+1] * z[i+1]);
        top += 0.5 + (temp1 - 0.5) / (temp2 * temp2);
    }
    temp1 = sin(sqrt(z[D-1] * z[D-1] + z[G[4]] * z[G[4]]));
	temp1 = temp1 * temp1;
    temp2 = 1.0 + 0.001 * (z[D-1] * z[D-1] + z[G[4]] * z[G[4]]);
    top += 0.5 + (temp1 - 0.5) / (temp2 * temp2);

	return top ;
}

double f20(double sol[D])				// hf10
{
	int i, j, tmp, num, cf_num = 6;
	char func_num[4] = "20";
    double sh_rate, top;
	int G[6], G_nx[6];
	double Gp[6] = {0.1, 0.1, 0.2, 0.2, 0.2, 0.2};
	lb = -100;	ub = 100;

	tmp = 0;
	for(i = 0; i < cf_num - 1 ; i++){
		G_nx[i] = ceil(Gp[i] * D);
		tmp += G_nx[i];
	}
	G_nx[cf_num-1] = D - tmp;

	G[0] = 0;
	for(i = 1; i < cf_num; i++){
		G[i] = G[i-1] + G_nx[i-1];
	}

	num = 3;
	if(fileflag == 0)
		read_cec2017_file( func_num, num );

	for(i = 0; i < D; i++)				//shift 
		z[i] = sol[i] - Os[i];

	rotatefunc (z, y);					//rotate
 
	for(i = 0; i < D; i++)				//shuffle
		z[i] = y[SS[i] - 1];


	double alpha, r2, sum_z;
	alpha = 1.0 / 4.0;	r2 = 0.0; 	sum_z = 0.0;
	sh_rate = 5.0 / 100.0;
	for(i = G[0]; i < G[1]; i++){				// hgbat_func (f14)
		z[i] = z[i] * sh_rate;
		z[i] = z[i] - 1.0;			//shift to orgin
        r2 += z[i] * z[i];
		sum_z += z[i];
    }
    top = pow(fabs(pow(r2, 2.0) - pow(sum_z, 2.0)), 2 * alpha) + (0.5 * r2 + sum_z) / G_nx[0] + 0.5;
	
	double temp, tmp1, tmp2, tmp3, top1 = 1.0;			/* Katsuura  */

	tmp3 = pow(1.0 * G_nx[1], 1.2);

	sh_rate = 5.0/100.0;					
	for(i = G[1]; i < G[2]; i++){
		z[i] = z[i] * sh_rate;
		temp = 0.0;	tmp1 = 1;
		for(j = 1; j <= 32; j++){
			tmp1 = tmp1 * 2.0;
			tmp2 = tmp1 * z[i];
			temp += fabs( tmp2 - floor(tmp2 + 0.5) ) / tmp1;
		}
		top1 *= pow(1.0 + (i + 1) * temp,10.0 / tmp3);
    }
	tmp1 = 10.0 / G_nx[1]/ G_nx[1];
    top1 = top1 * tmp1 - tmp1;
	top += top1;

    double sum1, sum2;						/* Ackley's  */
    sum1 = 0.0;    sum2 = 0.0;
	for(i = G[2]; i < G[3]; i++){
		sum1 += z[i] * z[i];
		sum2 += cos(2.0 * M_PI * z[i]);
	}
	sum1 = -0.2 * sqrt( sum1 / G_nx[2] );
	sum2 /= G_nx[2];

	top +=  exp(1.0) - 20.0 * exp(sum1) - exp(sum2) + 20.0; 
	
	sh_rate = 5.12 / 100.0;					// rastrigin_func (f5)
	for(i = G[3]; i < G[4]; i++){
		z[i] =  z[i] * sh_rate;
		top += (z[i] * z[i] - 10.0 * cos( 2.0 * M_PI * z[i]) + 10.0 );	
	}
	
	sh_rate = 1000.0 / 100.0;
	for(i = G[4]; i < G[5]; i++){				// schwefel_func (f10,11)
		z[i] = z[i] * sh_rate;
		z[i] += 4.209687462275036e+002;

		if(z[i] > 500){
			top -= (500.0 - fmod(z[i], 500)) * sin(pow(500.0 - fmod(z[i], 500), 0.5));
			tmp = (z[i] - 500.0) / 100;
			top += tmp * tmp / G_nx[4];
		} else if (z[i] < -500){
			top -= (-500.0 + fmod(fabs(z[i]), 500)) * sin(pow(500.0 - fmod(fabs(z[i]), 500), 0.5));
			tmp = (z[i] + 500.0) / 100;
			top += tmp * tmp / G_nx[4];
		}else
			top -= z[i] * sin(pow(fabs(z[i]), 0.5));
	}
	top += 4.189828872724338e+002 * G_nx[4];

	top1 = 0;						// schaffer_F7_func	
//	sh_rate = 1.0;
 	for(i = G[5]; i < D - 1; i++){
/*
		y[i] = pow(z[i] * z[i] + z[i+1] * z[i+1], 0.5);
		tmp = sin( 50.0 * pow(y[i], 0.2) );
		top1 += pow(z[i], 0.5) + pow(y[i], 0.5) * tmp * tmp ;
*/
		z[i] = pow(y[i] * y[i] + y[i+1] * y[i+1], 0.5);
		tmp = sin( 50.0 * pow(z[i], 0.2) );
		top1 += pow(z[i], 0.5) + pow(z[i], 0.5) * tmp * tmp ;

	}
	top1 = top1 * top1 / (G_nx[5] - 1) /(G_nx[5] - 1);
	top += top1;
	return top ;
}

double cf_Os[6][D], cf_Mr[6][D][D], cf_SS[6][D];	// cf need shift rotate and shuffle data

void read_cec2017_cf( char *funcnum, int num, int cf_num )		// for Composition Function
{ 
	char FName[70], tempc[5];
	int i, j, k;
	ifstream ffit;
	
	strcpy( FName, "../CEC2017_input_data/shift_data_" );
	strcat( FName, funcnum );
	strcat( FName, ".txt" );
	cout << FName << endl;

	ffit.open( FName );
	if( !ffit )
		cout << "Cannot open input file." << endl;

	for(k = 0; k < cf_num; k++)
		for(i = 0; i < D; i++)	
			ffit >> cf_Os[k][i];

	ffit.close( );

	strcpy( FName, "../CEC2017_input_data/M_" );
	strcat( FName, funcnum );
	strcat( FName, "_D" );
	itoa( D, tempc, 10 ); //cout<< tempc << endl;
	strcat( FName, tempc );
	strcat( FName, ".txt" );
	cout << FName << endl;

	ffit.open( FName );
	if( !ffit )
		cout << "Cannot open input file." << endl;

	for(k = 0; k < cf_num; k++)
		for(i = 0; i < D; i++)
			for (j = 0; j < D; j++)	
				ffit >> cf_Mr[k][i][j];
	ffit.close( );

	if( num == 3 ){
		strcpy( FName, "../CEC2017_input_data/shuffle_data_" );
		strcat( FName, funcnum );
		strcat( FName, "_D" );
		itoa( D, tempc, 10 );
		strcat( FName, tempc );
		strcat( FName, ".txt" );
		cout << FName << endl;

		ffit.open( FName );
		if( !ffit )
			cout << "Cannot open input file." << endl;

		for(k = 0; k < cf_num; k++)
			for(i = 0; i < D; i++)	
				ffit >> cf_SS[k][i];
		ffit.close( );
	}
	fileflag = 1;
}

void movedata_cf( int num, int cf_no, int r_flag )		// for Composition Function
{ 
	int i, j;

	for(i = 0; i < D; i++)	
		Os[i] = cf_Os[cf_no][i];

	if( r_flag == 1 )
		for(i = 0; i < D; i++)
			for (j = 0; j < D; j++)	
				Mr[i][j] = cf_Mr[cf_no][i][j];
	else
		for(i = 0; i < D; i++)
			for (j = 0; j < D; j++)	
				if( i == j )
					Mr[i][j] = 1;
				else
					Mr[i][j] = 0;

	if( num == 3 )
		for(i = 0; i < D; i++)	
			SS[i] = cf_SS[cf_no][i];
}

const double INF = 1.0e99;

double cf_cal(double *x, double * delta, double * bias, double * fit, int cf_num)
{
	int i, j;
	double w[6];
	double top, w_max = 0, w_sum = 0;

	for(i = 0; i < cf_num; i++){
		fit[i] += bias[i];
		w[i] = 0;
		for(j = 0; j < D; j++){
			w[i] += pow(x[j] -cf_Os[i][j],2.0);
		}
		if (w[i] != 0)
			w[i] = pow(1.0 / w[i], 0.5) * exp(-w[i] /2.0 / D / pow(delta[i], 2.0));
		else
			w[i] = INF;
		if (w[i] > w_max)
			w_max = w[i];
	}

	for (i = 0; i < cf_num; i++){
		w_sum = w_sum + w[i];
	}
	if(w_max == 0){
		for (i = 0; i < cf_num; i++)
			w[i] = 1;
		w_sum=cf_num;
	}
	top = 0.0;
    for(i = 0; i < cf_num; i++){
		top = top + w[i] / w_sum * fit[i];
    }

	return top;
}

//  ****** Basic Function ******

double bf00(double sol[D])				 /*  escaffer6_func  */
{
    int i;
    double sh_rate, top;
    double temp1, temp2;

	lb = -100;	ub = 100;

	sh_rate = 1.0;
	for(i = 0; i < D; i++)							//shift 
		y[i] = (sol[i] - Os[i]) * sh_rate;			//shrink to the original search range

	rotatefunc (y, z);			//rotate

	top = 0.0;
    for( i= 0; i < D - 1; i++){
        temp1 = sin( sqrt(z[i] * z[i] + z[i+1] * z[i+1]) );
		temp1 = temp1 * temp1;
        temp2 = 1.0 + 0.001 * (z[i] * z[i] + z[i+1] * z[i+1]);
        top += 0.5 + (temp1 - 0.5) / (temp2 * temp2);
    }
    temp1 = sin( sqrt( z[D-1] * z[D-1] + z[0] * z[0] ) );
	temp1 = temp1 * temp1;
    temp2 = 1.0 + 0.001 * (z[D-1] * z[D-1] + z[0] * z[0]);
    top += 0.5 + (temp1 - 0.5) / (temp2 * temp2);

	return top ;
}

double bf11(double sol[D])				 /*  ellips_func  */
{
    int i;
    double sh_rate, top;

	lb = -100;	ub = 100;

	sh_rate = 1.0;
	for(i = 0; i < D; i++)							//shift 
		y[i] = (sol[i] - Os[i]) * sh_rate;			//shrink to the original search range

	rotatefunc (y, z);			//rotate

	top = 0.0;
	for(i = 0; i < D; i++){
       top += pow(10.0, 6.0 * i / (D-1)) * z[i] * z[i];
	}

	return top ;
}

double bf12(double sol[D])				 /*  discus_func  */
{
    int i;
    double sh_rate, top;

	lb = -100;	ub = 100;

	sh_rate = 1.0;
	for(i = 0; i < D; i++)							//shift 
		y[i] = (sol[i] - Os[i]) * sh_rate;			//shrink to the original search range

	rotatefunc (y, z);			//rotate

	top = pow(10.0, 6.0) * z[0] * z[0];
	for(i = 1; i < D; i++){
		top += z[i] * z[i];		
	}

	return top ;
}

double bf13(double sol[D])				 /*  ackley_func  */
{
    int i;
    double sh_rate, top;
    double sum1, sum2;

	lb = -100;	ub = 100;
    sum1 = 0.0;    sum2 = 0.0;

	sh_rate = 1.0;
	for(i = 0; i < D; i++)							//shift 
		y[i] = (sol[i] - Os[i]) * sh_rate;			//shrink to the original search range

	rotatefunc (y, z);			//rotate

	top = 0.0;
	for(i = 0; i < D; i++){
		sum1 += z[i] * z[i];
		sum2 += cos(2.0 * M_PI * z[i]);
	}
	sum1 = -0.2 * sqrt(sum1/D);
	sum2 /= D;
	top =  exp(1.0) - 20.0 * exp(sum1) - exp(sum2) + 20.0;
	
	return top ;
}

double bf15(double sol[D])				 /*  griewank_func  */
{
    int i;
    double sh_rate, top;
    double s, p;

	lb = -100;	ub = 100;
    s = 0.0;    p = 1.0;

	sh_rate = 600.0 / 100.0;
	for(i = 0; i < D; i++)							//shift 
		y[i] = (sol[i] - Os[i]) * sh_rate;			//shrink to the original search range

	rotatefunc (y, z);			//rotate

	for(i = 0; i < D; i++){
		s += z[i] * z[i];
		p *= cos(z[i] / sqrt(1.0 + i));
	}
	top = 1.0 + s/4000.0 - p;		

	return top ;
}

double bf17(double sol[D])				 /*  happycat_func  */
{
    int i;
    double sh_rate, top;
	double alpha, r2, sum_z;

	lb = -100;	ub = 100;
	alpha = 1.0 / 8.0;

	sh_rate = 5.0 / 100.0;
	for(i = 0; i < D; i++)							//shift 
		y[i] = (sol[i] - Os[i]) * sh_rate;			//shrink to the original search range

	rotatefunc (y, z);			//rotate

	r2 = 0.0;	sum_z = 0.0;
	for(i = 0; i < D; i++){
		z[i] = z[i] - 1.0;		//shift to orgin
        r2 += z[i] * z[i];
		sum_z += z[i];
    }
    top = pow( fabs(r2 - D), 2 * alpha ) + (0.5 * r2 + sum_z ) / D + 0.5;

	return top ;
}

double bf18(double sol[D])				 /*  hgbat_func   */
{
     int i;
    double sh_rate, top;
	double alpha, r2, sum_z;

	lb = -100;	ub = 100;
	alpha = 1.0 / 4.0;

	sh_rate = 5.0 / 100.0;
	for(i = 0; i < D; i++)							//shift 
		y[i] = (sol[i] - Os[i]) * sh_rate;			//shrink to the original search range

	rotatefunc (y, z);			//rotate

	r2 = 0.0;	sum_z = 0.0;
	for(i = 0; i < D; i++){
		z[i] = z[i] - 1.0;		//shift to orgin
        r2 += z[i] * z[i];
		sum_z += z[i];
    }
	top = pow( fabs( pow(r2, 2.0) - pow(sum_z, 2.0) ), 2 * alpha ) + (0.5* r2 + sum_z) / D + 0.5;

	return top ;
}

//  ****** Basic Function END  ******


double f21(double sol[D])					// Composition Function 1
{
	int i, num, cf_num = 3;
	char func_num[4] = "21";
	double fit[3], top;
	double delta[3] = {10, 20, 30};
	double bias[3] = {0, 100, 200};

	num = 2;
	if(fileflag == 0)
		read_cec2017_cf( func_num, num, cf_num );

	i=0;
	movedata_cf( num, i, 0 );
	fit[i] = f4( sol );

	i=1;
	movedata_cf( num, i, 1 );
	fit[i] = bf11( sol );
	fit[i] = 10000 * fit[i] / 1e+10;

	i=2;
	movedata_cf( num, i, 1 );
	fit[i] = f5( sol );

	top = cf_cal( sol, delta, bias, fit, cf_num);

	return top;
}

double f22(double sol[D])					// Composition Function 2
{
	int i, num, cf_num = 3;
	char func_num[4] = "22";
	double fit[3], top;
	double delta[3] = {10, 20, 30};
	double bias[3] = {0, 100, 200};

	num = 2;
	if(fileflag == 0)
		read_cec2017_cf( func_num, num, cf_num );

	i=0;
	movedata_cf( num, i, 1 );
	fit[i] = f5( sol );

	i=1;
	movedata_cf( num, i, 1 );
	fit[i] = bf15( sol );
	fit[i] = 1000 * fit[i] / 100;

	i=2;
	movedata_cf( num, i, 1 );
	fit[i] = f10( sol );

	top = cf_cal( sol, delta, bias, fit, cf_num);

	return top;
}

double f23(double sol[D])					// Composition Function 3
{
	int i, num, cf_num = 4;
	char func_num[4] = "23";
	double fit[4], top;
	double delta[4] = {10, 20, 30, 40};
	double bias[4] = {0, 100, 200, 300};
	
	num = 2;
	if(fileflag == 0)
		read_cec2017_cf( func_num, num, cf_num );

	i=0;
	movedata_cf( num, i, 1 );
	fit[i] = f4( sol );
	
	i=1;
	movedata_cf( num, i, 1 );
	fit[i] = bf13( sol );
	fit[i] = 1000 * fit[i] / 100;

	i=2;
	movedata_cf( num, i, 1 );
	fit[i] = f10( sol );

	i=3;
	movedata_cf( num, i, 1 );
	fit[i] = f5( sol );
	
	top = cf_cal( sol, delta, bias, fit, cf_num);

	return top;
}

double f24(double sol[D])					// Composition Function 4
{
	int i, num, cf_num = 4;
	char func_num[4] = "24";
	double fit[4], top;
	double delta[4] = {10, 20, 30, 40};
	double bias[4] = {0, 100, 200, 300};
	
	num = 2;
	if(fileflag == 0)
		read_cec2017_cf( func_num, num, cf_num );

	i=0;
	movedata_cf( num, i, 1 );
	fit[i] = bf13( sol );
	fit[i] = 1000 * fit[i] / 100;
	
	i=1;
	movedata_cf( num, i, 1 );
	fit[i] = bf11( sol );
	fit[i] = 10000 * fit[i] / 1e+10;

	i=2;
	movedata_cf( num, i, 1 );
	fit[i] = bf15( sol );
	fit[i] = 1000 * fit[i] / 100;	

	i=3;
	movedata_cf( num, i, 1 );
	fit[i] = f5( sol );
	
	top = cf_cal( sol, delta, bias, fit, cf_num);

	return top;
}

double f25(double sol[D])					// Composition Function 5
{
	int i, num, cf_num = 5;
	char func_num[4] = "25";
	double fit[5], top;
	double delta[5] = {10, 20, 30, 40, 50};
	double bias[5] = {0, 100, 200, 300, 400};
	
	num = 2;
	if(fileflag == 0)
		read_cec2017_cf( func_num, num, cf_num );

	i=0;
	movedata_cf( num, i, 1 );
	fit[i] = f5( sol );
	fit[i] = 10000 * fit[i] / 1e+3;
	
	i=1;
	movedata_cf( num, i, 1 );
	fit[i] = bf17( sol );
	fit[i] = 1000 * fit[i] / 1e+3;

	i=2;
	movedata_cf( num, i, 1 );
	fit[i] = bf13( sol );
	fit[i] = 1000 * fit[i] / 100;	

	i=3;
	movedata_cf( num, i, 1 );
	fit[i] = bf12( sol );
	fit[i] = 10000 * fit[i] / 1e+10;
	
	i=4;
	movedata_cf( num, i, 0 );
	fit[i] = f4( sol );
	
	top = cf_cal( sol, delta, bias, fit, cf_num);

	return top;
}

double f26(double sol[D])					// Composition Function 6
{
	int i, num, cf_num = 5;
	char func_num[4] = "26";
	double fit[5], top;
	double delta[5] = {10, 20, 20, 30, 40};
	double bias[5] = {0, 100, 200, 300, 400};
	
	num = 2;
	if(fileflag == 0)
		read_cec2017_cf( func_num, num, cf_num );

	i=0;
	movedata_cf( num, i, 1 );
	fit[i] = bf00( sol );
	fit[i] = 10000 * fit[i] / 2e+7;
	
	i=1;
	movedata_cf( num, i, 1 );
	fit[i] = f10( sol );

	i=2;
	movedata_cf( num, i, 1 );
	fit[i] = bf15( sol );
	fit[i] = 1000 * fit[i] / 100;

	i=3;
	movedata_cf( num, i, 1 );
	fit[i] = f4( sol );
	
	i=4;
	movedata_cf( num, i, 1 );
	fit[i] = f5( sol );
	fit[i] = 10000 * fit[i] / 1e+3;
	
	top = cf_cal( sol, delta, bias, fit, cf_num);

	return top;
}

double f27(double sol[D])					// Composition Function 7
{
	int i, num, cf_num = 6;
	char func_num[4] = "27";
	double fit[6], top;
	double delta[6] = {10, 20, 30, 40, 50, 60};
	double bias[6] = {0, 100, 200, 300, 400, 500};
	
	num = 2;
	if(fileflag == 0)
		read_cec2017_cf( func_num, num, cf_num );

	i=0;
	movedata_cf( num, i, 1 );
	fit[i] = bf18( sol );
	fit[i] = 10000 * fit[i] / 1000;
	
	i=1;
	movedata_cf( num, i, 1 );
	fit[i] = f5( sol );
	fit[i] = 10000 * fit[i] / 1e+3;

	i=2;
	movedata_cf( num, i, 1 );
	fit[i] = f10( sol );
	fit[i] = 10000 * fit[i] / 4e+3;

	i=3;
	movedata_cf( num, i, 1 );
	fit[i] = f1( sol );
	fit[i] = 10000 * fit[i] / 1e+30;
	
	i=4;
	movedata_cf( num, i, 1 );
	fit[i] = bf11( sol );
	fit[i] = 10000 * fit[i] / 1e+10;

	i=5;
	movedata_cf( num, i, 1 );
	fit[i] = bf00( sol );
	fit[i] = 10000 * fit[i] / 2e+7;
		
	top = cf_cal( sol, delta, bias, fit, cf_num);

	return top;
}

double f28(double sol[D])					// Composition Function 8
{
	int i, num, cf_num = 6;
	char func_num[4] = "28";
	double fit[6], top;
	double delta[6] = {10, 20, 30, 40, 50, 60};
	double bias[6] = {0, 100, 200, 300, 400, 500};
	
	num = 2;
	if(fileflag == 0)
		read_cec2017_cf( func_num, num, cf_num );

	i=0;
	movedata_cf( num, i, 1 );
	fit[i] = bf13( sol );
	fit[i] = 1000 * fit[i] / 100;
	
	i=1;
	movedata_cf( num, i, 1 );
	fit[i] = bf15( sol );
	fit[i] = 1000 * fit[i] / 100;

	i=2;
	movedata_cf( num, i, 1 );
	fit[i] = bf12( sol );
	fit[i] = 10000 * fit[i] / 1e+10;

	i=3;
	movedata_cf( num, i, 1 );
	fit[i] = f4( sol );
	
	i=4;
	movedata_cf( num, i, 1 );
	fit[i] = bf17( sol );
	fit[i] = 1000 * fit[i] / 1e+3;
	
	i=5;
	movedata_cf( num, i, 1 );
	fit[i] = bf00( sol );
	fit[i] = 10000 * fit[i] / 2e+7;
	
	top = cf_cal( sol, delta, bias, fit, cf_num);

	return top;
}

double f29(double sol[D])					// Composition Function 9
{
	int i, num, cf_num = 3;
	char func_num[4] = "29";
	double fit[3], top;
	double delta[3] = {10, 30, 50};
	double bias[3] = {0, 100, 200};

	num = 3;
	if(fileflag == 0)
		read_cec2017_cf( func_num, num, cf_num );

	i=0;
	movedata_cf( num, i, 1 );
	fit[i] = f15( sol );

	i=1;
	movedata_cf( num, i, 1 );
	fit[i] = f16( sol );

	i=2;
	movedata_cf( num, i, 1 );
	fit[i] = f17( sol );

	top = cf_cal( sol, delta, bias, fit, cf_num);

	return top;
}

double f30(double sol[D])					// Composition Function 10
{
	int i, num, cf_num = 3;
	char func_num[4] = "30";
	double fit[3], top;
	double delta[3] = {10, 30, 50};
	double bias[3] = {0, 100, 200};

	num = 3;
	if(fileflag == 0)
		read_cec2017_cf( func_num, num, cf_num );

	i=0;
	movedata_cf( num, i, 1 );
	fit[i] = f15( sol );

	i=1;
	movedata_cf( num, i, 1 );
	fit[i] = f18( sol );

	i=2;
	movedata_cf( num, i, 1 );
	fit[i] = f19( sol );

	top = cf_cal( sol, delta, bias, fit, cf_num);

	return top;
}


// **********  END CEC2017 **********



