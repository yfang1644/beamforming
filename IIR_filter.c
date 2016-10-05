#include "beamer7.h"

/*
	Function: Calculate HPF/LPF coefficients, max. 6th order
	Import:	  type,Freq,type, freq in High Pass Filter/Low Pass Filter
	type:
		0: direct
		1: Butterworth, 2nd order(-12dB/oct)
		2: Butterworth, 3rd order(-18dB/oct)
		3: Butterworth, 4th order(-24dB/oct)
		4: Butterworth, 5th order(-30dB/oct)
		5: Butterworth, 6th order(-36dB/oct)
		6: Bessel, 2nd order
		7: Bessel, 3nd order
		8: Bessel, 4th order
		9: Bessel, 5th order
		10:Bessel, 6th order
		11:Linkwitz-Riley, 2nd order
		12:Linkwitz-Riley, 4th order

*/
#define	IIR_SECOND_ORDER(param)	\
	if(1 == low_pass)	rk = wo * wc;		\
	else	rk = wo / wc;					\
	temp = 1.0 / (rk * (rk + pk) + 1.0);	\
	param[0] = rk * rk * temp;				/*b2*/\
	param[1] = 2.0 * low_pass * param[0];	/*b1*/\
	param[2] = param[0];					/*b0*/\
	param[3] = (rk * (rk - pk) + 1.0) * temp;	/*a2*/\
	param[4] = 2.0 * low_pass * (rk * rk - 1) * temp;/*a1*/

#define	IIR_FIRST_ORDER(param)	\
	if(1 == low_pass)	rk = wo * wc;			\
	else	rk = wo / wc;						\
	temp = 1.0 / (rk + 1.0);					\
	param[2] = rk * temp;					/*b0*/\
	param[1] = low_pass * param[2];			/*b1*/\
	param[4] = low_pass * (rk - 1.0) * temp;/*a1*/

float coeff[2][5];

void filter_design(int low_pass, int freq, int type)
{
	float coeff_temp1[5];	//temp for exchanging..
	float coeff_temp2[5];	//temp for exchanging..
	// low_pass = 1 for low pass filter, low_pass = -1 for high pass filter
	float wc, wo, rk, temp, pk;

	//Direct refresh coeff_temp[], a0 = 1
	coeff_temp1[0] = 0;//b2
	coeff_temp1[1] = 0;//b1
	coeff_temp1[2] = 1;//b0
	coeff_temp1[3] = 0;//a2
	coeff_temp1[4] = 0;//a1

	coeff_temp2[0] = 0;//b2
	coeff_temp2[1] = 0;//b1
	coeff_temp2[2] = 1;//b0
	coeff_temp2[3] = 0;//a2
	coeff_temp2[4] = 0;//a1

	wc = tan(M_PI * freq / FS);
	switch (type)
	{
		case 1:	//Butterworth, 2nd order(-12dB/oct)
			wo = 1.0;
			pk = 1.414213562373095;		// 2 * sin(PI/4);
			IIR_SECOND_ORDER(coeff_temp1);
			break;

		case 2:	//Butterworth, 3rd order(-18dB/oct)
			wo = 1.0;
			pk = 1.0;				// 2 * sin(PI/6);
			IIR_SECOND_ORDER(coeff_temp1);

			IIR_FIRST_ORDER(coeff_temp2);
			break;

		case 3:	//Butterworth, 4th order(-24dB/oct)
			wo = 1.0;
			pk = 0.76536686473018;			// 2 * sin(PI/8);
			IIR_SECOND_ORDER(coeff_temp1);

			pk = 1.84775906502257351;		// 2 * sin(3*PI/8);
			IIR_SECOND_ORDER(coeff_temp2);
			break;

		case 6:	//Bessel, 2nd order(-12dB/oct)
			wo = 1.732050808;
			pk = 3.0 / wo;
			IIR_SECOND_ORDER(coeff_temp1);
			break;

		case 7:	//Bessel, 3rd order(-18dB/oct)
			wo = 2.541541401;
			pk = 3.6778146454 / wo;
			IIR_SECOND_ORDER(coeff_temp1);

			wo = 2.3221853546;
			IIR_FIRST_ORDER(coeff_temp2);
			break;

		case 8:	//Bessel, 4th order(-24dB/oct)
			wo = 3.389365793;
			pk = 4.2075787944 / wo;
			IIR_SECOND_ORDER(coeff_temp1);

			wo = 3.023264939;
			pk = 5.7924212056 / wo;
			IIR_SECOND_ORDER(coeff_temp2);
			break;

		case 11://Linkwitz-Riley, 2nd order(-12dB/oct)
			wo = 1.0;
			pk = 2.0;
			IIR_SECOND_ORDER(coeff_temp1);
			break;

		case 12://Linkwitz-Riley, 4th order(-24dB/oct)
			wo = 1;
			pk = 1.414213562373095;						//2*sin(PI/4)
			IIR_SECOND_ORDER(coeff_temp1);

			IIR_SECOND_ORDER(coeff_temp2);
			break;

		default:
			break;
	}

/*
	printf("b0, b1, b2 -- a0, a1, a2\n");
	printf("p1=[%4.8f %4.8f %4.8f]; q1=[1.0 %4.8f %4.8f];\n",
			coeff_temp1[2], coeff_temp1[1], coeff_temp1[0],
			coeff_temp1[4], coeff_temp1[3]);
	printf("p2=[%4.8f %4.8f %4.8f]; q2=[1.0 %4.8f %4.8f];\n",
			coeff_temp2[2], coeff_temp2[1], coeff_temp2[0],
			coeff_temp2[4], coeff_temp2[3]);
*/
	coeff[0][0] = coeff_temp1[0];
	coeff[0][1] = coeff_temp1[1];
	coeff[0][2] = coeff_temp1[2];
	coeff[0][3] = coeff_temp1[3];
	coeff[0][4] = coeff_temp1[4];
	
	coeff[1][0] = coeff_temp2[0];
	coeff[1][1] = coeff_temp2[1];
	coeff[1][2] = coeff_temp2[2];
	coeff[1][3] = coeff_temp2[3];
	coeff[1][4] = coeff_temp2[4];
}
/*
  y(n) = coeff2*x(n)+coeff1*x(n-1)+coeff0*x(n-2)-coeff4*y(n-1)-coeff3*y(n-2)


p1=[0.95363986 -1.90727973 0.95363986]; q1=[1.0 -1.90064657 0.91391283];
p2=[0.89892012 -1.79784024 0.89892012]; q2=[1.0 -1.79158771 0.80409276];

z=poly(0, 'z');
r = [z^2 z 1]';
hz =( p1*r )/(q1*r)*(p2*r)/(q2*r);

[hzm,fr]=frmag(hz,256);
*/

/*****************************************************************************/
/* DESCRIPTION                                                               */
/*   Infinite Impulse Response (IIR) filters fourth order type I and type II */
/*   Takes 3 numerator coefficients and 3 denominator coefficients.          */
/*                                                                           */
/*---------------------------------------------------------------------------*/
/*                                                                           */
/* An second infinite impulse response (IIR) filter can be represented by    */
/* the following equation:                                                   */
/*                                                                           */
/*        b0 + b1.z^-1 + b2.z^-2                                             */
/* H(z) = ----------------------                                             */
/*        a0 + a1.z^-1 + a2.z^-2                                             */
/*                                                                           */
/* where H(z) is the transfer function. a0 is always 1.000                   */
/*                                                                           */
/* To implement a fourth order filter, two of these stages are cascaded.     */
/*                                                                           */
/*****************************************************************************/

/* Numerator coefficients */
#define B0 2
#define B1 1
#define B2 0

/* Denominator coefficients */
#define A0 0
#define A1 4
#define A2 3

/*****************************************************************************/
/* IIR_direct_form_I()                                                       */
/*---------------------------------------------------------------------------*/
/*                                                                           */
/* Forth order direct form I IIR filter implemented by cascading 2 second    */
/* order filters.                                                            */
/*                                                                           */
/* This implementation uses two buffers, one for x[n] and the other for y[n] */
/*                                                                           */
/*****************************************************************************/

short IIR_direct_form_I(short input, struct FilterData *buf)
{
	float temp, y0;
	/* x(n), x(n-1), x(n-2), y(n), y(n-1), y(n-2). Must be static */
	unsigned int stages;

	temp = (float) input; /* Copy input to temp */

	for (stages = 0; stages < 2; stages++) {
		y0  = coeff[stages][B0] * temp;	            /* B0 * x(n)   */
		y0 += coeff[stages][B1] * buf->x1[stages];	/* B1 * x(n-1) */
		y0 += coeff[stages][B2] * buf->x2[stages];	/* B2 * x(n-2) */

		buf->x2[stages] = buf->x1[stages];          /* x(n-2) = x(n-1) */
		buf->x1[stages] = temp;                     /* x(n-1) = x(n)   */

		y0 -= coeff[stages][A1] * buf->y1[stages];	/* A1 * y(n-1) */
		y0 -= coeff[stages][A2] * buf->y2[stages];	/* A2 * y(n-2) */

		/* Shuffle values along one place for next time */

		buf->y2[stages] = buf->y1[stages];          /* y(n-2) = y(n-1) */
		buf->y1[stages] = y0;                       /* y(n-1) = y(n)   */

		/* temp is used as next section input*/
        temp =  y0;
	}

	return (short)temp;
}
/*
int main(int argc, char *argv[])
{
    struct FilterData fl[4];
    short inp, out;
    filter_design(-1, 300, 3);
    out  = IIR_direct_form_I(inp, &fl[2]);
}	
*/
