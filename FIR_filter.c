/************************************************************************
  File Name: FIR_filter.c
  Author: Fang Yuan
  Mail: yfang@nju.edu.cn
  Created Time: Thu 07 Jul 2016 10:18:25 PM CST
 ************************************************************************/

#include "beamer.h"

void Blackman_window(float f1, float f2, float *filter)
{
	float window, hd;
	float omega1, omega2;
	int i;

	omega1 = 2 * M_PI * f1 / FS;
	omega2 = 2 * M_PI * f2 / FS;

	for(i = 1; i <= (TAPS - 1) / 2; i++)
	{
		window = 0.42 + 0.5 * cos(2 * M_PI * i / (TAPS - 1))
			+ 0.08 * cos(4 * M_PI * i / (TAPS - 1));
		hd = (sin(omega2 * i) - sin(omega1 * i)) / M_PI / i;
		filter[(TAPS - 1)/2 + i] = hd * window;
		filter[(TAPS - 1)/2 - i] = filter[(TAPS - 1)/2 + i];
	}
	hd = (omega2 - omega1) / M_PI;
	filter[(TAPS - 1)/2] = hd;
}

// End of Filter.c
