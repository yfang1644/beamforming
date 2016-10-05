/************************************************************************
  File Name: filtertest.h
  Author: Fang Yuan
  Mail: yfang@nju.edu.cn
  Created Time: Thu 07 Jul 2016 03:39:58 PM CST
 ************************************************************************/

#ifndef _FILTERTEST_H
#define _FILTERTEST_H

#include <math.h>

#define LOWPASS    1
#define HIGHPASS   (-1)

#define FS      16000     /* sampling rate */

#define CHANNEL  (7*9)    /* Microphone matrix */
#define TAPS    67       /* FIRs TAPS */

struct FilterData {
    float x1[2], x2[2];
    float y1[2], y2[2];
} ;

void filter_design(int low_pass, int freq, int type);
short IIR_direct_form_I(short input, struct FilterData *buf);

void Blackman_window(float f1, float f2, float *filter);

void beamerinit(int argc);
void beamformer();

short noise_add(int channel);

#endif

