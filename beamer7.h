/************************************************************************
  File Name: beamer7.h
  Author: Fang Yuan
  Mail: yfang@nju.edu.cn
  Created Time: Fri 15 Jul 2016 04:31:23 PM CST
 ************************************************************************/

#ifndef _BEAMER7_H
#define _BEAMER7_H

#include <math.h>

#define LOWPASS    1
#define HIGHPASS   (-1)

#define FS      16000     /* sampling rate */

#define CHANNEL  7        /* Microphone matrix */
#define TAPS     67       /* FIRs TAPS */
#define BUFL     50      /* buffer length in sampling data,both signal and noise */
#define SAMPLE   (BUFL/2) /* current sample to be processed */
#define CS       34000.0  /* speed of sound, cm/s*/

struct FilterData {
    float x1[2], x2[2];
    float y1[2], y2[2];
} ;

void filter_design(int low_pass, int freq, int type);
short IIR_direct_form_I(short input, struct FilterData *buf);

void Blackman_window(float f1, float f2, float *filter);

void beamerinit(int argc);
void beamformer();
float beamformer_filter();
void savefile(int);

short noise_add(int channel);
void noise_gen();
void filter_init();

void init_delays();
float signal_resample(int channel);
void update_signal();

#endif	// _BEAMER7_H

