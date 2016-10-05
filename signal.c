/************************************************************************
  File Name: signal.c
  Author: Fang Yuan
  Mail: yfang@nju.edu.cn
  Created Time: Sun 17 Jul 2016 01:36:52 PM CST
 ************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "beamer7.h"

float   sig_dir = -108.0;              // 信号方向

// microphone positions in degree. R=4.25cm
float radius = 4.25;  /* cm */
float alpha[CHANNEL] = {270, 210, 150, 90, 30, -30, 0};
float delay[CHANNEL];   // 每路延迟，单位为采样点，浮点数

short source[BUFL][CHANNEL]; // 信号源

FILE *fp_signal;
void init_delays()
{
    int ch;
	printf("ok\n");
    for(ch = CHANNEL; --ch; )
        delay[ch] = FS * radius*cos((alpha[ch]-sig_dir)*M_PI/180.)/CS;

    fp_signal = fopen("/home/fang/Downloads/zero_degree.wav", "rb");
    fseek(fp_signal, 44, SEEK_SET);
}

void update_signal()
{
    int chn, tap;
    for(chn = CHANNEL; chn--; )
        for(tap = BUFL; --tap;)
            source[tap][chn] = source[tap-1][chn];

    chn = fread(source, 2, CHANNEL, fp_signal);
	if(chn < CHANNEL) {
		fseek(fp_signal, 44, SEEK_SET);
		printf("#########\n");
	}
}

float signal_resample(int ch)
{
    int n;
    float t, dnt, x1, x2;
    float sample = 0;

    if(ch == CHANNEL-1)
        sample = source[ch][SAMPLE];
    else {
        n = (int)delay[ch];
        dnt = delay[ch] - (float)n;
        x1 = source[ch][SAMPLE+n];
        x2 = source[ch][SAMPLE+n+1];
        sample = x1 + (x2 - x1)*dnt;
    }

    return sample;
}

