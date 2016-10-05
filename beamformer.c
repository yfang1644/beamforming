/************************************************************************
File Name: beamformer.c
Author: Fang Yuan
Mail: yfang@nju.edu.cn
Created Time: Sun 17 Jul 2016 12:59:23 PM CST
************************************************************************/

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <signal.h>
#include "beamer7.h"

float	mu_beamer = 2.e-11;        	// 阵列自适应步长
float	beamer[CHANNEL][TAPS];        // 滤波器系数(双字对齐)
short   input[CHANNEL][TAPS];        // 阵列输入序列,带训练噪声源叠加
short   mic[CHANNEL][TAPS];          // 阵列输入序列,不带训练噪声源
short   sum_input[TAPS];         	// signal tap sum

float	tap_sum[TAPS];                // filter tap sum
float   desired[TAPS];

float	r_chn = 1.0/CHANNEL;

void beamerinit(int argc)
{
    FILE *fp;
    int t, ch;
    pid_t pid;

    #ifdef FLAT
    for(t = TAPS; --t;)
       desired[t] = 0.0;
    desired[0] = 1.0;
    #else
    Blackman_window(300.0, FS/2.0, desired);
    #endif
    for(t = TAPS; t--;) {
        desired[t] *= r_chn;
        for(ch = CHANNEL; ch--;) {
            beamer[ch][t] = desired[t];
        }
    }

    if(argc == 2) {
        fp = fopen("filter.dat", "rb");
        fread(beamer, 4, CHANNEL*TAPS, fp);
        fclose(fp);
    }
}

//	y=sum(x .* w);
//	mu_y= mu*y;
// 	for j=1:J
//        w(:,j)=w(:,j) - mu_y*(x(:,j)-sum(x(:,j))/K) + (f(j)-sum(w(:,j)))/K;
//	end;

void beamformer()
{
    int	ch, t, n, i;
    int tmp;
    float	res;
    static int cnt = 0;
    static float energy = 0;
    static float ei = 0;

    res = 0;
    for(t = TAPS; t--;) {
        tap_sum[t] = 0.0;
        tmp = 0;
        for(ch = CHANNEL; ch--;) {
            tap_sum[t] += beamer[ch][t];
            tmp += input[ch][t];
            res += input[ch][t] * beamer[ch][t];
        }
        tap_sum[t] *= r_chn;
        sum_input[t] = (short)(tmp * r_chn);
    }

    energy += res*res;
    ei += (float)input[0][0]*input[0][0];
    cnt++;
    if(cnt == 5120) {
        printf("%8.1f   %8.1f\n", sqrt(energy), sqrt(ei));
        energy = 0;
        ei = 0;
        cnt = 0;
    }

    res = mu_beamer * res;

    for (ch = CHANNEL; ch--;) {
        for(t = TAPS; --t;) {
            beamer[ch][t] += res*(sum_input[t]-input[ch][t]) -
               (tap_sum[t] - desired[t]);
            input[ch][t] = input[ch][t-1];
            mic[ch][t] = mic[ch][t-1];
        }

        beamer[ch][0] += res*(sum_input[0]-input[ch][0]) -
            (tap_sum[0] - desired[0]);
        mic[ch][0] = signal_resample(ch);
        input[ch][0] = mic[ch][0] + noise_add(ch);
    }
}

/* only filter microphone signal */
float beamformer_filter()
{
    int	ch, t;
    float res = 0;

    for(t = TAPS; t--;) {
        for(ch = CHANNEL; ch--;) {
            res += mic[ch][t] * beamer[ch][t];
        }
    }

    return res;
}

void savefile(int num)
{
    FILE *fp;
    fp = fopen("filter.dat", "wb");
    fwrite(beamer, 4, CHANNEL*TAPS, fp);
    fclose(fp);

    exit(0);
}

