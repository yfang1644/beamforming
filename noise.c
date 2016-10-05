/************************************************************************
  File Name: noise.c
  Author: Fang Yuan
  Mail: yfang@nju.edu.cn
  Created Time: Sun 17 Jul 2016 01:24:28 PM CST
 ************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "beamer7.h"

short noise[4][BUFL];        // 4路噪声源，每路保留100点(左右两边最大延迟)
struct FilterData fl[4];
extern float radius;
extern float alpha[];       /* microphone positions defined in signal.c */
extern float delay[];       /* delay */
float noisedir[4] = {       // 噪声源方位, theta
    120, 
    -50,
    50, 
    -45,
    };

void filter_init()
{
    int ch, sec;
    for(ch = 0; ch < 4; ch++) {
		for(sec = 0; sec < 2; sec++) {
			fl[ch].x1[sec] = fl[ch].x2[sec] = 0.0;
			fl[ch].y1[sec] = fl[ch].y2[sec] = 0.0;
		}
	}

    filter_design(HIGHPASS, 300, 3);
}

void noise_gen()
{
    int i, j;
    int Ltmp;
    short tmp;
    static int cnt = 0;

    for(i = 0; i < 4; i++) {
        for(j = BUFL; --j;)    
            noise[i][j] = noise[i][j-1];
        Ltmp = rand() -(RAND_MAX/2);
        tmp = (short)(Ltmp >> 18);
//        noise[i][0] = IIR_direct_form_I(tmp, &fl[i]);
//        noise[i][0] = 5000*sin(2*M_PI*(i+1)*329/FS*cnt)
//                   + 3000*cos(2*M_PI*(i+1)*763/FS*cnt);
        noise[i][0] = tmp;
    }
    cnt++;
}

short noise_add(int ch)
{
    int  n, nt;
    float dnt, x1, x2, t;
    float sample = 0;

    if(ch == CHANNEL-1) {
        for(n = 0; n < 4; n++) {
            sample += noise[n][SAMPLE];
        }
    } else {
        for(n = 0; n < 2; n++) {
            dnt = FS * radius*cos((alpha[ch]-noisedir[n])*M_PI/180.)/CS;
            dnt = delay[ch] - dnt;
            nt = (int)dnt;
            dnt -= (float)nt;
            x1 = noise[n][SAMPLE+nt];
            x2 = noise[n][SAMPLE+nt+1];
            sample += x1 + (x2 - x1)*dnt;
        }
    }
    return (short)sample;
}

