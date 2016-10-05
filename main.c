/************************************************************************
  File Name: main.c
  Author: Fang Yuan
  Mail: yfang@nju.edu.cn
  Created Time: Fri 24 Jun 2016 12:35:15 PM CST
 ***********************************************************************/

/*
    o     o     o     o     o     o     o     o     o
 
    o     o     o     o     o     o     o     o     o
 
    o     o     o     o     o     o     o     o     o
 
    o     o     o     o     o     o     o     o     o
 
    o     o     o     o     o     o     o     o     o
 
    o     o     o     o     o     o     o     o     o
 
    o     o     o     o     o     o     o     o     o

    9x7阵列，Dx x Dy (cmxcm)
    水平方向行 dx=Dx/8(cm) 间距，垂直方向列 dy=Dy/6(cm) 间距 
    采样率 FS=48000Hz
    声速 c=34000cm/s
    水平方向角 theta, 垂直方向角 phi
    水平方向延迟 t(n)=dx*n*sin(theta)/c
    垂直方向延迟 t(m)=dy*m*sin(phi)/c

*/

#include <stdio.h>
#include <stdlib.h>
#include <signal.h>
#include "beamer.h"

float	mu_beamer = 1.e-12;        	// 阵列自适应步长
float	beamer[CHANNEL][TAPS];        // 滤波器系数(双字对齐)
short   input[CHANNEL][TAPS];        // 阵列输入序列
short   sum_input[TAPS];        	// signal tap sum

float	tap_sum[TAPS];                // filter tap sum
float   desired[TAPS];

float	r_chn = 1.0/CHANNEL;

float Dx = 35.0, dx = 35.0/8;
float Dy = 25.0, dy = 25.0/6;
float c = 34000.0;
float direction[CHANNEL] = {    // 每路延迟，单位为秒或毫秒
  0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0
};

short noise[4][100];        // 4路噪声源，每路保留100点(左右两边最大延迟)
struct FilterData fl[4];

#define SAMPLE  50          // 从中间取点，延迟有正负
float noisedir[4][2] = {    // 噪声源方位, theta, phi
    40, 50,
    -50, 30,
    50, -40,
    -45, -40
    };

FILE *fpn;

void noise_gen()
{
    int i, j;
    int Ltmp;
    short tmp;
    static int cnt = 0;

    for(i = 0; i < 4; i++) {
        for(j = 100; --j;)    
            noise[i][j] = noise[i][j-1];
        Ltmp = rand() -(RAND_MAX/2);
        tmp = (short)(Ltmp >> 18);
        noise[i][0] = IIR_direct_form_I(tmp, &fl[i]);
//        noise[i][0] = 5000*sin(2*M_PI*(i+1)*329/FS*cnt)
//                   + 3000*cos(2*M_PI*(i+1)*763/FS*cnt);
    }
    fwrite(&noise[0][0], 2, 1, fpn);
    cnt++;
}

short noise_add(int ch)
{
    int i, j, n, nt;
    float dnt, x1, x2, t;
    float sample = 0;

    i = ch % 9;
    j = ch / 9;
    for(n = 0; n < 4; n++) {
        t = (dx*i*sin(noisedir[n][0]) + dy*j*sin(noisedir[n][1]))/c;
        dnt = t * FS;
        nt = (int)dnt;
        dnt -= (float)nt;
        x1 = noise[n][SAMPLE+nt];
        x2 = noise[n][SAMPLE+nt+1];
        t = x1 + (x2 - x1)*dnt;
        sample += t;
    }
    return (short)sample;
}

void init_array(int argc)
{
    int i, j;
	int ch, t;
    float theta = 0, phi = 0;
    FILE *fp;

    for(j = 0; j < 7; j++)
        for(i = 0; i <9; i++)
            direction[j*9+i] = (dx*i*sin(theta) + dy*j*sin(phi))/c;

    for(i = 0; i < 4; i++) {        // degree -> arc
        noisedir[i][0] *= M_PI/180.0;
        noisedir[i][1] *= M_PI/180.0;
    }

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
			input[ch][t] = 0;
		}
	}

    if(argc == 2) {
        fp = fopen("filter.dat", "rb");
        fread(beamer, 4, CHANNEL*TAPS, fp);
        fclose(fp);
    }
}

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
        }

        beamer[ch][0] += res*(sum_input[0]-input[ch][0]) -
                            (tap_sum[0] - desired[0]);
        input[ch][0] = noise_add(ch);
    }
}

//	y=sum(x .* w);
//	mu_y= mu*y;
// 	for j=1:J
//        w(:,j)=w(:,j) - mu_y*(x(:,j)-sum(x(:,j))/K) + (f(j)-sum(w(:,j)))/K;
//	end;


void savefile()
{
    FILE *fp;
    fp = fopen("filter.dat", "wb");
    fwrite(beamer, 4, CHANNEL*TAPS, fp);
    fclose(fp);

    fclose(fpn);
    exit(0);
}

int main(int argc, char *argv)
{
    int ch, sec;
    signal(SIGINT, savefile);

    fpn = fopen("noise.dat", "wb");
    for(ch = 0; ch < 4; ch++) {
		for(sec = 0; sec < 2; sec++) {
			fl[ch].x1[sec] = fl[ch].x2[sec] = 0.0;
			fl[ch].y1[sec] = fl[ch].y2[sec] = 0.0;
		}
	}

    filter_design(HIGHPASS, 300, 3);
    init_array(argc);
    for(;;) {
        noise_gen();
        beamformer();
    }   
}
// main.c
