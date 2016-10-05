/************************************************************************
  File Name: channel7.c
  Author: Fang Yuan
  Mail: yfang@nju.edu.cn
  Created Time: Fri 15 Jul 2016 04:20:57 PM CST
 ***********************************************************************/

/*
             4
        3         5
             7
        2         6
             1

    话筒在半径4.25cm的圆周上均匀分布。中心一个(7号)。

    采样率 FS=16000Hz
    声速 c=34000cm/s
    话音来自话筒平面
    以中心为参照点，远场延迟计算方法：
        d = r*cos(theta - alpha_i)

*/

#include <stdio.h>
#include <signal.h>
#include "beamer7.h"

int main(int argc, char *argv)
{
    signal(SIGINT, savefile);

    init_delays();
    beamerinit(argc);
    filter_init();

    for(;;) {
        update_signal();
        noise_gen();
        beamformer();
        beamformer_filter();
    }   
}
// main.c
