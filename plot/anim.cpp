/*
 * Examples of gnuplot_i.c usage
 */

#include <stdio.h>
#include <stdlib.h>

#include "gnuplot_i.hpp"

#define SLEEP_LGTH  1
int main(int argc, char *argv[]) 
{
    Gnuplot gp = Gnuplot("lines");
    double phase;

    printf("*** example of gnuplot control through C ***\n") ;

    for (phase=0.1 ; phase<10 ; phase +=0.1) {
        gp.reset_plot();
        gp.cmd("plot sin(x+%g)", phase);
    }
    
    for (phase=10 ; phase>=0.1 ; phase -=0.1) {
        gp.reset_plot();
        gp.cmd("plot sin(x+%g)", phase);
    }
    
    delete &gp;  
    return 0 ;
}

