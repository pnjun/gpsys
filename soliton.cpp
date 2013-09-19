#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#include"gpsys.cpp"

#include"plot/gnuplot_i.hpp" //Plotting lib
#include"plot/gnuplot_i.cpp"

void plotState(Gnuplot* gp, State* state)
{
    FILE* fout = fopen("out.txt", "w");
    for (int i = 0; i < state->getStepNum(); i++)
        fprintf(fout, "%f %f\n", state->getDStart() + state->getStep()*i, 
                state->wf[2*i]*state->wf[2*i] + state->wf[2*i+1]*state->wf[2*i+1]);
    fclose(fout);
    system("sleep 0.1");
    gp->reset_plot();
    gp->cmd("plot 'out.txt' u 1:2 w l");
}

double getSigma(State* state)
{
    double sum = 0;
    double x = 0;
    for (int i = 0; i < state->getStepNum(); i++)
    {
        x = (i - state->getStepNum()/2)*state->getStep();
        sum += (state->wf[2*i]*state->wf[2*i] + state->wf[2*i+1]*state->wf[2*i+1])*
        x*x*state->getStep();
    }
    return sqrt(sum*2);
}

#define K 0.0
#define STARTPOS 0.0
#define MAXT 270.0
#define DSTART -800.0
#define DEND 800.0
#define SPACESTEPS 4096
#define TIMESTEP 0.02
#define HBAR 1
#define MASS 1
#define EFF_G_PARAM -0.03
#define ETA fabs(EFF_G_PARAM)*MASS/HBAR
void start_state(double x, double* real_out, double* imag_out)
{
    *real_out =  1.0/cosh( (x-STARTPOS)*ETA)*cos(K*x);
    *imag_out =  1.0/cosh( (x-STARTPOS)*ETA)*sin(K*x);
}

int main(int argc, char** argv)
{
    GPSys sys = GPSys(DSTART, DEND, SPACESTEPS, TIMESTEP, MASS, 
                      HBAR, EFF_G_PARAM);
    sys.state->setState(&start_state);
    sys.state->normalize();
    
    Gnuplot* gp = new Gnuplot();
    gp->cmd("set yr[0:0.02]");
    system("sleep 0.1");
    
    sys.printParams(stdout);
    FILE* fout = fopen("stddev.txt","w");
    fprintf(fout, "#t stddev\n");
    fprintf(fout, "%f %f\n",sys.getTime(), getSigma(sys.state));
    for (int j = 0; sys.getTime() < MAXT; j++)
    {
        printf("t:%.3f/%f\r", sys.getTime(), MAXT);
        if (j % 300 == 0) fprintf(fout, "%f %f\n",sys.getTime(), getSigma(sys.state));
        sys.next();
        if (j % 300 == 0) plotState(gp, sys.state);  //Uncomment to plot frames during computation
    }
    fclose(fout);
    return 0;
}