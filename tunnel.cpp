#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#include"gpsys.cpp"

#include"plot/gnuplot_i.hpp" //Plotting lib
#include"plot/gnuplot_i.cpp"

#define BARRIERHEIGHT 0.1
#define BARRIERWIDHT 35.

#define K 0.03
#define STARTPOS -250
#define THERESHOLD 0.005
#define DSTART -2000.0
#define DEND 2000.0
#define SPACESTEPS 8192 //4096
#define TIMESTEP 2
#define HBAR 1.
#define MASS 1.
#define EFF_G_PARAM -0.03
#define ETA fabs(EFF_G_PARAM)*MASS/HBAR

//NOTE: BARRIER HAS TO BE CENTERED IN 0, DOMAIN MUST THEREFORE INCLUDE 0

double pot(double x);

void getEnergies(State* state, double* kinetic, double* chemical, double* external)
{
    *kinetic = 0;
    *chemical = 0;
    *external = 0;
    
    //external energy integral. Btw, yes... I like the ? operator. 
    for (int i = 0; i < state->getStepNum(); i++)
    {
        *external +=  
            pot(DSTART + i*state->getStep())*
            (state->wf[2*i]*state->wf[2*i] + state->wf[2*i+1]*state->wf[2*i+1])*
            state->getStep();
    }
    //chemical energy
    for (int i = 0; i < state->getStepNum(); i++)
        *chemical += 
            state->getStep()*
            (state->wf[2*i]*state->wf[2*i] + state->wf[2*i+1]*state->wf[2*i+1])*
            (state->wf[2*i]*state->wf[2*i] + state->wf[2*i+1]*state->wf[2*i+1])*
            state->getEffGasParameter()*HBAR;
    
    //Kinetic
    for (int i = 0; i < state->getStepNum() - 1; i++)
    {
        double real_d = (state->wf[2*i + 2] - state->wf[2*i])/state->getStep();
        double imag_d = (state->wf[2*i + 3] - state->wf[2*i+1])/state->getStep();
        
        *kinetic += HBAR*HBAR/2/MASS * (real_d*real_d + imag_d*imag_d) * state->getStep();
    }
    
    
    //printf("kin chem ext tot %f %f %f %f\n", *kinetic, *chemical, *external, *kinetic + *chemical + *external);
}

void plotState(Gnuplot* gp, State* state)
{
    FILE* fout = fopen("state.txt", "w");   
    fprintf(fout, "#x |psi|^2");
    for (int i = 0; i < state->getStepNum(); i++)
    {
        fprintf(fout, "%f %f\n", state->getDStart() + state->getStep()*i, 
                state->wf[2*i]*state->wf[2*i] + state->wf[2*i+1]*state->wf[2*i+1]);
    }
    fclose(fout);
    system("sleep 0.1");
    gp->reset_plot();
    
    gp->cmd("plot 'state.txt' u 1:2 w l lw 2");  
    
    return;
    
    double kin, chem, ext;
    getEnergies(state, &kin,&chem,&ext);
    gp->cmd("plot 'state.txt' u 1:2 w l lw 2, %f ti 'kin', %f ti 'chem', %f ti 'ext', %f ti 'tot' ", 
            10*kin,10*chem,10*ext, 10*(kin + chem + ext));
}

double pot(double x)
{
    if (x > -BARRIERWIDHT/2 && x < BARRIERWIDHT/2) return BARRIERHEIGHT;
    else return 0;
}

void start_state(double x, double* real_out, double* imag_out)
{
    *real_out =  1.0/cosh( (x-STARTPOS)*ETA)*cos(K*x);
    *imag_out =  1.0/cosh( (x-STARTPOS)*ETA)*sin(K*x);
}

bool done(State* state)
{   
    if (fabs(state->wf[ 2*((int) ( (DEND-DSTART)*0.1/state->getStep())) ] ) > THERESHOLD)
        return true;
    if (fabs(state->wf[ 2*((int) ( (DEND-DSTART)*0.9/state->getStep())) ] ) > THERESHOLD) 
        return true;
    
    return false;
}


int main(int argc, char** argv)
{
    GPSys sys = GPSys(DSTART, DEND, SPACESTEPS, TIMESTEP, MASS, HBAR, 0);
    sys.setPotential(&pot);
    sys.state->setState(&start_state);
    sys.state->normalize();
    
    Gnuplot* gp = new Gnuplot();
    gp->cmd("set yr[-0.005:0.02]");  
    gp->cmd("set term aqua size 1024,768");
    sys.printParams(stdout);
    
    double energy = (sys.getHbar()*sys.getHbar()/2.0/sys.state->getMass())*  //Calculate packet energy
    (K*K + 1.0/3.0*(ETA*ETA));
    
    printf("Packet mean k: %.3f\n", K);
    printf("Packet energy: %f\n", energy);
    printf("Barrier height: %f\n", BARRIERHEIGHT);
    printf("Barrier width: %.3f\n", BARRIERWIDHT);

    printf("--BEGIN--\n");
    
    FILE* en_out = fopen("en_out.txt", "w");
    fprintf(en_out, "#t kin chem ext tot\n");
    
    for(int j = 0; !done(sys.state); j++)
    {
        printf("t:%.3f\r", sys.getTime());
        sys.next();
        
        if (j % 100 == 0) 
        {
            plotState(gp, sys.state);  //Uncomment to plot frames during computation
                                       
            //double kin, chem, ext;
            //getEnergies(sys.state, &kin,&chem,&ext);
            //fprintf(en_out, "%f %f %f %f %f\n", sys.getTime(), kin, chem, ext, kin + chem + ext);
        }
    }
    
    fclose(en_out);
    
    gp->reset_plot();
    gp->cmd("set autoscale");
    gp->cmd("plot 'en_out.txt' u 1:2 w l ti 'kin', 'en_out.txt' u 1:3 w l ti 'chem', 'en_out.txt' u 1:4 w l ti 'ext','en_out.txt' u 1:5 w l ti 'tot'");
    
    /// ***RESULTS*******
    double Rsum = 0; //Get experimental transmission and reflection coefficient
    for (int i = 0; i < fabs(DSTART)/sys.state->getStep(); i++)
        Rsum += (sys.state->wf[2*i]*sys.state->wf[2*i] +
                 sys.state->wf[2*i+1]*sys.state->wf[2*i+1])*sys.state->getStep();
    
    double Tsum = 0;
    for (int i = (int) fabs(DSTART)/sys.state->getStep(); i < sys.state->getStepNum(); i++)
        Tsum += (sys.state->wf[2*i]*sys.state->wf[2*i] +
                 sys.state->wf[2*i+1]*sys.state->wf[2*i+1])*sys.state->getStep();    
    

    double t,r, k1; //Get theoric trasmission and reflection coefficents
    if (energy >= BARRIERHEIGHT)
    {
        k1 = sqrt(2*sys.state->getMass()*(energy - BARRIERHEIGHT)/pow(sys.getHbar(),2));
        t = 1.0/(1.0+((pow(BARRIERHEIGHT,2)*pow(sin(k1*BARRIERWIDHT),2))/
                      (4*energy*(energy-BARRIERHEIGHT))));
    }
    else
    {
        k1 = sqrt(2*sys.state->getMass()*(BARRIERHEIGHT - energy)/pow(sys.getHbar(),2));
        t = 1.0/(1.0+((pow(BARRIERHEIGHT,2)*pow(sinh(k1*BARRIERWIDHT),2))/
                      (4*energy*(BARRIERHEIGHT-energy))));
    }
    r = 1 - t;
        
    //Print results
    printf("\n--RESULTS--\n");
    printf("Reflection coefficent: %.3f (plane wave: %.3f)\n", Rsum, r);
    printf("Transmission coefficent: %.3f (plane wave: %.3f)\n", Tsum, t);
    printf("T + R: %.3f\n", Tsum + Rsum);
    return 0;
}
