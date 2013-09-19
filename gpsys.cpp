/*
 * Gross-Pitaevskij Equation time-splitting integrator in one dimension by Fabiano Lever, 
 * released under CreativeCommons Attribution Licence.
 */


#include<math.h>
#include<float.h>
#include <gsl/gsl_fft_complex.h>

/*
 * Represents a wawefunction in a finite discretized space
 */
class State 
{
public:
    State(double start, double end, int slices, double mass_in,
          double eff_gas_parameter_in);
    ~State();
    double *wf;             //Pointers to real and imaginary parts of the wawefunction
                            //even indexes are real part, odd imaginary
    void setState(void (*state)(double, double*, double*)); //Set a state from a given function (see below)
    void setState(void (*state)(double, double*, double*, void*), void*); //As above but parametric (see below)

    void normalize();                                       //Normalize the wavefunction
    
    inline int getStepNum() {return num;}                   //Used to make variables "read only"
    inline double getMass() {return mass;} 
    inline double getStep() {return dx;}
    inline double getDEnd() {return domain_end;}
    inline double getDStart() {return domain_start;}
    inline double getEffGasParameter() {return eff_gas_parameter;}

private:
    int num;                 //Number of slices
    double domain_start;     //Domain of the wawefunction starting and ending points
    double domain_end;
    double eff_gas_parameter;    //Effective gas parameter: number of bosons in the condensate
                                 //times the scattering lenght times the transverse trapping frequency
    double dx;               //Spatial step between points.
    double mass;             //Mass of the particles
};

State::State(double start, double end, int slices, double mass_in, 
             double eff_gas_parameter_in)
{
    num = slices;
    eff_gas_parameter = eff_gas_parameter_in;
    domain_start = start;
    domain_end = end;
    wf = new double[2*num];
    mass = mass_in;
    dx = (end - start)/slices;
}

State::~State()
{
    delete[] wf;
}

//first parameter of the given function is a point in space
//sencond and third paramenters are pointers to respectively 
//the real and  imaginary values of the wavefunction at that 
//point that are to be filled by the given function 
void State::setState(void (*state)(double, double*, double*)) 
{
    if (state)
        for (int i = 0; i < num; i++)
            state(domain_start + dx*i, &wf[2*i], &wf[2*i+1]);
}

//As the above function, but it allows to pass a parameter to 
//the state function
void State::setState(void (*state)(double, double*, double*, void*), void* param) 
{
    if (state)
        for (int i = 0; i < num; i++)
            state(domain_start + dx*i, &wf[2*i], &wf[2*i+1], param);
}


void State::normalize()
{
    double sum = 0;
    for (int i = 0; i < num; i++)
        sum += (wf[2*i]*wf[2*i] + wf[2*i + 1]*wf[2*i + 1])*dx;
    
    sum = sqrt(sum);
    for (int i = 0; i < num; i++)
    {
        wf[2*i] /= sum;
        wf[2*i + 1] /= sum;
    }
}

/*This class manages the simulation*/
class GPSys
{
public:
    //Constructor, see below for parameters meaning.
    GPSys(double start, double end, int slices, double dt_in, double mass, 
          double hbar, double eff_gas_parameter);
    ~GPSys();
    
    State* state;  //State of the system
    
    void next();   //Calls the update equations, moving state dt forward in time

    void setPotential(double (*pot)(double));       //Sets potential, pot must accept values between start and end.
                                                    //The time step gest automatically resized if needed.
    void setPotential(double (*pot)(double,void*),  //As above but parametric
                                    void*);      

    
    inline double getTime() {return time;}          //Used to make variables "read only"
    inline double getHbar() {return hbar;}    
    
    void printParams(FILE*);
    
private:    
    double *pot_phase;                 //Array that get filled with precalulated phase-shift 
                                       //values due to potential-particle interaction
    
    double *momentum_phase_real;       //Arrays of precalculated sines and cosines of
    double *momentum_phase_imag;       //momentum-space phase-shifts
    
    double dt;                         //Time step
    double hbar;                       //Reduced plank's constant
    double time;                       //Time elapsed since beginning of the simulation
};

/* GPSys Constructor:
 * <start> and <end> are the domain boundaries, <slices> is the number of space 
 * point to use, <dt_in> is the time step to be used, <mass>,  and <hbar> are what 
 * their name suggests; <eff_gas_parameter> (effective gas parameter) is the signed 
 * gas parameter (number of bosons in the condensate times the scattering lenght) 
 * times the transverse trapping frequency.
 */
GPSys::GPSys(double start, double end, int slices, double dt_in, double mass, 
             double hbar_in, double eff_gas_parameter)
{
    state = new State(start, end, slices, mass, eff_gas_parameter);
    hbar = hbar_in;
    time = 0;
    dt = dt_in;
    
    pot_phase = new double[state->getStepNum()];  
    setPotential(NULL);
    
    momentum_phase_real = new double[state->getStepNum()];   
    momentum_phase_imag = new double[state->getStepNum()];   
    for (int k = 0; k < state->getStepNum(); k++)
    {
        double theta;
        if (k <= state->getStepNum()/2)
            theta = -dt*hbar/2/state->getMass()*(k*2*M_PI/(state->getDEnd()-state->getDStart()))
            *(k*2*M_PI/(state->getDEnd()-state->getDStart()));
        else 
            theta = -dt*hbar/2/state->getMass()*((k-state->getStepNum())*2*M_PI/(state->getDEnd()-state->getDStart()))
            *((k-state->getStepNum())*2*M_PI/(state->getDEnd()-state->getDStart()));

        momentum_phase_real[k] = cos(theta);
        momentum_phase_imag[k] = sin(theta);
    }
}

GPSys::~GPSys()
{
    delete state;
    delete[] pot_phase;
    delete[] momentum_phase_real;
    delete[] momentum_phase_imag;
}

//Precaculate array of phase-shift values due to potential-particle interation
void GPSys::setPotential(double (*potential)(double)) 
{
    if (potential)
    {
        for (int i = 0; i < state->getStepNum(); i++)
            pot_phase[i] = -dt*potential(state->getDStart() + state->getStep()*i)/hbar;
    }
    else 
    {
        for (int i = 0; i < state->getStepNum(); i++)
        {
            pot_phase[i] = 0;
        }
    }
}

//****TODO: change everything*****
//Aa s above but allow a parameter to be passed to potential function
void GPSys::setPotential(double (*potential)(double, void*), void* param) 
{
    if (potential)
    {
        for (int i = 0; i < state->getStepNum(); i++)
            pot_phase[i] = -dt*potential(state->getDStart() + state->getStep()*i, param)/hbar;
    }
    else 
    {
        for (int i = 0; i < state->getStepNum(); i++)
        {
            pot_phase[i] = 0;
        }
    }
}


//Here lies the actual physics
void GPSys::next()
{    
    gsl_fft_complex_radix2_forward(state->wf, 1, state->getStepNum());
    for (int k = 0; k< state->getStepNum(); k++)
    {    
        double temp = state->wf[2*k];
        state->wf[2*k] = momentum_phase_real[k]*state->wf[2*k] 
                       - momentum_phase_imag[k]*state->wf[2*k + 1];
        
        state->wf[2*k+1] = momentum_phase_real[k]*state->wf[2*k + 1] 
                         + momentum_phase_imag[k]*temp;
    }
    gsl_fft_complex_radix2_inverse(state->wf, 1, state->getStepNum());
    
    for (int i = 0; i < state->getStepNum(); i++)
    {
        double temp = state->wf[2*i];
        double gp_phase = -dt*2*state->getEffGasParameter()*
                          (state->wf[2*i]*state->wf[2*i] + state->wf[2*i+1]*state->wf[2*i+1]);
        
        state->wf[2*i] =  cos(pot_phase[i] + gp_phase)*state->wf[2*i] - 
                          sin(pot_phase[i] + gp_phase)*state->wf[2*i+1];
        
        state->wf[2*i+1] =  cos(pot_phase[i] + gp_phase)*state->wf[2*i+1] + 
                          sin(pot_phase[i] + gp_phase)*temp;
    }
    time += dt;
}


void GPSys::printParams(FILE* stream)
{
    fprintf(stream, "#Simulation parameters at time t=%.3f\n", time);
    fprintf(stream, "#Function domain: [%.2f,%.2f]\n", state->getDStart(), state->getDEnd());
    fprintf(stream, "#Number of space cells: %i\n", state->getStepNum());
    fprintf(stream, "#Spatial step: %.4f\n", state->getStep());
    fprintf(stream, "#Particles mass: %.2f\n", state->getMass());
    fprintf(stream, "#Effective Gas parameter: %.4f\n", state->getEffGasParameter());
    fprintf(stream, "#Reduced Plank's constant: %.4f\n", hbar);
}