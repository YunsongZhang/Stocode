#ifndef RACIPELIB_h
#define RACIPELIB_h

#include <stdio.h>

struct topo {
	int  numG;
	int  numR;
	int  *SourceG;
	int  *TargetG;
	int  *TypeR;
	char **Gname;
	int  *ParasPos;
	double **prsrandrange;
    int  *numover;
    int  *numdown;
};

struct opts {
	int    flag; 			  //										  - Only produce .cfg file or not
	int    KDID;			  // 										  - The gene or link to knockdown	
	char   *distname;
	int    dist;      		  // "Distribution"            			      - Distribution used for randomization
	double SF;    			  //                           			      - Scale factor
	int    num_findT;         // "NumberOfSimulationToFindThreshold"      - Number of simulaitons used to estimate threshold;
	int    num_paras;      	  // "NumberOfRACIPEModels"    			      - Number of RACIPE models to generate
	int    num_ode;      	  // "NumberOfRIVs"            			      - Number of Random initial values to solve ODEs
	int    num_stability;     // "NumberOfStatesToStore"   			      - Number of stable states to count
	double thrd;    		  // "ThresholdForConvergence" 			      - Cutoff for convergence of states
	int    Toggle_f_p;        // "ToggleOfSaveParameters"  			      - Toggle to save parameter or not
	double stepsize;    	  // "Stepsize"                			      - Stepsize for solving ODE;
	int    maxiters;      	  // "MaximumOfIterations"     			      - Maximum of Iteration for solving ODE at each RIVs
	int    Toggle_T_test;     // "TestThreshold"           			      - Toggle to test threshold assumption or not
};

extern double Hillshift (double x, double x0, double nx, double lamda);
extern double randu(double minV, double maxV);
extern double randg(double m, double stdvalue);
extern double randpg(double m, double stdvalue);
extern double randexp(double m);
extern double randfd(double m, double stdvalue, int dist);
extern double randd(int dist);
extern double median(double *x, int n);
extern double sumdelta(double *y, double *ytmp, int NEQN);

extern void   initial_simuopts (char *fileextension, int argc, char **argv, struct opts *simu_opts);
extern void   Model_generate (char *topofile, struct topo *topoinfo, struct opts *simu_opts);
extern void   check_topo(FILE *f_in, FILE *f_out, struct topo *topoinfo);
extern void   generate_random_range(FILE *f_paras, struct topo *topoinfo, int num, int dist, double SF);
extern void   estimate_threshold(int num, int ID, double minP, double maxP, double minK, double maxK, double minN, double maxN, double minF, double maxF, double *minT, double *maxT, struct topo *topoinfo, int dist, double SF);
extern void   generate_ODE(FILE *f_model, struct topo *topoinfo);
extern void   generate_riv(FILE *f_rivs, struct topo *topoinfo);

extern void   read_cfg(char *configfile, struct topo *topoinfo, struct opts *simu_opts);
extern void   save_model_paras(char *inputfile, struct opts *simu_opts, struct topo *topoinfo, int num, int cnt_loop, int cnt, double *paras);
extern void   save_model_solns(char *inputfile, struct opts *simu_opts, struct topo *topoinfo, int num, int cnt_loop, int cnt, double *soln);
extern void   model_ODE(double t,  double *ytmp, double *yp, double *p, struct topo *topoinfo);
extern void   RIVs(double *y, double *ytmp, double *p, struct topo *topoinfo);
extern int    solve_ODE (double *y_store, int j, double *paras, struct opts *simu_opts, struct topo *topoinfo);
extern void   set_parameters ( double *paras, struct opts *simu_opts, struct topo *topoinfo);
extern void   T_test(char *inputfile, struct opts *simu_opts, struct topo *topoinfo, int num, int cnt, double *soln, double *paras);
extern int    count_state (double *y_store, double *soln, struct opts *simu_opts, struct topo *topoinfo);
extern void   release_memory(char *fileextension, struct topo *topoinfo, int *cnt_store, double *y_store, double *soln, double *paras);

#endif
