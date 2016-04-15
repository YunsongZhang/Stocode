# include <stdlib.h>
# include <stdio.h>
# include <stdint.h>
# include <math.h>
# include <time.h>
# include <string.h>

# include "RACIPELIB.h"
# include "pcg_basic.h"

int main (int argc, char **argv)
{  
  clock_t begin, end;
  double  time_spent;

  begin = clock();

  struct  topo   topoinfo;
  struct  opts   simu_opts;

  int  i        = 0;
  int  j        = 0;
  int  cnt      = 0;
  int  tmp_loop = 0;
  int  cnt_loop = 0;
  char *token;
  char *modelname;
  char *fileextension;
  char *inputfile;

  double *y_store;        // y_store to store the solutions for each RIVs
  double *soln;           // soln to store the solutions for each RACIPE model
  int    *cnt_store;      // cnt_store to count the stability
  double *paras;          // paras to store the parameters for each RACIPE model

  if (argc == 1) {
    printf("Missing input file!\n");
    exit(3); 
  }

  inputfile     = strdup(argv[1]);
  token         = strtok(argv[1], ".");
  modelname     = strdup(token);
  token         = strtok(NULL, ".");
  fileextension = strdup(token);

  if (strcmp(fileextension, "topo") == 0) {
    initial_simuopts (fileextension, argc, argv, &simu_opts);
    Model_generate (inputfile, &topoinfo, &simu_opts);
  }
  else if ((strcmp(fileextension, "cfg") == 0)) {
    read_cfg(inputfile, &topoinfo, &simu_opts);
    initial_simuopts (fileextension, argc, argv, &simu_opts);
  }
  else {
    printf("Input file is not recognized. It should be .topo or .cfg file!\n");
    exit(4);
  }

  cnt_store       = (int *)    calloc(simu_opts.num_stability,                sizeof(int));
  y_store         = (double *) calloc(simu_opts.num_ode*topoinfo.numG,        sizeof(double));
  soln            = (double *) calloc(simu_opts.num_stability*topoinfo.numG,  sizeof(double));
  paras           = (double *) calloc(3*topoinfo.numR+2*topoinfo.numG,        sizeof(double));
  
  // Seed the random value generator
  pcg32_srandom(time(NULL) ^ (intptr_t)&printf, 54u);
  
  // RACIPE main function
  for (i = 0; i < simu_opts.num_paras; i++){
    cnt_loop = 0;

    set_parameters(paras, &simu_opts, &topoinfo);
    
    for (j = 0; j < simu_opts.num_ode; j++){
        tmp_loop = solve_ODE(y_store, j, paras, &simu_opts, &topoinfo);
        cnt_loop = cnt_loop + tmp_loop;
    }

    cnt_loop = cnt_loop/simu_opts.num_ode;

    cnt  = count_state (y_store, soln, &simu_opts, &topoinfo);
    
    cnt_store[cnt-1] = cnt_store[cnt-1] + 1;

    save_model_paras(inputfile, &simu_opts, &topoinfo, i, cnt_loop, cnt, paras);
    save_model_solns(inputfile, &simu_opts, &topoinfo, i, cnt_loop, cnt, soln);
    T_test          (inputfile, &simu_opts, &topoinfo, i, cnt, soln, paras);
  }

  // Screen printout
  printf("\n-----------------Stability-----------------\n");
  printf("#states  Count\n");
  for (i = 0; i < simu_opts.num_stability; i++){
     printf("%d\t%d\n", i+1, cnt_store[i]);
  }

// Release the memory  ---- 
  release_memory(fileextension, &topoinfo, cnt_store, y_store, soln, paras);

  end = clock();
  time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
  printf("---> Run time : %f seconds ( %f hours)\n", time_spent, time_spent/3600.0);

  return 0;
}


