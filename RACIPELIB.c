# include <stdlib.h>
# include <stdint.h>
# include <stdio.h>
# include <math.h>
# include <time.h>
# include <string.h>

# include "RACIPELIB.h"
# include "pcg_basic.h"

/*********Shared Functions*********/
// Generate random value in (minV, maxV) following uniform distribution
double randu(double minV, double maxV)
{
  double u;
  
  do {
    u = ((double)pcg32_random()/UINT32_MAX);
  } while (u==0);

  return (minV + (maxV - minV)*u); 
}

// Generate random value following standard guassian distribtuion
double randg(double m, double stdvalue)
{
  double u1 = 0.0, u2 = 0.0, rsq = 0.0, fac = 0.0;
  static double z1 = 0.0, z2 = 0.0;   
  static int call = 0;

  if (call == 1){
    call = !call;
    return z2*stdvalue + m;
  }

  do {
    u1 = 2.0*randu(0, 1) - 1;
    u2 = 2.0*randu(0, 1) - 1;
    rsq = u1*u1 + u2*u2;
  } while (rsq >= 1.0 || rsq == 0.0);
  
  fac = sqrt(-2.0*log(rsq)/rsq);
  
  z1 = u1*fac;
  z2 = u2*fac;

  call = !call; 

  return z1*stdvalue + m;  
}

// Generate random value following non-negative standard guassian distribtuion
double randpg(double m, double stdvalue)
{
  double u = 0.0;
  do {
    u = randg(m, stdvalue);
  } while (u <= 0);

  return u;
}

// Generate random value following exponential distribution
double randexp(double m)
{
  double z = 0.0;
  double exp_value = 0.0;

  z = randu(0, 1);

  exp_value = -m * log(z); // Compute exponential random variable using inversion method

  return(exp_value);
}

// Generate random value for fold change parameter which need to be greater than or equal to 1
double randfd(double m, double stdvalue, int dist)
{
  double u = 0.0;

  switch (dist) {
    case 2:  // Guassian distribution
      do {
        u = randg(m, stdvalue);
      } while (u < 1);
      break;
    case 3:  // Exponential distribution
      do {
        u = randexp(m);
      } while (u < 1);
      break;
    }

    return u;
}

// Generate discrete integer between 1 and 6 following uniform distribution
double randd(int dist)
{
  double u = 0.0;
  double z = 0.0;
  int cnt  = 0;

  do {
    if (dist == 1) {
      u = randu(0, 7);
    }
    else if (dist == 2){
      u = randpg(3.5, 3);
    }
    else if (dist == 3){
      u = randexp(3.5);
    }

    if (u >= 0.5 && u < 1.5){
        z = 1.0;
        cnt = 1;
    }
    else if ( u >= 1.5 && u < 2.5){
        z = 2.0;
        cnt = 1;
    }
    else if ( u >= 2.5 && u < 3.5){
        z = 3.0;
        cnt = 1;
    }
    else if ( u >= 3.5 && u < 4.5){
        z = 4.0;
        cnt = 1;
    } 
    else if ( u >= 4.5 && u < 5.5){
        z = 5.0;
        cnt = 1;
    }
    else if ( u >= 5.5 && u < 6.5){
        z = 6.0;
        cnt = 1;
    }
  } while (u == 6.5 || cnt == 0);

  return z; 
}

// Shifted Hill Function
double Hillshift (double x, double x0, double nx, double lamda)
{
  double out;
  out = lamda + (1.0 - lamda) * (1.0/(1.0 + pow((x/x0), nx)));
  return out;
}

// Find median of an array
double median(double *x, int n)
{
  double  tmp;
  int     i = 0;
  int     j = 0;

  // Sorting
  for (i = 0; i < n-1; i++) {
    for (j = i+1; j < n; j++){
      if (x[j] < x[i]) {
        tmp  = x[i];
        x[i] = x[j];
        x[j] = tmp;
      }
    }
  }
  
  if (n%2==0) {
    return ((x[n/2] + x[n/2 - 1]) / 2.0);
  }
  else {
    return x[n/2];
  }
}

double sumdelta (double *y, double *ytmp, int NEQN)
{
    int      i = 0;
    double out = 0.0;
    
    for (i = 0; i < NEQN; i++){
        out = out + (y[i] - ytmp[i])*(y[i] - ytmp[i]);
    }

    return sqrt(out);
}

/*********Model_generate Functions*********/
void initial_simuopts (char *fileextension, int argc, char **argv, struct opts *simu_opts)
{

  int i = 0;

  // Default setting
  if (strcmp(fileextension, "topo") == 0){
    simu_opts->KDID          = 0;
    simu_opts->flag          = 0;
    simu_opts->distname      = strdup("Uniform");
    simu_opts->dist          = 1;
    simu_opts->SF            = 1.0;
    simu_opts->num_findT     = 1000;
    simu_opts->num_paras     = 100;
    simu_opts->num_ode       = 100; 
    simu_opts->num_stability = 10;
    simu_opts->thrd          = 1.0;
    simu_opts->Toggle_f_p    = 1;
    simu_opts->stepsize      = 0.1;
    simu_opts->maxiters      = 20;
    simu_opts->Toggle_T_test = 1;

    if (argc == 2){
      printf("Uniform distribution is used.\n");
    }
  }

  if (argc >= 3){
    if ((argc - 2)%2 != 0){
      if (argc == 3 && strcmp(argv[2], "-h") == 0){
        printf("Available options:\n");
        printf("-h             : Show all available options\n");
        printf("-dist          : Distribution used for randomization (only use with .topo file)\n");
        printf("                 1 ---> Uniform Distribution (Default)\n");
        printf("                 2 ---> Guassian Distribution\n");
        printf("                 3 ---> Exponential Distribution\n");
        printf("-SF            : Scale the distribution range, usually smaller than 1 (Default 1) (only use with .topo file)\n");
        printf("-num_findT     : Number of simulaitons used to estimate threshold (Default 1000) (only use with .topo file)\n");
        printf("-flag          : Only produce .cfg file or not (Default 0, not only produce .cfg file) (only use with .topo file);\n");
        printf("-num_paras     : Number of RACIPE models to generate (Default 100)\n");
        printf("-num_ode       : Number of Random initial values to solve ODEs (Default 100)\n");
        printf("-num_stability : Maximum number of stable states to save for one RACIPE model (Default 10)\n"); 
        printf("-thrd          : Cutoff for convergence of steady states for numerically solving ODEs (Default 1.0)\n");
        printf("-Toggle_f_p    : Save parameters of each RACIPE model or not (Default 1 (yes))\n");
        printf("-stepsize      : Stepsize for solving ODE (Default 0.1)\n");
        printf("-maxiters      : Maximum of Iteration for solving ODE at each RIVs (Default 20)\n");
        printf("-Toggle_T_test : Test threshold assumption or not (Default 1 (yes))\n");
        printf("-KDID          : Gene or link (See their ID in .cfg file) to be knockdown (Default 0)\n");
        exit(6);
      }
      else {
        printf("Input arguments are not enough!\n");
        exit(8);
      }
    }

    for (i = 2; i < argc; i = i + 2){

      if (strcmp(argv[i], "-dist") == 0){

        if (strcmp(fileextension, "cfg") == 0){
          printf("Please use .topo file instead of .cfg file.\n");
          exit(7);
        }

        switch (atoi(argv[i+1])) {
          case 1 :
            printf("Uniform distribution is used.\n");
            simu_opts->distname = strdup("Uniform");
            simu_opts->dist = 1;
            break;
          case 2 :
            printf("Guassion distribution is used.\n");
            simu_opts->distname = strdup("Guassian");
            simu_opts->dist = 2;
            break;
          case 3 :
            printf("Exponential distribution is used.\n");
            simu_opts->distname = strdup("Exponential");
            simu_opts->dist = 3;
            break;
          default:
            printf("The distribution you selected is not recognized! Please select the follwing distribution:\n");
            printf("1 ---> Uniform Distribution\n");
            printf("2 ---> Guassian Distribution\n");
            printf("3 ---> Exponential Distribution\n");
            exit(3);
        }
      }
      else if (strcmp(argv[i], "-SF") == 0){

        if (strcmp(fileextension, "cfg") == 0){
          printf("Please use .topo file instead of .cfg file.\n");
          exit(7);
        }

        simu_opts->SF = atof(argv[i+1]);
      }
      else if (strcmp(argv[i], "-num_findT") == 0){

        if (strcmp(fileextension, "cfg") == 0){
          printf("Please use .topo file instead of .cfg file.\n");
          exit(7);
        }

        simu_opts->num_findT = atoi(argv[i+1]);
      }
      else if (strcmp(argv[i], "-num_paras") == 0){
        simu_opts->num_paras = atoi(argv[i+1]);
      }
      else if (strcmp(argv[i], "-num_ode") == 0){
        simu_opts->num_ode = atoi(argv[i+1]);
      }
      else if (strcmp(argv[i], "-num_stability") == 0){
        simu_opts->num_stability = atoi(argv[i+1]);
      }
      else if (strcmp(argv[i], "-thrd") == 0){
        simu_opts->thrd = atof(argv[i+1]);
      }
      else if (strcmp(argv[i], "-Toggle_f_p") == 0){
        simu_opts->Toggle_f_p = atoi(argv[i+1]);
      }
      else if (strcmp(argv[i], "-stepsize") == 0){
        simu_opts->stepsize = atof(argv[i+1]);
      }
      else if (strcmp(argv[i], "-maxiters") == 0){
        simu_opts->maxiters = atoi(argv[i+1]);
      }
      else if (strcmp(argv[i], "-Toggle_T_test") == 0){
        simu_opts->Toggle_T_test = atoi(argv[i+1]);
      }
      else if (strcmp(argv[i], "-KDID") == 0){
        simu_opts->KDID = atoi(argv[i+1]);
      }
      else if (strcmp(argv[i], "-flag") == 0){

        if (strcmp(fileextension, "cfg") == 0){
          printf(".cfg file is already produced.\n");
          exit(7);
        }

        simu_opts->flag = atoi(argv[i+1]);
      }
      else{
        printf("Can not recognize the input arguments. Please use one of the follwing:\n");
        printf("-h             : Show all available options\n");
        printf("-dist          : Distribution used for randomization (only use with .topo file)\n");
        printf("                 1 ---> Uniform Distribution (Default)\n");
        printf("                 2 ---> Guassian Distribution\n");
        printf("                 3 ---> Exponential Distribution\n");
        printf("-SF            : Scale the distribution range, usually smaller than 1 (Default 1) (only use with .topo file)\n");
        printf("-num_findT     : Number of simulaitons used to estimate threshold (Default 1000) (only use with .topo file)\n");
        printf("-flag          : Only produce .cfg file or not (Default 0, not only produce .cfg file) (only use with .topo file);\n");
        printf("-num_paras     : Number of RACIPE models to generate (Default 100)\n");
        printf("-num_ode       : Number of Random initial values to solve ODEs (Default 100)\n");
        printf("-num_stability : Maximum number of stable states to save for one RACIPE model (Default 10)\n"); 
        printf("-thrd          : Cutoff for convergence of steady states for numerically solving ODEs (Default 1.0)\n");
        printf("-Toggle_f_p    : Save parameters of each RACIPE model or not (Default 1 (yes))\n");
        printf("-stepsize      : Stepsize for solving ODE (Default 0.1)\n");
        printf("-maxiters      : Maximum of Iteration for solving ODE at each RIVs (Default 20)\n");
        printf("-Toggle_T_test : Test threshold assumption or not (Default 1 (yes))\n");
        printf("-KDID          : Gene or link (See their ID in .cfg file) to be knockdown (Default 0)\n");
        exit(6);
      }
    }
  }
}

void Model_generate (char *topofile, struct topo *topoinfo, struct opts *simu_opts)
{

  // Section 1 parameters
  FILE *fin;
  FILE *fout;
  char *modelname;
  char configname[100] = "";
  
  // Section 2 parameters
  FILE *fparas; 
  char parasname[100] = "";
  
  // Section 1 --- Printout the number of genes involved and transform the topo file
  fin  = fopen(topofile, "r");

  if (fin == NULL){
    printf("No topology files are found!\n");
    exit(2);
  }

  modelname = strtok(topofile, ".");
  strcpy(configname, modelname);
  strcat(configname, ".cfg");
  fout = fopen(configname, "w");

  // Parameters for the simulations, stored in .cfg file
  fprintf(fout, "Distribution\t%d\t%s\t%f\n",              simu_opts->dist, simu_opts->distname, simu_opts->SF);
  fprintf(fout, "NumberOfSimulationToFindThreshold\t%d\n", simu_opts->num_findT);
  fprintf(fout, "NumberOfRACIPEModels\t%d\n",              simu_opts->num_paras);
  fprintf(fout, "NumberOfRIVs\t%d\n",                      simu_opts->num_ode);
  fprintf(fout, "NumberOfStatesToStore\t%d\n",             simu_opts->num_stability);
  fprintf(fout, "ThresholdForConvergence\t%f\n",           simu_opts->thrd);
  fprintf(fout, "ToggleOfSaveParameters\t%d\n",            simu_opts->Toggle_f_p);
  fprintf(fout, "Stepsize\t%f\n",                          simu_opts->stepsize);
  fprintf(fout, "MaximumOfIterations\t%d\n",               simu_opts->maxiters);
  fprintf(fout, "TestThreshold\t%d\n",                     simu_opts->Toggle_T_test);

  check_topo(fin, fout, topoinfo);
  fclose(fin);
  fclose(fout); 
  
  // S2 --- Generate the randomization range
  strcpy(parasname, modelname);
  strcat(parasname, ".prs");
  fparas = fopen(parasname, "w+"); 
  generate_random_range(fparas, topoinfo, simu_opts->num_findT, simu_opts->dist, simu_opts->SF);
  fclose(fparas);

  fin    = NULL;
  fout   = NULL;
  fparas = NULL;

  if (simu_opts->flag == 1){
    exit(10);
  }

}

void check_topo(FILE *f_in, FILE *f_out, struct topo *topoinfo)
{ 
  int  tmpnumG  = 0;
  int  tmpnumR  = 0; 
  int  RT       = 0;
  int  count    = 0;
  int  i        = 0;
  int  j        = 0;
  int  cntP     = 0;
  int  matchnum = 0;
  char G1[100]  = "";
  char G2[100]  = "";

  topoinfo->Gname = (char**)malloc(sizeof(char *));

  if (f_in == NULL){
    printf("Topology file is missing!\n");
    exit(1);
  }

  rewind(f_in);
  
  fscanf(f_in, "%*[^\n]\n", NULL); //skip first line of topo file
  
  while (fscanf(f_in, "%s\t%s\t%d\n", G1, G2, &RT) == 3) {
    tmpnumR++;
                
    // Check the types of regualtions: 1 --> Activation; 2 --> Inhibition;                
    if (RT != 1 && RT != 2){
      printf("ERROR: Number %d regulation is not recognized\n", tmpnumR);
      exit(1);
    } 
    
    if (count == 0) { // second line of topo file
      
      if (strcmp(G1, G2) == 0) {
        count = count + 1;
        topoinfo->Gname = (char**)realloc(topoinfo->Gname, (count)*sizeof(char *));
        topoinfo->Gname[count-1] = (char*)malloc(sizeof(G1));
        strcpy(topoinfo->Gname[count-1], G1);
        tmpnumG = tmpnumG + 1;
      }
      else {
        count = count + 2;
        topoinfo->Gname = (char**)realloc(topoinfo->Gname, (count)*sizeof(char *));
        topoinfo->Gname[count-2] = (char*)malloc(sizeof(G1));
        topoinfo->Gname[count-1] = (char*)malloc(sizeof(G2));
        strcpy(topoinfo->Gname[count-2], G1);
        strcpy(topoinfo->Gname[count-1], G2);
        tmpnumG = tmpnumG + 2;
      }
    }
    else {
      matchnum = 0;
      // compare G1 first
      for (i = 0; i < count; i++) {
        if (strcmp(G1, topoinfo->Gname[i]) == 0) {
          matchnum = matchnum + 1;
        }
      }
      
      if (matchnum == 0) {
        count = count + 1;
        topoinfo->Gname = (char**)realloc(topoinfo->Gname, (count)*sizeof(char *));
        topoinfo->Gname[count-1] = (char*)malloc(sizeof(G1));
        strcpy(topoinfo->Gname[count-1], G1);
        tmpnumG = tmpnumG + 1;
      }
      
      if (strcmp(G1, G2) != 0) {
        matchnum = 0;
        for (i = 0; i < count; i++) {
          if (strcmp(G2, topoinfo->Gname[i]) == 0) {
            matchnum = matchnum + 1;
          }
        }
      
        if (matchnum == 0) {
          count = count + 1;
          topoinfo->Gname = (char**)realloc(topoinfo->Gname, (count)*sizeof(char *));
          topoinfo->Gname[count-1] = (char*)malloc(sizeof(G2));
          strcpy(topoinfo->Gname[count-1], G2);
          tmpnumG = tmpnumG + 1;
        }
      }
    }

  }
  
  printf("\n-------------------------------------------\n");
  printf("The topo file contains the following genes:\n");
  printf("Gene_ID\tGene_Name\n");
  for (i = 0; i < count; i++) {
    printf("%d -- %s\n", i, topoinfo->Gname[i]);
  } 
  printf("-------------------------------------------\n");
  
  topoinfo->numR = tmpnumR;
  topoinfo->numG = tmpnumG;
  
  printf("Total number of genes = %d\n",       tmpnumG);
  printf("Total number of regulations = %d\n", tmpnumR);

  printf("-------------------------------------------\n");
  
  // transformation of the topo file  
  topoinfo->SourceG   = (int *)calloc(tmpnumR, sizeof(int));
  topoinfo->TargetG   = (int *)calloc(tmpnumR, sizeof(int));
  topoinfo->TypeR     = (int *)calloc(tmpnumR, sizeof(int));
  topoinfo->ParasPos  = (int *)calloc(tmpnumR, sizeof(int));
  topoinfo->numover   = (int *)calloc(tmpnumG, sizeof(int));
  topoinfo->numdown   = (int *)calloc(tmpnumG, sizeof(int));
  count   = 0;
  
  rewind(f_in);
  fscanf(f_in, "%*[^\n]\n", NULL); //skip first line of topo file
  while (fscanf(f_in, "%s\t%s\t%d\n", G1, G2, &RT) == 3) {
    
    count++;
    
    for (i = 0; i < tmpnumG; i++) {
      if (strcmp(G1, topoinfo->Gname[i]) == 0) {
        topoinfo->SourceG[count-1] = i;
      }
      
      if (strcmp(G2, topoinfo->Gname[i]) == 0) {
        topoinfo->TargetG[count-1] = i;
      }
    }
    
    topoinfo->TypeR[count-1] = RT;
  }
  
  fprintf(f_out, "NumberOfRegulations\t%d\n",       topoinfo->numR);
  fprintf(f_out, "NumberOfGenes\t%d\n",             topoinfo->numG);
    
  for (i = 0; i < topoinfo->numG; i++) {
    fprintf(f_out, "%d\t%s\n", i+1, topoinfo->Gname[i]);
  } 

  for (i = 0; i < topoinfo->numR; i++){
    fprintf(f_out, "%d\t%d\t%d\t%d\n", topoinfo->numG+i+1, topoinfo->SourceG[i], topoinfo->TargetG[i], topoinfo->TypeR[i]);
  }

  for(i = 0; i < topoinfo->numG; i++){
    for (j = 0; j < topoinfo->numR; j++){
      if (topoinfo->TargetG[j] == i){ 
        topoinfo->ParasPos[j] = 3*cntP + 2*topoinfo->numG;  // Position of threshold parameters for each regulation.
        cntP = cntP + 1;
      }
    }
  }

  // for (i = 0; i < topoinfo->numR; i++){
  //   printf("%d\n", topoinfo->ParasPos[i]);
  // }
}

void generate_random_range(FILE *f_paras, struct topo *topoinfo, int num, int dist, double SF)
{
  char   tmpparasname[100] = "";
  double tmpmin = 0.0;
  double tmpmax = 0.0;
  int    typeR  = 0;

  /**** Default ranges for each class of parameters ****/
  double minP_d = 1.0;
  double maxP_d = 100.0;
  double minK_d = 0.1;
  double maxK_d = 1.0;
  double minN_d = 1.0;
  double maxN_d = 6.0;
  double minF_d = 1.0;
  double maxF_d = 100.0;
  /**** Default ranges for each class of parameters ****/

  double minP = 0.0;
  double maxP = 0.0;
  double minK = 0.0;
  double maxK = 0.0;
  double minN = 0.0;
  double maxN = 0.0;
  double minF = 0.0;
  double maxF = 0.0;

  double meanP = 0.0;
  double stdP  = 0.0;
  double meanK = 0.0;
  double stdK  = 0.0;
  double meanN = 0.0;
  double stdN  = 0.0;
  double meanF = 0.0;
  double stdF  = 0.0;

  int i = 0;
  int j = 0;

  double *minT;
  double *maxT;
  minT = (double *)calloc(topoinfo->numG, sizeof(double));
  maxT = (double *)calloc(topoinfo->numG, sizeof(double));

  double *meanT;
  double *stdT;
  meanT = (double *)calloc(topoinfo->numG, sizeof(double));
  stdT  = (double *)calloc(topoinfo->numG, sizeof(double));


  switch (dist) {
    case 1 :
      minP = ((minP_d + maxP_d) - SF * (maxP_d - minP_d))/2.0;
      maxP = ((minP_d + maxP_d) + SF * (maxP_d - minP_d))/2.0;
      minK = ((minK_d + maxK_d) - SF * (maxK_d - minK_d))/2.0;
      maxK = ((minK_d + maxK_d) + SF * (maxK_d - minK_d))/2.0;
      minN = minN_d;
      maxN = maxN_d;
      minF = ((minF_d + maxF_d) - SF * (maxF_d - minF_d))/2.0;
      maxF = ((minF_d + maxF_d) + SF * (maxF_d - minF_d))/2.0;

      fprintf(f_paras, "Parameter\tMinimum_value\tMaximum_Value\tRegulation_type\n");

      // Format of parameter files: name minV maxV
      // production rate  
      for (i = 0; i < topoinfo->numG; i++){
        fprintf(f_paras, "Prod_of_%s\t%f\t%f\t%d\n", topoinfo->Gname[i], minP, maxP, 0);
      }
      
      // degradation rate
      for (i = 0; i < topoinfo->numG; i++){
        fprintf(f_paras, "Deg_of_%s\t%f\t%f\t%d\n", topoinfo->Gname[i], minK, maxK, 0);
      }

      for(i = 0; i < topoinfo->numG; i++){
        // estimate threshold and printout
        estimate_threshold(num, i, minP, maxP, minK, maxK, minN, maxN, minF, maxF, minT, maxT, topoinfo, dist, SF);
        // printf("%s = %f\n", topoinfo->Gname[i], (minT+maxT)/2.0);
      }

      for(i = 0; i < topoinfo->numG; i++){
        for (j = 0; j < topoinfo->numR; j++){
          if (topoinfo->TargetG[j] == i){ 
            // Threshold
            fprintf(f_paras,   "Trd_of_%sTo%s\t%f\t%f\t%d\n", topoinfo->Gname[topoinfo->SourceG[j]], topoinfo->Gname[topoinfo->TargetG[j]], minT[topoinfo->SourceG[j]], maxT[topoinfo->SourceG[j]], 0);
            // Number of binding sites
            fprintf(f_paras,   "Num_of_%sTo%s\t%f\t%f\t%d\n", topoinfo->Gname[topoinfo->SourceG[j]], topoinfo->Gname[topoinfo->TargetG[j]], minN, maxN, 0);
            // Fold change of a regulation
            if (topoinfo->TypeR[j] == 1) {  //Activation
              fprintf(f_paras, "Act_of_%sTo%s\t%f\t%f\t%d\n", topoinfo->Gname[topoinfo->SourceG[j]], topoinfo->Gname[topoinfo->TargetG[j]], minF, maxF, topoinfo->TypeR[j]);
            }
            else if (topoinfo->TypeR[j] == 2) { //Inhibition
              fprintf(f_paras, "Inh_of_%sTo%s\t%f\t%f\t%d\n", topoinfo->Gname[topoinfo->SourceG[j]], topoinfo->Gname[topoinfo->TargetG[j]], minF, maxF, topoinfo->TypeR[j]);
            }
          }
        }
      }

      free(minT);
      free(maxT);
      free(meanT);
      free(stdT);
      break;
    case 2 :
      // Format of parameter files: name meanV stdV
      meanP = (minP_d+maxP_d)/2.0;
      stdP  = (maxP_d-minP_d)*SF/2.0;
      meanK = (minK_d+maxK_d)/2.0;
      stdK  = (maxK_d-minK_d)*SF/2.0;
      meanN = 3.5;
      stdN  = 3;
      meanF = (minF_d+maxF_d)/2.0;
      stdF  = (maxF_d-minF_d)*SF/2.0;

      fprintf(f_paras, "Parameter\tMean\tStandard_deviation\tRegulation_type\n");

      // production rate  
      for (i = 0; i < topoinfo->numG; i++){
        fprintf(f_paras, "Prod_of_%s\t%f\t%f\t%d\n", topoinfo->Gname[i], meanP, stdP, 0);
      }
      
      // degradation rate
      for (i = 0; i < topoinfo->numG; i++){
        fprintf(f_paras, "Deg_of_%s\t%f\t%f\t%d\n", topoinfo->Gname[i], meanK, stdK, 0);
      }

      for(i = 0; i < topoinfo->numG; i++){
        // estimate threshold and printout
        estimate_threshold(num, i, meanP, stdP, meanK, stdK, meanN, stdN, meanF, stdF, meanT, stdT, topoinfo, dist, SF);
        // printf("%s = %f\n", topoinfo->Gname[i], (minT+maxT)/2.0);
      }

      for(i = 0; i < topoinfo->numG; i++){
        for (j = 0; j < topoinfo->numR; j++){
          if (topoinfo->TargetG[j] == i){ 
            // Threshold
            fprintf(f_paras,   "Trd_of_%sTo%s\t%f\t%f\t%d\n", topoinfo->Gname[topoinfo->SourceG[j]], topoinfo->Gname[topoinfo->TargetG[j]], meanT[topoinfo->SourceG[j]], stdT[topoinfo->SourceG[j]], 0);
            // Number of binding sites
            fprintf(f_paras,   "Num_of_%sTo%s\t%f\t%f\t%d\n", topoinfo->Gname[topoinfo->SourceG[j]], topoinfo->Gname[topoinfo->TargetG[j]], meanN, stdN, 0);
            // Fold change of a regulation
            if (topoinfo->TypeR[j] == 1) {  //Activation
              fprintf(f_paras, "Act_of_%sTo%s\t%f\t%f\t%d\n", topoinfo->Gname[topoinfo->SourceG[j]], topoinfo->Gname[topoinfo->TargetG[j]], meanF, stdF, topoinfo->TypeR[j]);
            }
            else if (topoinfo->TypeR[j] == 2) { //Inhibition
              fprintf(f_paras, "Inh_of_%sTo%s\t%f\t%f\t%d\n", topoinfo->Gname[topoinfo->SourceG[j]], topoinfo->Gname[topoinfo->TargetG[j]], meanF, stdF, topoinfo->TypeR[j]);
            }
          }
        }
      }

      free(minT);
      free(maxT);
      free(meanT);
      free(stdT);
      break;
    case 3 :
      // Format of parameter files: name meanV stdV
      meanP = (minP_d+maxP_d)*SF/2.0;
      meanK = (minK_d+maxK_d)*SF/2.0;
      meanN = (minN_d+maxN_d)/2.0;
      meanF = (minF_d+maxF_d)*SF/2.0;

      fprintf(f_paras, "Parameter\tMean\tNo_sense\tRegulation_type\n");

      // production rate  
      for (i = 0; i < topoinfo->numG; i++){
        fprintf(f_paras, "Prod_of_%s\t%f\t%f\t%d\n", topoinfo->Gname[i], meanP, stdP, 0);
      }
      
      // degradation rate
      for (i = 0; i < topoinfo->numG; i++){
        fprintf(f_paras, "Deg_of_%s\t%f\t%f\t%d\n", topoinfo->Gname[i], meanK, stdK, 0);
      }

      for(i = 0; i < topoinfo->numG; i++){
        // estimate threshold and printout
        estimate_threshold(num, i, meanP, stdP, meanK, stdK, meanN, stdN, meanF, stdF, meanT, stdT, topoinfo, dist, SF);
        // printf("%s = %f\n", topoinfo->Gname[i], (minT+maxT)/2.0);
      }

      for(i = 0; i < topoinfo->numG; i++){
        for (j = 0; j < topoinfo->numR; j++){
          if (topoinfo->TargetG[j] == i){ 
            // Threshold
            fprintf(f_paras,   "Trd_of_%sTo%s\t%f\t%f\t%d\n", topoinfo->Gname[topoinfo->SourceG[j]], topoinfo->Gname[topoinfo->TargetG[j]], meanT[topoinfo->SourceG[j]], stdT[topoinfo->SourceG[j]], 0);
            // Number of binding sites
            fprintf(f_paras,   "Num_of_%sTo%s\t%f\t%f\t%d\n", topoinfo->Gname[topoinfo->SourceG[j]], topoinfo->Gname[topoinfo->TargetG[j]], meanN, stdN, 0);
            // Fold change of a regulation
            if (topoinfo->TypeR[j] == 1) {  //Activation
              fprintf(f_paras, "Act_of_%sTo%s\t%f\t%f\t%d\n", topoinfo->Gname[topoinfo->SourceG[j]], topoinfo->Gname[topoinfo->TargetG[j]], meanF, stdF, topoinfo->TypeR[j]);
            }
            else if (topoinfo->TypeR[j] == 2) { //Inhibition
              fprintf(f_paras, "Inh_of_%sTo%s\t%f\t%f\t%d\n", topoinfo->Gname[topoinfo->SourceG[j]], topoinfo->Gname[topoinfo->TargetG[j]], meanF, stdF, topoinfo->TypeR[j]);
            }
          }
        }
      }

      free(minT);
      free(maxT);
      free(meanT);
      free(stdT);
      break;
  }

  topoinfo->prsrandrange    = (double **)calloc(3,                                        sizeof(double *));
  topoinfo->prsrandrange[0] = (double *) calloc(3*topoinfo->numR+2*topoinfo->numG,        sizeof(double));
  topoinfo->prsrandrange[1] = (double *) calloc(3*topoinfo->numR+2*topoinfo->numG,        sizeof(double));
  topoinfo->prsrandrange[2] = (double *) calloc(3*topoinfo->numR+2*topoinfo->numG,        sizeof(double));

  rewind(f_paras);
  fscanf(f_paras, "%*[^\n]\n", NULL);
  i = 0;
  while (fscanf(f_paras, "%s\t%lf\t%lf\t%d", tmpparasname, &tmpmin, &tmpmax, &typeR) == 4){
    topoinfo->prsrandrange[0][i] = tmpmin;
    topoinfo->prsrandrange[1][i] = tmpmax;
    topoinfo->prsrandrange[2][i] = (double)typeR;
    // printf("%f\t%f\t%f\n", topoinfo->prsrandrange[0][i], topoinfo->prsrandrange[1][i], topoinfo->prsrandrange[2][i]);
    i++;
  }
}

void estimate_threshold(int num, int ID, double minP, double maxP, double minK, double maxK, double minN, double maxN, double minF, double maxF, double *minT, double *maxT, struct topo *topoinfo, int dist, double SF)
{
  int    i    = 0;
  int    j    = 0;
  int    numA = 0;
  int    numI = 0;

  double g      = 0.0;
  double k      = 0.0;
  double T      = 0.0;
  double n      = 0.0;
  double lambda = 0.0;
  double MA     = 0.0;
  double MB     = 0.0;
  double f1     = 0.0;
  double f2     = 0.0;

  double *A; // A is a standalone gene
  double *B;

  A  = (double *)calloc(num, sizeof(double));
  B  = (double *)calloc(num, sizeof(double));

  for (i = 0; i < topoinfo->numR; i++){
    if (topoinfo->TargetG[i] == ID){
      if (topoinfo->TypeR[i] == 1){ //Activation
        numA = numA + 1;        
      }
      else if (topoinfo->TypeR[i] == 2){ //Inhibition
        numI = numI + 1;
      }
    }
  }

  f1 = (2.0 - SF*1.96)/2.0;
  f2 = (2.0 + SF*1.96)/2.0;

  switch (dist){
    case 1:    // Uniform distribution
      for (i = 0; i < num; i++){
        g    = randu(minP, maxP);
        k    = randu(minK, maxK);
        A[i] = g/k;
      }   

      MA = median(A, num); 

      for (i = 0; i < num; i++){
        g    = randu(minP, maxP);
        k    = randu(minK, maxK);
        B[i] = g/k;
        
        if (numA != 0){
          for (j = 0; j < numA; j++){
            g      = randu(minP, maxP);
            k      = randu(minK, maxK);
            n      = randd(dist);
            T      = randu(MA*f1, MA*f2);
            lambda = randu(minF, maxF);

            B[i] = B[i]*Hillshift(g/k, T, n, lambda)/lambda;
          }
        }
        
        if (numI != 0){
          for (j = 0; j < numI; j++){
            g      = randu(minP, maxP);
            k      = randu(minK, maxK);
            n      = randd(dist);
            T      = randu(MA*f1, MA*f2);
            lambda = 1.0/randu(minF, maxF);

            B[i] = B[i]*Hillshift(g/k, T, n, lambda);
          }
        }
      }

      MB = median(B, num);
      minT[ID] = MB*f1;
      maxT[ID] = MB*f2;
      free(A);
      free(B);
      break;
    case 2:    // non-negative Guassian distribution
      for (i = 0; i < num; i++){
          g    = randpg(minP, maxP);
          k    = randpg(minK, maxK);
          A[i] = g/k;
        }   

        MA = median(A, num); 

        for (i = 0; i < num; i++){
          g    = randpg(minP, maxP);
          k    = randpg(minK, maxK);
          B[i] = g/k;
          
          if (numA != 0){
            for (j = 0; j < numA; j++){
              g      = randpg(minP, maxP);
              k      = randpg(minK, maxK);
              n      = randd(dist);
              T      = randpg(MA, (MA*f2-MA*f1)/2.0);
              lambda = randfd(minF, maxF, dist);

              B[i] = B[i]*Hillshift(g/k, T, n, lambda)/lambda;
            }
          }
          
          if (numI != 0){
            for (j = 0; j < numI; j++){
              g      = randpg(minP, maxP);
              k      = randpg(minK, maxK);
              n      = randd(dist);
              T      = randpg(MA, (MA*f2-MA*f1)/2.0);
              lambda = 1.0/randfd(minF, maxF, dist);

              B[i] = B[i]*Hillshift(g/k, T, n, lambda);
            }
          }
        }

        MB = median(B, num);
        minT[ID] = MB;
        maxT[ID] = (MB*f2-MB*f1)/2.0;
        free(A);
        free(B);
        break;
    case 3:    // Exponential distribution
      for (i = 0; i < num; i++){
        g    = randexp(minP);
        k    = randexp(minK);
        A[i] = g/k;
      }   

      MA = median(A, num); 

      for (i = 0; i < num; i++){
        g    = randexp(minP);
        k    = randexp(minK);
        B[i] = g/k;
        
        if (numA != 0){
          for (j = 0; j < numA; j++){
            g      = randexp(minP);
            k      = randexp(minK);
            n      = randd(dist);
            T      = randexp(MA);
            lambda = randfd(minF, 0, dist);

            B[i] = B[i]*Hillshift(g/k, T, n, lambda)/lambda;
          }
        }
        
        if (numI != 0){
          for (j = 0; j < numI; j++){
            g      = randexp(minP);
            k      = randexp(minK);
            n      = randd(dist);
            T      = randexp(MA);
            lambda = 1.0/randfd(minF, 0, dist);

            B[i] = B[i]*Hillshift(g/k, T, n, lambda);
          }
        }
      }

      MB = median(B, num);
      minT[ID] = MB;
      maxT[ID] = 0.0;
      free(A);
      free(B);
      break;
    }
}

/*********RACIPE Functions*********/
void read_cfg(char *inputfile, struct topo *topoinfo, struct opts *simu_opts)
{
  int tmpID = 0;
  int cntP  = 0;
  int i     = 0;
  int j     = 0;

  FILE   *fprs;
  char   *modelname;
  char   configname[100] = "";
  char   prsname[100] = "";
  double tmpmin = 0.0;
  double tmpmax = 0.0;
  int    typeR  = 0;

  modelname = strtok(inputfile, ".");
  strcpy(configname, modelname);
  strcat(configname, ".cfg");

  FILE *fcfg;
  char tmpparasname[1000]  = "";
  char tmpparasname2[1000]  = "";

  fcfg = fopen(configname, "r");

  if (fcfg == NULL){
    printf("No configure file provided!\n");
    exit(2);
  }

  rewind(fcfg);

  fscanf(fcfg, "%s\t%d\t%s\t%lf\n",  tmpparasname,  &simu_opts->dist, tmpparasname2, &simu_opts->SF);
  fscanf(fcfg, "%s\t%d\n",           tmpparasname,  &simu_opts->num_findT);
  fscanf(fcfg, "%s\t%d\n",           tmpparasname,  &simu_opts->num_paras);
  fscanf(fcfg, "%s\t%d\n",           tmpparasname,  &simu_opts->num_ode);
  fscanf(fcfg, "%s\t%d\n",           tmpparasname,  &simu_opts->num_stability);
  fscanf(fcfg, "%s\t%lf\n",          tmpparasname,  &simu_opts->thrd);
  fscanf(fcfg, "%s\t%d\n",           tmpparasname,  &simu_opts->Toggle_f_p);
  fscanf(fcfg, "%s\t%lf\n",          tmpparasname,  &simu_opts->stepsize);
  fscanf(fcfg, "%s\t%d\n",           tmpparasname,  &simu_opts->maxiters);
  fscanf(fcfg, "%s\t%d\n",           tmpparasname,  &simu_opts->Toggle_T_test);
  fscanf(fcfg, "%s\t%d\n",           tmpparasname,  &topoinfo->numR);
  fscanf(fcfg, "%s\t%d\n",           tmpparasname,  &topoinfo->numG);

  simu_opts->distname = strdup(tmpparasname2);
  simu_opts->KDID     = 0;
  simu_opts->flag     = 0;

  switch (simu_opts->dist) {
    case 1 :
      printf("Uniform distribution is used.\n");
      break;
    case 2 :
      printf("Guassion distribution is used.\n");
      break;
    case 3 :
      printf("Exponential distribution is used.\n");
      break;
    default:
      printf("The distribution you used is not recognized! Please check .cfg file:\n");
      printf("1 ---> Uniform Distribution\n");
      printf("2 ---> Guassian Distribution\n");
      printf("3 ---> Exponential Distribution\n");
      exit(3);
  }

  // printf("%s\t%d\t%s\t%f\n",   tmpparasname,  simu_opts->dist,  simu_opts->distname, simu_opts->SF);
  // printf("%s\t%d\n",           tmpparasname,  simu_opts->num_findT);
  // printf("%s\t%d\n",           tmpparasname,  simu_opts->num_paras);
  // printf("%s\t%d\n",           tmpparasname,  simu_opts->num_ode);
  // printf("%s\t%d\n",           tmpparasname,  simu_opts->num_stability);
  // printf("%s\t%lf\n",          tmpparasname,  simu_opts->thrd);
  // printf("%s\t%d\n",           tmpparasname,  simu_opts->Toggle_f_p);
  // printf("%s\t%lf\n",          tmpparasname,  simu_opts->stepsize);
  // printf("%s\t%d\n",           tmpparasname,  simu_opts->maxiters);
  // printf("%s\t%d\n",           tmpparasname,  simu_opts->Toggle_T_test);
  // printf("%s\t%d\n",           tmpparasname,  topoinfo->numR);
  // printf("%s\t%d\n",           tmpparasname,  topoinfo->numG);

  topoinfo->SourceG         = (int *)    calloc(topoinfo->numR,                           sizeof(int));
  topoinfo->TargetG         = (int *)    calloc(topoinfo->numR,                           sizeof(int));
  topoinfo->TypeR           = (int *)    calloc(topoinfo->numR,                           sizeof(int));
  topoinfo->ParasPos        = (int *)    calloc(topoinfo->numR,                           sizeof(int));
  topoinfo->numover         = (int *)    calloc(topoinfo->numG,                           sizeof(int));
  topoinfo->numdown         = (int *)    calloc(topoinfo->numG,                           sizeof(int));
  topoinfo->prsrandrange    = (double **)calloc(3,                                        sizeof(double *));
  topoinfo->prsrandrange[0] = (double *) calloc(3*topoinfo->numR+2*topoinfo->numG,        sizeof(double));
  topoinfo->prsrandrange[1] = (double *) calloc(3*topoinfo->numR+2*topoinfo->numG,        sizeof(double));
  topoinfo->prsrandrange[2] = (double *) calloc(3*topoinfo->numR+2*topoinfo->numG,        sizeof(double));

  for (i = 0; i < topoinfo->numG; i++) {
    fscanf(fcfg, "%*[^\n]\n", NULL);
  }

  for (i = 0; i < topoinfo->numR; i++) {
    fscanf(fcfg, "%d\t%d\t%d\t%d\n",  &tmpID, topoinfo->SourceG + i, topoinfo->TargetG + i, topoinfo->TypeR + i);
    // printf("%d\t%d\t%d\t%d\n",        tmpID, topoinfo->SourceG[i], topoinfo->TargetG[i], topoinfo->TypeR[i]);
  }

  for(i = 0; i < topoinfo->numG; i++){
    for (j = 0; j < topoinfo->numR; j++){
      if (topoinfo->TargetG[j] == i){ 
        topoinfo->ParasPos[j] = 3*cntP + 2*topoinfo->numG;  // Position of threshold parameters for each regulation.
        cntP = cntP + 1;
      }
    }
  }

  // for (i = 0; i < topoinfo->numR; i++){
  //   printf("%d\n", topoinfo->ParasPos[i]);
  // }

  strcpy(prsname, modelname);
  strcat(prsname, ".prs");

  fprs = fopen(prsname, "r");

  if (fprs == NULL){
    printf("No .prs file found!\n");
    exit(2);
  }

  rewind(fprs);
  fscanf(fprs, "%*[^\n]\n", NULL); //skip first line of topo file
  i = 0;
  while (fscanf(fprs, "%s\t%lf\t%lf\t%d", tmpparasname, &tmpmin, &tmpmax, &typeR) == 4){
    topoinfo->prsrandrange[0][i] = tmpmin;
    topoinfo->prsrandrange[1][i] = tmpmax;
    topoinfo->prsrandrange[2][i] = (double)typeR;
    // printf("%f\t%f\t%f\n", topoinfo->prsrandrange[0][i], topoinfo->prsrandrange[1][i], topoinfo->prsrandrange[2][i]);
    i++;
  }
}

void save_model_paras(char *inputfile, struct opts *simu_opts, struct topo *topoinfo, int num, int cnt_loop, int cnt, double *paras)
{
  static FILE *f_p = NULL;                  
  char fpname [100] = "";
  char *modelname;

  int i = 0;

  if (simu_opts->Toggle_f_p == 1) {
    if (f_p == NULL) {
      modelname = strtok(inputfile, ".");
      strcpy(fpname, modelname);
      strcat(fpname, "_parameters");
      strcat(fpname, ".dat");
      f_p   = fopen(fpname,"w");
    }

    fprintf(f_p, "%d\t%d\t%d", num+1, cnt_loop, cnt);
    for (i = 0; i < 3*topoinfo->numR+2*topoinfo->numG; i++){
      fprintf(f_p, "\t%f", paras[i]);
    }
    fprintf(f_p, "\n");

    if (simu_opts->num_paras == num + 1){
      fclose(f_p);
      f_p = NULL;
    }
  }
}

void save_model_solns(char *inputfile, struct opts *simu_opts, struct topo *topoinfo, int num, int cnt_loop, int cnt, double *soln)
{

  static FILE **f_s;

  if (num == 0) {
    f_s = (FILE **) calloc (simu_opts->num_stability, sizeof(FILE *));
  }
  
  char fsname[simu_opts->num_stability][100];
  char *modelname;
  char tmpparasname[100] = "";
  modelname = strtok(inputfile, ".");

  int i = 0;
  int h = 0;
  int h2= 0;

  if (f_s[cnt-1] == NULL) {
    sprintf(tmpparasname, "%d", cnt);
    strcpy (fsname[cnt-1], modelname);
    strcat (fsname[cnt-1], "_solution_");
    strcat (fsname[cnt-1], tmpparasname);
    strcat (fsname[cnt-1], ".dat");
    f_s[cnt-1] = fopen(fsname[cnt-1],"w");
  }

  fprintf(f_s[cnt-1], "%d\t%d\t%d", num+1, cnt_loop, cnt);
    for (h = 0; h < cnt; h++){
        h2 = 1;
        while (h2 <= topoinfo->numG) {
            fprintf(f_s[cnt-1], "\t%f", log2(soln[topoinfo->numG*h + h2 - 1]));
            h2++;
        }
    }
  fprintf(f_s[cnt-1], "\n");

  if (simu_opts->num_paras == num + 1) {
    for (i = 0; i < simu_opts->num_stability; i++){
      if (f_s[i] != NULL){
        fclose(f_s[i]);
      }
    }
    free(f_s);
  }
}

void set_parameters ( double *paras, struct opts *simu_opts, struct topo *topoinfo)
{ 

  int i = 0;

  switch (simu_opts->dist) {
    case 1:
      // Production rate
      for (i = 0; i < topoinfo->numG; i++){
          paras[i] = randu(topoinfo->prsrandrange[0][i], topoinfo->prsrandrange[1][i]);
      }

      // Degradation rate
      for (i = 0; i < topoinfo->numG; i++){
          paras[i + topoinfo->numG] = randu(topoinfo->prsrandrange[0][i + topoinfo->numG], topoinfo->prsrandrange[1][i + topoinfo->numG]);
      }

      // Threshold
      for (i = 0; i < topoinfo->numR; i++){
          paras[3*i + 2*topoinfo->numG] = randu(topoinfo->prsrandrange[0][3*i + 2*topoinfo->numG], topoinfo->prsrandrange[1][3*i + 2*topoinfo->numG]);
      }

      // Coefficient
      for (i = 0; i < topoinfo->numR; i++){
          paras[3*i + 1 + 2*topoinfo->numG] = randd(simu_opts->dist);
      }

      // lambda
      for (i = 0; i < topoinfo->numR; i++){
          if (topoinfo->prsrandrange[2][3*i + 2 + 2*topoinfo->numG] == 1) {      // Activation
            paras[3*i + 2 + 2*topoinfo->numG] =     randu(topoinfo->prsrandrange[0][3*i + 2 + 2*topoinfo->numG], topoinfo->prsrandrange[1][3*i + 2 + 2*topoinfo->numG]);
          }
          else if (topoinfo->prsrandrange[2][3*i + 2 + 2*topoinfo->numG] == 2) { // Inhibition
            paras[3*i + 2 + 2*topoinfo->numG] = 1.0/randu(topoinfo->prsrandrange[0][3*i + 2 + 2*topoinfo->numG], topoinfo->prsrandrange[1][3*i + 2 + 2*topoinfo->numG]);
          }
      }
      break;
    case 2:
      // Production rate
      for (i = 0; i < topoinfo->numG; i++){
          paras[i] = randpg(topoinfo->prsrandrange[0][i], topoinfo->prsrandrange[1][i]);
      }

      // Degradation rate
      for (i = 0; i < topoinfo->numG; i++){
          paras[i + topoinfo->numG] = randpg(topoinfo->prsrandrange[0][i + topoinfo->numG], topoinfo->prsrandrange[1][i + topoinfo->numG]);
      }

      // Threshold
      for (i = 0; i < topoinfo->numR; i++){
          paras[3*i + 2*topoinfo->numG] = randpg(topoinfo->prsrandrange[0][3*i + 2*topoinfo->numG], topoinfo->prsrandrange[1][3*i + 2*topoinfo->numG]);
      }

      // Coefficient
      for (i = 0; i < topoinfo->numR; i++){
          paras[3*i + 1 + 2*topoinfo->numG] = randd(simu_opts->dist);
      }

      // lambda
      for (i = 0; i < topoinfo->numR; i++){
          
          if (topoinfo->prsrandrange[2][3*i + 2 + 2*topoinfo->numG] == 1) {      // Activation
            paras[3*i + 2 + 2*topoinfo->numG] =     randfd(topoinfo->prsrandrange[0][3*i + 2 + 2*topoinfo->numG], topoinfo->prsrandrange[1][3*i + 2 + 2*topoinfo->numG], simu_opts->dist);
          }
          else if (topoinfo->prsrandrange[2][3*i + 2 + 2*topoinfo->numG] == 2) { // Inhibition
            paras[3*i + 2 + 2*topoinfo->numG] = 1.0/randfd(topoinfo->prsrandrange[0][3*i + 2 + 2*topoinfo->numG], topoinfo->prsrandrange[1][3*i + 2 + 2*topoinfo->numG], simu_opts->dist);
          }
      }
      break;
    case 3:
      // Production rate
      for (i = 0; i < topoinfo->numG; i++){
          paras[i] = randexp(topoinfo->prsrandrange[0][i]);
      }

      // Degradation rate
      for (i = 0; i < topoinfo->numG; i++){
          paras[i + topoinfo->numG] = randexp(topoinfo->prsrandrange[0][i + topoinfo->numG]);
      }

      // Threshold
      for (i = 0; i < topoinfo->numR; i++){
          paras[3*i + 2*topoinfo->numG] = randexp(topoinfo->prsrandrange[0][3*i + 2*topoinfo->numG]);
      }

      // Coefficient
      for (i = 0; i < topoinfo->numR; i++){
          paras[3*i + 1 + 2*topoinfo->numG] = randd(simu_opts->dist);
      }

      // lambda
      for (i = 0; i < topoinfo->numR; i++){
          
          if (topoinfo->prsrandrange[2][3*i + 2 + 2*topoinfo->numG] == 1) {      // Activation
            paras[3*i + 2 + 2*topoinfo->numG] =     randfd(topoinfo->prsrandrange[0][3*i + 2 + 2*topoinfo->numG], 0, simu_opts->dist);
          }
          else if (topoinfo->prsrandrange[2][3*i + 2 + 2*topoinfo->numG] == 2) { // Inhibition
            paras[3*i + 2 + 2*topoinfo->numG] = 1.0/randfd(topoinfo->prsrandrange[0][3*i + 2 + 2*topoinfo->numG], 0, simu_opts->dist);
          }
      }
      break;
  }
  
  if (simu_opts->KDID != 0) {
    if (simu_opts->KDID <= topoinfo->numG) {
      paras[simu_opts->KDID - 1] = 0.0;
    }
    else {
      paras[topoinfo->ParasPos[simu_opts->KDID - topoinfo->numG - 1] + 2] = 1.0;
    }
  }
}

void model_ODE(double t,  double *ytmp, double *yp, double *p, struct topo *topoinfo)
{
  int i     = 0;
  int j     = 0;
  int count = 1;

  for(i = 0; i < topoinfo->numG; i++){
    yp[i] = p[i]; 

    for (j = 0; j < topoinfo->numR; j++){
      if (topoinfo->TargetG[j] == i){
        if (topoinfo->TypeR[j] == 1){ //Activation 
          yp[i] = yp[i]*(Hillshift(ytmp[topoinfo->SourceG[j]], p[2*topoinfo->numG+3*(count-1)], p[2*topoinfo->numG+3*(count-1)+1], p[2*topoinfo->numG+3*(count-1)+2])/p[2*topoinfo->numG+3*(count-1)+2]);
        }
        else if (topoinfo->TypeR[j] == 2) { //Inhibition
          yp[i] = yp[i]*Hillshift(ytmp[topoinfo->SourceG[j]], p[2*topoinfo->numG+3*(count-1)], p[2*topoinfo->numG+3*(count-1)+1], p[2*topoinfo->numG+3*(count-1)+2]);
        }
        count++;
      }
      
    }
    yp[i] = yp[i] - p[i + topoinfo->numG]*ytmp[i];
  }
}

void RIVs(double *y, double *ytmp, double *p, struct topo *topoinfo)
{
  int    i         = 0;
  int    j         = 0;
  int    count     = 1;
  double minV      = 0;
  double maxV      = 0;

  for(i = 0; i < topoinfo->numG; i++){

    if (p[i] != 0){
      minV = p[i];
      maxV = p[i];

      for (j = 0; j < topoinfo->numR; j++){
        if (topoinfo->TargetG[j] == i){
          if (topoinfo->TypeR[j] == 1){ //Activation  
              minV = minV*(1.0/p[2*topoinfo->numG+3*(count-1)+2]); 
          }
          else if (topoinfo->TypeR[j] == 2) { //Inhibition
              minV = minV*p[2*topoinfo->numG+3*(count-1)+2];
          }       
          count++;
        }
      }

      minV = minV/p[topoinfo->numG + i];

      maxV = maxV/p[topoinfo->numG + i];

      y[i] = exp2(randu(log2(minV), log2(maxV)));
    }
    else {
      y[i]    = 0.0;
      ytmp[i] = 0.0;
    }
  }
}

int solve_ODE (double *y_store, int j, double *paras, struct opts *simu_opts, struct topo *topoinfo)
{

  int    n_step     = 1000;
  int    i_step     = 1;
  int    i          = 0;
  double testdelta  = 0.0;
  double t          = 0.0;
  double t_start    = 0.0;
  double t_stop     = 0.0;
  double *y;
  double *yp;
  double *ytmp;

  y     = (double *)calloc(topoinfo->numG,      sizeof(double));
  yp    = (double *)calloc(topoinfo->numG,      sizeof(double));
  ytmp  = (double *)calloc(topoinfo->numG,      sizeof(double));

  for (i = 0; i < topoinfo->numG; i++){
    ytmp[i] = 2000.0;
  }
 
  int cnt_loop = 0;
  
  RIVs(y, ytmp, paras, topoinfo);
    
  testdelta = sumdelta(y, ytmp, topoinfo->numG);
  
  while (testdelta != 0 && cnt_loop < 20) {
    t_start = t_stop;
    t_stop  = t_stop + 100;
    
    cnt_loop = cnt_loop + 1;
  
    for ( i_step = 1; i_step <= n_step; i_step++ )
    { 
      for (i = 0; i < topoinfo->numG; i++){
          ytmp[i] = y[i];
      }
        
      model_ODE ( t, ytmp, yp, paras, topoinfo );
      
      for (i = 0; i < topoinfo->numG; i++){
              y[i] = ytmp[i] + yp[i]*simu_opts->stepsize;
          }
          
          t = t + simu_opts->stepsize;
      }
      
      testdelta = sumdelta(y, ytmp, topoinfo->numG);
    }
  
  for (i = 0; i < topoinfo->numG; i++){
    y_store[topoinfo->numG*j + i] = y[i];
  }
  
  free(y);
  free(yp);
  free(ytmp);

  return cnt_loop;
}

int count_state (double *y_store, double *soln, struct opts *simu_opts, struct topo *topoinfo)
{
    int i         = 0;
    int j         = 0;
    int h         = 0;
    int count     = 0;
    int cnt       = 1;
    double delta  = 0.0;
    double sumpow = 0.0;
    
    for (h = 1; h <= topoinfo->numG; h++){
        soln[h - 1] = y_store[h - 1]; 
    }
    
    for (i = 2; i <= simu_opts->num_ode; i++){
        count = 0;
        
        for (j = 1; j <= cnt; j++){
            h = 1;
            sumpow = 0.0;
            while (h <= topoinfo->numG) {
                sumpow = sumpow + pow((y_store[topoinfo->numG*(i-1) + h - 1] - soln[topoinfo->numG*(j-1) + h - 1]), 2);
                h++;
            }
        
            delta = sqrt(sumpow);
            
            if (delta > simu_opts->thrd){
                count = count + 1;
            } 
        }    
        
        if (count == cnt){
            cnt = cnt + 1;
            
            if (cnt <= simu_opts->num_stability) {
                for (h = 1; h <= topoinfo->numG; h++){
                    soln[(cnt-1)*topoinfo->numG + h - 1] = y_store[topoinfo->numG*(i-1) + h - 1]; 
                }
            }
            else{
                cnt = simu_opts->num_stability;
                break;
            }
        }
    }
    
    return cnt;
}

void T_test(char *inputfile, struct opts *simu_opts, struct topo *topoinfo, int num, int cnt, double *soln, double *paras)
{
  static FILE *f_test = NULL;
  char ftname [100] = "";
  char *modelname;

  int i = 0;
  int j = 0;
  int h = 0;
  double tmp = 0.0;
  int *localnumover;
  int *localnumdown;

  modelname = strtok(inputfile, ".");

  if (simu_opts->Toggle_T_test == 1) {
    if (f_test == NULL) {
      strcpy(ftname, modelname);
      strcat(ftname, "_T_test");
      strcat(ftname, ".dat");
      f_test   = fopen(ftname, "w");
    }

    localnumover = (int *) calloc(topoinfo->numG, sizeof(int));
    localnumdown = (int *) calloc(topoinfo->numG, sizeof(int));

    for (i = 0; i < cnt; i++){
      for (j = 0; j < topoinfo->numG; j++){
        for (h = 0; h < topoinfo->numR; h++){
          if (topoinfo->SourceG[h] == j){
            tmp = soln[topoinfo->numG*i + j]/paras[topoinfo->ParasPos[h]];
            if (tmp >= 1){
              localnumover[j] = localnumover[j] + 1;
            }
            else {
              localnumdown[j] = localnumdown[j] + 1;
            }
          }
        }
      }
    }

    fprintf(f_test, "%d", num + 1);
    for (j = 0; j < topoinfo->numG; j++){
      fprintf(f_test, "\t%d\t%d", localnumover[j], localnumdown[j]);
      topoinfo->numover[j] = topoinfo->numover[j] + localnumover[j];
      topoinfo->numdown[j] = topoinfo->numdown[j] + localnumdown[j];
    }
    fprintf(f_test, "\n");

    free(localnumover);
    free(localnumdown);
  
    if (simu_opts->num_paras == num + 1){
      fclose(f_test);
      f_test = NULL;

      printf("\n-------------------T_test------------------\n");
      printf("Gene_ID\tProbs_over_T\n");
      for (i = 0; i < topoinfo->numG; i++){
         printf("%d\t%f\n", i, (double)topoinfo->numover[i]/(double)(topoinfo->numover[i]+topoinfo->numdown[i]));
      }
      printf("-------------------------------------------\n");
    }
  
  }
}

void release_memory(char *fileextension, struct topo *topoinfo, int *cnt_store, double *y_store, double *soln, double *paras)
{
  int i = 0;

  free(topoinfo->SourceG);
  free(topoinfo->TargetG);
  free(topoinfo->TypeR);
  free(topoinfo->ParasPos);
  free(topoinfo->numover);
  free(topoinfo->numdown);

  if (strcmp(fileextension, "topo") == 0){
    for (i = 0; i < topoinfo->numG; i++) {
      free(topoinfo->Gname[i]);
    }
    free(topoinfo->Gname);
  }

  free(topoinfo->prsrandrange[0]);
  free(topoinfo->prsrandrange[1]);
  free(topoinfo->prsrandrange[2]);
  free(topoinfo->prsrandrange);

  free(cnt_store);
  free(y_store);
  free(soln);
  free(paras);
}







