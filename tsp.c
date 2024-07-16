/*****************************************************************************
  A template program for 2-dimensional euclidean symmetric TSP solver. 
  Subroutines to read instance data and compute the objective value of a given 
  tour (solution) are included. The one to output the computed tour in the
  TSPLIB format is also included. 

  The URL of TSPLIB is:
       http://www.iwr.uni-heidelberg.de/groups/comopt/software/TSPLIB95/

  NOTE: Indices of nodes range from 0 to n-1 in this program,
        while it does from 1 to n in "README.eng" and the data files of
        instances and of tours. 

  If you would like to use various parameters, it might be useful to modify
  the definition of struct "Param" and mimic the way the default value of
  "timelim" is given and how its value is input from the command line.
******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <math.h>
#include "cpu_time.c"
#include <stdbool.h>


/***** constants *************************************************************/
#define MAX_STR    1024

/***** macros ****************************************************************/
#define dist(k,l) ( (int)( sqrt( (tspdata->x[k]-tspdata->x[l])*(tspdata->x[k]-tspdata->x[l]) + (tspdata->y[k]-tspdata->y[l])*(tspdata->y[k]-tspdata->y[l]) ) + 0.5 ) )

/***** default values of parameters ******************************************/
#define TIMELIM    300 /* the time limit for the algorithm in seconds */
#define GIVESOL    0   /* 1: input a solution; 0: do not give a solution */
#define OUTFORMAT  2   /* 0: do not output the tour;
			  1: output the computed tour in TSPLIB format;
			  2: output the computed tour in TSP_VIEW format */
#define TOURFILE   "result.tour"
                       /* the output file of computed tour */

typedef struct {
  int    timelim;              /* the time limit for the algorithm in secs. */
  int    givesol;              /* give a solution (1) or not (0) */
  int    outformat;            /* 1: output the computed tour in TSPLIB format;
				  0: do not output it */
  char   tourfile[MAX_STR];    /* the output file of computed tour */
  /* NEVER MODIFY THE ABOVE VARIABLES.  */
  /* You can add more components below. */

} Param;                /* parameters */


typedef struct{
  char     name[MAX_STR];         /* name of the instance */
  int      n;                     /* number of nodes */
  double   *x;                    /* x-coordinates of nodes */
  double   *y;                    /* y-coordinates of nodes */
  int      min_node_num;          /* minimum number of nodes the solution contains */
} TSPdata;              /* data of TSP instance */

typedef struct {
  double        timebrid;       /* the time before reading the instance data */
  double        starttime;      /* the time the search started */
  double        endtime;        /* the time the search ended */
  int           *bestsol;       /* the best solution found so far */
  /* NEVER MODIFY THE ABOVE FOUR VARIABLES. */
  /* You can add more components below. */

} Vdata;                /* various data often necessary during the search */

// 解の表現を変更：訪問する都市とその順序を保持
typedef struct {
    int *cities;  // 訪問する都市のリスト
    int size;     // 実際に訪問する都市の数
} Solution;

/************************ declaration of functions ***************************/
FILE *open_file( char *fname, char *mode );
void *malloc_e( size_t size );

void copy_parameters( int argc, char *argv[], Param *param );
void prepare_memory( TSPdata *tspdata, Vdata *vdata );
void read_header( FILE *in, TSPdata *tspdata );
void read_tspfile( FILE *in, TSPdata *tspdata, Vdata *vdata );
void read_tourfile( FILE *in, TSPdata *tspdata, int *tour );
void output_tour( FILE *out, TSPdata *tspdata, int *tour );
void output_tour_for_tsp_view( FILE *out, TSPdata *tspdata, int *tour );
void recompute_obj( Param *param, TSPdata *tspdata, Vdata *vdata );
void my_algorithm( Param *param, TSPdata *tspdata, Vdata *vdata );

int compute_cost( TSPdata *tspdata, int *tour );
int is_feasible( TSPdata *tspdata, int *tour );

/***** open the file with given mode *****************************************/
FILE *open_file( char *fname, char *mode ){
  FILE *fp;
  fp=fopen(fname,mode);
  if(fp==NULL){
    fprintf(stderr,"file not found: %s\n",fname);
    exit(EXIT_FAILURE);
  }
  return fp;
}

/***** malloc with error check ***********************************************/
void *malloc_e( size_t size ){
  void *s;
  if ( (s=malloc(size)) == NULL ) {
    fprintf( stderr, "malloc : not enough memory.\n" );
    exit(EXIT_FAILURE);
  }
  return s;
}


/***** copy and read the parameters ******************************************/
/***** Feel free to modify this subroutine. **********************************/
void copy_parameters( int argc, char *argv[], Param *param ){
    
  /**** copy the default parameters ****/
  param->timelim    = TIMELIM;
  param->givesol    = GIVESOL;
  param->outformat  = OUTFORMAT;
  strcpy(param->tourfile,TOURFILE);
  
  /**** read the parameters ****/
  if(argc>0 && (argc % 2)==0){
    printf("USAGE: ./%s [param_name, param_value] [name, value]...\n",argv[0]);
    exit(EXIT_FAILURE);
  }
  else{
    int i;
    for(i=1; i<argc; i+=2){
      if(strcmp(argv[i],"timelim")==0)    param->timelim    = atoi(argv[i+1]);
      if(strcmp(argv[i],"givesol")==0)    param->givesol    = atoi(argv[i+1]);
      if(strcmp(argv[i],"outformat")==0)  param->outformat  = atoi(argv[i+1]);
      if(strcmp(argv[i],"tourfile")==0)   strcpy(param->tourfile,argv[i+1]);
    }
  }
}


/***** prepare memory space **************************************************/
/***** Feel free to modify this subroutine. **********************************/
void prepare_memory( TSPdata *tspdata, Vdata *vdata ){
  int k,n;
  n=tspdata->n;
  tspdata->x       = (double*)malloc_e(n*sizeof(double));
  tspdata->y       = (double*)malloc_e(n*sizeof(double));
  vdata->bestsol   = (int*)malloc_e(n*sizeof(int));
  /* the next line is just to give an initial solution */
  for(k=0;k<n;k++)
    vdata->bestsol[k]=k;
}

/***** reading the header of a file in TSPLIB format *************************/
/***** NEVER MODIFY THIS SUBROUTINE! *****************************************/
void read_header( FILE *in, TSPdata *tspdata ){
  char str[MAX_STR],name[MAX_STR],dim[MAX_STR],type[MAX_STR],edge[MAX_STR],min[MAX_STR];
  int flag=0;

  for(;;){
    char *w,*u;
    /* error */
    if(fgets(str,MAX_STR,in)==NULL){
      fprintf(stderr,"error: invalid data input.\n");
      exit(EXIT_FAILURE);
    }
    /* halt condition */
    if(strcmp(str,"NODE_COORD_SECTION\n")==0){ break; }
    if(strcmp(str,"TOUR_SECTION\n")==0){ flag=1; break; }
    /* data input */
    w = strtok(str," :\n");
    u = strtok(NULL," :\n");
    if(w==NULL || u==NULL) continue;
    if(strcmp("NAME",w)==0)                  strcpy(name,u);
    if(strcmp("DIMENSION",w)==0)             strcpy(dim,u);
    if(strcmp("TYPE",w)==0)                  strcpy(type,u);
    if(strcmp("EDGE_WEIGHT_TYPE",w)==0)      strcpy(edge,u);
    if(strcmp("MIN_NODE_NUM",w)==0)          strcpy(min,u);
  }

  /* read a TSP instance */
  if(flag==0){
    strcpy(tspdata->name,name);
    tspdata->min_node_num=atoi(min);
    tspdata->n=atoi(dim);
    if(strcmp("TSP",type)!=0 || strcmp("EUC_2D",edge)!=0 ){
      fprintf(stderr,"error: invalid instance.\n");
      exit(EXIT_FAILURE);
    }
  }
  /* read a tour */
  else{
    if(strcmp("TOUR",type)!=0){
      fprintf(stderr,"error: invalid tour.\n");
      exit(EXIT_FAILURE);
    }
  }
}

/***** reading the file of TSP instance **************************************/
/***** NEVER MODIFY THIS SUBROUTINE! *****************************************/
void read_tspfile( FILE *in, TSPdata *tspdata, Vdata *vdata ){
  char str[MAX_STR];
  int k;

  /* reading the instance */
  read_header(in,tspdata);
  prepare_memory(tspdata,vdata);
  for(k=0;k<tspdata->n;k++){
    int dummy;
    if(fgets(str,MAX_STR,in)==NULL) break;
    if(strcmp(str,"EOF\n")==0) break;
    sscanf(str,"%d%lf%lf",&dummy, &(tspdata->x[k]),&(tspdata->y[k]));
  }
  if(k!=tspdata->n){
    fprintf(stderr,"error: invalid instance.\n");
    exit(EXIT_FAILURE);
  }
}

/***** read the tour in the TSPLIB format with feasibility check *************/
/***** NEVER MODIFY THIS SUBROUTINE! *****************************************/
void read_tourfile( FILE *in, TSPdata *tspdata, int *tour ){
  int k;

  read_header(in,tspdata);
  for(k=0;k<tspdata->n;k++){
    int val;
    if(fscanf(in,"%d",&val)==EOF) break;
    if(val==-1) break;
    tour[k]=val-1;
  }
  if(k!=tspdata->n){
    fprintf(stderr,"error: invalid tour.\n");
    exit(EXIT_FAILURE);
  }
}


/***** output the tour in the TSPLIB format **********************************/
/***** note: the output tour starts from the node "1" ************************/
void output_tour( FILE *out, TSPdata *tspdata, int *tour ){
  int k,idx=0;

  fprintf(out,"NAME : %s\n",tspdata->name);
  fprintf(out,"COMMENT : tour_length=%d\n",compute_cost(tspdata,tour));
  fprintf(out,"TYPE : TOUR\n");
  fprintf(out,"DIMENSION : %d\n",tspdata->n);
  fprintf(out,"TOUR_SECTION\n");
  for(k=0;k<tspdata->n;k++)
    if(tour[k]==0)
      { idx=k; break; }
  for(k=idx;k<tspdata->n;k++){
    if(tour[k]<0) break;
    fprintf(out,"%d\n",tour[k]+1);
  }
  for(k=0;k<idx;k++)
    fprintf(out,"%d\n",tour[k]+1);
  fprintf(out,"-1\n");
  fprintf(out,"EOF\n");
}

/***** output the tour in the TSP_VIEW format **********************************/
/***** note: the indices of the tour starts from "0" ***************************/
void output_tour_for_tsp_view( FILE *out, TSPdata *tspdata, int *tour ){
  int k;

  fprintf(out, "%d\n", tspdata->n);
  for(k=0; k<tspdata->n; k++){
    fprintf(out, "%g %g\n", tspdata->x[k], tspdata->y[k]);
  }
  for(k=0; k<tspdata->n; k++){
    if(tour[k]<0) break;
    fprintf(out,"%d\n", tour[k]);
  }
}

/***** check the feasibility and recompute the cost **************************/
/***** NEVER MODIFY THIS SUBROUTINE! *****************************************/
void recompute_obj( Param *param, TSPdata *tspdata, Vdata *vdata ){
  if(!is_feasible(tspdata,vdata->bestsol)){
    fprintf(stderr,"error: the computed tour is not feasible.\n");
    exit(EXIT_FAILURE);
  }
  printf("recomputed tour length = %d\n",
	 compute_cost(tspdata,vdata->bestsol));
  printf("time for the search:   %7.2f seconds\n",
         vdata->endtime - vdata->starttime);
  printf("time to read the instance: %7.2f seconds\n",
         vdata->starttime - vdata->timebrid);
}

/***** cost of the tour ******************************************************/
/***** NEVER MODIFY THIS SUBROUTINE! *****************************************/
int compute_cost( TSPdata *tspdata, int *tour ){
  int k,cost=0,n;
  n=tspdata->n;

  for(k=0;k<n-1;k++){
    if(tour[k+1]<0) break;
    cost += dist(tour[k],tour[k+1]);
  }
  cost += dist(tour[k],tour[0]);
  return cost;
}

/***** check the feasibility of the tour *************************************/
/***** NEVER MODIFY THIS SUBROUTINE! *****************************************/
int is_feasible( TSPdata *tspdata, int *tour ){
  int k,n,*visited,flag=1,num_visited=0;
  n=tspdata->n;
  visited=(int*)malloc_e(n*sizeof(int));

  for(k=0;k<n;k++)
    visited[k]=0;
  for(k=0;k<n;k++){
    if(tour[k]<0)
      { break; }
    if(tour[k]>=n)
      { flag=0; break; }
    if(visited[tour[k]])
      { flag=0; break; }
    else{
      visited[tour[k]]=1;
      num_visited++;
    }
  }
  if(num_visited < tspdata->min_node_num) flag=0;

  free(visited);
  /* if tour is feasible (resp., not feasible), then flag=1 (resp., 0) */
  return flag;
}

/***** Simulated Annealing algorithm **************************************/

// 部分巡回路の妥当性チェック
bool is_valid_partial_tour(Solution *sol, TSPdata *tspdata) {
  if (sol->size < tspdata->min_node_num || sol->size > tspdata->n) {
      return false;
  }
  
  bool *visited = calloc(tspdata->n, sizeof(bool));
  for (int i = 0; i < sol->size; i++) {
    if (sol->cities[i] < 0 || sol->cities[i] >= tspdata->n || visited[sol->cities[i]]) {
        free(visited);
        return false;
    }
    visited[sol->cities[i]] = true;
  }
  free(visited);
  return true;
}

// 部分巡回路のコスト計算 (修正版)
int compute_partial_cost(TSPdata *tspdata, Solution *sol) {
    if (!is_valid_partial_tour(sol, tspdata)) {
      return INT_MAX;
    }
    
    int cost = 0;
    for (int i = 0; i < sol->size - 1; i++) {
      cost += dist(sol->cities[i], sol->cities[i+1]);
    }
    cost += dist(sol->cities[sol->size-1], sol->cities[0]);
    return cost;
}

// 初期解の生成
void generate_initial_solution(TSPdata *tspdata, Solution *sol) {
    sol->size = tspdata->min_node_num;
    sol->cities = (int*)malloc_e(sol->size * sizeof(int));
    
    // ランダムに都市を選択
    bool *selected = (bool*)calloc(tspdata->n, sizeof(bool));
    for (int i = 0; i < sol->size; i++) {
        int city;
        do {
            city = rand() % tspdata->n;
        } while (selected[city]);
        sol->cities[i] = city;
        selected[city] = true;
    }
    free(selected);
}

// 2-opt操作（部分巡回路用）
void two_opt_partial(Solution *sol, int i, int j) {
  int k;
  for (k = 0; k < (j-i+1)/2; k++) {
      int temp = sol->cities[i+k];
      sol->cities[i+k] = sol->cities[j-k];
      sol->cities[j-k] = temp;
  }
}

// 都市の挿入操作 (修正版)
bool insert_city(Solution *sol, TSPdata *tspdata) {
  if (sol->size >= tspdata->n) return false;
  
  int new_city;
  bool *in_tour = calloc(tspdata->n, sizeof(bool));
  for (int i = 0; i < sol->size; i++) {
    in_tour[sol->cities[i]] = true;
  }
  
  do {
    new_city = rand() % tspdata->n;
  } while (in_tour[new_city]);
  
  free(in_tour);
  
  sol->size++;
  sol->cities = realloc(sol->cities, sol->size * sizeof(int));
  int insert_pos = rand() % sol->size;
  
  for (int i = sol->size - 1; i > insert_pos; i--) {
    sol->cities[i] = sol->cities[i-1];
  }
  sol->cities[insert_pos] = new_city;
  
  return true;
}

// 都市の削除操作 (修正版)
bool remove_city(Solution *sol, TSPdata *tspdata) {
  if (sol->size <= tspdata->min_node_num) return false;
  
  int remove_pos = rand() % sol->size;
  for (int i = remove_pos; i < sol->size - 1; i++) {
    sol->cities[i] = sol->cities[i+1];
  }
  sol->size--;
  
  return true;
}

// ソリューションの深いコピーを行う関数
void copy_solution(Solution *dest, const Solution *src) {
  dest->size = src->size;
  dest->cities = (int*)malloc_e(dest->size * sizeof(int));
  memcpy(dest->cities, src->cities, dest->size * sizeof(int));
}

// ソリューションのメモリを解放する関数
void free_solution(Solution *sol) {
  if (sol->cities != NULL) {
      free(sol->cities);
      sol->cities = NULL;
  }
  sol->size = 0;
}

//焼きなまし法
void simulated_annealing_partial_tsp(Param *param, TSPdata *tspdata, Vdata *vdata) {
  Solution current_sol, new_sol, best_sol;
  generate_initial_solution(tspdata, &current_sol);
  copy_solution(&best_sol, &current_sol);
  int current_cost = compute_partial_cost(tspdata, &current_sol);
  int best_cost = current_cost;
  double temp = 10000.0, cooling_rate;

  if (tspdata->n < 500) {
    cooling_rate = 0.9999999;
  } else if (tspdata->n < 1000) {
    cooling_rate = 0.99999977;
  } else if (tspdata->n < 5000) {
    cooling_rate = 0.99999975;
  } else if (tspdata->n < 10000) {
    cooling_rate = 0.9999995;
  } else if (tspdata->n < 15000) {
    cooling_rate = 0.9999992;
  } else {
    cooling_rate = 0.999998;
  }

  while (cpu_time() - vdata->starttime < param->timelim && temp > 1.0) {
    copy_solution(&new_sol, &current_sol);
    // 近傍操作の選択
    int op = rand() % 4;  // 4つの基本操作に戻す
    bool operation_successful = false;
      
    switch (op) {
      case 0: // 2-opt
        {
          int i = rand() % new_sol.size;
          int j = (i + 1 + rand() % (new_sol.size - 1)) % new_sol.size;
          two_opt_partial(&new_sol, i, j);
          operation_successful = true;
        }
        break;
      case 1: // 都市の挿入
        operation_successful = insert_city(&new_sol, tspdata);
        break;
      case 2: // 都市の削除
        operation_successful = remove_city(&new_sol, tspdata);
        break;
      case 3: // 都市の交換
        {
          int i = rand() % new_sol.size;
          int j = rand() % new_sol.size;
          int temp = new_sol.cities[i];
          new_sol.cities[i] = new_sol.cities[j];
          new_sol.cities[j] = temp;
          operation_successful = true;
        }
        break;
    }
    if (!operation_successful) {
      free_solution(&new_sol);
      continue;
    }
    int new_cost = compute_partial_cost(tspdata, &new_sol);
    if (new_cost == INT_MAX) {
      free_solution(&new_sol);
      continue;
    }
    double acceptance_prob = exp((current_cost - new_cost) / temp);
    
    if (new_cost < current_cost || ((double)rand() / RAND_MAX) < acceptance_prob) {
      free_solution(&current_sol);
      current_sol = new_sol;
      new_sol.cities = NULL;
      current_cost = new_cost;
      if (current_cost < best_cost) {
        free_solution(&best_sol);
        copy_solution(&best_sol, &current_sol);
        best_cost = current_cost;
      }
    } else {
      free_solution(&new_sol);
    }
      
    temp *= cooling_rate;
  }
  // 最良解を vdata->bestsol にコピー
  for (int i = 0; i < tspdata->n; i++) {
    vdata->bestsol[i] = -1;  // 初期化：すべての都市を未訪問に設定
  }
  for (int i = 0; i < best_sol.size; i++) {
    vdata->bestsol[i] = best_sol.cities[i];
  }
  free_solution(&current_sol);
  free_solution(&best_sol);
}

/***** main ******************************************************************/
int main(int argc, char *argv[]){

  Param     param;     /* parameters */
  TSPdata   tspdata;   /* data of TSP instance */
  Vdata     vdata;     /* various data often needed during search */

  vdata.timebrid = cpu_time();
  copy_parameters(argc, argv, &param);
  read_tspfile(stdin,&tspdata,&vdata);
  if(param.givesol==1) read_tourfile(stdin,&tspdata,vdata.bestsol);
  vdata.starttime = cpu_time();

  simulated_annealing_partial_tsp(&param,&tspdata,&vdata);

  vdata.endtime = cpu_time();
  recompute_obj(&param,&tspdata,&vdata);
  if(param.outformat==1){
    output_tour(open_file(param.tourfile,"w"),&tspdata,vdata.bestsol);}
  else if(param.outformat==2){
    output_tour_for_tsp_view(open_file(param.tourfile,"w"),&tspdata,vdata.bestsol);}

  return EXIT_SUCCESS;
}
