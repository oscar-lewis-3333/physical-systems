//
//  kinetic_fokker_planck_solver.c
//  Kinetic Fokker-Planck solver
//
//  Created by Oscar Lewis on 15/01/2026.
//

#include <stdlib.h>
#include <stdio.h>
#include <math.h> /* we take absolute values of doubles later in the code, and use exp.*/
#include <lapacke.h>

/* We reuse the code given in assignment 4 "bandutility.c"*/

struct band_mat{
  long ncol;        /* Number of columns in band matrix */
  long nbrows;      /* Number of rows (bands in original matrix) */
  long nbands_up;   /* Number of bands above diagonal */
  long nbands_low;  /* Number of bands below diagonal */
  double *array;    /* Storage for the matrix in banded format */
  /* Internal temporary storage for solving inverse problem */
  long nbrows_inv;  /* Number of rows of inverse matrix */
  double *array_inv;/* Store the inverse if this is generated */
  int *ipiv;        /* Additional inverse information */
};
/* Define a new type band_mat */
typedef struct band_mat band_mat;

/* Initialise a band matrix of a certain size, allocate memory,
   and set the parameters.  */
int init_band_mat(band_mat *bmat, long nbands_lower, long nbands_upper, long n_columns) {
  bmat->nbrows = nbands_lower + nbands_upper + 1;
  bmat->ncol   = n_columns;
  bmat->nbands_up = nbands_upper;
  bmat->nbands_low= nbands_lower;
  bmat->array      = (double *) malloc(sizeof(double)*bmat->nbrows*bmat->ncol);
  bmat->nbrows_inv = bmat->nbands_up*2 + bmat->nbands_low + 1;
  bmat->array_inv  = (double *) malloc(sizeof(double)*(bmat->nbrows+bmat->nbands_low)*bmat->ncol);
  bmat->ipiv       = (int *) malloc(sizeof(int)*bmat->ncol);
  if (bmat->array==NULL||bmat->array_inv==NULL) {
    return 0;
  }
  /* Initialise array to zero */
  long i;
  for (i=0;i<bmat->nbrows*bmat->ncol;i++) {
    bmat->array[i] = 0.0;
  }
  return 1;
}

/* Finalise function: should free memory as required */
void finalise_band_mat(band_mat *bmat) {
  free(bmat->array);
  free(bmat->array_inv);
  free(bmat->ipiv);
}

/* Get a pointer to a location in the band matrix, using
   the row and column indexes of the full matrix.           */
double *getp(band_mat *bmat, long row, long column) {
  int bandno = bmat->nbands_up + row - column;
  if(row<0 || column<0 || row>=bmat->ncol || column>=bmat->ncol ) {
    printf("Indexes out of bounds in getp: %ld %ld %ld \n",row,column,bmat->ncol);
    exit(1);
  }
  return &bmat->array[bmat->nbrows*column + bandno];
}

/* Retrun the value of a location in the band matrix, using
   the row and column indexes of the full matrix.           */
double getv(band_mat *bmat, long row, long column) {
  return *getp(bmat,row,column);
}

/* Set an element of a band matrix to a desired value based on the pointer
   to a location in the band matrix, using the row and column indexes
   of the full matrix.           */
double setv(band_mat *bmat, long row, long column, double val) {
  *getp(bmat,row,column) = val;
  return val;
}

/* Solve the equation Ax = b for a matrix a stored in band format
   and x and b real arrays                                          */
int solve_Ax_eq_b(band_mat *bmat, double *x, double *b) {
  /* Copy bmat array into the temporary store */
  int i,bandno;
  for(i=0;i<bmat->ncol;i++) {
    for (bandno=0;bandno<bmat->nbrows;bandno++) {
      bmat->array_inv[bmat->nbrows_inv*i+(bandno+bmat->nbands_low)] = bmat->array[bmat->nbrows*i+bandno];
    }
    x[i] = b[i];
  }

  long nrhs = 1;
  long ldab = bmat->nbands_low*2 + bmat->nbands_up + 1;
  int info = LAPACKE_dgbsv( LAPACK_COL_MAJOR, bmat->ncol, bmat->nbands_low, bmat->nbands_up, nrhs, bmat->array_inv, ldab, bmat->ipiv, x, bmat->ncol);
  return info;
}

int printmat(band_mat *bmat) {
  long i,j;
  for(i=0; i<bmat->ncol;i++) {
    for(j=0; j<bmat->nbrows; j++) {
      printf("%ld %ld %g \n",i,j,bmat->array[bmat->nbrows*i + j]);
    }
  }
  return 0;
}


int main(void) {
/* Initialising variables ready to be read from input.txt*/
    
    double t_f, L, v_m, C;
    long N_x, N_v, I_min;
    
    
    
    
/* Opening given input.txt file, reading if the constants have the correct type, and returning and error if they do not
     */
    FILE *finput = fopen("input.txt", "r");
    if (finput == NULL) {
        printf("Error: Input File is empty\n");
        return 1;
    }
    if (fscanf(finput, "%ld", &N_x) != 1 ||
        fscanf(finput, "%ld", &N_v) != 1||
        fscanf(finput, "%lf", &t_f) != 1 ||
        fscanf(finput, "%lf", &L) != 1 ||
        fscanf(finput, "%lf", &v_m) != 1 ||
        fscanf(finput, "%lf", &C) != 1 ||
        fscanf(finput, "%ld", &I_min) != 1)
    {
        printf("Error: Wrong formatting of input file\n");
        return 1;
    }
    fclose(finput);
    
/* Allocating memory for F in order to read it's coefficients */
    
    double *F = malloc(N_x*sizeof(double));
/* Checking successful allocation returning an error if not successful*/
    if (F == NULL) {
        printf("Error allocating memory for coefficient function\n");
        return 1;
    }
    
    FILE *fcoeff = fopen("coefficients.txt", "r");
    if (fcoeff == NULL) {
        printf("Coefficients file is empty\n");
        return 1;
    }
    for (long i = 0; i < N_x; i++) {
        if (fscanf(fcoeff, "%lf", &F[i]) != 1) {
            printf("Wrong formatting of coefficient file\n");
            return 1;
        }
    }
    fclose(fcoeff);
    
/* Defining dx, dv to be used later, giving x_j = j * dx, v_j = -v_m + j * dv as desired in spec */
    if (v_m == 0) {
        printf("Error: v_m = 0 causes an ill-posed, unsolvable PDE, as there is no way to satisfy all boundary conditions\n");
        free(F);
        return 1;
    }
/* For the above statement, if v_m = 0, we have v = 0 everywhere, hence f = 0 everywhere. But at x = 0, we must have f(0,0,t) = exp(0) = 1, and f(L,0,t) = exp(0) = 1, so we cannot satisfy both cases.*/
    
    double dx = L/(N_x - 1);
    double dv = 2* v_m / (N_v - 1);

/* We need to use upwinding scheme to solve, and we need to satisfy CFL conditions for x and v. As dv/dt is defined to be our coefficients function. Rerrange to get dt < dx / v, but we aim to maximise denominator. This gives the denominators to be v_m and the max absolute value of F*/
    
    double max_F = 0.0;
    for (long i = 0; i < N_x; i++) {
        if (max_F < fabs(F[i])) {
            max_F = fabs(F[i]);
        }
    }
    
    double dt_x = dx/fabs(v_m);
    double dt_v;
    if (max_F > 0) {
        dt_v = dv / max_F;
    }
    else  {dt_v = 1e10;}
    
/* Imposing CFL condtion for both advection terms. from individual online research, -lm is required to compile functions ceil and fmin on vonneuman (i am coding on my own system and cannot get vonneuman to work), so we use loops to do the same thing, because, as said in the spec, it must compile with "-Wall -std=c99"  */
    
    double dt_required = dt_x;
    if (dt_v < dt_required) {
        dt_required = dt_v;
    }
    
    
/* we need dt < dt_required, so we originally choose a dt slightly less than dt_required */
    
    double dt_1 = 0.9 * dt_required;
    
/* We need to find the optimal N_t for the requirements set out in the specification. It must be at least I_min, but it could be more if I_min violates CFL. We want to use the ceiling function as we want the smallest integer N_t such that N_t > t_f / dt, but for reasons voiced above, we use the following loop to imitate it. If that N_t is still less than I_min, we overwrite N_t to be I_min.*/
    
    double ratio = t_f/dt_1;
    long N_t = ratio;
    if (ratio > (double) N_t) {
        N_t = N_t + 1;
    }
    if (N_t < I_min) {
        N_t = I_min;
    }
/* we reset dt after finding N_t, else we may end up with an uneven time loop if N_t < I_min*/
    double dt = t_f/N_t;
    
    
    /* creating memory to store x, v grids*/
    double *x = malloc(N_x * sizeof(double));
    double *v = malloc(N_v * sizeof(double));
    
/* Checking successful allocation again. If something goes wrong, we exit immediately, so we have to free any memory previously correctly allocated to avoid a memory leak. This pattern will follow below and not be commented on again since it is the same pattern. The if (.) free (.) pattern is because only one may fail, so the other must also be freed. Same pattern will be repeated and not commented on again for the same reasons. */
    if(v == NULL || x == NULL ) {
        printf("Grid Memory allocation failed\n");
        free(F);
        if (x) free(x);
        if (v) free(v);
        return 1;
    }
    
    for (long i = 0; i < N_x; i++) {
        x[i] = i * dx;
    }
    
    for (long j = 0; j < N_v; j++) {
        v[j] = -v_m + j*dv;
    }
/* We create matrices for use later, we need tridiagonal matrices as we have terms involving i-1, i, i+1, (0 < i < N - 1) from finite difference methods, and then handling cases i = 0, i = N-1 with given boundary conditions. We also ensure nothing goes wrong simultaenously, with "init_band_mat" in the if loop also initialising as wanted. */
    
    band_mat xmat;
    band_mat vmat;
    
    if (!init_band_mat(&xmat, 1, 1, N_x)) {
        printf("Error creating x band matrix\n");
        free(x);
        free(v);
        free(F);
        return 1;
    }
    if (!init_band_mat(&vmat, 1, 1, N_v)) {
        printf("Error creating v band matrix\n");
        finalise_band_mat(&xmat);
        free(x);
        free(v);
        free(F);
        return 1;
    }
    
/* allocating memory for the function of interest. Calloc is given the size of what 'array' you want, and what types are going to be placed into your array, and returns that array but with all zeroes instead of not specifying the values. This gives us a function starting at zero, satisfying the initial condition f(x, v, 0) = 0 here, meaning we do not have to worry about it later. */
        
    double *f  = calloc(N_x * N_v, sizeof(double));
    if (f == NULL) {
        printf("Error allocating memory to solver function\n");
        free(F);
        free(x);
        free(v);
        finalise_band_mat(&xmat);
        finalise_band_mat(&vmat);
        return 1;
    }

/* before starting our time loop, we introduce a temporary solution as it neccessary later, and the memory allocation for the RHS solution of each linear system as we iterate over time*/
    double *f_temp = malloc(N_x * N_v * sizeof(double));
    double *rhs_x = malloc(N_x * sizeof(double));
    double *rhs_v = malloc(N_v * sizeof(double));
    double *sol_x = malloc(N_x * sizeof(double));
    double *sol_v = malloc(N_v * sizeof(double));
    if (sol_v == NULL||sol_x == NULL||rhs_v == NULL || rhs_x == NULL || f_temp == NULL) {
        printf("Error allocating memory just before time loop\n");
        free(f);
        free(F);
        free(x);
        free(v);
        finalise_band_mat(&vmat);
        finalise_band_mat(&xmat);
        if (f_temp) free(f_temp);
        if (rhs_x) free (rhs_x);
        if (rhs_v) free (rhs_v);
        if (sol_x) free (sol_x);
        if (sol_v) free (sol_v);
        return 1;
    }
    
/* We start the time loop, and write the rest of the code inside there*/
    for (long n = 0; n < N_t; n++) {
        
        /* We solve implicitly in each direction, using the independence of the x, v coordinate system. Since they are independent, it allows us to say that dv/dx = 0, and dF(x)/ dv = 0 since F is strictly a function of x. We can use the product rule and previous 2 relations to rewrite our PDE as (where here d/dx,d/dv, d/dt represent partial not full derivatives and we use conditions (2) dx/dt = v, dv/dt = F(x)): df/dt = - v * df/dx - F(x) * df/dv + C * d^2(f)/dv^2. We split into 2 PDE's, solve one first,then store that solution, and use it to (independently of the first solution) solve the second PDE, hence having solved the original PDE (obviously applying the given boundary conditions where neccessary). Also similar to approaches ive seen in many maths PDE modules ive taken. */
        
        /* We start by solving df/dt = -v * df/dx. At the start of every loop, we want our band matrix to be zero to avoid any coefficients being left over from previous loops, hence we explicitly set every thing to be zero.*/
        for (long j = 0; j < N_v; j++) {
            for(long i = 0; i < N_x; i++) {
                setv(&xmat, i, i, 0);
                if (i > 0) {
                    setv(&xmat, i, i-1, 0);
                }
                if (i < N_x - 1) {
                    setv(&xmat, i, i + 1, 0);
                }
                rhs_x[i] = f[i + j*N_x];
            }
            /* If v = 0, then df/dt = 0 and nothing changes in that time step, hence we just copy f into f_temp, and continue with our loop*/
            if (v[j] == 0.0) {
                for (long i = 0;i< N_x; i++) {
                    f_temp[i + j*N_x] = f[i + j* N_x];
                }
                continue;
            }
            /* Due to use of upwinding scheme, we need to consider the cases of v>0, v<0, v = 0, and using the slides from week 6 and rearranging, we get the matrix coefficients described. */
            
            else if (v[j] > 0.0) {
                /* In top row, we need f_temp[0] = exp(-v[j]^2) to enforce BC at each time step, so we need to set matrix coefficients to be one and zero respectively*/
                setv(&xmat, 0, 0, 1);
                setv(&xmat, 0, 1, 0);
                rhs_x[0] = exp(-v[j]*v[j]);
                
                /* Using last slide of week 6 presentation to get coefficients (where f_temp := f^(j+1)_{k})*/
                for (long i = 1; i < N_x; i++) {
                    setv(&xmat, i, i, 1.0 + (v[j] * dt/dx));
                    setv(&xmat, i, i - 1, -v[j] * dt/dx);
                }
            }
            
            
            /*continuing analogously to above but using the other equation from week 6 presentation, and setting the boundary condtion analagously by defining the matrix to be one on the diagonal and zero elsewhere*/
            
            else { /* this is equivalent to v[j] < 0.0*/
                rhs_x[N_x - 1] = exp(-v[j]*v[j]);
                setv(&xmat, N_x - 1, N_x - 1, 1);
                setv(&xmat, N_x - 1, N_x - 2, 0);
                for (long i = 0; i < N_x - 1; i++) {
                    setv(&xmat, i, i, 1 - v[j]*dt/dx);
                    setv(&xmat, i, i + 1, v[j] * dt/dx);
                }
            }
            /* we solve the equation, and store the solution as f_temp for each i. We then sum over the j's*/
            int info = solve_Ax_eq_b(&xmat, sol_x, rhs_x);
            if (info != 0) {
                printf("Error when solving system in x-direction at time step %ld, x_index %ld\n", n, j);
            }
            for (long i = 0; i < N_x; i++) {
                f_temp[i + N_x*j] = sol_x[i];
            }
        }
        /* we now repeat the same process, but instead of solving df/dt = -v*df/dx, we are solving df/dt = C*d^2 f/dv^2 - F(x) df/dv. We set vmat to be zero everywhere at the start of the loop to avoid coefficients carrying over from previous loops. We also build rhs_v from f_temp, seeing the connection between the solution of both pde's which was discussed at the start of the loop. To make sure that we retain correct boundary conditions for x throughout the v-loop, we also strictly enforce them at the start. If they were to be enforced at the end, then we would see a gradual 'diffusion' of this error in our results, but not at the boundary. We avoid the values for which x-boundary conditions apply in our later v-loop. We now also see the use of f_temp, if it was f used, we would lose our initial value of f, and hence mess up a future loop at some stage. */
        for(long i = 0; i < N_x; i++) {
            for(long j = 0; j < N_v; j++) {
                setv(&vmat, j, j, 0);
                if (j > 0) {
                    setv(&vmat, j, j - 1, 0);
                }
                if (j < N_v - 1) {
                    setv(&vmat, j, j + 1, 0);
                }
                rhs_v[j] = f_temp[i + j*N_x];
            }
            
            if (i == 0) {
                for(long j = 0; j < N_v; j++){
                    if (v[j] >= 0) {
                        setv(&vmat, j, j, 1);
                        rhs_v[j] = exp(-v[j] * v[j]);
                    }
                }
            }
            if (i == N_x - 1) {
                for (long j = 0; j < N_v; j++) {
                    if (v[j] <= 0) {
                        setv(&vmat, j, j, 1);
                        rhs_v[j] = exp(-v[j]*v[j]);
                    }
                }
            }
            /* we apply the boundary condition for v to our new loop here. We applied this same BC at the start to f, but it also need to apply to sol_v, vmat and rhs_v. This same trick was used above (set RHS to what you want, and set LHS to be 1 on diag, 0 elsewhere in that row). It also means that we only have to deal with the interior components of any matrix below.*/
            setv(&vmat, 0, 0, 1);
            setv(&vmat, N_v - 1, N_v - 1, 1);
            setv(&vmat, 0, 1, 0);
            setv(&vmat, N_v - 1, N_v - 2, 0);
            rhs_v[0] = 0;
            rhs_v[N_v - 1] = 0;
/* We continue analogously to the other PDE. Using the fact that v_dot is defined as F, we basically swap x for v, and then v for F in week 6 last slide. We use standard finite difference method to deal with diffusion term, and rearrange to get whats seen below. We also avoid the terms for which the old boundary conditions apply, to avoid violating the boundary conditions for x in the v-loop.*/
            for (long j = 1; j < N_v - 1; j++) {
/*skipping terms for reasons previously mentioned*/
                if((i == 0 && v[j] >=0) || (i == N_x - 1 && v[j] <=0)) {
                    continue;
                }
                if (F[i] == 0.0) {
                    setv(&vmat, j, j, 1 + (2*C *dt)/(dv*dv));
                    setv(&vmat, j, j + 1, -(C * dt)/(dv*dv));
                    setv(&vmat, j, j - 1, -(C*dt)/(dv*dv));
                }
                if (F[i] > 0.0) {
                    setv(&vmat, j, j, 1 + (2*C * dt)/(dv*dv) + (F[i]*dt)/dv);
                    setv(&vmat, j, j + 1, -(C * dt)/(dv*dv));
                    setv(&vmat, j, j - 1, -(C*dt) /(dv*dv) - (F[i]*dt)/dv);
                    
                }
                if (F[i] < 0.0) {
                    setv(&vmat, j, j, 1 + (2*C*dt/(dv*dv)) - (F[i] * dt)/dv);
                    setv(&vmat, j, j+1, -(C*dt)/(dv*dv) + (F[i]*dt)/dv);
                    setv(&vmat, j, j-1, -(C*dt)/(dv*dv));
                }
            }
                /* We now solve for sol_v, which contains all info about the other PDE for this timestep, and hence our actual solution, f, is just the solution of the second PDE, given that it contains the solution of the first PDE in its (i) index.*/
            int info2 = solve_Ax_eq_b(&vmat, sol_v, rhs_v);
            if (info2 != 0) {
                printf("Error when solving system in v-direction at time step %ld, x_index %ld\n", n, i);
            }
            for (long j = 0; j < N_v; j++) {
                f[i + N_x*j] = sol_v[j];
                }
/* We do not need to apply boundary conditions at the end of our loop since they were applied and unaltered throughout both x,v loops*/
        }
    }
/*Note above that t was never explicitly used in the time loop. This is because t never appears in any finite difference approximation, or any boundary condition. Even though only dt is used in the time loop, the system of (basically) solve PDE1 -> solve PDE2 -> get solution to both is iterated N_t times, and hence we are still 'advancing in time' in a physical sense. */
    
/* opening output.txt and checking nothing goes wrong, if it does, freeing all allocated memory as done numerous times above.*/
    FILE *foutput = fopen("output.txt", "w");
    if (foutput == NULL) {
        printf("Error when opening output file\n");
        finalise_band_mat(&xmat);
        finalise_band_mat(&vmat);
        free(f);
        free(x);
        free(v);
        free(rhs_v);
        free(rhs_x);
        free(sol_x);
        free(sol_v);
        free(f_temp);
        free(F);
        return 1;
    }
    
/* write output.txt (successfully opened) in correct order (we need every i in chronological order for each j first as shown in spec, hence the i loop must sit inside the j loop). We use a 12 decimal approximation for our answer, which perhaps could be overkill, but just in case.*/
    for (long j = 0; j < N_v; j++) {
        for (long i = 0; i < N_x; i++) {
            fprintf(foutput,"%.12e\n", f[i + j*N_x]);
        }
    }
    fclose(foutput);
    
/* Finalising band matrices and freeing all previously allocated memory */
    finalise_band_mat(&xmat);
    finalise_band_mat(&vmat);
    free(f);
    free(x);
    free(v);
    free(rhs_v);
    free(rhs_x);
    free(sol_x);
    free(sol_v);
    free(f_temp);
    free(F);
    
    return(0);
}
    
    

