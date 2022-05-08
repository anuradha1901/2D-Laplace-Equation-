// Name: Anuradha Agarwal
// Title: Solving LAPLACE 2D Numerically
// Date: 12/14/2021

// Required Headers
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#define PI 3.14159265359 // defining value of Pi
#define E 2.718281828459 // defining value of e

/*
 This code solves the numerical solution of Laplace's equation using the
 Finite Difference Method (FDM). This method converts the PDE into a set of linear
 algebraic equations. This code solves this system of equations using the Conjugate
 Gradient (CG) Method.
*/
// All the functions that are used
void LAPLACEWCG(int n, int m, double a, double b, double hx, double hy, double tol, double max1, double** U); // finite difference method
void CGUPDATE(int n, int m, double** UIJ, double** APK, double** R0, double** P0); // iterative solution
double ERROR_METRIC(int n, int m, double* A, int N, int option); // finds the norm
void Real_function(int n, int m, double hx, double hy, double** W); // analytical solution
void error_matrix(int n, int m, double** U, double** W, double* F); // error between the analytical and approximation
double BDYVAL(int option, double w); // boundary value
double f1(double x); // boundary 1
double f2(double x); // boundary 2
double g1(double x); // boumdary 3
double g2(double x); // boundary 4

int main(void) // main function
{
  // initializing all the variables
  // i and j are for iterations
  // n and m are the number of steps in x and y directions, respectively
  // hx and hy are step sizes; in this project, hx = hy for the sake of simplicity
  // tol is the tolerance
  // max1 is the maximum error
  // W is the matrix with the exact solution (2 - Dimensional)
  // U is the matrix with the approximate solution (2 - Dimensional)
  // F is the error matrix (1 - Dimensional) - to calculate the error
  double time_spent = 0.0;
  clock_t begin = clock(); // beginning the clock
  int i, j, n, m; // declaring the iterators i and j; initializing number of rows and columns
  int a = 1; // domain
  int b = 1; // domain
  double hx = 0.005; // step-size
  double hy = 0.005; // step-size
  double tol = 0.000000001; // tolerance of 1e-9
  double max1 = 100000000;  // max number of iterations
  double** W; // declaring the matrix for the analytical solution
  double** U; // declaring the matrix for the approximate solution
  double* F; // declaring the vector to compute the error
  FILE *fptr; // declaring the output file

  // Opening the output file
  fptr = fopen("approximation.txt","w");
  // Calculating the number of steps
  n = (a/hy) + 1;
  m = (b/hx) + 1;

  // Allocating space for matrix U which is a matrix of size n by m
  U = (double **) malloc(n * sizeof(double*)); // n rows
  for (i = 0; i < n; i++) {
	   U[i] = (double *) malloc(m * sizeof(double)); // each row has m columns
  }
  if (!U)
  {
	// Memory allocation failed
     printf("Memory allocation of U array failed \n");
	   exit(EXIT_FAILURE);
  }

  // allocating space for matrix W  which is a matrix of size n by m
  W = (double **) malloc(n * sizeof(double*)); // n rows
  for (i = 0; i < n; i++) {
	   W[i] = (double *) malloc(m * sizeof(double)); // each row has m columns
  }
  if (!W)
  {
	// Memory allocation failed
	  printf("Memory allocation of W array failed \n");
	  exit(EXIT_FAILURE);
  }

  // allocating space for vector F  which is of size m*n
  F = (double *) malloc((m*n) * sizeof(double)); // n*m rows
  if (!F)
  {
	// Memory allocation failed
	printf("Memory allocation of F vector failed \n");
	exit(EXIT_FAILURE);
  }

  // calling functions
  // LAPLACEWCG to output the matrix approximation
  // Real_function to output the exact solution
  // error_matrix to output the error
  LAPLACEWCG(n, m, a, b, hx, hy, tol, max1, U); // The appromiate solution gets saved in the U matrix
  Real_function(n,m,hx,hy,W); // The analytical solution gets saves in the W matrix
  error_matrix(n, m, U, W, F); // Uses the previous two matrices (analytical and approximation) to find the error
  // The vector of differences gets saved in the F vector
  // prints the error in the command window

  // Printing out the approximation to an output file
  for (i = 0; i < n; i++){ // Iterating over rows
      for (j = 0; j < m; j ++){ // iterating over columns
        fprintf(fptr, "%0.15lf\t", U[i][j]); // printing every element
      }
      fprintf(fptr, "\n");
    }

  // Closing the file
  fclose(fptr);
  // Freeing the space
  free(W);
  free(F);
  free(U);

  clock_t end = clock(); // end clock
  time_spent += (double)(end - begin) / CLOCKS_PER_SEC; // calculating run time
  printf("The elapsed time is %lf seconds", time_spent);
  return 0;
}

// FUNCTION LAPLACEWG
void LAPLACEWCG(int n, int m, double a, double b, double hx, double hy, double tol, double max1, double** U){
  // This funtion takes in n, m which are the number of rows and columns, length of the domain a,b , the
  // step size in both the x and y direction hx, hy; error tolerance tol; maximum number of iterations max1
  // as the inputs. This procedure applied the Forward difference method to the Laplaces equation and Uses
  // the method of Conjugate gradient method to find the problem solution.
  // Matrix U here is the required approximation.
  int i,j,cnt; // declaring the iteratives and the counter
  double err,rx,ry,ave;
  double* X; // declaring vector X
  double* Y; // declaring  vector Y
  double** R0; // declaring matrix R0
  double** P0; // declaring matrix P0
  double** APK; // declaring matrix APK
  double* RV; // declaring vector RV
  cnt = 0; // initializing the counter

  // allocating space for matrix R0 which is a matrix of size n by m
  R0 = (double **) malloc(n * sizeof(double*));
  for (i = 0; i < n; i++) {
    R0[i] = (double *) malloc(m * sizeof(double));
  }
  if (!R0){
  // Memory allocation failed
    printf("Memory allocation of R0 array failed \n");
  exit(EXIT_FAILURE);
  }

// allocating space for matrix P0 of size n by m
  P0 = (double **) malloc(n * sizeof(double*));
  for (i = 0; i < n; i++) {
      P0[i] = (double *) malloc(m * sizeof(double));
  }
  if (!P0){
  // Memory allocation failed
    printf("Memory allocation of P0 array failed \n");
  exit(EXIT_FAILURE);
  }

// allocating space for matrix APK which is a matrix of size n by m
  APK = (double **) malloc(n * sizeof(double*));
  for (i = 0; i < n; i++) {
    APK[i] = (double *) malloc(m * sizeof(double));
  }
  if (!APK){
  // Memory allocation failed
      printf("Memory allocation of APK array failed \n");
  exit(EXIT_FAILURE);
  }

// allocating space for vector X of size m
  X = (double *) malloc((m) * sizeof(double));
  if (!X){
	// Memory allocation failed
	   printf("Memory allocation of X vector failed \n");
	exit(EXIT_FAILURE);
  }

// allocating space for vector Y of size n
  Y = (double *) malloc((n) * sizeof(double));
  if (!Y){
  // Memory allocation failed
    printf("Memory allocation of Y vector failed \n");
  exit(EXIT_FAILURE);
  }

// allocating space for vector RV of size n*m
  RV = (double *) malloc((m*n) * sizeof(double));
  if (!RV){
  // Memory allocation failed
    printf("Memory allocation of RV vector failed \n");
  exit(EXIT_FAILURE);
  }


  for (i = 0; i < n; i++){ // iterating over rows
    for (j = 0; j < m; j ++){ // iterating over columns
      U[i][j] = 1.0; // intializing our matrix with ones
    }
  }

  for (j = 0; j < m; j++){ // iterating over rows
    X[j] = j*hx; // The vector x which forms the x axis coordinates
  }

  for (j = 0; j < n; j++){ // iterating over rows
    Y[j] = b - j*hy; // The vector y which forms the y axis coordinates
  }

  for (i = 0; i < n; i++){ // iterating over rows
    for (j = 0; j < m; j ++){ // iterating over columns
      R0[i][j] = 0.0; // initializing R0 to be a matrix of size m*n with zeroes
      RV[(n*i)+j] = R0[i][j]; // converting the matrix into a vector
    }
  }

  for (i = 0; i < n; i++){ // iterating over rows
    for (j = 0; j < m; j ++){ // iterating over columns
      P0[i][j] = 0.0; // initializing P0 to be a matrix of size m*n with zeroes
    }
  }

  for (i = 0; i < n; i++){ // iterating over rows
    for (j = 0; j < m; j ++){ // iterating over columns
      APK[i][j] = 0.0; // initializing APK to be a matrix of size m*n with zeroes
    }
  }

  rx = 1/(hx*hx); // assigning rx
  ry = 1/(hy*hy); // assigning ry
  ave = (a*(BDYVAL(1,0)+BDYVAL(2,0)) + b*(BDYVAL(3,0) + BDYVAL(4,0)))/(2*a + 2*b); // initial condition

  for(i = 0; i < n; i++){ // iterating over rows
    for(j = 0; j < m; j++){ // iterating over columns
      U[i][j] = ave * U[i][j];  // updating the matrix that stores the approximation
      // getting the initial condition
      }
    }

// This loop updates the boundary conditions on left and right sides of the matrix
  for (i = 0; i < n; i++){ // iterating over rows
    U[i][0] = BDYVAL(3, Y[i]); // Boundary 1
    U[i][m-1] = BDYVAL(4, Y[i]); // Boundary 2
  }

// This loop updates the boundary conditions on the top and bottom side of the matrix
  for (j = 0; j < m; j++){ // iterating over columns
    U[0][j] = BDYVAL(1, X[j]); // Boundary 3
    U[n-1][j] = BDYVAL(2, X[j]); // Boundary 4
  }

  U[0][0] = (U[0][1] + U[1][0])/2; // Averaging the top-left corner value
  U[0][m-1] = (U[0][m-2] + U[1][m-1])/2; // Averaging the top-right corner value
  U[n-1][0] = (U[n-2][0] + U[n-1][1])/2; // Averaging the bottom left corner value
  U[n-1][m-1] = (U[n-2][m-1] + U[n-1][m-2])/2; // Averaging the bottom-right corner value

  for (j = 1; j < m-1; j++){ // iterating over the rows
    for (i = 1; i < n-1; i++){ // iterating over the columns
      // going through every element of R0 and updating it
      R0[i][j] = ((rx*U[i][j+1]) + (rx*U[i][j-1]) + (ry*U[i+1][j]) + (ry*U[i-1][j]) - (2*(rx+ry)*U[i][j])); // updating the R0 matrix
      RV[(n*i)+j] = R0[i][j]; // converting the matrix into a vector
    }
  }

  for (j = 0; j < m; j++){ // iterating over rows
    for (i = 0; i < n; i++){ // iterating over columns
      P0[i][j] = R0[i][j]; // Assigning the values of P0 as the elements of R0
    }
  }

 err = ERROR_METRIC(n, m, RV, m*n, 2); // finding the norm of the vector
 // this vector is the same as the matrix R0
 //printf("The norm of the matrix is: %lf \n", err);
 // The following part or code runs when the error is more than the tolerance
 // and the our iterations that is kept track by cnt (counter) does not exceed the
 // maximum number of iterations
  while ((err > tol) && (cnt <= max1)){ // condition
    for ( j = 1; j < m-1; j++){ // iterating over rows
      for (i = 1; i < n-1; i++){ // iterating over columns
        if ( j == 1){ // going step by step - if column 2
          if (i == 1){ // if row 2
            APK[i][j] = (-rx *P0[i][j+1]) - (ry*P0[i+1][j]) + 2*(rx + ry)*(P0[i][j]);
          } // updating APK accordingly
          else if (i == n - 2){ // if column - 1
            APK[i][j] = (-rx *P0[i][j+1]) - (ry*P0[i-1][j])+ 2*(rx + ry)*(P0[i][j]);
          } // updating APK accordingly
          else {
            APK[i][j] = (-rx *P0[i][j+1]) - (ry*P0[i+1][j]) - (ry*P0[i-1][j]) + 2*(rx + ry)*(P0[i][j]);
          } // updating APK accordingly
        }

        else if (j == m-2){ // j not 1 but now m - 2
          if (i == 1){ // if row 2
            APK[i][j] = (-rx *P0[i][j-1]) - (ry*P0[i+1][j]) + 2*(rx + ry)*(P0[i][j]);
          } // updating APK accordingly
          else if (i == n-2){ // if row n -1
            APK[i][j] = (-rx *P0[i][j-1]) - (ry*P0[i-1][j]) + 2*(rx + ry)*(P0[i][j]);
          } // updating APK accordingly
          else {
            APK[i][j] = (-rx *P0[i][j-1]) - (ry*P0[i+1][j]) - (ry*P0[i-1][j]) + 2*(rx + ry)*(P0[i][j]);
          } // updating APK accordingly
        }

        else if (i == n - 2){ // if column m -1 and row n -1
          APK[i][j] = (-rx *P0[i][j+1]) - (ry*P0[i][j-1]) - (ry*P0[i-1][j]) + 2*(rx + ry)*(P0[i][j]);
        } // updating APK accordingly

        else if (i == 1){ // if column m -1 and row 2
          APK[i][j] = (-rx *P0[i][j+1]) - (ry*P0[i][j-1]) - (ry*P0[i+1][j]) + 2*(rx + ry)*(P0[i][j]);
        } // updating APK accordingly

        else {
          APK[i][j] = (-rx *P0[i][j+1]) - (rx*P0[i][j-1]) - (ry*P0[i+1][j])  - (ry*P0[i-1][j]) + 2*(rx + ry)*(P0[i][j]);
        } // updating APK accordingly
      } // end of for i (iterations)
    } // end of for j (iterations)
    CGUPDATE(n,m, U, APK, R0, P0); // updating all the matrices for every iteration
    for (i = 0; i < n; i++){ // iterating over rows
      for (j = 0; j < m; j++){ // iterating over columns
        RV[(n*i)+j] = R0[i][j]; // converting the updated R0 matrix to vector RV, to calculate the norm again
      }
    }
    err = ERROR_METRIC(n , m, RV, m*n, 2); // finding the norm and looping around
	//printf("The norm of the matrix is: %lf \n", err);
    cnt = cnt + 1; // updating the counter for the number of iterations
  } // end of while
  if (cnt >= max1){ // The procedure prints out an error message when the counter exceeds the max number of iterations
    printf("Maximum number of iterations reached");
  }
  printf("The number of iterations is: %d \n", cnt); // printing the final number of iterations

// The bottom print statements can be used to verify the solution
  /*
  printf("\n");
  printf("The approximate solution is below \n");
  for (i = 0; i < n; i++){
    for (j = 0; j < m; j ++){
      printf("%lf\t", U[i][j]);
    }
    printf("\n");

  }
  */

  //printf("The norm is %lf \n", err);
  // freeing all the matrices
  free(X);
  free(Y);
  free(R0);
  free(RV);
  free(P0);
  free(APK);
} // end of the function

// FUNCTION CGUPDATE
void CGUPDATE(int n, int m, double** UIJ, double** APK,double** R0,double** P0){
  // This procedure is used by LAPLACEWG to update the iterative solution of the
  // PDE using the CG method.
  /*
  This function takes n, m which are the number of rows and columns of the matrices;
  4 matrices of the size n*m and updates them.
  */
  int i,j; // declaring the iterators
  double nak, dak, ak, nbk, bk; // declaring constant
  double** RK; // declaring a matrix of size n*m
  double** PK; // declaring a matrix of size m*n

  // allocating space for matrix RK which is a matrix of size n by m
  RK = (double **) malloc(n * sizeof(double*)); // rows
  for (i = 0; i < n; i++) { // columns
    RK[i] = (double *) malloc(m * sizeof(double));
  }
  if (!RK)
  {
  // Memory allocation failed
    printf("Memory allocation of W array failed \n");
    exit(EXIT_FAILURE);
  }

  // allocating space for vector PK which is a matirx of size n by m
  PK = (double **) malloc(n * sizeof(double*)); // rows
  for (i = 0; i < n; i++) {
    PK[i] = (double *) malloc(m * sizeof(double)); // columns
  }
  if (!PK)
  {
  // Memory allocation failed
    printf("Memory allocation of W array failed \n");
    exit(EXIT_FAILURE);
  }

  nak = 0.0; // initializing n alpha k to 0
  for (i = 0; i < n; i ++){ // iterating over rows
    for (j = 0; j < m; j ++){ // iterating over columns
      nak = nak + ((R0[i][j])*(R0[i][j])); // summation of dot product of R0 * R0
    } // assigning the scalar product to n alpha k
  }
  dak = 0.0; // initializing d alpha k to 0
  for (i = 0; i < n; i ++){ // iterating over rows
    for (j = 0; j < m; j ++){ // iterating over columns
      dak = dak + ((P0[i][j])*(APK[i][j])); // summation of dot product of P0 * APK
    } // assigning the scalar product to d alpha k
  }
  ak = nak/dak; // dividing nak and dak; assigning to ak

    for (i = 0; i < n; i++){ // iterating over rows
      for (j = 0; j < m; j++) { // iterating over columns
        UIJ[i][j] = UIJ[i][j] + ak * P0[i][j]; // updating every element of matrix U
        // with U = U + alpha k * P0
      }
    }

    for (i = 0; i < n; i++){ // iterating over rows
      for (j = 0; j < m; j++) { // iterating over columns
        RK[i][j] = R0[i][j] - ak * APK[i][j]; // updating every element of matrix Rk
        // with Rk = R0 - alpha k * APK
      }
    }

    nbk = 0.0; // initializing n beta k to 0

    for (i = 0; i < n; i ++){ // iterating over rows
      for (j = 0; j < m; j ++){ // iterating over columns
        nbk = nbk + ((RK[i][j])*(RK[i][j])); // dot product of Rk with itself
      } // assigning the dot product to n beta k
    }

    bk = nbk/nak; // assigning beta k as n beta k / n alpha k

    for (i = 0; i < n; i++){ // iterating over rows
      for (j = 0; j < m; j++) { // iterating over columns
        PK[i][j] = RK[i][j] + bk * P0[i][j]; // updating every element of matrix PK
        // with RK + beta k * P0
      }
    }

    for (i = 0; i < n; i++){ // iterating over rows
      for (j = 0; j < m; j++) { // iterating over columns
        P0[i][j] = PK[i][j]; // assigning every element of PK to P0
      } // This loop updates P0; P0 is used in calculations in LAPLACEWG
    }

    for (i = 0; i < n; i++){ // iterating over rows
      for (j = 0; j < m; j++) { // iterating over columns
        R0[i][j] = RK[i][j]; // assigning every element of RK to R0
      }
    }
  // freeing the space
  free(RK);
  free(PK);
}

// FUNCTION ERROR_METRIC
double ERROR_METRIC(int n, int m, double* A, int N, int option)
/*
This function takes in n, m which are the number of rows and columns of the
matrix; a vector A which is of size m*n, which is denoted as N and the option
for different types of norms as the inputs. The output of this function is
a double depicting the norm of the vector corresponding to the option.

This function essentially finds the norm of the vector A of size N (m*n in our case)
provided by the input array. This function is used by LAPLACEWCG.

L1 norm is the Taxi cab norm
L2 norm is the euclidean norm
Linf norm is the infinity norm

*/
{
	int i, j; // declaring the iterators
	double norm; // declaring variable norm which saves out result

	norm = 0; // initializing norm to 0
	if (option == 1) {  // option 1: L1 norm

		for (i = 0; i < N; i++){ // iterating over rows
			norm = norm + fabs(A[i]); // summation of the absolute values
		}
	}
	else if (option == 2) { // option 2: L2 norm

		for (i = 0; i < N; i++){ // iterating over rows
			norm = norm + (A[i] * A[i]); // summation of the squares of the norm
		}
		norm = sqrt(norm); // square root of the sum of squares
	}
	else if (option == 3) { // option 3: L inf norm

		for (i = 0; i < N; i++) { // iterating over rows
			if (norm < fabs(A[i]))
				norm = fabs(A[i]); // absolute maximum
		}
	}

	return norm; // this function returns the value of the norm - Type double
}

// FUNCTION Real_function
void Real_function(int n, int m, double hx, double hy, double** W)
// Anuradha Agarwal
// Date: 12/14/2021
/*
This function takes in n, m which are the size of the matrix, hx, hy which are
are the step sizes as inputs. It outputs a matrix of size W which is the matrix
of the exact or analytical solution
*/
{
    int i, j;  // initializing the iterators
    for (i = 0; i < n; i++){ // iterating over every row
      for (j = 0; j < m; j++){ // interating over columns
        W[i][j] = (pow(E, (PI*hy*j)))*cos(PI*hx*(n-1-i)); // exact solution of the laplace
      }
    }
/* Printing results for verification
    printf("The real solution is below \n");
    for (i = 0; i < n; i++){
      for (j = 0; j < m; j++){
        printf("%lf\t", W[i][j]);
      }
      printf("\n");
    }
*/
}

// FUNCTION error_matrix
void error_matrix(int n, int m, double** U, double** W, double* F)
// Anuradha Agarwal
// Date: 12/14/2021

/*
This function finds the difference between the analytical solution and the
approximation. After finding the difference this function transforms the matrix
which is of size [n]*[m] to a vector of size n*m.
Inside this function another function ERROR_METRIC is called which finds the Norms
of the matrix, with three options : L1 norm, L2 norm, L inf norm.
Here the input parameters are n and m which are the size of the matrix, a double
matrix U which is the analytical solution, a double matrix W which is the approximate
solution.
*/
{
  int i,j;  // declaring i and j
  double** G;

  // allocating space for matrix G which is a matrix of size n by m
  // Matrix G finds the difference between the analytical and approximate solutions
  G = (double **) malloc((n-2) * sizeof(double*));
  for (i = 0; i < n; i++) {
    G[i] = (double *) malloc((m-2) * sizeof(double));
  }
  if (!G){
  // Memory allocation failed
    printf("Memory allocation of W array failed \n");
    exit(EXIT_FAILURE);
  }

  int k;       // declaring the iterator
  k = 0;       // initializing the iterator
  // calculating the error only with the interior coordinates
  for (i = 1; i < n-1; i++){ // looping through every row
    for (j = 1; j < m-1; j++){ // looping through every column
      G[i][j] = U[i][j] - W[i][j];   // finding the difference between two matrices
      F[k] = G[i][j];   // converting the matrix into a vector
      k++;
    }
  }

// Printing results for verification
/*
    for (i = 1; i < n-1; i++){
      for (j = 1; j < m-1; j++){
        printf("%0.10lf\t", G[i][j]);
      }
      printf("\n");
    }
  printf("Norms \n");
	for (i = 0; i < (n-2)*(m-2); i++){
        printf("%lf\n", F[i]);
      }
*/
    double g; // initializing the variable g of type double
    g = ERROR_METRIC(n, m, F, (m-1)*(n-1), 3); // calling the ERROR_METRIC to calculate the L inf norm
    printf("The error is %lf \n",g);
    free(G); // freeing te space
}

// FUNCTION BDYVAL
double BDYVAL(int option, double w)
// Anuradha Agarwal
// Date: 12/14/2021
// The inputs of the function are the option which denotes the function and double w
// which is a constant double
// This procedure is used to provide the algorithm with the boundary values. It has
// input parameters option and w. Parameter w is the grid point coordinate value x
// or y chosen by LAPLACEWCG. Parameter option is used to choose the boundary value
// function
{
    double func;
    if (option == 1){
      // option 1 returns f2(w)
      func = f2(w);
    }
    else if (option == 2){
      // option 2 returns f1(w)
      func = f1(w);
    }
    else if (option == 3){
      // option 3 returns g1(w)
      func = g1(w);
    }
    else if (option ==4){
      // option 4 returns g2(w)
      func = g2(w);
    }
    return func;
}

// All the functions below are the function depicting the boundary value conditions
double f1(double x)
// Anuradha Agarwal
// Date: 12/14/2021
// This function is the boundary value condition for y = 0 and x in (0,1)
// This takes in double x which is the variable and returns a double value
{
	double result = pow(E, PI*x);	// function
  return result;
}
double f2(double x)
// Anuradha Agarwal
// Date: 12/14/2021
// This function is the boundary value condition for y = 1 and x in (0,1)
// This takes in double x which is the variable and returns a double value
{
  double result = -1*pow(E, PI*x);	// function
  return result;
}
double g1(double y)
// Anuradha Agarwal
// Date: 12/14/2021
// This function is the boundary value condition for x = 0 and y in (0,1)
// This takes in double y which is the variable and returns a double value
{
  double result = cos(PI*y);	// function
  return result;
}
double g2(double y)
// Anuradha Agarwal
// Date: 12/14/2021
// This function is the boundary value condition for x = 1 and y in (0,1)
// This takes in double y which is the variable and returns a double value
{
  double result = (pow(E, PI))*cos(PI*y);	 // function
  return result;
}
