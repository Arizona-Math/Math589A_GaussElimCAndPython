/*----------------------------------------------------------------
* File:     gauss_solve.c
*----------------------------------------------------------------
*
* Author:   Marek Rychlik (rychlik@arizona.edu)
* Date:     Sun Sep 22 15:40:29 2024
* Copying:  (C) Marek Rychlik, 2020. All rights reserved.
*
*----------------------------------------------------------------*/

#include <math.h>
#include <assert.h>
#include <stdio.h>
#include <string.h>

#include "gauss_solve.h"
#include "helpers.h"


void gauss_solve_in_place(const int n, double A[n][n], double b[n])
{
  for(int k = 0; k < n; ++k) {
    for(int i = k+1; i < n; ++i) {
      /* Store the multiplier into A[i][k] as it would become 0 and be
	 useless */
      A[i][k] /= A[k][k];
      for( int j = k+1; j < n; ++j) {
	A[i][j] -= A[i][k] * A[k][j];
      }
      b[i] -= A[i][k] * b[k];
    }
  } /* End of Gaussian elimination, start back-substitution. */
  for(int i = n-1; i >= 0; --i) {
    for(int j = i+1; j<n; ++j) {
      b[i] -= A[i][j] * b[j];
    }
    b[i] /= A[i][i];
  } /* End of back-substitution. */
}

void lu_in_place(const int n, double A[n][n])
{
  for(int k = 0; k < n; ++k) {
    for(int i = k; i < n; ++i) {
      for(int j=0; j<k; ++j) {
	/* U[k][i] -= L[k][j] * U[j][i] */
	A[k][i] -=  A[k][j] * A[j][i]; 
      }
    }
    for(int i = k+1; i<n; ++i) {
      for(int j=0; j<k; ++j) {
	/* L[i][k] -= A[i][k] * U[j][k] */
	A[i][k] -= A[i][j]*A[j][k]; 
      }
      /* L[i][k] /= U[k][k] */
      A[i][k] /= A[k][k];	
    }
  }
}

/* Simply finds A = L.U from the packed form A = L+U */
void lu_in_place_reconstruct(int n, double A[n][n])
{
  for(int k = n-1; k >= 0; --k) {
    for(int i = k+1; i<n; ++i) {
      A[i][k] *= A[k][k];
      for(int j=0; j<k; ++j) {
	A[i][k] += A[i][j]*A[j][k];
      }
    }
    for(int i = k; i < n; ++i) {
      for(int j=0; j<k; ++j) {
	A[k][i] +=  A[k][j] * A[j][i];
      }
    }
  }
}

void plu_in_place(const int n, double A[n][n], int P[n])
{
  /* (Re)Initialize P */
  for(int k = 0; k < n; ++k) P[k] = k;

  for(int k = 0; k < n; ++k) {
    /* Find the pivot element */
    double pivot = A[k][k]; int r = k;

    for(int i = k + 1; i < n; ++i) {
      if( fabs(pivot) < fabs(A[i][k]) ) {
	pivot = A[i][k];
	r = i;
      }
    }

    /* row exchange */
    if(r != k)  {
      SWAP(P[k], P[r], int);

      /* Swap rows k and r in L, only cols 1 to k-1 */
      for(int j = 0; j < k; ++j) {
	// SWAP(L[k][j], L[r][j], double);
	SWAP(A[k][j], A[r][j], double);
      }
      /* In this version, we swap the rest of the row, too */
      for(int j = k; j < n; ++j) {
	// SWAP(U[k][j], U[r][j], double);
	SWAP(A[k][j], A[r][j], double);
      }
    }

    for(int i = k+1; i < n; ++i) {
      // L[i][k] = A[i][k] / A[k][k];
      A[i][k] /= A[k][k];
      for(int j = k+1; j < n; ++j) {
	// A[i][j] -=  L[i][k] * A[k][j]; 
	A[i][j] -=  A[i][k] * A[k][j]; 
      }
    }
  }
}
