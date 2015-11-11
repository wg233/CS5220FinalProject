#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define INPUT "model_def.txt" // Map of path to input file
#define OUTPUT "results.txt" // Map of path to output file

/*
 3D_geom_nonlin_truss.c is a gemetrically nonlinear FE code for analyzing truss structures>
 Copyright (C) 2009  Christopher J. Earls, Cornell University
 
 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.
 
CEE 7790 Truss Code (8/11/2009)
Finite element analysis for space trusses
Geometrically-nonlinear; Newton-Raphson method
Materially-linear
Total Lagrangian formulation

Input data:
    units: kip, inch, radian
    enter number of elements, number of joints (in main) - ne, nj
    enter member incidences (in struc) - minc[1,i], minc[2,i]; i = 1 to ne
    enter joint constraint(s) (in struc) - jnum, jdir, end = 0,0
    enter joint coordinates (in prop) - x[1,j], x[2,j], x[3,j]; j = 1 to nj
    enter member properties (in prop) - area[i], emod[i]; i = 1 to ne
    enter joint load(s) (in load) - jnum, jdir, force, end = 0,0,0
    enter load proportionality parameters:
        max lambda, initial lambda, increment of lambda (in main) - qimax, qi, dqi
    enter max number of iterations (in main) - itemax
    enter tolerances on out-of-balance forces and energy (in main) - tolfor, tolener
*/

// Define primary variables - number of elements, joints, and equations.
int ne, nj, neq;

// Define file pointers.
FILE *ifp; // Input file
FILE *ofp; // Output file

/* Define function prototypes - convention: "p" precedes the name of the variable in
   main being pointed to in the function */

// Model definition functions:
int struc (int *pjcode, int *pminc);
void codes (int *pmcode, int *pjcode, int *pminc);
void prop (double *px, double *parea, double *pemod, double *peleng, double *pc1,
    double *pc2, double *pc3, int *pminc);
int load (double *pq, int *pjcode);

// Stiffness functions:
void skylin (int *pkht, int *pmaxa, int *pmcode, int *plss);
void stiff (double *pss, double *parea, double *pemod, double *peleng, double *pc1,
    double *pc2, double *pc3, double *pelong, int *pmaxa, int *pmcode, int *plss);

// Internal force vector function:
void forces (double *pf, double *parea, double *pemod, double *pc1, double *pc2,
    double *pc3, double *pelong, double *peleng, int *pmcode);

// Solver functions:
int solve (double *pss, double *pr, double *pdd, int *pmaxa);
void test (double *pf, double *pfp, double *pqtot, double *pdd, double *pfpi,
    double *pintener1, int *pinconv, int *pneq, double *ptolfor, double *ptolener);

// Miscellaneous functions:
void updatc (double *px, double *pdd, double *pc1, double *pc2, double *pc3,
    double *pelong, double *peleng, double *pdeflen, int *pminc, int *pjcode);

int main (void)
{
    // Open I/O for business!
    ifp = fopen(INPUT, "r"); // Open input file for reading
    ofp = fopen(OUTPUT, "w"); // Open output file for writing

    // Verify that I/O files are open
    if (ifp != 0) {
         printf("Input file is open\n");
    } else {
        printf("***ERROR*** Unable to open the input file\n");
        getchar();
        return 0;
    }
    if (ofp != 0) {
        printf("Output file is open\n");
    } else {
        printf("***ERROR*** Unable to open the output file\n");
        getchar();
        return 0;
    }

    // Read in number of elements, joints from input file
    fscanf(ifp, "%d,%d\n", &ne, &nj);
    fprintf(ofp, "Control Variables:\n\tNumber of Elements: %d\n", ne);
    fprintf(ofp, "\tNumber of Joints: %d\n\n", nj);

    /*
    Associate pointers with base address of arrays in function calls: parea = area;
        Note: in this context "&" can be omitted to obtain "&area[0]"
        Note: previous remark only works for 1D arrays: px=&x[0][0]
    */

    // Define secondary variables which DO NOT depend upon neq
    double x[3][nj];							// Current joint coordinates
    int mcode[6][ne], jcode[3][nj], minc[2][ne];// Member incidence and constraint data
    double emod[ne];							// Element material properties
	double eleng[ne], deflen[ne], elong[ne];	// Element length, deformed length, and elongation
    double area[ne];							// Element cross-sectional properties
    double c1[ne], c2[ne], c3[ne];				// Direction cosines
    double qi, dqi;								// Current and incremental load proportionality factor
    double qimax;								// Maximum allowable load proportionality factor
    
	// Convergence parameters
    double tolfor, tolener;						// Tolerances on force and energy
    double intener1;							// Internal energy from first equilibrium iteration
    int inconv;									// Flag for convergence test
    int itecnt, itemax;							// Iteration counter and max number of iterations
    int errchk;									// Error check variable on user defined functions
    int i;										// Counter variable

    // Pass control to struc function
    errchk = struc(&jcode[0][0], &minc[0][0]);

    // Terminate program if errors encountered
    if (errchk == 1) {
        fprintf(ofp, "\n\nSolution failed\n");
        printf("Solution failed, see output file\n");

        // Close the I/O
        if (fclose(ifp) != 0) {
            printf("***ERROR*** Unable to close the input file\n");
        } else {
            printf("Input file is closed\n");
        }
        if (fclose(ofp) != 0) {
            printf("***ERROR*** Unable to close the output file\n");
        } else {
            printf("Output file is closed\n");
        }
        getchar();
        return 0;
    }

    // Pass control to codes function
    codes (&mcode[0][0], &jcode[0][0], &minc[0][0]);

    // Define secondary variables which DO depend upon neq
    double q[neq];						// Reference load vector
    double qtot[neq];					// Total load vector, i.e. qi * q[neq]
    double d[neq], dd[neq];				// Total and incremental nodal displacement vectors
    double f[neq];						// Internal force vector
    double fp[neq];						// Internal force vector from previous load increment
    double fpi[neq];					// Internal force vector from previous iteration
    double r[neq];						// Residual force vector, i.e. qtot - f
    int maxa[neq + 1], kht[neq], lss;	// Skyline storage parameters for stiffness matrix

    // Pass control to skylin function
    skylin (kht, maxa, &mcode[0][0], &lss);

    // Define secondary variable which depends upon lss (the size of the stiffness marix as a 1D array)
    double ss[lss];						// Tangent stiffness matrix stored as an array

    // Pass control to prop function
    prop (&x[0][0], area, emod, eleng, c1, c2, c3, &minc[0][0]);

    // Pass control to load function
    errchk = load (q, &jcode[0][0]);

    // Terminate program if errors encountered
    if (errchk == 1) {
        fprintf(ofp, "\n\nSolution failed\n");
        printf("Solution failed, see output file\n");

        // Close the I/O
        if (fclose(ifp) != 0) {
            printf("***ERROR*** Unable to close the input file\n");
        } else {
            printf("Input file is closed\n");
        }
        if (fclose(ofp) != 0) {
            printf("***ERROR*** Unable to close the output file\n");
        } else {
            printf("Output file is closed\n");
        }
        getchar();
        return 0;
    }

    // Read in solver parameters from input file
    fscanf(ifp, "%lf,%lf,%lf\n", &qimax, &qi, &dqi);
    fscanf(ifp, "%d\n", &itemax);
    fscanf(ifp, "%lf,%lf\n", &tolfor, &tolener);

    // Print layout for output of results
    fprintf(ofp, "\nNonlinear Equilibrium Path:\n\tLambda\t");
    for (i = 0; i <= neq - 1; ++i) {
        fprintf(ofp, "\tDOF %d\t", i + 1);
    }

    // Initialize generalized nodal displacement and internal force vectors
    for (i = 0; i <= neq - 1; ++i) {
        d[i] = 0;
        f[i] = 0;
    }

    // Initialize element length and elongation variables
    for (i = 0; i <= ne - 1; ++i) {
        deflen[i] = eleng[i];
        elong[i] = 0;
    }

    /* Begin load incrementation; load will be incremented until load proportionality
       factor is equal to user specified maximum */
    do {
        /* Compute the generalized joint total load vector, store generalized internal
           force vector from previous configuration and re-initialize */
        for (i = 0; i <= neq - 1; ++i) {
            qtot[i] = q[i] * qi;
            fp[i] = f[i];
            f[i] = 0;
        }

        // Pass control to forces function
        forces (f, area, emod, c1, c2, c3, elong, eleng, &mcode[0][0]);

        itecnt = 1;  // Re-initialize iteration counter at the start of each increment

        /* Start of each equilibrium iteration within load increment; iterations will
           continue until convergence is reached or iteration count exceeds user
           specified maximum */
        do {
            // Compute residual force vector for use in evaluating displacement increment
            for (i = 0; i <= neq - 1; ++i) {
                r[i] = qtot[i] - f[i];
            }

            // Pass control to stiff function
            stiff (ss, area, emod, eleng, c1, c2, c3, elong, maxa, &mcode[0][0], &lss);

            // Solve the system for incremental displacements
            if (lss == 1) {
                // Carry out computation of incremental displacement directly for lss = 1
                dd[0] = r[0] / ss[0];
            } else {
                // Pass control to solve function
                errchk = solve (ss, r, dd, maxa);

                // Terminate program if errors encountered
                if (errchk == 1) {
                    fprintf(ofp, "\n\nSolution failed\n");
                    printf("Solution failed, see output file\n");

                    // Close the I/O
                    if (fclose(ifp) != 0) {
                        printf("***ERROR*** Unable to close the input file\n");
                    } else {
                        printf("Input file is closed\n");
                    }
                    if (fclose(ofp) != 0) {
                        printf("***ERROR*** Unable to close the output file\n");
                    } else {
                        printf("Output file is closed\n");
                    }
                    getchar();
                    return 0;
                }
            }

            /* Update generalized nodal displacement vector, store generalized internal
               force vector from previous iteration, and re-initialize generalized
               internal force vector */
            for (i = 0; i <= neq - 1; ++i) {
                d[i] += dd[i];
                fpi[i] = f[i];
                f[i] = 0;
            }

            // Pass control to forces function
            updatc (&x[0][0], dd, c1, c2, c3, elong, eleng, deflen, &minc[0][0],
                &jcode[0][0]);

            // Pass control to forces function
            forces (f, area, emod, c1, c2, c3, elong, eleng, &mcode[0][0]);

            if (itecnt == 1) {
                // Compute internal energy from first iteration
                intener1 = 0;
            	for (i = 0; i <= neq - 1; ++i) {
            		intener1 += dd[i] * (qtot[i] - fp[i]);
            	}
            }

            // Pass control to test function
            test (f, fp, qtot, dd, fpi, &intener1, &inconv, &neq, &tolfor, &tolener);

            itecnt ++; // Advance solution counter

        } while (inconv != 0 && itecnt <= itemax);

        // Store generalized internal force vector from current configuration
        for (i = 0; i <= neq - 1; ++i) {
            fp[i] = f[i];
        }

        if (inconv != 0) {
            fprintf(ofp, "\n\n***ERROR*** Convergence not reached within specified");
            fprintf(ofp, " number of iterations (INCONV = %d)", inconv);
            fprintf(ofp, "\n\nSolution failed\n");
            printf("Solution failed, see output file\n");

            // Close the I/O
            if (fclose(ifp) != 0) {
                printf("***ERROR*** Unable to close the input file\n");
            } else {
                printf("Input file is closed\n");
            }
            if (fclose(ofp) != 0) {
                printf("***ERROR*** Unable to close the output file\n");
            } else {
                printf("Output file is closed\n");
            }
            getchar();
            return 0;
        }

        // Output model response
        fprintf(ofp, "\n\t%lf", qi);
        for (i = 0; i <= neq - 1; ++i) {
            fprintf(ofp, "\t%lf", d[i]);
        }

        qi += dqi; // Increment load proportionality factor

    } while (qi <= qimax);

    if (qi >= qimax && inconv == 0) {
        fprintf(ofp, "\n\nSolution successful\n");
        printf("Solution successful\n");
    }

    // Close the I/O
    if (fclose(ifp) != 0) {
        printf("***ERROR*** Unable to close the input file\n");
    } else {
        printf("Input file is closed\n");
    }
    if (fclose(ofp) != 0) {
        printf("***ERROR*** Unable to close the output file\n");
    } else {
        printf("Output file is closed\n");
    }
    getchar();
    return 0;
}

/* This function reads in member incidences and joint constraints and checks for errors
   in the joint constraint input */
int struc (int *pjcode, int *pminc)
{
    int i, j, k; // Initialize function variables

    // Establish member incidences
    fprintf(ofp, "Member Incidences:\n\tMember\tA-End\tB-End\n");
    for (i = 0; i <= 2 * ne - 1; ++i) {
        *(pminc+i) = 0; // Zero-out minc
    }
    for (i = 0; i <= ne - 1; ++i) {
        fscanf(ifp, "%d,%d\n", &j, &k); // Read in A- and B-Ends from input file
        *(pminc+i) = j;
        *(pminc+i+ne) = k;
        fprintf(ofp, "\t%d\t%d\t%d\n", i + 1, *(pminc+i), *(pminc+i+ne));
    }

    // Establish joint constraints
    fprintf(ofp, "\nJoint Constraints:\n\tJoint\tDirection\n");
    for (i = 0; i <= 3 * nj - 1; ++i) {
        *(pjcode+i) = 1; // Initialize jcode with ones
    }
    // Potentially all joints are constrained, early exit likely
    for (i = 0; i <= 3 * nj - 1; ++i) {
        // Read in joint number and constraint direction from input file
        fscanf(ifp, "%d,%d\n", &j, &k);
        if ((j != 0) && (k != 0)) {
            switch (k) {
                case 1:
                    *(pjcode+j-1) = 0;
                    fprintf(ofp, "\t%d\t%d\n", j, k);
                    break;
                case 2:
                    *(pjcode+j+nj-1) = 0;
                    fprintf(ofp, "\t%d\t%d\n", j, k);
                    break;
                case 3:
                    *(pjcode+j+2*nj-1) = 0;
                    fprintf(ofp, "\t%d\t%d\n", j, k);
                    break;
                default:
                    fprintf(ofp, "\n***ERROR*** Joint constraint input not recognized");
                    return 1;
                    break;
            }
        } else {
            break;
        }
    }
    return 0;
}

/* This function generates joint code, jcode, by assigning integers in sequence, by
   columns, to all nonzero elements of jcode from 1 to neq; generate the member code,
   mcode, by transferring, via minc, columns of jcode into columns of mcode.  Follows
   procedure in Chapter 6 of Bathe and Wilson as well as Chapter 6 of Holzer */
void codes (int *pmcode, int *pjcode, int *pminc)
{
    int i, j, k, m; // Initialize function variables

    // Generate jcode
    neq = 0;
    for (i = 0; i <= nj - 1; ++i) {
        for(j = 0; j <= 2; ++j) {
            if(*(pjcode+i+j*nj) != 0) {
                neq ++;
                *(pjcode+i+j*nj) = neq;
            }
        }
    }

    fprintf(ofp, "\nNumber of equations (system DOFs): %d\n\n", neq);

    // Generate mcode from jcode using member incidence matrix
    for (i = 0; i <= ne - 1; ++i) {
        j = *(pminc+i); // Establish A-End joint number
        k = *(pminc+i+ne); // Establish B-end joint number
        // Increment of DOFs for each end
        for (m = 0; m <= 2; ++m) {
            *(pmcode+i+m*ne) = *(pjcode+j-1+m*nj); // Assign mcode entry using jcode
            *(pmcode+i+3*ne+m*ne) = *(pjcode+k-1+m*nj); // As above, for DOFs 4 thru 6
        }
    }
}

/* This function reads and echoes joint coordinates, element properties, as well as
   computing element lengths and direction cosines */
void prop (double *px, double *parea, double *pemod, double *peleng, double *pc1,
    double *pc2, double *pc3, int *pminc)
{
    // Initialize function variables
    int i, j, k;
    double x1, x2, x3, el1, el2, el3, xarea, mod;

    // Initialize joint coordinates
    for (i = 0; i <= nj - 1; ++i) {
        *(px+i) = 0;
        *(px+i+nj) = 0;
        *(px+i+2*nj) = 0;
    }

    fprintf(ofp, "Joint Coordinates:\n\tJoint\tDirection-1\tDirection-2\tDirection-3\n");

    for (i = 0; i <= nj - 1; ++i) {
        // Read in joint coordinates from input file
        fscanf(ifp, "%lf,%lf,%lf\n", &x1, &x2, &x3);
        *(px+i) = x1;
        *(px+i+nj) = x2;
        *(px+i+2*nj) = x3;
        fprintf(ofp, "\t%d\t%lf\t%lf\t%lf\n", i + 1, *(px+i), *(px+i+nj), *(px+i+2*nj));
    }

    fprintf(ofp, "\nElement Properties:\n");
    fprintf(ofp, "\tElement\t\tArea\t\tModulus\t\tLength\n");

    // Compute element lengths and direction cosines
    for (i = 0; i <= ne - 1; ++i) {
        j = *(pminc+i);
        k = *(pminc+i+ne);
        el1 = *(px+k-1) - *(px+j-1);
        el2 = *(px+k-1+nj) - *(px+j-1+nj);
        el3 = *(px+k-1+2*nj) - *(px+j-1+2*nj);
        *(peleng+i) = sqrt(el1*el1 + el2*el2 + el3*el3);
        *(pc1+i) = el1 / (*(peleng+i));
        *(pc2+i) = el2 / (*(peleng+i));
        *(pc3+i) = el3 / (*(peleng+i));
        // Read in element properties from input file
        fscanf(ifp, "%lf,%lf\n", &xarea, &mod);
        *(parea+i) = xarea;
        *(pemod+i) = mod;
        fprintf(ofp, "\t%d\t\t%lf\t%lf\t%lf\n",
            i + 1, *(parea+i), *(pemod+i), *(peleng+i));
    }
}

/* This function reads in the joint number, jnum, the joint direction, jdir, and the
   applied force, force */
int load (double *pq, int *pjcode)
{
    // Initialize function variables
    int i, k, jnum, jdir;
    double force;

    // Zero-out generalized load vector
    for (i = 0; i <= neq - 1; ++i) {
      *(pq+i) = 0;
    }
    fscanf(ifp, "%d,%d,%lf\n", &jnum, &jdir, &force);
    fprintf(ofp, "\nJoint Loads:\n\tGlobal Joint\tDirection\tForce\n");
    if (jnum != 0) { // Check for joint loading
        while (jnum != 0) { // Check for last joint load
            fprintf(ofp, "\t%d\t\t%d\t\t%lf\n", jnum, jdir, force);
            k = *(pjcode+jnum-1+(jdir-1)*nj); // Scan and load jcode
            // Store only joint loads corresponding to active global DOFs
            switch (k) {
                case (0):
                    break; // Do not store loads at supports
                default:
                    *(pq+k-1) = force;
                    break;
            }
            fscanf(ifp, "%d,%d,%lf\n", &jnum, &jdir, &force);
        }
    } else {
        fprintf(ofp, "\n***ERROR*** No joint forces present in input file");
        return 1;
    }
    return 0;
}

// This function determines kht using mcode, and determines maxa from kht
void skylin (int *pkht, int *pmaxa, int *pmcode, int *plss)
{
    int i, j, k, min; // Initialize function variables

    // Zero-out kht
    for (i = 0; i <= neq-1; ++i) {
        *(pkht+i) = 0;
    }

    /* Define column height array, kht.  Each address in kht corresponds to a column in
       the stiffness matrix; the value of the address defines the skyline height above
       the diagonal entry. */
    // Iterate over elements to span columns in mcode (LM array using Bathe's notation)
    for (i = 0; i <= ne - 1; ++i) {
        min = neq; // Guess at a reasonably large "minimum" to get things started
        for (j = 0; j <= 5; ++j) {
            // Does mcode entry correspond to a global DOF? Smaller than min?
            if ((*(pmcode+i+j*ne) > 0) && (*(pmcode+i+j*ne) < min)) {
                min = *(pmcode+i+j*ne); // New min...
            }
        }
        for (j = 0; j <= 5; ++j) {
            k = *(pmcode+i+j*ne);
            // Does the mcode entry correspond to a global DOFs?
            if (k != 0) {
                /* Use the mcode to discern column height kht.  The maximum difference
                   between non-zero mcode entries, (k - min), corresponding to given
                   element, defines the column height in the stiffness matrix,
                   corresponding to the degree of freedom "k". */
                if ((k - min) > *(pkht+k-1)) {
                    *(pkht+k-1) = k - min;
                }
            }
        }
    }

    /* Generate maxa array - provides the addresses in stiffness array for the diagonal
       elements of original stiffness matrix */
    *pmaxa = 1;
    for (i = 0; i <= neq - 1; ++i) { // Set counter to number of diagonal elements
        /* Last diagonal term + column height + one equals new diagonal address
           Last array element is used to define length of stiffness array */
        *(pmaxa+i+1) = *(pmaxa+i) + *(pkht+i) + 1;
    }
    /* Length of stiffness matrix.  This is true since last element of maxa is the
       address of the final skyline element, plus 1 */
    *plss = *(pmaxa+neq)-1;
    fprintf(ofp, "Length of stiffness array: %d\n\n", *plss);
}

/* This function computes the generalized tangent stiffness matrix and stores it as an
   array */
void stiff (double *pss, double *parea, double *pemod, double *peleng, double *pc1,
    double *pc2, double *pc3, double *pelong, int *pmaxa, int *pmcode, int *plss)
{
    // Initialize function variables
    int i, n, je, j, ie, k, L;
    // Matrix of constants that guides element stiffness assembly from G's and H's
    int index[6][6];
    double gamma; // Axial stiffness of truss element
    double strain; // Axial strain in truss element
    double G[7]; // Non-zero elements of linear stiffness matrix; upper left quadrant
    /* Elements of non-linear stiffness matrix coinciding with non-zero elements G;
       upper left quadrant */
    double H[7];

    index[0][0] = index[3][3] = 1;
    index[0][1] = index[3][4] = 3;
    index[0][2] = index[3][5] = 6;
    index[1][0] = index[4][3] = 3;
    index[1][1] = index[4][4] = 2;
    index[1][2] = index[4][5] = 5;
    index[2][0] = index[5][3] = 6;
    index[2][1] = index[5][4] = 5;
    index[2][2] = index[5][5] = 4;

    index[3][0] = index[0][3] = -1;
    index[3][1] = index[0][4] = -3;
    index[3][2] = index[0][5] = -6;
    index[4][0] = index[1][3] = -3;
    index[4][1] = index[1][4] = -2;
    index[4][2] = index[1][5] = -5;
    index[5][0] = index[2][3] = -6;
    index[5][1] = index[2][4] = -5;
    index[5][2] = index[2][5] = -4;

    // Initialize stiffness matrix to zero before each element call
    for (i = 0; i <= (*plss) - 1; ++i) {
        *(pss+i) = 0;
    }

    for (n = 0; n <= ne - 1; ++n) {
        gamma = (*(parea+n)) * (*(pemod+n)) / (*(peleng+n));
        strain = (*(pelong+n)) / (*(peleng+n));
        G[1] = gamma * (*(pc1+n)) * (*(pc1+n));
        G[2] = gamma * (*(pc2+n)) * (*(pc2+n));
        G[3] = gamma * (*(pc1+n)) * (*(pc2+n));
        G[4] = gamma * (*(pc3+n)) * (*(pc3+n));
        G[5] = gamma * (*(pc2+n)) * (*(pc3+n));
        G[6] = gamma * (*(pc1+n)) * (*(pc3+n));
        H[1] = H[2] = H[4] = (*(parea+n)) * (*(pemod+n)) * strain /
            ((*(peleng+n)) + (*(pelong+n)));
        H[3] = H[5] = H[6] = 0;

        /* Initialize index and then assign stiffness coefficients, G & H, of element n
           to the stiffness matrix by index, mcode, and maxa */
        for (je = 0; je <= 5; ++je) {
            j = *(pmcode+n+je*ne);
            if (j != 0) {
                // Check mcode above current entry to find rank of "j"
                for (ie = 0; ie <= je; ++ie) {
                    i = *(pmcode+n+ie*ne);
                    if (i != 0) {
                        if (i > j) { // Find element address as diagonal address + delta
                            k = *(pmaxa+i-1) + (i - j);
                        } else {
                            k = *(pmaxa+j-1) + (j - i);
                        }
                        L = index[ie][je];
                        /* Add current element stiffness to previous elements'
                           contributions to the given DOFs */
                        if (L > 0) {
                            *(pss+k-1) += G[L] + H[L];
                        } else {
                            *(pss+k-1) -= G[-1 * L] + H[-1 * L];
                        }
                    }
                }
            }
        }
    }
}

// This function computes the generalized internal force vector
void forces (double *pf, double *parea, double *pemod, double *pc1, double *pc2,
    double *pc3, double *pelong, double *peleng, int *pmcode)
{
    // Initialize function variables
    int i, j, k;

    for (i = 0; i <= ne - 1; ++i) {
        for (j = 0; j <= 5; ++j) {
            k = *(pmcode+i+j*ne);
            if (k != 0) {
                if (j == 0) {
                    *(pf+k-1) -= (*(parea+i)) * (*(pemod+i)) *
                    ((*(pelong+i)) / (*(peleng+i)) +
                    0.5 * (*(pelong+i)) / (*(peleng+i)) * (*(pelong+i)) /
                    (*(peleng+i))) * (((*(peleng+i)) + (*(pelong+i))) / (*(peleng+i))) *
                    (*(pc1+i));
                } else if (j == 1) {
                    *(pf+k-1) -= (*(parea+i)) * (*(pemod+i)) *
                    ((*(pelong+i)) / (*(peleng+i)) +
                    0.5 * (*(pelong+i)) / (*(peleng+i)) * (*(pelong+i)) /
                    (*(peleng+i))) * (((*(peleng+i)) + (*(pelong+i))) / (*(peleng+i))) *
                    (*(pc2+i));
                } else if (j == 2) {
                    *(pf+k-1) -= (*(parea+i)) * (*(pemod+i)) *
                    ((*(pelong+i)) / (*(peleng+i)) +
                    0.5 * (*(pelong+i)) / (*(peleng+i)) * (*(pelong+i)) /
                    (*(peleng+i))) * (((*(peleng+i)) + (*(pelong+i))) / (*(peleng+i))) *
                    (*(pc3+i));
                } else if (j == 3) {
                    *(pf+k-1) += (*(parea+i)) * (*(pemod+i)) *
                    ((*(pelong+i)) / (*(peleng+i)) +
                    0.5 * (*(pelong+i)) / (*(peleng+i)) * (*(pelong+i)) /
                    (*(peleng+i))) * (((*(peleng+i)) + (*(pelong+i))) / (*(peleng+i))) *
                    (*(pc1+i));
                } else if (j == 4) {
                    *(pf+k-1) += (*(parea+i)) * (*(pemod+i)) *
                    ((*(pelong+i)) / (*(peleng+i)) +
                    0.5 * (*(pelong+i)) / (*(peleng+i)) * (*(pelong+i)) /
                    (*(peleng+i))) * (((*(peleng+i)) + (*(pelong+i))) / (*(peleng+i))) *
                    (*(pc2+i));
                } else if (j == 5) {
                    *(pf+k-1) += (*(parea+i)) * (*(pemod+i)) *
                    ((*(pelong+i)) / (*(peleng+i)) +
                    0.5 * (*(pelong+i)) / (*(peleng+i)) * (*(pelong+i)) /
                    (*(peleng+i))) * (((*(peleng+i)) + (*(pelong+i))) / (*(peleng+i))) *
                    (*(pc3+i));
                }
            }
        }
    }
}

/* This function determines the solution within the increment using the compact Gaussian
   elimination strategy (Holzer pp. 290, 296, 307) and based on subroutine COLSOL -
   "active column solution or the skyline (or column) reduction method"
   (Bathe 1976, p. 257) */
int solve (double *pss, double *pr, double *pdd, int *pmaxa)
{
    // Initialize function variables
    int i, n, kn, kl, ku, kh, k, ic, klt, j, ki, nd, kk, l;
    double b, c, q[neq];

    /* Pass residual to local array variable "q" which will later be used to pass the
       incremental displacements back out to main */
    for (i = 1; i <= neq; ++i) {
        q[i - 1] = *(pr+i-1);
    }

    // Perform LDL^t factorization of the stiffness matrix
    for (n = 1; n <= neq; ++n) {
        kn = *(pmaxa+n-1);
        kl = kn + 1;
        ku = *(pmaxa+n) - 1;
        kh = ku - kl;
        if (kh < 0) {
            if (*(pss+kn-1) <= 0) {
                fprintf(ofp, "\n\n***ERROR*** Non-positive definite stiffness matrix");
                return 1;
            }
        } else if (kh > 0) {
            k = n - kh;
            ic = 0;
            klt = ku;
            for (j = 1; j <= kh; ++j) {
                ic ++;
                klt -= 1;
                ki = *(pmaxa+k-1);
                nd = *(pmaxa+k) - ki - 1;
                if (nd > 0) {
                    if (nd < ic) {
                        kk = nd;
                    } else {
                        kk = ic;
                    }
                    c = 0;
                    for (l = 1; l <= kk; ++l) {
                        c += *(pss+ki-1+l) * (*(pss+klt-1+l));
                    }
                    *(pss+klt-1) -= c;
                }
                k ++;
            }
            k = n;
            b = 0;
            for (kk = kl; kk <= ku; ++kk) {
                k -= 1;
                ki = *(pmaxa+k-1);
                c = *(pss+kk-1) / (*(pss+ki-1));
                b += c * (*(pss+kk-1));
                *(pss+kk-1) = c;
            }
            *(pss+kn-1) -= b;
            if (*(pss+kn-1) <= 0) {
                fprintf(ofp, "\n\n***ERROR*** Non-positive definite stiffness matrix");
                return 1;
            }
        } else if (kh == 0) {
            k = n;
            b = 0;
            for (kk = kl; kk <= ku; ++kk) {
                k -= 1;
                ki = *(pmaxa+k-1);
                c = *(pss+kk-1) / (*(pss+ki-1));
                b += c * (*(pss+kk-1));
                *(pss+kk-1) = c;
            }
            *(pss+kn-1) -= b;
            if (*(pss+kn-1) <= 0) {
                fprintf(ofp, "\n\n***ERROR*** Non-positive definite stiffness matrix");
                return 1;
            }
        }
    }

    // Reduce right-hand-side load vector
    for (n = 1; n <= neq; ++n) {
        kl = *(pmaxa+n-1) + 1;
        ku = *(pmaxa+n) - 1;
        kh = ku - kl;
        if (kh >= 0) {
            k = n;
            c = 0;
            for (kk = kl; kk <= ku; ++kk) {
                k -= 1;
                c += *(pss+kk-1) * q[k - 1];
            }
            q[n - 1] -= c;
        }
    }

    // Back-substitute
    for (n = 1; n <= neq; ++n) {
        k = *(pmaxa+n-1);
        q[n - 1] = q[n - 1] / (*(pss+k-1));
    }
    if (neq == 1) {
        return 0;
    }
    n = neq;
    for (i = 2; i <= neq; ++i) {
        kl = *(pmaxa+n-1) + 1;
        ku = *(pmaxa+n) - 1;
        kh = ku - kl;
        if (kh >= 0) {
            k = n;
            for (kk = kl; kk <= ku; ++kk) {
                k -= 1;
                q[k - 1] -= *(pss+kk-1) * q[n - 1];
            }
        }
        n -= 1;
    }

    /* Pass the local array variable "q" to the incremental displacements for use in
       main */
    for (i = 1; i <= neq; ++i) {
        *(pdd+i-1) = q[i - 1];
    }
    return 0;
}

/* This function performs convergence check; assume that convergence is reached,
   i.e. inconv = 0, until proven otherwise */
void test (double *pf, double *pfp, double *pqtot, double *pdd, double *pfpi,
    double *pintener1, int *pinconv, int *pneq, double *ptolfor, double *ptolener)
{
    // Initialize function variables
    int i;
    double unbfi = 0;
    double unbfp = 0;
    double inteneri = 0;
    double c = 0;
    *pinconv = 0;

    if (*ptolfor < 1) {
        // Compute Euclidean vector norm of unbalanced forces
        for (i = 0; i <= neq - 1; ++i) {
            unbfi += (*(pqtot+i) - (*(pf+i))) * (*(pqtot+i) - (*(pf+i)));
            unbfp += (*(pqtot+i) - (*(pfp+i))) * (*(pqtot+i) - (*(pfp+i)));
        }
        // Check against force tolerance
        if (unbfp != 0) {
            c = sqrt(unbfi) / sqrt(unbfp);
            if (c > *ptolfor) {
               *pinconv += 10;
            }
        }
    }

    if (*ptolener < 1) {
        // Compute incremental internal energy during each iteration
        for (i = 0; i <= neq - 1; ++i) {
            inteneri += (*(pdd+i) * (*(pqtot+i) - (*(pfpi+i))));
        }
        // Check against tolerances
        if (*pintener1 != 0) {
            c = fabs(inteneri / *pintener1);
            if (c > *ptolener) {
                *pinconv += 100;
            }
        }
    }
}

// This function updates member geometry
void updatc (double *px, double *pdd, double *pc1, double *pc2, double *pc3,
    double *pelong, double *peleng, double *pdeflen, int *pminc, int *pjcode)
{
    // Initialize functin variables
    int i, j, k;
    double el1, el2, el3;

    // Update nodal coordinates during each increment
    for (i = 0; i <= nj - 1; ++i) {
        for (j = 0; j <= 2; ++j) {
            k = *(pjcode+i+j*nj);
            if (k != 0) {
                *(px+i+j*nj) += *(pdd+k-1);
            }
        }
    }

    // Update member lengths and direction cosines during each increment
    for (i = 0; i <= ne - 1; ++i) {
        j = *(pminc+i);
        k = *(pminc+i+ne);
        el1 = *(px+k-1) - *(px+j-1);
        el2 = *(px+k-1+nj) - *(px+j-1+nj);
        el3 = *(px+k-1+2*nj) - *(px+j-1+2*nj);
        *(pdeflen+i) = sqrt(el1*el1 + el2*el2 + el3*el3);
        *(pelong+i) = *(pdeflen+i) - *(peleng+i);
        *(pc1+i) = el1 / (*(pdeflen+i));
        *(pc2+i) = el2 / (*(pdeflen+i));
        *(pc3+i) = el3 / (*(pdeflen+i));
    }
}
