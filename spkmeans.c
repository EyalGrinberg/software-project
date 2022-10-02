#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>


/*declarations*/

/* general auxilary functions */
void print2Darray(double**, int, int); /* print in requested output format */
double** mem2D(int, int); /* allocate memory for 2D array */
void free2D(double **); /* free memory for 2D array */
double norm(double*, double*, int); /* calculate norm */
double dist(double *, double *); /* calculate distance between 2 points */
double** extractFromFile(char*); /* retrieve data points from file */

/* wam */
double** wam(double**); /* create weights adjacency matrix */
/* ddg */
double** ddg(double**); /* create D matrix */
/* Lnorm functions */
double** leftDiagMult(double**, double**); /* matrix multiplication when left is a diagonal matrix */
double** rightDiagMult(double**, double**); /* matrix multiplication when right is a diagonal matrix */
double** invertRoot(double**); /* diagonal matrix to the power of -0.5 */
double** lnorm(double**); /* create Lnorm matrix */

/* jacobi functions */
double** newP(int, int, double, double); /* update P rotation matrix */
double* obtain(double, double, double); /* obtain s,c values for jacobi algorithm */
int convergance(double**, double**); /* check stop condition */
double** updateVects(double**, double**); /* matrix multipication */
double* find_I_J(double**); /* find largest non diagonal value in a matrix */
double** jacobi(double**); /* jacobi algorithm main function */

/* spk functions */
double** sortMat(double**); /* sort matrix's columns by first row */
void swap(int*, int*); /* auxiliary function for bubble sort */
int eigenGap(double**); /* find k by heuristic */
double** trans(double**, int , int); /* transpose matrix */

/* C-api functions */
double **kmeans(double **, int *, int, int); /* kmeans algorithm */
double **mainFuncCapi(double **, int , int , int ); /* main function called from C-api module */

/*global variables*/
int row, col;


/*functions*/
void print2Darray(double** a, int x, int y){ /* print in requested output format */
    int rows, columns;
 	for(rows = 0; rows < x; rows++)
  	{
  		for(columns = 0; columns < y; columns++)
  		{
  			if (columns == y - 1){
                printf("%.4f \n" , a[rows][columns]);
            } else{
                printf("%.4f,", a[rows][columns]);
            }
		}
  	}
    printf("\n");   	
}

double** mem2D(int x, int y){ /* allocate memory for 2D array */
    double** mat;
    int i;
    mat = (double**)calloc(x, sizeof(double*));
    assert(mat != NULL && "An Error Has Occured");
    for(i=0; i<x; i++){ 
        mat[i] = (double*)calloc(y, sizeof(double));
        assert(mat[i] != NULL && "An Error Has Occured");
    }
    return mat;
}

void free2D(double **mat){ /* free memory for 2D array */
    free(*mat);
    free(mat);
}

double norm(double *point1, double *point2 , int x){ /* calculate norm */
    double res = 0;
    int i = 0;
    for (i=0; i<x; i++){
         res += pow(*(point1+i) - *(point2+i) , 2);
    }
    res = pow(res , 0.5);
    return res;
}

double ** extractFromFile(char* fileName){ /* retrieve data points from file */

    int i, j;    
    char *cont, *ptr;
    char line[1000]; 
    char lineToSplit[1000];
    double **pointsArray;
    int maxCol = 0;
    int currCol = 0;
    FILE* myFile1 = NULL;
    FILE* myFile = NULL;
    myFile1 = fopen(fileName, "r");
    if (myFile1 == NULL){
        printf("Invalid Input!3");
        exit(1);
    }
    while (fgets(line , sizeof(line) , myFile1)){ /*get dims*/
        row++;
        cont = strtok(line , ",");
        if ((maxCol != currCol) && (maxCol != 0)){ /*check dims validity*/
                printf("Invalid Input!");
                exit(1);
            }
        currCol = 0;
        while(cont != NULL){
            cont = strtok(NULL , ",");
            currCol++;
        }
        if ((maxCol < currCol) && (maxCol == 0)){
            maxCol = currCol;
        }
    }
    col = maxCol;
    fclose(myFile1);
    /*dynamic mat alloc*/ 
    pointsArray = mem2D(row, col);
    myFile = fopen(fileName, "r");
    assert(myFile != NULL && "Invalid Input!5");
    i=0;
    while (fgets(line , sizeof(lineToSplit) , myFile)){ /*initial points matrix*/
         j = 0;
         cont = strtok(line , ","); /* read until next value */    
         while (cont != NULL){
             pointsArray[i][j] = strtod(cont, &ptr); /* convert string to double */
             j++;
             cont = strtok(NULL , ",");
         }
         i++;
    }
    fclose(myFile);
    return pointsArray;
}

double** wam(double** dataPoints){ /* create weights adjacency matrix */

    double weight;
    int i,j;
    double** weightsMat;

    weightsMat = mem2D(row , row);

    for (i = 0; i < row; i++){
        for (j = i; j < row; j++){
            if (i == j) {
                weightsMat[i][j] = 0;
            }
            else {
                weight = exp(norm(dataPoints[i] , dataPoints[j], col) / -2);  /* determine weight by requsted function */
                weightsMat[i][j] = weight;
                weightsMat[j][i] = weight; /* create a symetric weights matrix */
            }
        }
    }
    return weightsMat;
}

double** ddg(double** dataPoints){ /* create D matrix */

    double** diagMat;
    double** weightsMat;
    int i, j;
    double sum = 0;

    weightsMat = mem2D(row, row);
    diagMat = mem2D(row, row);
    weightsMat = wam(dataPoints);

    for(i = 0; i < row; i++){
        sum = 0;
        for (j = 0; j < row; j++){ /* sum weights of W[i] and insert to D[i][i] */
            sum += weightsMat[i][j];
        }
        diagMat[i][i] = sum;
    }
    return diagMat;
}

double** leftDiagMult(double** diagMat, double** otherMat){ /* other[i][j] = other[i][j] * diag[i][i] */

    int i, j;
    double** res;

    res = mem2D(row, row);

    for (i = 0; i < row; i++){
        for (j = 0; j < row; j++){
            res[i][j] = otherMat[i][j] * diagMat[i][i];
        }
    }
    return res;
}

double** rightDiagMult(double** otherMat, double** diagMat){ /* other[i][j] = other[i][j] * diag[j][j] */
    int i, j;
    double** res;

    res = mem2D(row, row);

    for (i = 0; i < row; i++){
        for (j = 0; j < row; j++){
            res[i][j] = otherMat[i][j] * diagMat[j][j];
        }
    }
    return res;
}

double** invertRoot(double** diagMat){ /* diagonal matrix to the power of -0.5 */

    double** invSqrtMat;
    int i;

    invSqrtMat = mem2D(row, row);

    for (i = 0; i < row; i++){
        invSqrtMat[i][i] = pow(diagMat[i][i] , -0.5);
    }
    return invSqrtMat;
}

double** lnorm(double** datapoints){ /* create Lnorm matrix */

    double** diagMat;
    double** mat_L;
    double** weightsMat;
    double** invSqrtMat;
    double** res;
    int i, j;

    diagMat = mem2D(row, row);
    mat_L = mem2D(row, row);
    weightsMat = mem2D(row, row);
    invSqrtMat = mem2D(row, row);
    res = mem2D(row, row);
    /* call previous steps in the algorithm */
    weightsMat = wam(datapoints);
    diagMat = ddg(datapoints);

    for (i = 0; i < row; i++){ /*create L matrix*/
        for (j = 0; j< row; j++){
            mat_L[i][j] = diagMat[i][j] - weightsMat[i][j];
        }
    }
    for (i = 0; i < row; i++){ /* create D^-0.5 matrix */
        invSqrtMat[i][i] = pow(diagMat[i][i] , -0.5);
    }
    res = rightDiagMult( leftDiagMult(invSqrtMat , mat_L) , invSqrtMat); /*triple matrix multiplication*/
    return res;
}

double* obtain(double aij, double aii, double ajj){ /* get s and c for jacobi */

    double* result;
    int sign;
    double theta, t;
    
    if (aij == 0.0){ /* if largest non diagonal is 0 ---> theta is INFINITY */
        theta = strtod("Inf", NULL);
    } else {
        theta = (ajj-aii)/(2*aij);
    }
    sign = theta>=0 ? 1 : -1;
    t = sign/(fabs(theta) + pow((pow(theta, 2) + 1), 0.5));
    result = (double*)calloc(2, sizeof(double)); /* return an array of c,s */
    assert(result != NULL && "An Error Has Occured");
    result[0] = pow(pow(t,2) + 1, -0.5); /* c */
    result[1] = t*result[0]; /* s */
    return result;
}

double** newP(int iInd, int jInd, double c, double s){ /* create rotation matrix P */

    double** pMat;
    int i;

    pMat = mem2D(row,row);

    for(i=0; i<row; i++){ /* start from a unit matrix */
        pMat[i][i] = 1;
    }
    /* update with c,s */
    pMat[iInd][iInd] = c;
    pMat[jInd][jInd] = c;
    pMat[iInd][jInd] = s;
    pMat[jInd][iInd] = -s;

    return pMat;
}

double* find_I_J(double** aMat){ /* get aij - largest non diagonal value and relevant indices */

    double aij, iRes, jRes;
    double* vals;
    int i, j, iInd, jInd;
    aij = 0.0;
    iRes = 0.0; 
    jRes = 0.0;

    for(i = 0; i<row; i++){
        for(j=i; j<row; j++){
            if ( (fabs(aMat[i][j])>fabs(aij)) && (i!=j) ){
                aij = aMat[i][j];
                iRes = i;
                jRes = j;
            }
        }
    }
    vals = (double*)calloc(3, sizeof(double)); /* return an array of max aij, i, j */
    assert(vals != NULL && "An Error Has Occured");
    iInd = (int)iRes;
    jInd = (int)jRes; /* cast to int for comparison and return */
    if ((iInd == jInd) && (row != 1)){ /* in case of I matrix */
        iRes = 0.0;
        jRes = 1.0;
    }
    vals[0]= aij;
    vals[1] = iRes;
    vals[2] = jRes;

    return vals;
}

int convergance(double** aMat, double** prevA){ /* check stop condition */
    int i, j;
    double offA = 0;
    double offPrev = 0;
    double eps = 1.0 * pow(10, -5);
    for(i=0; i<row; i++){
        for(j=0; j<row; j++){
            if (i!=j){
                offA += pow(aMat[i][j], 2);
                offPrev += pow(prevA[i][j], 2);
            }
        }
    }
    if ((offPrev - offA) <= eps){
        return 1;
    }
    return 0;

}

double** updateVects(double** tempV, double** pMat){ /* matrix multiplication */

    double** res;
    int i,j,k;

    res = mem2D(row, row);

    for(i=0; i<row; i++){
        for(j=0; j<row; j++){
            for(k=0; k<row; k++){
                res[i][j] += tempV[i][k]*pMat[k][j];
            }
        }
    }
    return res;
}

double** trans(double** mat,int x ,int y){ /* transpose matrix */

    double** res;
    int i, j;

    res = mem2D(x, y);

    for(i = 0; i < x; i++){
        for (j = 0; j < y; j++){
            res[i][j] = mat[j][i];    
        }
    }
    return res;
}

void swap(int *a, int *b){ /* auxiliary function for bubble sort */
    int tmp;
    tmp = *a;
    *a = *b;
    *b = tmp;
}

double** sortMat(double** eigenMat){ /* returns the matrix of eigenvals/vecs sorted decreasingly */

    double **tmp1Mat, **tmp2Mat, **res;
    int* indices;
    int i, j;

    indices = (int*)calloc(row, sizeof(int));
    assert(indices!= NULL && "An Error Has Occured");
    for(i=0; i<row; i++){ /* create indices array */
        indices[i] = i;
    }
    tmp1Mat = mem2D(row , row);
    for(i=0; i<row - 1; i++){/*bubble sort indices*/
        for(j=0; j<(row - i - 1); j++){
            if (eigenMat[0][indices[j]] > eigenMat[0][indices[j+1]]){
                swap(&indices[j], &indices[j+1]);               
            } 
        }
    }
    for(i=0; i<row; i++){/*create EigenMat transpose, vectors are rows*/
        for(j = 1; j<row + 1; j++){
            tmp1Mat[i][j-1] = eigenMat[j][i];
        }
    }
    tmp2Mat = mem2D(row, row);
    for (i = 0; i < row; i++){ /* sort vectors rows by eigenvalues according to sorted indices decreasingly */
        for (j = 0; j < row; j++){
            tmp2Mat[i][j] = tmp1Mat[indices[row - i - 1]][j];
        }
    }
    tmp2Mat = trans(tmp2Mat, row, row); /* transpose again, vectors are columns */
    res = mem2D(row + 1, row);
    for (i = 0; i < row; i++){
        res[0][i] = eigenMat[0][indices[row - i - 1]]; /* value for each column */
    }
    for (i = 1; i < row + 1; i++){ /* column vector for each value */
        for (j = 0; j < row; j++){
            res[i][j] = tmp2Mat[i - 1][j];
        }
    }
    return res;
}

int eigenGap(double** sortedVmat){ /* find k by heuristic */

    int i, k;
    double maxGap = 0;

    double *gaps = (double*)calloc(row / 2 , sizeof(double)); /* eigen gaps array */
    k = 0; 
    for(i = 0; i < row / 2; i++){
        gaps[i] = fabs(sortedVmat[0][i] - sortedVmat[0][i+1]);
    }
    for (i = 0; i < row / 2; i++){ /* find max gap ---> wanted k */
        if (gaps[i] > maxGap){
            maxGap = gaps[i];
            k = i;
        }
    }
    return k + 1;
}

double** jacobi(double** aMat){ /* jacobi algorithm main function */

    double **vResult, **pMat, **tempV, **prevA;
    double *parameters, *c_s;
    double aij, c, s;
    int i,j, iRes, jRes;
    int flag = 0;
    int cnt = 0;

    parameters = (double*)calloc(3, sizeof(double)); /* an array of aij, i, j */
    assert(parameters != NULL && "An Error Has Occured");
    c_s = (double*)calloc(2, sizeof(double)); /* an array of c, s */
    assert(c_s != NULL && "An Error Has Occured");

    tempV = mem2D(row,row);
    pMat = mem2D(row,row);
    prevA = mem2D(row,row);

    for(i=0; i<row; i++){/* init V as unit matrix to start matrix multiplication */
        tempV[i][i] = 1;
    }

    while ((flag == 0) && (cnt<100)) { /*Jacobi main loop*/
        parameters = find_I_J(aMat); /* obtain parameters aij, i, j */
        aij = parameters[0];
        iRes = (int)parameters[1];
        jRes = (int)parameters[2];
        if (iRes == jRes){ /* case of a single entry matrix */
            vResult = mem2D(2, 1);
            vResult[0][0] = aMat[0][0];
            vResult[1][0] = 1.0;
            return vResult;
        }
        c_s = obtain(aij, aMat[iRes][iRes], aMat[jRes][jRes]); /* obtain c, s */
        c = c_s[0];
        s = c_s[1];
        pMat = newP(iRes, jRes, c, s); /* update P matrix */
        for(i=0; i<row; i++){/*save current A values for convergence checking */
            for(j=0; j<row; j++){
                prevA[i][j] = aMat[i][j];
            }
        }
        aMat = updateVects(trans(pMat, row, row), updateVects(prevA, pMat)); /* create A' */
        tempV = updateVects(tempV, pMat); /*   multiply V * Pi   */
        flag = convergance(aMat, prevA); /* check condintion */
        cnt+=1;
    }
    vResult = mem2D((row+1), row); 
    for(i=0; i<row; i++){ /* insert eigenvalues to V */
        vResult[0][i] = aMat[i][i];
    }
    for(i=1; i<=row; i++){ /* insert eigenvectors to V */
        for(j=0; j<row; j++){
            vResult[i][j] = tempV[i-1][j];
        }
    }
    return vResult;
}

double dist(double *point, double *centroid){ /* calculate distance between two points */
    double res = 0;
    int i = 0;
    for (i=0; i<col; i++){
         res += pow(*(point+i) - *(centroid+i) , 2);
    }
    return res;
}

double **kmeans(double **pointsArray, int *initCentroids, int rows, int col){ /* kmeans algorithm */

    int i, j, l, k, flag;
    int *clusters, *denominators;
    double **sums, **prevAvgs, **currAvgs;
    int cnt = 0;
    int max_iter = 300; 
    i = 0;
    flag = 1;
    k = col;
    row = rows;

    /* memory allocations */
    denominators = (int *)calloc(k , sizeof(int));
    clusters = (int *)calloc(row , sizeof(int));
    sums = mem2D(k, col);
    prevAvgs = mem2D(k, col);
    currAvgs = mem2D(k, col);

    for (i=0; i<k; i++){ /*init centroids according to indices received from k++ */
            for (j=0; j<col; j++){
                currAvgs[i][j] = pointsArray[initCentroids[i]][j];
            }
    }
    while ((flag == 1) && (cnt < max_iter)){ /* main kmeans loop */
        for (i = 0; i < k; i++){ /*init denom to zeroes*/
            denominators[i] = 0;
        }
        for (i=0; i<k; i++){ /*init sums to zeroes*/
            for (j=0; j<col; j++){
                sums[i][j] = 0;
            }
        } 
        for (i=0; i<row; i++){
            double currDist = __DBL_MAX__;
            double minDist = __DBL_MAX__;
            for (j=0; j<k; j++){ /*find closest centroid*/
                currDist = dist((double *)pointsArray[i] , (double *)currAvgs[j]);
                if (currDist < minDist){
                    minDist = currDist;
                    clusters[i] = j;
                }
            }/*end of inner loop for centroids*/
            for (l=0; l<col; l++){ /*adding point's value in each dimension to relevant centroid*/
                sums[clusters[i]][l] += pointsArray[i][l]; 
            }
            denominators[clusters[i]] += 1; /* increasing divisor */
        }/*end of points loop*/

        for (i=0; i<k; i++){ /*update avgs*/
            for (j=0; j<col; j++){
                currAvgs[i][j] = (sums[i][j])/(denominators[i]);
            }
        }        
        flag = 0;
        for (i = 0; i < k; i++) /*check eps condition*/
        {
            double currNorm = dist((double *)currAvgs[i] , (double *)prevAvgs[i]);
            if (currNorm > 0){ 
                flag = 1;
            }
        }
        for (i = 0; i < k; i++){ /*prevAvgs <-- currAvgs*/
            for (j=0; j<col; j++){
                prevAvgs[i][j] = currAvgs[i][j];
            }
        }
        cnt++; /*count iters*/
    }

    free(clusters);
    free(denominators);
    free2D(sums);
    free2D(prevAvgs);
    free2D(pointsArray);

    return currAvgs;

}

double **mainFuncCapi(double **inputMat, int goal, int x, int y){ /* input from module, output to module */

    double **res, **lnormOutput;
    row = x;
    col = y;
    if (goal == 1){ /* W matrix */
        res = mem2D(x, x);
        res = wam(inputMat);
        return res;
    }
    if (goal == 2){ /* D matrix */
        res = mem2D(x, x);
        res = ddg(inputMat);
        return res;
    }
    if (goal == 3){ /* Lnorm matrix */
        res = mem2D(x, x);
        res = lnorm(inputMat);
        return res;
    }
    if (goal == 4){ /* jacobi */
        res = mem2D(x + 1, x); /* jacobi matrix is n+1 X n */
        res = jacobi(inputMat); /* this inputMat is a symetric Lnorm and not a regular data set */
        return res; 
    }
    if (goal == 5){ /* spk - step 1 receive data points and return jacobi */
        lnormOutput = mem2D(x, x); 
        res = mem2D(x + 1, x);
        lnormOutput = lnorm(inputMat); /* get Lnorm to pass to jacobi */
        res = jacobi(lnormOutput); /* get jacobi */
        return res;
    }
    return inputMat; /* should not get here */
}

int main(int argc, char *argv[]){ /* input from CMD, print wam|ddg|lnorm|jacobi matrix */ 

    double **dataPoints, **res;

    char *goal = argv[1];
    char *fileName = argv[2];

    if (argc != 3){ 
        printf("Invalid Input!1");
        exit(1);
    }
    if ( (strcmp(goal, "wam") != 0) && (strcmp(goal, "ddg") != 0) && (strcmp(goal, "lnorm") != 0) && (strcmp(goal, "jacobi") != 0) ){  /*check goal validity*/ 
        printf("Invalid Input!2");
        exit(1);
    }

    dataPoints = mem2D(row, row);
    dataPoints = extractFromFile(fileName); /* recieve data points */

    if (strcmp(goal, "jacobi") == 0){ /* in case of jacobi: output dimensions are n+1 X n */
        res = mem2D(row + 1, row);
        res = jacobi(dataPoints); 
        print2Darray(res, row + 1, row);
        free2D(dataPoints);
        free2D(res);
        return 0;
    }

    res = mem2D(row, row); /* all other goals: output dimensions are n X n */
    if (strcmp(goal, "wam") == 0){
        res = wam(dataPoints);
    }
    if (strcmp(goal, "ddg") == 0){
        res = ddg(dataPoints);
    }
    if (strcmp(goal, "lnorm") == 0){
        res = lnorm(dataPoints);
    }
    print2Darray(res, row, row);
    free2D(dataPoints);
    free2D(res);
    return 0; 

}
