#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include "spkmeans.h"

double **pyToCmat(PyObject *mat, int x, int y){ /* convert python 2D list to C 2D array */
    
    double **res = mem2D(x, y);
    int i, j;
    
    for (i = 0; i < x; i++){
        for (j = 0; j < y; j++){
            res[i][j] = PyFloat_AsDouble(PyList_GetItem(PyList_GetItem(mat, i), j)); /* retrieve value and change type */
        }
    }
    
    return res;

}

PyObject *cToPyMat(double **mat, int x, int y){ /* convert C 2D array to python 2D list */

    PyObject *res = PyList_New(x);
    int i, j;
    
    for (i = 0; i < x; i++){
        PyObject *row = PyList_New(y);
        for (j = 0; j < y; j++){
            PyList_SetItem(row, j, PyFloat_FromDouble(mat[i][j])); /* change type and insert */
        }
        PyList_SetItem(res, i, row);
    }
    
    return res;

}


static PyObject *pyToCtoPy(PyObject *self, PyObject *args){ /* all goals except spk - one back and forth */

    double **outputC, **pointsMat;
    int row, col;
    int goal;
    PyObject *dataPoints, *resPy;

    if (!PyArg_ParseTuple(args, "Oi", &dataPoints, &goal)){ /* assign args */
        printf("Invalid Input C-api1!");
        exit(1);
    }
    if ( (!PyList_Check(dataPoints)) || (!PyList_Check(PyList_GetItem(dataPoints, 0))) ){/* check args */
        printf("Invalid Input C-api2!");
        exit(1);
    }
    
    row = (int)PyList_Size(dataPoints); /* get dimensions */
    col = (int)PyList_Size(PyList_GetItem(dataPoints, 0));
    pointsMat = pyToCmat(dataPoints, row, col); /* convert Py data points to C data points */
    outputC = mainFuncCapi(pointsMat, goal, row, col); /* send to main func in C and receive requested output */
    if (goal == 4){ /* 4 means "jacobi" */
        resPy = cToPyMat(outputC, row + 1, row); /* +1 row for eigenvalues */
    }
    else{
        resPy = cToPyMat(outputC, row, row);
    }

    free2D(pointsMat);
    free2D(outputC);

    return resPy;
}

static PyObject *spkCapi(PyObject *self, PyObject *args){ /* first back and forth for spk: data points ---> V matrix from jacobi */
    
    double **outputVmatC, **pointsMat;
    int row, col;
    PyObject *dataPoints, *resPyVmat;

    if (!PyArg_ParseTuple(args, "O", &dataPoints)){ /* assign args */
        printf("Invalid Input spk C-api1!");
        exit(1);
    }
    if ( (!PyList_Check(dataPoints)) || (!PyList_Check(PyList_GetItem(dataPoints, 0))) ){ /* check args */
        printf("Invalid Input  spk C-api2!");
        exit(1);
    }
    
    row = (int)PyList_Size(dataPoints); /* get dimensions */
    col = (int)PyList_Size(PyList_GetItem(dataPoints, 0));
    pointsMat = pyToCmat(dataPoints, row, col); /* convert Py data points to C data points */
    outputVmatC = mainFuncCapi(pointsMat, 5, row, col); /* send to main func in C and receive V matrix (goal 5 is "spk") */
    resPyVmat = cToPyMat(outputVmatC, row + 1, row); /* convert back to python */

    free2D(pointsMat);
    free2D(outputVmatC);

    return resPyVmat; /* send V matrix to Py - next module func recieves T matrix and k indices from kmeans++ */
}

static PyObject *eigenGapCapi(PyObject *self, PyObject *args){ /* returns k (eigen heuristic) + sorted V matrix */

    double **outputSortedVmatC, **inputC_V_mat;
    int row, col, k;
    PyObject *inputPy_V_mat, *resPySortedVmat;

    if (!PyArg_ParseTuple(args, "O", &inputPy_V_mat)){ /* assign args */
        printf("Invalid Input eigenGap C-api1!");
        exit(1);
    }
    if ( (!PyList_Check(inputPy_V_mat)) || (!PyList_Check(PyList_GetItem(inputPy_V_mat, 0))) ){ /* check args */
        printf("Invalid Input eigenGap C-api2!");
        exit(1);
    }

    row = (int)PyList_Size(inputPy_V_mat); /* get dimensions */
    col = (int)PyList_Size(PyList_GetItem(inputPy_V_mat, 0));
    inputC_V_mat = pyToCmat(inputPy_V_mat, row, col); /* convert V matrix to C */
    outputSortedVmatC = sortMat(inputC_V_mat); /* sort V matrix by eigen values */
    resPySortedVmat = cToPyMat(outputSortedVmatC, row, col); /* convert sorted V to Py */
    k = eigenGap(outputSortedVmatC); /* obtain k by hueristic */

    free2D(inputC_V_mat);
    free2D(outputSortedVmatC);

    return Py_BuildValue("(iO)" , k , resPySortedVmat); /* return k , sorted V matrix */

}

static PyObject *kmeansCapi(PyObject *self, PyObject *args){ /* second back and forth - T mat + k indices ---> k clusters */
    
    double **tMatC, **outputClustersC;
    int *centroidsArrayC;
    int row, col, k, i;
    PyObject *tMatPy, *initCentroidsPy, *resClustersPy;

    if (!PyArg_ParseTuple(args, "OOi", &tMatPy, &initCentroidsPy, &k)){ /* assign args */
        printf("Invalid Input kmeans C-api1!");
        exit(1);
    }
    if ( (!PyList_Check(tMatPy)) || (!PyList_Check(PyList_GetItem(tMatPy, 0))) || (!PyList_Check(initCentroidsPy)) ){ /* check args */
        printf("Invalid Input  kmeans C-api2!");
        exit(1);
    }
    row = (int)PyList_Size(tMatPy);
    col = (int)PyList_Size(PyList_GetItem(tMatPy, 0));
    tMatC = pyToCmat(tMatPy, row, col); /* convert T mat from Py to C */
    centroidsArrayC = (int*)calloc(k, sizeof(int));
    for (i = 0; i < k; i++){
	centroidsArrayC[i] = (int)PyFloat_AsDouble(PyList_GetItem(initCentroidsPy, i));
    }
    outputClustersC = kmeans(tMatC, centroidsArrayC , row, col); /* send T matrix to main func to calculate k clusters in kmeans algo */
    resClustersPy = cToPyMat(outputClustersC, k, col); /* receive clusters from C */

    free2D(outputClustersC);
    return resClustersPy;
}

static PyMethodDef capiMethods[] = { /* an array of the module methods */ 
    {"pyToCtoPy", (PyCFunction)pyToCtoPy, METH_VARARGS, PyDoc_STR("receive input from Py and outputs C to Py")},
    {"spkCapi", (PyCFunction)spkCapi, METH_VARARGS, PyDoc_STR("first back and forth - data points ---> V matrix")}, 
    {"eigenGapCapi", (PyCFunction)eigenGapCapi, METH_VARARGS, PyDoc_STR("returns k (eigen heuristic) + sorted V matrix")},
    {"kmeansCapi", (PyCFunction)kmeansCapi, METH_VARARGS, PyDoc_STR("second back and forth - T mat + k indices ---> k clusters")},
    {NULL, NULL, 0, NULL} 
};

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT, "spkmeansmodule", "spkmeans module", -1, capiMethods
};

PyMODINIT_FUNC PyInit_spkmeansmodule(void)
{
    PyObject *m;

    m = PyModule_Create(&moduledef);
    if (!m)
    {
        return NULL;
    }
    return m;
}
