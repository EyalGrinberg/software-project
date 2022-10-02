# software-project
spectral k-means clustering
This is a project from the course Software Project from my third year in Tel Aviv University. This code performs a spectral k-means clustering algorithm to a list of coordinates:

image

In order to compile and run the code, run the following commands after cloning the repository:

gcc -ansi -Wall -Wextra -Werror -pedantic-errors spkmeans.c -lm -o spkmeans
python setup.py build_ext --inplace
If you want to run the code, create a txt/csv file with points in equal amount of coordinates in each line, separated by commas only, and run the following: The program (both on C and python) gets the following arguments:

k - number of clusters (If k = 0, will run the eigengap heuristics algorithm to find it)
goal - one of the following
spk - run full spectral clustering
wam - calculate the weighted adjacency matrix and print it
ddg - calculate the diagonal degree matrix and print it
lnorm - calculate lnorm matrix and print it
jacobi - in this case, the file should be a matrix to calculate its eigenvalues and eigenvectors using Jacobi's algorithm
file_name - the file where the data is in
Running on Python, use spkmeans.py, and running on C, use ./spkmeans
