import numpy as np
import matplotlib.pyplot as plt


'''

Im gonna have next object in my program:

	Different methods of solving system of linear equations

	Information that characterize concrete СЛАУ

class LinearEqSystemSolver:

	+ __init__ (A, b)

	+ SolveWithGaussian

	+ SolveWithLU

	+ SolveWithZeidel

	+ SolveWithJacoby

	+ SolveWithUpperRelaxation

'''


# Create matrix from the task: г (гэ)

def construct_task_matrix(N):

	A = np.triu(np.ones((N, N),  dtype = np.float64), -1)
	A = A + np.identity(N, dtype = np.float64) * 9

	A[0][0] = A[0][0] - 9
	A[N-1][N-1] = A[N-1][N-1] - 9

	b = np.arange(100, 0, -1, dtype = np.float64)

	return A, b

if __name__ == '__main__':

	A, b = construct_task_matrix(100)

	print(A)
	print(b)
