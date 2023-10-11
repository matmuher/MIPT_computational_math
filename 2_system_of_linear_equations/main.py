import numpy as np
import matplotlib.pyplot as plt
from abc import ABC, abstractmethod

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

class LinearEqSystemSolver(ABC):

	@abstractmethod
	def solve(self, A, b):
		pass

# class GaussianSovler(LinearEqSystemSolver):

# 	def solve(slef, A, b):



# Create matrix from the task: y

def construct_task_matrix(N):

	A = np.identity(N, dtype = np.float64)

	for i in range(1, N+1):

		for j in range(1, N+1):

			if i != j:
				
				A[i-1][j-1] = 1 / (i*i + j) 

	b = np.arange(1, N+1, dtype = np.float64)
	
	for i in range(1, N+1):

		b[i-1] = 1/(i)

	return A, b

if __name__ == '__main__':

	A, b = construct_task_matrix(12)

	print(A)
	print(b)

