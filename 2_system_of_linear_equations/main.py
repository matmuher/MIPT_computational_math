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

'''

1. Exclude x1 from down equations. By adding 1st equation muliplied by -a11/ai1.
2. Then exlide x2 in the same way.
...
3. We get triangle matrix. Answer can be got by walking from down to up.

'''

class GaussianSolver(LinearEqSystemSolver):

	def solve(self, A, b):

		equs = np.column_stack((A, b))

		N = b.size

		for exluder_equ_index in range(0, N-1):

			exluder_equ = equs[exluder_equ_index]

			for current_equ_index in range(exluder_equ_index + 1,  N):

				current_equ = equs[current_equ_index] 

				if current_equ[exluder_equ_index] == 0:

					continue

				mutiplier = - current_equ[exluder_equ_index] / exluder_equ[exluder_equ_index]
				current_equ += exluder_equ * mutiplier

		print_matrix(equs)

def print_matrix(matrix):

	for i in range(0, np.shape(matrix)[0]):

		for j in range(0, np.shape(matrix)[1]):

			print(f'\t{matrix[i][j]:.2}', end = ' ')

		print('\n')

# Create matrix from the task: y

def construct_task_matrix():

	N = 12

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

	A, b = construct_task_matrix()

	solver = GaussianSolver()
	solver.solve(A, b)
