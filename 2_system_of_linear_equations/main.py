import numpy as np
from numpy import linalg as LA

import matplotlib.pyplot as plt
from abc import ABC, abstractmethod

class LinearSolver(ABC):


	def __init__(self, name, A, b):

		self.name = name
		self.A = A
		self.b = b

	def concat(A, b):

		return np.column_stack((A, b))


	def unwind_U(equs):

		rowN = np.shape(equs)[0]
		columnN = np.shape(equs)[1]

		FillerValue = -7
		x = np.full(rowN, FillerValue, dtype = np.float64)

		for equ_index in range(rowN-1, -1, -1):

			consider_computed_x = 0

			for x_index in range(equ_index + 1, rowN):

				consider_computed_x += x[x_index] * equs[equ_index][x_index]

			x[equ_index] = (equs[equ_index][-1] - consider_computed_x) / equs[equ_index][equ_index]

		return x


	def unwind_L(equs):

		rowN = np.shape(equs)[0]
		columnN = np.shape(equs)[1]

		FillerValue = -7
		x = np.full(rowN, FillerValue, dtype = np.float64)

		for equ_index in range(0, rowN, 1):

			consider_computed_x = 0

			for x_index in range(0, equ_index):

				consider_computed_x += x[x_index] * equs[equ_index][x_index]

			x[equ_index] = (equs[equ_index][-1] - consider_computed_x) / equs[equ_index][equ_index]

		return x

	def swap_rows(equs, i, j):

		equs[[i, j]] = equs[[j, i]]


	def print_matrix(matrix):

		for i in range(0, np.shape(matrix)[0]):

			for j in range(0, np.shape(matrix)[1]):

				print(f'\t{matrix[i][j]:.2}', end = ' ')

			print('\n')

	@abstractmethod
	def solve(self):
		pass

	def compute_residual(x, A, b):

		return LA.norm(A.dot(x) - b)

	def print_residual(self, x):

		residual = LinearSolver.compute_residual(x, self.A, self.b)
		print(f'[{self.name}] residual is {residual}')

'''

1. Exclude x1 from down equations. By adding 1st equation muliplied by -a11/ai1.
2. Then exlide x2 in the same way.
...
3. We get triangle matrix. Answer can be got by walking from down to up.

'''

class GaussianSolver(LinearSolver):

	def __init__(self, A, b):

		LinearSolver.__init__(self, 'Gauss', A, b)

	def solve(self):

		equs = LinearSolver.concat(self.A, self.b)
		
		N = np.shape(equs)[0]

		# Прямой ход алгоритма Гаусса

		for exluder_equ_index in range(0, N-1):

			exluder_equ = equs[exluder_equ_index]

			# find element to maxmize value in denominator

			for current_equ_index in range(exluder_equ_index,  N):

				current_equ = equs[current_equ_index] 

				if current_equ[exluder_equ_index] > exluder_equ[exluder_equ_index]:

					swap_rows(equs, exluder_equ_index, current_equ_index)

			for current_equ_index in range(exluder_equ_index + 1,  N):

				current_equ = equs[current_equ_index] 

				if current_equ[exluder_equ_index] == 0:

					continue

				mutiplier = - current_equ[exluder_equ_index] / exluder_equ[exluder_equ_index]
				current_equ += exluder_equ * mutiplier

		# Обратный ход

		x = LinearSolver.unwind_U(equs)

		return x

'''

Ax = LUx = b

Ux = y

Ly = b -> unwind_L

Ux = y -> unwind_U

'''

class LUSolver(LinearSolver):

	def __init__(self, A, b):

		LinearSolver.__init__(self, 'LU decomposition', A, b)

	def LU_decomposition(self):

		N = np.shape(self.A)[0]

		L = np.copy(self.A)
		U = np.copy(self.A)

		for i in range(0, N):

			for j in range(0, N):

				L[i][j] = 0
				U[i][j] = 0

				if i == j:

					L[i][j] = 1


		for i in range(0, N):

			for j in range(0, N):

				if i <= j:

					U[i][j] = self.A[i][j] - np.dot(L[i], U.T[j])

				else:

					L[i][j] = (self.A[i][j] - np.dot(L[i], U.T[j])) / U[j][j]

		return L, U

	def solve(self):

		L, U = LUSolver.LU_decomposition(self)

		y = LinearSolver.unwind_L(LinearSolver.concat(L, self.b))
		x = LinearSolver.unwind_U(LinearSolver.concat(U, y))

		return x


class IterativeSolver(LinearSolver):

	count = 0
	ax = plt.gca()


	ax.set_yscale('log')

	ax.set_ylabel('Residual')
	ax.set_xlabel('Iteration number')

	ax.set_title(f'Residuals for iterative methods of solving linear systems')


	def __init__(self, name, A, b, B, F):

		LinearSolver.__init__(self, name, A, b)

		self.B = B
		self.F = F
		self.Id = IterativeSolver.count

		IterativeSolver.count += 1

	def solve(self, x0, N):

		x_prev = x0

		residual_list = []

		for i in range(N):

			x = self.F + self.B.dot(x_prev)
			x_prev = x

			residual_list.append(LinearSolver.compute_residual(x, self.A, self.b))

		iteration_list = np.arange(N)

		IterativeSolver.ax.plot(iteration_list, residual_list, label = f'{self.name} method')

'''

https://en.wikipedia.org/wiki/Jacobi_method

'''

class JacobiMethod(IterativeSolver):

	def __init__(self, A, b):

		U = np.triu(A, k = 1)
		L = np.tril(A, k = -1)
		D = np.diag(np.diag(A))

		D_inv = LA.inv(D)

		B = -np.matmul(D_inv, L + U)
		F = D_inv.dot(b)

		IterativeSolver.__init__(self, 'Jacobi', A, b, B, F)

'''

https://en.wikipedia.org/wiki/Gauss%E2%80%93Seidel_method

'''

class SeidelMethod(IterativeSolver):

	def __init__(self, A, b):

		U = np.triu(A, k = 1)
		LD = np.tril(A, k = 0) # L + D

		LD_inv = LA.inv(LD)

		B = -np.matmul(LD_inv, U)
		F = LD_inv.dot(b)

		IterativeSolver.__init__(self, 'Seidel', A, b, B, F)


'''

https://en.wikipedia.org/wiki/Successive_over-relaxation

'''

class SORMethod(IterativeSolver):

	def __init__(self, A, b, w):

		U = np.triu(A, k = 1)
		L = np.tril(A, k = -1)
		D = np.diag(np.diag(A))

		B = -np.matmul(LA.inv(D + w*L), (w-1)*D + w*U)
		F = w * LA.inv(D + L*w).dot(b)

		IterativeSolver.__init__(self, f'SOR {w}', A, b, B, F)



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

	# print(f'det(A) = {np.linalg.det(A)}')

	return A, b

if __name__ == '__main__':

	A, b = construct_task_matrix()

	print('task matrix:')
	LinearSolver.print_matrix(LinearSolver.concat(A, b))

	straight_solvers = [LUSolver(A, b), GaussianSolver(A, b)] 

	for solver in straight_solvers:
		
		solver.print_residual(solver.solve())

	iterative_solvers = [
						SeidelMethod(A, b),
						JacobiMethod(A, b),
						SORMethod(A, b, 1.2),
						SORMethod(A, b, 1.5)
						]

	N = 80

	for solver in iterative_solvers:

		solver.solve(b, N)

	IterativeSolver.ax.legend()

	plt.savefig(f'methods_cmp.svg')