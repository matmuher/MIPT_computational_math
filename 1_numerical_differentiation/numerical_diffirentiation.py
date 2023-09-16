"""
IN:
	list of functions
	numerical derivative formulas

OUT:
	for every function in list of functions
	log plot of increase of error

Plan:

	1. Get Machine Error
	2. func(src_func, num_diff_func, src_diff_func)
	3. plot it
	4. make list of funcs
"""

import math as m
import matplotlib.pyplot as plt

#----------------------------------------------
# [UTILITYS]

def GetMachineZero(x):

	machineZero = 1.;

	lastNonZeroVal = machineZero

	while x != x + machineZero:

		lastNonZeroVal = machineZero
		machineZero /= 10

	return lastNonZeroVal

def GetNumDiffError(srcFunc, numDiffFunc, h, analytDiffFunc, x):
	
	computError = abs(numDiffFunc(srcFunc, h, x) - analytDiffFunc(x))

	return computError

def getH(n):

	return 2 ** (1 - n)

#----------------------------------------------
# [SRC FUNCTIONS] 

#---1---

def srcFunc1(x):

	return m.sin(x ** 2)

def analyticalDiffFunc1(x):
	return m.cos(x ** 2) * 2 * x

#---2---

def srcFunc2(x):

	return m.cos(m.sin(x))

def analyticalDiffFunc2(x):

	return -m.sin(m.sin(x)) * m.cos(x)

#---3---

def srcFunc3(x):

	return m.exp(m.sin(m.cos(x)))

def analyticalDiffFunc3(x):

	return srcFunc3(x) * m.cos(m.cos(x)) * (-m.sin(x))

#---4---

def srcFunc4(x):

	return m.log(x + 3)

def analyticalDiffFunc4(x):

	return 1 / (x + 3)

#---5---

def srcFunc5(x):

	return (x + 3) ** 0.5

def analyticalDiffFunc5(x):

	return 0.5 * (x + 3) ** (-0.5)

#----------------------------------------------

class Function:
	
	def __init__(self, srcFunc, srcAnalyticalDiffFunc, name):
	
		self.srcFunc = srcFunc
		self.srcAnalyticalDiffFunc = srcAnalyticalDiffFunc
		self.name = name

	def __str__(self):
		return self.name

#----------------------------------------------
# [DIF FUNCTIONS]

def DiffFunc1(srcFunc, h, x):

	return ( srcFunc(x + h) - srcFunc (x) ) / h

def DiffFunc2(srcFunc, h, x):
	
	return ( srcFunc(x) - srcFunc(x - h) ) / h

def DiffFunc3(srcFunc, h, x):

	return ( srcFunc(x + h) - srcFunc(x - h) ) / (2 * h) 

def DiffFunc4(srcFunc, h, x):

	return 3 ** (-1) * (
						4 * DiffFunc3(srcFunc, h, x) -
						 	DiffFunc3(srcFunc, 2 * h, x)
						)

def DiffFunc5(srcFunc, h, x):

	return 10 ** (-1) * (
						15 * DiffFunc3(srcFunc, h, x) 	  -
						6  * DiffFunc3(srcFunc, 2 * h, x) +
							 DiffFunc3(srcFunc, 3 * h, x)
						)
#----------------------------------------------
# [INITS]

def initNumDiffFuncsList():

	numDiffFuncs = []
	
	# 1
	numDiffFuncs.append(DiffFunc1)

	# 2
	numDiffFuncs.append(DiffFunc2)

	# 3
	numDiffFuncs.append(DiffFunc3)

	# 4
	numDiffFuncs.append(DiffFunc4)

	# 5
	numDiffFuncs.append(DiffFunc5)

	return numDiffFuncs

def initFuncs():

	funcs = []

	funcs.append(Function(srcFunc1, analyticalDiffFunc1, "sin(x^2)"))

	funcs.append(Function(srcFunc2, analyticalDiffFunc2, "cos(sin(x))"))

	funcs.append(Function(srcFunc3, analyticalDiffFunc3, "exp(sin(cos(x)))"))

	funcs.append(Function(srcFunc4, analyticalDiffFunc4, "ln(x+3)"))

	funcs.append(Function(srcFunc5, analyticalDiffFunc5, "(x+3)**0.5"))

	return funcs

def initSteps(N):

	steps = [getH(i) for i in range(1, 1 + N)]

	return steps
#----------------------------------------------

def plotComputeError(func, numDiffFuncs, steps, x):

	fg = plt.figure()
	ax = plt.gca()

	fig_height = 15
	fig_width = 15
	fg.set_size_inches(fig_width, fig_height)	

	for diffFuncIndex in range(len(numDiffFuncs)):

		computeErrorList = []

		for step in steps:

			computeError = GetNumDiffError(	func.srcFunc,
											numDiffFuncs[diffFuncIndex],
											step,
											func.srcAnalyticalDiffFunc,
											x)
			
			computeErrorList.append(computeError)

		ax.tick_params(labelsize = 20)

		ax.plot(steps, computeErrorList, marker = '_', label = f'method #{diffFuncIndex}')
		
		ax.set_ylabel('computational error', fontsize = 20)
		ax.set_xlabel('step size', fontsize = 20)

		ax.set_xscale('log', base = 2)
		ax.set_yscale('log', base = 10)

		ax.set_title(f'${func}$', fontsize = 20)
		ax.legend(fontsize = 20)

	plt.savefig(f'plots/{func}.svg')


if __name__ == '__main__':

	#configs

	N = 20
	x = 1.5

	steps = initSteps(N)

	numDiffFuncs = initNumDiffFuncsList()

	funcs = initFuncs()

	print(f'Machine zero for x = {x} is {GetMachineZero(x)}')

	for func in funcs:

		plotComputeError(func, numDiffFuncs, steps, x)