# Numerical Differentiation

## Description

Here we research computational errors of different numerical differentiation methods (for function's 1st derivative).

In the lab we consider next functions:

0) $sin(x^2)$
1) $cos(sin(x)$
2) $exp(sin(cos(x)$
3) $ln(x + 3)$
4) $(x + 3)^{0.5}$

And approximate their 1st derivatives with methods:

0) $\frac{f(x+h) - f(x)}{h}$
1) $\frac{f(x) - f(x-h)}{h}$
2) $\frac{f(x+h) - f(x-h)}{2h}$
3) $\frac{4}{3}\frac{f(x+h) - f(x-h)}{2h} -\frac{1}{3}\frac{f(x+2h) - f(x-2h)}{4h}$
4) $\frac{3}{2}\frac{f(x+h) - f(x-h)}{2h} -\frac{3}{5}\frac{f(x+2h) - f(x-2h)}{4h} + \frac{1}{10}\frac{f(x+3h) - f(x-3h)}{6h}$

## Structure of computational error

Computational error ~ `$h^{n} + \frac{Epsilon}{h}$`

* $h^(n)$ part in general can be decreased by choosing methods that use more points (like method 4)

* $\frac{Epsilon}{h}$ part is initially small, but becomes valuable when h gets close to Epsilon (machine zero)


## Observations

For mentioned functions we compute difference between analytically computed derivative and numerical computed derivative.

![](plots/cos(sin(x)).svg)
![](plots/exp(sin(cos(x))).svg)
![](plots/ln(x+3).svg)
![](plots/sin(x^2).svg)
![](plots/(x+3)**0.5.svg)

 In these plots we observe similar trends:

+ methods that use more points give more accurate results
+ at first error decreases with decrease of step size (h)
+ then error increases from some point in methods that use more points

It happens because `$h^{n}$` part becomes small. So $\frac{Epsilon}{h}$ part becomes significant. And this part increases with decrease of step size.

That's why we see this extremes in methods 1-4. In method 0 we don't see it as error of $h^{n}$ is much bigger than $\frac{Epsilon}{h}$ part.

