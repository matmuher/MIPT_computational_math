# Numerical Differentiation

## Description

In this lab we research computational errors of different formulas for function's 1st derivative.

We consider next functions:

0) $sin(x^2)$
1) $ cos(sin(x) $
2) $ exp(sin(cos(x) $
3) $ ln(x + 3) $ 
4) $ (x + 3)^0.5 $

And approximate their 1st derivatives with formulas:

0) 
1)
2)
3)
4)

## Structure of computational error

Computational error ~ `h^(n) + Epsilon / h`

* `h^(n)` part in general can be decreased by choosing methods that uses more points (like method 4)

* `Epsilon / h` part is valuable when h is close to or less than Epsilon (machine zero)


## Observations

![](plots/cos(sin(x)).svg)
![](plots/exp(sin(cos(x))).svg)
![](plots/ln(x+3).svg)
![](plots/sin(x^2).svg)
![](plots/(x+3)**0.5.svg)

In all plots we observe the similar trends:

+ method that uses the most amount of function points gives the least computation error
+ with decrease of step (h) computation error at first decreases
+ in methods that use more points computation error increases from some point

It happens because `h^(n)` part becomes small enough to `Epsilon / h` part become valuable.

That's why we see this extremes in methods 1-4. In method 0 we don't see it as error of method is 

