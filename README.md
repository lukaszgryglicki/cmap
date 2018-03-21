# cmap
C++ program for generatin complex map in EPS format

# description

This is code based on **Kouznetsov**: conto.cin, ado.cin.

Basically You compile it via c++ like (in Makefile):
c++ -o cmap cmap.cc libfparser.cc -lm

It generates EPS file (like PDF) with complex map (constant levels, real=red, imag=blue, modulo=green) for given function.

You pass real range, imag range, number of grid lines, levels to draw etc.
You also give function name (last argument) as literal (put it in ' ' or " ").
There is a function parser included (all functions are compelx).
So you can pass for example:
`'sin(x)' or 'exp(exp(x*2)-log(x)' or 'i^x' or '(3+2*i)^(1/x^3)'` etc.

There are example scripts how to call program (see files `*.sh`).
For example:
`./cmap 'out.eps' 500 500 -4 4 -4 4 8 8 -10 -10 0 10 10 20 -2 2 .1 -2 2 .1 0 4 .1 1 1 1 'log(log(x))'`.

Also there are all options commented in cmap.cc in `main()` function and `eps_cont()` function.
libfparser.h and libfparser.cc are my function parser files, see libfparser.h to see which functions are supported (by sin(x) I mean complex sine function of course).

Files ado.cc and conto.cc are just renamed ado.cin and conto.cin from TORI webpage (**by Kouznetsov**).

You may find it useful (I've used ado.cin and conto.cin - because the're available for all and autor allows to use them).
I've added function parser and added all parameters to control what area to draw, what levels, cut off values etc - see details in cmap.cc.

I hope someone find it usefull - examples in `*.sh` files are generating some nice fractals.

