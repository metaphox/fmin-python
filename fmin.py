"""
DO WHAT THE FUCK YOU WANT TO PUBLIC LICENSE 
                    Version 2, December 2004 

 Copyright (C) 2004 Sam Hocevar <sam@hocevar.net> 

 Everyone is permitted to copy and distribute verbatim or modified 
 copies of this license document, and changing it is allowed as long 
 as the name is changed. 

            DO WHAT THE FUCK YOU WANT TO PUBLIC LICENSE 
   TERMS AND CONDITIONS FOR COPYING, DISTRIBUTION AND MODIFICATION 

  0. You just DO WHAT THE FUCK YOU WANT TO.
"""

from math import sqrt

def fmin(a, b, func, tol):
    """
        /*
         *
         * Here is a copy of the Netlib documentation:
         *
         * c c An approximation x to the point where f attains a minimum on c
         * the interval (ax,bx) is determined. c c input.. c c ax left endpoint
         * of initial interval c bx right endpoint of initial interval c f
         * function subprogram which evaluates f(x) for any x c in the interval
         * (ax,bx) c tol desired length of the interval of uncertainty of the
         * final c result (.ge.0.) c c output.. c c fmin abcissa approximating
         * the point where f attains a c minimum c c The method used is a
         * combination of golden section search and c successive parabolic
         * interpolation. Convergence is never much slower c than that for a
         * Fibonacci search. If f has a continuous second c derivative which is
         * positive at the minimum (which is not at ax or c bx), then
         * convergence is superlinear, and usually of the order of c about
         * 1.324. c The function f is never evaluated at two points closer
         * together c than eps*abs(fmin)+(tol/3), where eps is approximately the
         * square c root of the relative machine precision. If f is a unimodal c
         * function and the computed values of f are always unimodal when c
         * separated by at least eps*abs(x)+(tol/3), then fmin approximates c
         * the abcissa of the global minimum of f on the interval (ax,bx) with c
         * an error less than 3*eps*abs(fmin)+tol. If f is not unimodal, c then
         * fmin may approximate a local, but perhaps non-global, minimum to c
         * the same accuracy. c This function subprogram is a slightly modified
         * version of the c Algol 60 procedure localmin given in Richard Brent,
         * Algorithms For c Minimization Without Derivatives, Prentice-Hall,
         * Inc. (1973). c
         */
    """

    c = .5 * (3.0 - sqrt(5.0))
    d = 0.0

    # 1.1102e-16 is machine precision
    eps = 1.2e-16
    eps = sqrt(eps)

    v = a + c * (b - a)
    w = v
    x = v
    e = 0.0
    fx = func(x)
    fv = fx
    fw = fx
    tol3 = tol / 3.0

    xm = .5 * (a + b)
    tol1 = eps * abs(x) + tol3
    t2 = 2.0 * tol1

    # main loop
    while abs(x - xm) > (t2 - .5 * (b - a)):
        p = q = r = 0.0
        if abs(e) > tol1:
            # fit the parabola
            r = (x - w) * (fx - fv)
            q = (x - v) * (fx - fw)
            p = (x - v) * q - (x - w) * r
            q = 2.0 * (q - r)

            if (q > 0.0):
                p = -p
            else:
                q = -q

            r = e
            e = d

        if abs(p) < abs(.5 * q * r) and p > q * (a - x) and p < q * (b - x):
            # a parabolic interpolation step
            d = p / q
            u = x + d
            # f must not be evaluated too close to a or b
            if (u - a) < t2 or (b - u) < t2:
                d = tol1
                if x >= xm:
                    d = -d
        else:
            # a golden-section step
            if x < xm:
                e = b - x
            else:
                e = a - x
            d = c * e

        # f must not be evaluated too close to x
        if abs(d) >= tol1:
            u = x + d
        else:
            if d > 0.0:
                u = x + tol1
            else:
                u = x - tol1

        fu = func(u)

        # Update a, b, v, w, and x
        if fx <= fu:
            if u < x:
                a = u
            else:
                b = u

        if fu <= fx:
            if u < x:
                b = x
            else:
                a = x

            v = w
            fv = fw
            w = x
            fw = fx
            x = u
            fx = fu

            xm = .5 * (a + b)
            tol1 = eps * abs(x) + tol3
            t2 = 2.0 * tol1
        else:
            if fu <= fw or w == x:
                v = w
                fv = fw
                w = u
                fw = fu

                xm = .5 * (a + b)
                tol1 = eps * abs(x) + tol3
                t2 = 2.0 * tol1
            elif fu > fv and v != x and v != w:
                xm = .5 * (a + b)
                tol1 = eps * abs(x) + tol3
                t2 = 2.0 * tol1
            else:
                v = u
                fv = fu

                xm = .5 * (a + b)
                tol1 = eps * abs(x) + tol3
                t2 = 2.0 * tol1
    return x
