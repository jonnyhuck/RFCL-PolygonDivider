#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# file pyroots/ridder.py
#
#############################################################################
# Copyright (c) 2013 by Panagiotis Mavrogiorgos
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# * Redistributions of source code must retain the above copyright notice,
#   this list of conditions and the following disclaimer.
# * Redistributions in binary form must reproduce the above copyright notice,
#   this list of conditions and the following disclaimer in the documentation
#   and/or other materials provided with the distribution.
# * Neither the name(s) of the copyright holders nor the names of its
#   contributors may be used to endorse or promote products derived from this
#   software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AS IS AND ANY EXPRESS OR
# IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
# MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO
# EVENT SHALL THE COPYRIGHT HOLDERS BE LIABLE FOR ANY DIRECT, INDIRECT,
# INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA,
# OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
# LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
# EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#############################################################################
#
# @license: http://opensource.org/licenses/BSD-3-Clause
# @authors: see AUTHORS.txt


"""
Ridder's algorithm for root finding.
"""

from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals
from __future__ import absolute_import

from math import sqrt

from .utils import EPS, nearly_equal
from .base import BaseSolver


class Ridder(BaseSolver):
    """
    Defines a Solver for the equation `f(x) = 0` in the interval `[xa, xb]` using Ridder's Method.

    Function `f` must be solvable in `[xa, xb]`. Also `f(xa)` and `f(xb)` must
    have different signs.

    """

    def __init__(self, epsilon=1e-6, xtol=EPS, max_iter=500, raise_on_fail=True, debug_precision=10):
        super(Ridder, self).__init__(
            epsilon=epsilon,
            xtol=xtol,
            max_iter=max_iter,
            raise_on_fail=raise_on_fail,
            debug_precision=debug_precision,
            solver_name="Ridder"
        )

    def _solve(self, f, xa, xb, *args, **kwargs):
        """ Ridder implementation.  """
        # local names
        xtol = self.xtol

        # initialize counters
        i = 0
        x_steps = []
        fx_steps = []

        #check that the bracket's interval is sufficiently big.
        if nearly_equal(xa, xb, xtol):
            return self._return_result(None, None, i, x_steps, fx_steps, False, "small bracket")

        # check lower bound
        fa = f(xa, *args, **kwargs)               # First function call
        x_steps.append(xa)
        fx_steps.append(fa)
        if self.is_root(fa):
            return self._return_result(xa, fa, i, x_steps, fx_steps, True, "lower bracket")

        # check upper bound
        fb = f(xb, *args, **kwargs)               # Second function call
        x_steps.append(xb)
        fx_steps.append(fb)
        self._debug(i, len(fx_steps), xa, xb, fa, fb)
        if self.is_root(fb):
            return self._return_result(xb, fb, i, x_steps, fx_steps, True, "upper bracket")

        # check if the root is bracketed.
        if fa * fb > 0.0:
            return self._return_result(None, None, i, x_steps, fx_steps, False, "no bracket")

        # start iterations
        for i in range(1, self.max_iter + 1):
            # Bisect the bracket and calculate the new function value.
            xm = 0.5 * (xa + xb)
            fm = f(xm, *args, **kwargs)           # New function call.
            x_steps.append(xm)
            fx_steps.append(fm)
            self._debug(i, len(fx_steps), xa, xm, fa, fm)

            # check for convergence.
            if self.is_root(fm):
                return self._return_result(xm, fm, i, x_steps, fx_steps, True, "convergence")

            # `t` is the denominator followingly
            # if `t == 0` then the ridder's method cannot be applied due to a
            # ZeroDivisionError. Not sure when this can happen though...
            # Perhaps on an almost vertical line.
            # TODO write a test for this case.
            # NOTE 1: Scipy's ridder method doesn't check for t == 0.
            # NOTE 2: `t` is always positive, so we do not need a try/except block
            #         to check if the value is positive.
            t = sqrt(fm ** 2 - fa * fb)
            #if t == 0.0:
                #return Result(xm, fm, i + 1, 3 + 2 * i, False, "Solution is not possible.")

            # calculate the improved x from Ridder's formula and
            # NOTE: Scipy's version is a little bit different. Didn't find any
            # reference though.
            sign = -1 if fa < fb else 1
            xs = xm + (xm - xa) * sign * fm / t
            fs = f(xs, **kwargs)
            x_steps.append(xs)
            fx_steps.append(fs)
            self._debug(i, len(fx_steps), xa, xs, fa, fs)

            if self.is_root(fs):
                return self._return_result(xs, fs, i, x_steps, fx_steps, True, "convergence")

            # When ftol is very small (e.g. 1e-15) then there are cases that the
            # method can't converge in a reasonable amount of iterations.
            # One such case is f(x) = 24.5 * x**5 - 92.2 * x**3 + 23 * x - 12
            # In order to prevent the method from becoming stagnant we keep record
            # of the old values of xm and xs and we check if they change values
            # during the iterations.
            # NOTE: Perhaps this check is not very robust.
            if i > 1 and abs(xs - xs_old) < xtol and abs(xm - xm_old) < xtol:
                result = self._return_result(xs, fs, i, x_steps, fx_steps, False, "Precision not achieved. Iteration stagnant.")
                self.logger.debug(result)
                return result


            # Re-bracket the root as tightly as possible
            if fm * fs > 0.0:
                if fa * fs < 0.0:
                    xb, fb = xs, fs
                else:
                    xa, fa = xs, fs
            else:
                xa, fa = xm, fm
                xb, fb = xs, fs

            # Ensure that xa < xb
            # NOTE: This is not necessary but it is useful for debugging.
            if xa > xb:
                xa, xb = xb, xa
                fa, fb = fb, fa

            # check for the new bracket size.
            #print(abs(max(xa, xb)) * xtol)
            #if abs(xb - xa) < abs(max(xa, xb)) * xtol:
            if nearly_equal(xa, xb, xtol):
                return self._return_result(xs, fs, i, x_steps, fx_steps, False, "small bracket")

            # Store values of the previous iteration.
            xm_old = xm
            xs_old = xs

        return self._return_result(xm, fm, i, x_steps, fx_steps, False, "iterations")
