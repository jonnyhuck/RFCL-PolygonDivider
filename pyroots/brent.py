#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# file pyroots/brent.py
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
Brent's algorithm for root finding.
"""

from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals
from __future__ import absolute_import

from .utils import EPS, nearly_equal
from .base import BaseSolver


class _Brent(BaseSolver):

    def _extrapolation(self, *args, **kwargs):
        raise NotImplementedError

    def _solve(self, f, xa, xb, *args, **kwargs):
        # local names
        xtol = self.xtol
        _extrapolate = self._extrapolate

        # initialize counters
        i = 0
        x_steps = []
        fx_steps = []

        # rename variables in order to be consistent with scipy's code.
        xpre, xcur = xa, xb
        xblk, fblk, spre, scur = 0, 0, 0, 0

        #check that the bracket's interval is sufficiently big.
        if nearly_equal(xa, xb, xtol):
            return self._return_result(None, None, i, x_steps, fx_steps, False, "small bracket")

        # check lower bound
        fpre = f(xpre, *args, **kwargs)             # First function call
        x_steps.append(xpre)
        fx_steps.append(fpre)
        if self.is_root(fpre):
            return self._return_result(xpre, fpre, i, x_steps, fx_steps, True, "lower bracket")

        # check upper bound
        fcur = f(xcur, *args, **kwargs)             # Second function call
        x_steps.append(xcur)
        fx_steps.append(fcur)
        self._debug(i, len(fx_steps), xpre, xcur, fpre, fcur)
        if self.is_root(fcur):
            return self._return_result(xcur, fcur, i, x_steps, fx_steps, True, "upper bracket")

        # check if the root is bracketed.
        if fpre * fcur > 0.0:
            return self._return_result(None, None, i, x_steps, fx_steps, False, "no bracket")

        # start iterations
        for i in range(self.max_iter):
            if (fpre*fcur < 0):
                xblk = xpre
                fblk = fpre
                spre = scur = xcur - xpre

            if (abs(fblk) < abs(fcur)):
                xpre = xcur
                xcur = xblk
                xblk = xpre
                fpre = fcur
                fcur = fblk
                fblk = fpre

            # check for convergence
            #if self.is_root(fcur):
                #return self._return_result(xcur, fcur, i + 1, x_steps, fx_steps, True, "convergence")

            # check bracket
            sbis = (xblk - xcur) / 2;
            if abs(sbis) < xtol:
                return self._return_result(xcur, fcur, i + 1, x_steps, fx_steps, False, "small bracket")

            # calculate short step
            #self.logger.debug("spre %f; fcur %f; fpre %f; xblk %f; sbis %f", spre, fcur, fpre, xblk, sbis)
            if abs(spre) > xtol and abs(fcur) < abs(fpre):
                if xpre == xblk:
                    # interpolate
                    stry = -fcur * (xcur - xpre) / (fcur - fpre)
                    #self.logger.debug("Interpolate: stry %f", stry)
                else:
                    # extrapolate
                    dpre = (fpre - fcur) / (xpre - xcur)
                    dblk = (fblk - fcur) / (xblk - xcur)
                    stry = _extrapolate(fcur, fpre, fblk, dpre, dblk)
                    #self.logger.debug("Extrapolate: stry %f", stry)

                # check short step
                if (2 * abs(stry) < min(abs(spre), 3 * abs(sbis) - xtol)):
                    # good short step
                    spre = scur
                    scur = stry
                    #self.logger.debug("Good step")
                else:
                    # bisect
                    spre = sbis
                    scur = sbis
                    #self.logger.debug("Bad step, bisecting")
            else:
                # bisect
                spre = sbis
                scur = sbis

            xpre = xcur;
            fpre = fcur;
            if (abs(scur) > xtol):
                xcur += scur
            else:
                xcur += xtol if (sbis > 0) else -xtol

            fcur = f(xcur, *args, **kwargs)     # function evaluation
            x_steps.append(xcur)
            fx_steps.append(fcur)
            self._debug(i + 1, len(fx_steps), xpre, xcur, fpre, fcur)
            if self.is_root(fcur):
                return self._return_result(xcur, fcur, i, x_steps, fx_steps, True, "convergence")

        return self._return_result(xcur, fcur, i + 1, x_steps, fx_steps, False, "iterations")


class Brentq(_Brent):
    """
    Defines a Solver for the equation `f(x) = 0` in the interval `[xa, xb]` using Brent's Method.

    Function `f` must be solvable in `[xa, xb]`. Also `f(xa)` and `f(xb)` must
    have different signs.

    """

    def __init__(self, epsilon=1e-6, xtol=EPS, max_iter=500, raise_on_fail=True, debug_precision=10):
        super(Brentq, self).__init__(
            epsilon=epsilon,
            xtol=xtol,
            max_iter=max_iter,
            raise_on_fail=raise_on_fail,
            debug_precision=debug_precision,
            solver_name="Brentq"
        )

    def _extrapolate(self, fcur, fpre, fblk, dpre, dblk):
        return -fcur * (fblk * dblk - fpre * dpre) / (dblk * dpre * (fblk - fpre))


class Brenth(_Brent):
    """
    Defines a Solver for the equation `f(x) = 0` in the interval `[xa, xb]` using Brent's Method.

    Function `f` must be solvable in `[xa, xb]`. Also `f(xa)` and `f(xb)` must
    have different signs.

    """

    def __init__(self, epsilon=1e-6, xtol=EPS, max_iter=500, raise_on_fail=True, debug_precision=10):
        super(Brenth, self).__init__(
            epsilon=epsilon,
            xtol=xtol,
            max_iter=max_iter,
            raise_on_fail=raise_on_fail,
            debug_precision=debug_precision,
            solver_name="Brenth"
        )

    def _extrapolate(self, fcur, fpre, fblk, dpre, dblk):
        return -fcur * (fblk - fpre) / (fblk * dpre - fpre * dblk)
