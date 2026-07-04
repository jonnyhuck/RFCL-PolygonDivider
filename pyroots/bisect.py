#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# file pyroots/bisect.py
#
#############################################################################
# Copyright (c) 2014 by Panagiotis Mavrogiorgos
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
bisect method.
"""

from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals
from __future__ import absolute_import

from math import copysign

from .utils import EPS, nearly_equal
from .base import BaseSolver


class Bisect(BaseSolver):
    """
    Defines a Solver for the equation `f(x) = 0` in the interval `[xa, xb]` using the Bisection Method.

    Function `f` must be solvable in `[xa, xb]`. Also `f(xa)` and `f(xb)` must
    have different signs.

    """

    def __init__(self, epsilon=1e-6, xtol=EPS, max_iter=500, raise_on_fail=True, debug_precision=10):
        super(Bisect, self).__init__(
            epsilon=epsilon,
            xtol=xtol,
            max_iter=max_iter,
            raise_on_fail=raise_on_fail,
            debug_precision=debug_precision,
            solver_name="Bisect"
        )

    def _solve(self, f, xa, xb, *args, **kwargs):
        """ Bisect implementation.  """
        # local names
        xtol = self.xtol

        # initialize counters
        i = 0
        x_steps = []
        fx_steps = []

        # check that the bracket's interval is sufficiently big.
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

            # close the bracket
            if copysign(1, fm) == copysign(1, fa):
                xa = xm
                fa = fm
            else:
                xb = xm
                fb = fm
            self._debug(i, len(fx_steps), xa, xb, fa, fb)

            # check for convergence.
            if self.is_root(fm):
                return self._return_result(xm, fm, i, x_steps, fx_steps, True, "convergence")

            # check for the new bracket size.
            if nearly_equal(xa, xb, xtol):
                return self._return_result(xm, fm, i, x_steps, fx_steps, False, "small bracket")

        return self._return_result(xm, fm, i, x_steps, fx_steps, False, "iterations")
