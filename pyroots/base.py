#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# file pyroots/base.py
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
BaseSolver class
"""

from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals
from __future__ import absolute_import

import abc
import logging

from .utils import Result, ConvergenceError, LOG_MSG, EPS, nearly_equal


class BaseSolver(object):
    """ Base Solver for pyroots. """

    messages = {
        "small bracket": "Bracket is smaller than tolerance.",
        "lower bracket": "Root is equal to the lower bracket",
        "upper bracket": "Root is equal to the upper bracket",
        "no bracket": "Root is not bracketed.",
        "convergence": "Solution converged.",
        "iterations": "Exceeded max iterations.",
    }

    def __init__(self, epsilon=1e-6, xtol=EPS, max_iter=500, raise_on_fail=True, solver_name="BaseSolver", debug_precision=10):
        """
        Parameters
        ----------
        :param float ftol:
            The required accuracy for `f`.
        :param float xtol:
            Equals machine accuracy.
        :param int max_inter:
            The maximum allowed number of iterations.

        """
        # sanity check
        if xtol < EPS:
            raise ArithmeticError("'xtol' can't be smaller than EPS (xtol=%.15f, EPS=%.15f)" % (xtol, EPS))
        if epsilon < EPS:
            raise ArithmeticError("'epsilon' can't be smaller than EPS (epsilon=%.15f, EPS=%.15f)" % (xtol, EPS))
        if (not isinstance(max_iter, int)) or max_iter < 0:
            raise ArithmeticError("max_iter must be a positive integer, not: %r <%r>" % (max_iter, type(max_iter)))

        self.log_msg = LOG_MSG.format(precision=debug_precision)
        self.xtol = xtol
        self.epsilon = epsilon
        self.max_iter = max_iter
        self.raise_on_fail = raise_on_fail
        self.solver_name = solver_name
        self.logger = logging.getLogger("pyroots.{solver_name}".format(solver_name=solver_name))

    def __repr__(self):
        return """ `Pyroots.{solver_name}`:\n
              xtol : {xtol}
           epsilon : {epsilon}
        iterations : {max_iter}
             raise : {raise_on_fail}
        """.format(**self.__dict__)

    def _return_result(self, x0, fx0, iterations, x_steps, fx_steps, converged, condition):
        msg = self.messages[condition]
        result = Result(x0, fx0, iterations, converged, self.xtol, self.epsilon, x_steps, fx_steps, msg)
        if not result.converged and self.raise_on_fail:
            self.logger.info("Solution did not converge: %r", result)
            raise ConvergenceError(msg)
        else:
            self.logger.info("Solution converged: %r", result)
            return result

    def _debug(self, i, fcalls, xa, xb, fa, fb):
        self.logger.debug(self.log_msg, i, fcalls, xa, xb, xb - xa, fa, fb)

    def __call__(self, f, xa, xb, *args, **kwargs):
        """
        Parameters
        ----------
        :param function f:
            A single-variable function that we are searching for it's root within
            the given interval. It must take at least one argument and any number
            of positional or keyword arguments.
        :param float xa:
            The lower bound of the interval.
        :param float xb:
            The upper bound of the interval.
        :param tuple *args:
            Function's `f` positional arguments.
        :param dict kwargs:
            Function's `f` keyword arguments.

        Returns
        -------

        :returns: `Result`'s instance.

        :raises: `ConvergenceError` if `raise_on_fail` is True and the function fails to
            converge. If `raise_on_fail` is `False` then it returns a `Result` instance.

        """
        return self._solve(f, xa, xb, *args, **kwargs)

    def is_root(self, root):
        """
        Return True if the root is sufficiently close to 0.

        Just a wrapper around `nearly_equal()`.
        """
        return nearly_equal(0, root, self.epsilon)

    @abc.abstractmethod
    def _solve(self, f, xa, xb, *args, **kwargs):
        """ Return a result object or raise a ConvergenceError. """

