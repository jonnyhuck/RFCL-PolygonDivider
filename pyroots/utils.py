#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# file pyroots/utils.py
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
pyroots/utils.py
"""

from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals
from __future__ import absolute_import

import sys

# Constants
EPS = sys.float_info.epsilon
LOG_MSG = "Iter: %3d; fcall: %3d; x=[% .{precision}f, % .{precision}f]; Δx=% .{precision}f; f=[% .{precision}f, %+.{precision}f]"


class Result(object):
    """ Solver's result. Used for providing a summary. """

    _result_representation = """
 converged : {converged}
   message : {msg}
iterations : {iterations:3d}
func calls : {func_calls:3d}
        x0 : {x0: 22.16f}
      xtol : {xtol: 22.16f}
     f(x0) : {fx0: 22.16f}
   epsilon : {epsilon: 22.16f}
   x_steps : {x_steps}
  fx_steps : {fx_steps}
""".rstrip()

    _result_representation_x0_none = """
 converged : {converged}
   message : {msg}
iterations : {iterations:3d}
func calls : {func_calls:3d}
        x0 : None
      xtol : {xtol}
     f(x0) : None
   epsilon : {epsilon: 22.16f}
   x_steps : {x_steps}
  fx_steps : {fx_steps}
""".rstrip()

    def __init__(self, x0, fx0, iterations, converged, xtol, epsilon, x_steps, fx_steps, msg=""):
        self.x0 = x0
        self.fx0 = fx0
        self.iterations = iterations
        self.func_calls = len(fx_steps)
        self.converged = converged
        self.msg = msg
        self.xtol = xtol
        self.epsilon = epsilon
        self.x_steps = x_steps
        self.fx_steps = fx_steps

    def __repr__(self):
        if self.x0 is None:
            representation = self._result_representation_x0_none
        else:
            representation = self._result_representation
        return representation.format(**self.__dict__)


def nearly_equal(a, b, epsilon):
    """
    Return `True` if the "difference" between `a` and `b` is smaller than `epsilon`, `False` otherwise.

    This function tries to solve the problem of float comparison in a way that handles all cases
    (i.e. comparing floats with `inf`, `nan` etc).

    Do take note that when the difference between `a` and `b` is "equal" to `epsilon` then the
    results may not be what you expect them to be...  For example::

        a = 1
        b = a + 1e-7
        epsilon = 1e-7
        assert nearly_equal(a, b, epsilon)          # will raise an exception!

    In order to understand why this happens, we have to print lots of decimal digits::

        f = "{:.30f}".format
        print(f(a))                 # 1.000000000000000000000000000000
        print(f(b))                 # 1.000000100000000058386717682879
        print(f(epsilon))           # 0.000000099999999999999995474811

    As we can see, the float that represents the Real number "1e-7" is actually smaller than "1e-7".
    Similarly, the float that represents the Real number "1 + 1e-7" is actually greater than
    "1 + 1e-7".  So, when you try to compare "1" with "1 + 1e-7" using an "epsilon" value of
    "1e-7" the function will return `False` even though you would expect it to return `True`.

    Interestingly, if you reverse the sign of `a` this no longer happens::

        a = -1
        b = a + 1e-7
        epsilon = 1e-7
        assert nearly_equal(a, b, epsilon)          # No exceptions will be raised!

    Another thing that should be noted is that you shouldn't set epsilon to values lower than
    `sys.float_info.epsilon` (it's a value close to 2e-16).  If you do, you may run into things like
    these::

        import sys
        epsilon = sys.float_info.epsilon / 10    #
        a = 2
        b = 2 + 1e-16
        assert nearly_equal(a, b, epsilon)          # Will not raise an exception!!!

    """
    if a == b:
        return True                         # shortcut. Handles infinities etc
    diff = abs(a - b)
    max_ab = max(abs(a), abs(b), 1)
    if max_ab >= diff or max_ab > 1:
        return diff <= epsilon              # absolute error
    else:
        return diff < epsilon * max_ab      # relative  error


class PyRootsError(Exception):
    pass


class ConvergenceError(PyRootsError):
    pass
