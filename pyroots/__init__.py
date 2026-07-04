#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# file pyroots/__init__.py
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

""" A python-only implementation of single-variable root finding methods.  """

from __future__ import division
from __future__ import print_function
from __future__ import absolute_import

# package information
__package_name__ = "pyroots"
__version__ = "0.5.0"
__license__ = "BSD"
__description__ = __doc__.split(".")[0]
__url__ = "http://github.com/pmav99/%s" % __package_name__
__download_url__ = "http://github.com/pmav99/%s/archive/master.zip" % __package_name__
__author__ = "Panagiotis Mavrogiorgos"
__author_email__ = "gmail pmav99"

# set up logging
import logging
logging.getLogger('pyroots').addHandler(logging.NullHandler())

# Package imports
from .utils import ConvergenceError
from .bisect import Bisect
from .ridder import Ridder
from .brent import Brentq, Brenth

__all__ = ["Bisect", "Ridder", "Brenth", "Brentq", "ConvergenceError"]
