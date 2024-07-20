"""
The MIT License (MIT)

Copyright (c) 2014 Daniel Smith

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.


This code is obtained from: https://github.com/dgasmith/opt_einsum/blob/master/opt_einsum/backends/object_arrays.py

This file is part of the ``opt_einsum`` package and is distributed under the MIT License.
We add it here to avoid the dependency on the package.

Functions for performing contractions with array elements which are objects.
"""

import functools
import operator

import numpy as np


def object_einsum(eq, *arrays):
    """A ``einsum`` implementation for ``numpy`` arrays with object dtype.
    The loop is performed in python, meaning the objects themselves need
    only to implement ``__mul__`` and ``__add__`` for the contraction to be
    computed. This may be useful when, for example, computing expressions of
    tensors with symbolic elements, but note it will be very slow when compared
    to ``numpy.einsum`` and numeric data types!

    Parameters
    ----------
    eq : str
        The contraction string, should specify output.
    arrays : sequence of arrays
        These can be any indexable arrays as long as addition and
        multiplication is defined on the elements.

    Returns
    -------
    out : numpy.ndarray
        The output tensor, with ``dtype=object``.
    """

    # when called by ``opt_einsum`` we will always be given a full eq
    lhs, output = eq.split("->")
    inputs = lhs.split(",")

    sizes = {}
    for term, array in zip(inputs, arrays):
        for k, d in zip(term, array.shape):
            sizes[k] = d

    out_size = tuple(sizes[k] for k in output)
    out = np.empty(out_size, dtype=object)

    inner = tuple(k for k in sizes if k not in output)
    inner_size = tuple(sizes[k] for k in inner)

    for coo_o in np.ndindex(*out_size):

        coord = dict(zip(output, coo_o))

        def gen_inner_sum():
            for coo_i in np.ndindex(*inner_size):
                coord.update(dict(zip(inner, coo_i)))
                locs = (tuple(coord[k] for k in term) for term in inputs)
                elements = (array[loc] for array, loc in zip(arrays, locs))
                yield functools.reduce(operator.mul, elements)

        out[coo_o] = functools.reduce(operator.add, gen_inner_sum())

    return out
