#!/usr/bin/python2

import numpy as np
import sqlite3 as sql
from itertools import groupby
import re

class SummaryParser(object):
    """
    A parser for the hop summary files. The following must be fulfilled:
      * The columns MUST be 20 bytes (chars) long
      * The first line of the file must include the version of the hopping
        simulation
      * the second line includes the column titles
      * the third line includes the variable types

    """

    dtype = None
    data = None
    working_data = None
    expression_matcher = None
    
    def __init__(self, filename=None):
        self.version = 2.0
        self.col_width = 20

        if filename:
            self.read(filename)
        
    def read(self, filename='summary.dat', cols=31):
        f = open(filename, 'r')

        # some helper function
        def chunks(l, n):
            for i in xrange(0, len(l), n):
                yield l[i:i+n].strip('# ')
        
        # check the version
        version = f.readline().strip().split()[1]
        if float(version) != self.version:
            raise Exception("Parser and summary file versions don't match!")

        # column titles
        col_line = f.readline().strip()
        columns = list(chunks(col_line, self.col_width))
        
        # types
        type_line = f.readline().strip()
        types = list(chunks(type_line, self.col_width))
        
        # build the dtypes
        self.dtype = np.dtype(
            [(c, self._get_type(t)) for c, t in zip(columns, types)]
            )
        f.close()
        
        # read the data
        self.data = np.genfromtxt(filename, dtype=self.dtype, 
            delimiter=self.col_width, names=None, autostrip=True,
            usecols = range(len(self.dtype)))

        # copy the data
        self.reset()

        return self

    def print_column_names(self):
        print "Available Column keys:"
        print "   ", "\n    ".join(self.dtype.names)
    
    def reset(self):
        self.working_data = self.data
        return self

    def sort(self, column, desc=False):
        """
        Sort by column and reverse it if desc is True
        """

        if not column in self.dtype.names:
            raise Exception("Unknown column")

        self.working_data = np.sort(self.working_data, order=column)

        if desc:
            self.working_data = self.working_data[::-1]

        return self

    def filter(self, left, right, operator="=="):
        """
        This function filters the data. You must provide a left and right
        operand, and may provide an operator (>, >=, <, <=, ==, !=). The 
        default is ==. Operands can be any number that is convertable into a
        float or one of the column names as a string.

        Example: a.filter("dos_exponent", 2.0, "<")
        """
        operators = {
            ">": np.greater,
            ">=": np.greater_equal,
            "<": np.less,
            "<=": np.less_equal,
            "==": np.equal,
            "!=": np.not_equal,
            "str==": np.core.defchararray.equal,
            "str!=": np.core.defchararray.not_equal
        }

        # test the operator validity
        if not operator in operators:
            raise Exception("Unknown operator.")

        #test the left side
        if left in self.dtype.names:
            left = self.working_data[left]
        else:
            try:
                float(left)
            except:
                pass # seems to be astring...

        # test the right side
        if right in self.dtype.names:
            right = self.working_data[right]
        else:
            try:
                float(right)
            except:
                pass # seems to be a string...
                
        self.working_data = self.working_data[operators[operator](left, right)]

        return self

    def get(self, cols=None):
        """
        Return the data as numpy array. If cols is given, return only those 
        columns. cols must be a list!
        """

        if not cols:
            ret = self.working_data
        else:
            for col in cols:
                if not col in self.dtype.names:
                    raise Exception("Unknown column name!")

            ret = self.working_data[cols]

        self.reset()
        return ret

    def _get_type(self, typestr):
        types = {
            'int' : '<i4',
            'float' : '<f8',
            'long' : '<i8',
            'str' : '|S{}'.format(self.col_width),
            'datetime' : 'datetime64[s]'
            }
        return types[typestr]
