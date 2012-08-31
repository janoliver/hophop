#!/usr/bin/python2

import numpy as np
import sqlite3 as sql
from itertools import groupby

class SummaryParser(object):
    """
    A parser for the hop summary files. The following must be fulfilled:
      * The columns MUST be 20 bytes (chars) long
      * The first line of the file must include the version of the hopping
    simulation
      * the second line includes the column titles
      * the third line includes the variable types

    """
    
    def __init__(self):
        self.version = 2.0
        self.col_width = 20
        
    def read(self, filename='summary.dat', cols=31):
        f = open(filename, 'r')

        # some helper function
        def chunks(l, n):
            for i in xrange(0, len(l), n):
                yield l[i:i+n].strip('# ')
        
        # check the version
        version = f.readline().strip().split()[1]
        if float(version) != self.version:
            print version
            raise Exception('Wrong version!')

        # column titles
        col_line = f.readline().strip()
        columns = list(chunks(col_line, self.col_width))
        
        # types
        type_line = f.readline().strip()
        types = list(chunks(type_line, self.col_width))
        
        # build the dtypes
        dt = np.dtype([(c, self.get_type(t)) for c, t in zip(columns, types)])
        print dt
        f.close()
        
        # read the data
        data = np.genfromtxt(filename, dtype=dt, delimiter=self.col_width,
                             names=None, autostrip=True,
                             usecols = range(len(dt)))

        print data
        return self
    
    def get_type(self, typestr):
        types = {
            'int' : '<i4',
            'float' : '<f8',
            'long' : '<i8',
            'str' : '|S{}'.format(self.col_width),
            'datetime' : 'datetime64'
            }
        return types[typestr]
