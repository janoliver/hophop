#!/usr/bin/python2

import numpy as np
import sqlite3 as sql
from itertools import groupby

class SummaryParser(object):
    """
    A parser for the hop summary files. The following must be fulfilled:
      * The columns MUST be 18 bytes (chars) long
      * The first line of the file must include the column names
      * there should not be empty fields in the columns that are read in

      The number of columns to be read in can be specified using the parameter
      cols of the read function. A typical usage example would be something
      like this:

        >>> from summaryparser import SummaryParser
        >>> s = SummaryParser()
        >>> s.read('summaryfile.dat', cols=30)
        >>> print s.get_columns()
        >>> print s.q('SELECT col1, col2 FROM data WHERE col3="val"').get_groups()
        >>> s.close()
    """
    
    def __init__(self):
        self.con = sql.connect(":memory:")
        self.raw = None
        self.cols = None
        self.data = None
        return self.
        
    def read(self, filename='summary.dat', cols=31):
        self.raw = np.genfromtxt('compare.dat', dtype=None, delimiter=18,
                                 autostrip=True, names=True,
                                 usecols=range(cols),
                                 comments=None)
        self.cols = self.raw.dtype.names

        self._insert()
        return self

    def q(self, query):
        self.data = np.array(self.con.execute(query).fetchall())
        return self

    def get_data(self):
        return self.data
        
    def get_groups(self, col=0):
        groups = {}
        for k, g in groupby(self.data.tolist(), lambda x: x[col]):
            groups[k] = np.delete(np.array(list(g)), col, 1).T
        return groups

    def get_columns(self):
        return self.cols

    def get_raw(self):
        return self.raw
    
    def _insert(self):
        self.con.execute("create table data({0})".format(",".join(self.cols)))
        self.con.executemany("insert into data values (%s)" %
                             ",".join(['?']*len(self.cols)), self.raw.tolist())

    def close(self):
        self.con.close()
    
