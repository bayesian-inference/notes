#!/usr/bin/env python

from cliqueTree import *
from numpy import *
from factor import *

def createCliqueTree():
    f_01 = Factor([0,1], [2,2], array(range(4)))
    f_12 = Factor([1,2], [2,2], array(range(4)) + 1)
    f_23 = Factor([2,3], [2,2], array(range(4)) * 2)
    f_24 = Factor([2,4], [2,2], array(range(4)) * 10)
    factors = [[f_01], [f_12], [f_23], [f_24]]
    cliques = [[0,1], [1,2], [2,3], [2,4]]
    edges = [[0,1], [1,2], [1,3]]
    ct = CliqueTree(cliques, edges, factors)
    return ct

ct = createCliqueTree()
ct.beliefPropagation()

for i in range(len(ct.beliefs)):
    print "CLIQUE " + str(i)
    print ct.beliefs[i]
    print
