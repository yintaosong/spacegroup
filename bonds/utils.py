'''
Created on Sep 11, 2014

@author: Yintao Song
'''
from __future__ import print_function, division, absolute_import
from .globals import *
import numpy as np
from math import cos, sin, radians

def lp2B(a, b, c, beta):
    """lattice parameters to base, for monoclinic"""
    B = np.array([[a,0,0],[0,b,0],[c*np.cos(radians(beta)), 0, c*np.sin(radians(beta))]]).T
    return B

@jit('boolean(f8[:], f8[:, :], i4, f8)')
def in_list(ary, arys, start, tol):
    assert start < len(arys)
    for i in range(start, len(arys)):
        b = arys[i]
        if max(np.abs(ary - b)) < tol:
            return True
    return False

def rm_dup(arys, tol=1e-7):
    """remove duplicated arrays"""
    dup = []
    for i in range(len(arys) - 1):
        a = arys[i]
        if in_list(a, arys, i+1, tol):
            dup.append(i)
    new_arys = [arys[i] for i in range(len(arys)) if not i in dup]
    return new_arys

@jit("f8[:](f8[:])")
def mv2cell(p):
    for i in range(len(p)):
        if p[i] < 0: p[i] += 1
        if p[i] >= 1: p[i] -= 1
    return p