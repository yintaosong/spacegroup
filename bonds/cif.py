'''
Created on Sep 11, 2014
 
@author: Yintao Song
'''
from __future__ import print_function, division, absolute_import
from .globals import *
from .site import Site
 
def parse_site(line):
    site, occpy, x, y, z, _, _, atom = line.split()
    return Site(atom, [float(x), float(y), float(z)], site, float(occpy))
 
def read_cif(filename):
    f = open(filename, 'r')
     
    syms = []
    sites = []
     
    in_pos = False
    in_sym = False
    for line in f:
        # lattice parameters
        if "_cell_length_a" in line: a = float(line.split()[1])
        if "_cell_length_b" in line: b = float(line.split()[1])
        if "_cell_length_c" in line: c = float(line.split()[1])
        if "_cell_angle_beta" in line: beta = float(line.split()[1])
         
        # symmetry
        if in_sym:
            if "'" in line:
                syms.append("".join(line.split())[1:-1])
            else:
                in_sym = False
             
        if "_symmetry_equiv_pos_as_xyz" in line: in_sym = True
         
        # sites
        if in_pos: sites.append(parse_site(line))            
        if "_atom_site_type_symbol" in line: in_pos = True
             
    return (a, b, c, beta), syms, sites
