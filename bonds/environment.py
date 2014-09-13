'''
Created on Sep 11, 2014

@author: Yintao Song
'''
from __future__ import print_function, division, absolute_import
from .globals import *
from .site import Site
from itertools import groupby


class Environment(object):
    '''
    An environment of an atom site is defined by
    the bonds between itself and its nearest and second-nearest neighbors.
    '''

    def __init__(self, center, ngbrs, B, coord={}, ngbr_type='global'):
        '''
        Constructor.
        '''
        assert isinstance(center, Site)
        self._center = center
        self._centerpos = center.pos()
        
        
        try:
            B = np.array(B)
        except Exception:
            raise TypeError("the input base cannot be parsed to a numpy.ndarray")
        assert B.shape == (3,3)        
        self._B = B
        
        assert isinstance(ngbrs, list) or isinstance(ngbrs, tuple) 
        for ngbr in ngbrs:
            assert isinstance(ngbr, Site)
        
        assert isinstance(ngbr_type, str)
        if ngbr_type == 'global':
            self._ngbrs = ngbrs
            self._bonds = self.glb2loc(ngbrs)
        elif ngbr_type == 'local':
            self._bonds = ngbrs
            self._ngbrs = self.loc2glb(ngbrs)
        else:
            raise ValueError("unknown neighborhood type: {:s}".format(ngbr_type))
            
        assert isinstance(coord, dict)
        self._coord = coord
    
    def glb2loc(self, sites):
        local_sites = []
        for site in sites:
            newpos = np.dot(self._B, site.pos() - self._centerpos)
            newsite = site.move_to(newpos)
            local_sites.append(newsite)
        return local_sites
    
    def loc2glb(self, sites):
        global_sites = []
        for site in sites:
            newpos = la.inv(self._B).dot(site.pos()) + self._centerpos
            newsite = site.move_to(newpos)
            global_sites.append(newsite)
        return global_sites
    
    def base(self):
        return self._B.copy()
        
    def center(self):
        return self._center
    
    def ngbrs(self, limit=-1):
        """
        Get neighbor sites up to the limit,
        positions in the global coordinates.
        When limit is -1, get all. 
        """
        limit = min(limit, len(self._ngbrs))
        return self._ngbrs[:limit]
        
    def bonds(self, limit=-1):
        """
        Get neighbor sites up to the limit,
        positions in the global coordinates.
        When limit is -1, get all. 
        """
        limit = min(limit, len(self._bonds))
        return self._bonds[:limit]
        
    def bonds_by_coord(self):
        """
        output bonds according to the assigned coordination numbers
        """
        bonds = dict()
        # group by atom type at end point
        for atom, itr in groupby(sorted(self.bonds(), key=lambda s: s.atom()), key=lambda s: s.atom()):
            limit = self._coord[(self._center.atom(), atom)]  \
                if (self._center.atom(), atom) in self._coord \
                else -1
            bonds[atom] = sorted(itr, key=lambda s: la.norm(s.pos()))
            limit = min(limit, len(bonds[atom]))
            bonds[atom] = bonds[atom][:limit]            
        return bonds
    
    def tostr(self, show_vec=False):
        label = self._center.label()
        
        msg = "\"{:s}\" site at ".format(label)
        msg += "({:>6.3f}, {:>6.3f}, {:>6.3f}) :".format(*self._centerpos)
        output = msg + "\n" + "=" * len(msg)
        
        for atom, bonds in self.bonds_by_coord().items():
            msg = "{:s}-{:s} bonds:".format(self._center.atom(), atom)
            output += "\n" + msg + "\n"
            output += "-"*len(msg) + "\n"
            for bond in bonds:
                output += "  {:s}-{:s}: ".format(label, bond.label())
                output += "{:>9.6f}".format(la.norm(bond.pos()))
                if show_vec:
                    output += ",  ({:>6.3f}, {:>6.3f}, {:>6.3f})".format(*bond.pos())
                output += "\n"
        return output
    
    def __str__(self):
        return self.tostr(True)
        
        
        
        
        
    
    