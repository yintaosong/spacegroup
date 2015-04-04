'''
Created on Sep 11, 2014

@author: Yintao Song
'''
from __future__ import print_function, division, absolute_import
from .globals import *

class Site(object):
    '''
    A site is a distinct atom position up to symmetry in a cell.
    To define a site, we need to know site position as linear combinations
    of base vectors, atom type, equivalent positions and give it a label.
    
    This classe is an immutable.
    '''
    
    def __init__(self, atom, pos, label=None, ocpy=1.0):
        '''
        Constructor
        '''
        assert isinstance(atom, str)
        self._atom = atom
        
        assert isinstance(pos, list) or isinstance(pos, tuple) or isinstance(pos, np.ndarray) 
        assert len(pos) == 3
        self._pos = np.array(pos)
                
        if label is not None:
            assert isinstance(label, str)
            self._label = label
        else:
            self._label = atom
        
        assert ocpy >= 0.0 and ocpy <= 1.0
        self._ocpy = ocpy
    
    def atom(self):
        return self._atom
    def pos(self):
        return self._pos.copy()
    def label(self):
        return self._label
    def ocpy(self):
        return self._ocpy
    
    def move_to(self, c):
        """
        move the position to c
        :return: a new Site object
        """
        return Site(self.atom(), c, self.label(), self.ocpy())
    
    def copy(self):
        return self.move_to(self.pos())  
    
    def __repr__(self):
        return '{:s}({:s}) @ ({:f}, {:f}, {:f})'.format(self._atom, self._label, *self._pos)
    
    def __str__(self):
        return '{:6s} {:>9.6f} {:>9.6f} {:>9.6f} {:3s}'.format(self._label, self._pos[0], self._pos[1], self._pos[2], self._atom)  
    
    def cif_record(self):
        output = '{:6s} {:>5.2f}'.format(self._label, self._ocpy)
        output += ' {:>9.6f} {:>9.6f} {:>9.6f} '.format(*self._pos)
        output += 'Biso 1.000000 '
        output += '{:3s}'.format(self._atom)
        return output

    def cri_record(self):
        output = '{:2s}{:3s}'.format(self._atom, self._label)
        output += ' {:>9.6f} {:>9.6f} {:9.6f}'.format(*self._pos)
        output += '  {:5.5f}'.format(self._ocpy)
        output += '  {:3s}'.format(self._atom)
        return output
    