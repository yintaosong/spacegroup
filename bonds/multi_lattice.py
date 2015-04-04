'''
Created on Sep 11, 2014

@author: Yintao
'''
from __future__ import print_function, division, absolute_import
from .globals import *
from itertools import groupby

from .site import Site
from .environment import Environment
from .utils import rm_dup, mv2cell

@jit()
def _all_in_nb_jit(site, syms, rng):
    """
    all copies of sites in neighboring cells
    :return: positions
    :rtype: :py:class:`list` of :py:class:`numpy.ndarray`
    """
    cell = np.array(MultiLattice.all_in_cell(site, syms))
    ngbr = []
    for i in range(-rng, rng + 1):
        for j in range(-rng, rng + 1):
            for k in range(-rng, rng + 1):
                c = np.array([i, j, k])
                for s in cell:                    
                    ngbr.append(s + c)
    return ngbr 

@jit()
def _rm_dup_by_sym(sites, syms):
    # remove duplicates
    uniquesites = []
    while len(sites) > 0:
        # generate equivalent positions
        equiv = MultiLattice.all_in_cell(sites[0], syms)
        equiv = sorted(equiv, key=la.norm)
        # if there are remained
        rm_idx = [0]
        for i in range(len(sites)):
            remain = sites[i]
            diff = la.norm(equiv - remain.pos(), axis=1) 
            if min(diff) < 1e-7:
                rm_idx.append(i)
        uniquesites.append(sites[0].move_to(equiv[0]))
        sites = [sites[i] for i in range(len(sites)) if not i in rm_idx]                     
    return uniquesites
            
_TRIVAL_SYMS=('x,y,z',)

class MultiLattice(object):
    '''
    A multilattice is a lattice with more than one sites and general equivalent positions.
    '''
    
    @staticmethod
    def all_in_cell(site, syms=_TRIVAL_SYMS):
        """
        all distinct copies of the site in the cell.
        :return: positions of all the copies
        :rtype: :py:class:`list` of :py:class:`numpy.ndarray`
        """
        x, y, z = site.pos()
        positions = []
        for sym in syms:
            positions.append(mv2cell(np.array(eval(sym))))
        return rm_dup(positions)
    
    @staticmethod
    def all_in_nb(site, syms=_TRIVAL_SYMS, rng=1):
        """
        all copies of sites in neighboring cells
        :return: positions
        :rtype: :py:class:`list` of :py:class:`numpy.ndarray`
        """
        return _all_in_nb_jit(site, syms, rng)
    
    @staticmethod
    def one_to_one_bonds(s1, s2, B, syms=_TRIVAL_SYMS):
        """
        find all bonds from s1 to any copy of s2 in neighboring cells
        :return: a list of (site1, site2, bond_length)
        """
        # position vectors
        p1 = s1.pos()
        p2_cpys = MultiLattice.all_in_nb(s2, syms)
        bnds = []
        
        for p2 in p2_cpys:
            bond = p2 - p1
            new_s2 = s2.move_to(p2)
            # site1, site2, bond_length
            record = (s1, new_s2, la.norm(B.dot(bond)))
            bnds.append(record)
        bnds = sorted(bnds, key=lambda x: x[2])
        
        return bnds[1:] if abs(bnds[0][2]) < 1e-7 else bnds
    
    @staticmethod
    def one_to_many_bonds(s1, s2_list, B, syms=_TRIVAL_SYMS):
        """find all bonds from s1 to any copy of each site in s2_list"""
        bnds = []
        for s2 in s2_list:
            bnds.extend(MultiLattice.one_to_one_bonds(s1, s2, B, syms))
        return sorted(bnds, key=lambda x: x[2])
    
    def __init__(self, B, sites, syms=_TRIVAL_SYMS):
        '''
        Constructor
        '''
        assert np.array(B).shape == (3, 3)
        self._B = B
        
        assert isinstance(sites, list) or isinstance(sites, tuple)
        for site in sites:
            assert isinstance(site, Site)
        self._sites = sites
        self._allSites = None
        
        assert isinstance(syms, list) or isinstance(syms, tuple)
        for sym in syms:
            assert isinstance(sym, str)
        self._syms = tuple(syms)
        
    def base(self):
        return self._B.copy()
    def sites(self):
        return self._sites
    def syms(self):
        return self._syms
    
    def atoms(self):
        """ all different atoms """
        atoms = []
        for k, _ in groupby(self.all_sites(), key=lambda s: s.atom()):
            atoms.append(k)
        return tuple(set(atoms))
    
    def site_labels(self, gry_by_atom=False):
        if gry_by_atom:
            pairs = [(s.atom(), s.label()) for s in self.all_sites()]
            pairs = list(set(pairs))
            pairs = sorted(pairs, key=lambda p: p[0])            
            labels = {}
            for k, v in groupby(pairs, key=lambda s: s[0]):
                labels[k] = tuple(p[1] for p in v)
            return labels
        else:       
            labels = []
            for k, _ in groupby(self.all_sites(), key=lambda s: s.label()):
                labels.append(k)
            return tuple(set(labels))
        
    def rebase(self, L, neworigin, newsyms=None):
        """
        rewrite the multilattice in a new base and at a new origin, with a new set of symmetries
        """
        newsyms = self._syms if newsyms is None else newsyms
        neworigin = np.array(neworigin)
        rng = max(int(np.max(L)), 1)        
        Linv = la.inv(L)
        # get all atoms in the new cell
        newsites = []
        for site in self.sites():
            candidates = np.array(MultiLattice.all_in_nb(site, self._syms, rng))
            candidates = Linv.dot(candidates.T).T - neworigin
            for cand in candidates:
                if min(cand) >= 0.0 and max(cand) < 1.0 - 1e-10:
                    newsites.append(site.move_to(cand))
        
        newsites = _rm_dup_by_sym(newsites, newsyms)
        
        # rename labels
        atoms = {}
        finalsites = []
        for site in newsites:
            if site.atom() not in atoms:
                atoms[site.atom()] = 1
            else:
                atoms[site.atom()] += 1
            finalsites.append(Site(site.atom(), site.pos(), site.atom()+str(atoms[site.atom()]), site.ocpy()))
        return MultiLattice(self.base().dot(L), finalsites, newsyms)                    
                      
        
    def all_sites(self):
        """
        :return: all atom positions after applying the symmetry
        :rtype: :py:class:`list` of :py:class:`Site` objects
        """
        if self._allSites is None:
            all_sites = []
            for site in self._sites:
                for pos in MultiLattice.all_in_cell(site, self._syms):
                    all_sites.append(Site(site.atom(), pos, site.label(), site.ocpy()))
            self._allSites = all_sites        
        return self._allSites
    
    def all_sites_by_atom(self):
        """
        :return: all sites grouped by atom type
        :rtype: :py:class:`dict` of <atom, sites> pairs
        """
        grps = dict()
        for k, v in groupby(self.all_sites(), key=lambda s: s.atom()):
            grps[k] = list(v)
        return grps
    
    def sites_by_atom(self):
        """
        :return: all sites grouped by atom type
        :rtype: :py:class:`dict` of <atom, sites> pairs
        """
        grps = dict()
        for k, v in groupby(self.sites(), key=lambda s: s.atom()):
            grps[k] = list(v)
        return grps
    
    def local_env(self, site, coord={}):
        """
        get the local environments for a particular site, 
        using assigned coordination numbers
        """
        bonds = MultiLattice.one_to_many_bonds(site, self.sites(), self._B, self._syms)
        ngbrs = tuple(bond[1] for bond in bonds)
        return Environment(site, ngbrs, self.base(), coord)
    
    def all_local_envs(self, coord={}):
        return tuple(self.local_env(site, coord) for site in self.sites())
    
    def sites_tostr(self, grp_by_atom=False):
        grps = self.sites_by_atom()
        output = ""
        for k, v in grps.items():
            for site in v:
                output += str(site) + "\n"
        return output
    
    def to_xmas_cri(self, phase_name='New structure', Int_Tbl=1):
        
        output = '{:s}\n'.format(phase_name)
        
        M = self.base().T.dot(self.base())
        a, b, c = np.sqrt([M[0, 0], M[1, 1], M[2, 2]])
        alpha, beta, gamma = np.degrees(np.arccos([M[1, 2]/(b*c), M[0, 2]/(a*c), M[0, 1]/(a*b)]))

        output += '{:d}\n'.format(Int_Tbl)

        output += '{:.5f}   {:.5f}   {:.5f}   {:.5f}   {:.5f}   {:.5f}\n'.format(a, b, c, alpha, beta, gamma)

        output += '{:2d}\n'.format(len(self.all_sites()))

        for atom, sites in self.sites_by_atom().items():
            for site in sites:
                output += site.cri_record() + "\n"
        return output

    def to_vesta_cif(self, phase_name='New structure', space_group="P 1", Int_Tbl=1):

        output = \
        "#======================================================================\n\n" + \
        "# CRYSTAL DATA \n\n" + \
        "#----------------------------------------------------------------------\n\n" + \
        "data_VESTA_phase_1\n\n" + \
        "_pd_phase_name                         '{:s}'\n".format(phase_name)

        M = self.base().T.dot(self.base())
        a, b, c = np.sqrt([M[0, 0], M[1, 1], M[2, 2]])
        alpha, beta, gamma = np.degrees(np.arccos([M[1, 2]/(b*c), M[0, 2]/(a*c), M[0, 1]/(a*b)]))

        output += "_cell_length_a                         {:.5g}\n".format(a)
        output += "_cell_length_b                         {:.5g}\n".format(b)
        output += "_cell_length_c                         {:.5g}\n".format(c)
        output += "_cell_angle_alpha                      {:.5g}\n".format(alpha)
        output += "_cell_angle_beta                       {:.5g}\n".format(beta)
        output += "_cell_angle_gamma                      {:.5g}\n".format(gamma)

        output += "_symmetry_space_group_name_H-M         '{:s}'\n".format(space_group)
        output += "_symmetry_Int_Tables_number            {:d}\n\n".format(Int_Tbl)
#
#
#         _cell_length_a                         5.83514
#         _cell_length_b                         4.04300
#         _cell_length_c                         9.28064
#         _cell_angle_alpha                      90
#         _cell_angle_beta                       97.83205
#         _cell_angle_gamma                      90
#         output = "loop_\n_symmetry_equiv_pos_as_xyz\n"



        output += "loop_\n_symmetry_equiv_pos_as_xyz\n"
        for sym in self._syms:
            output += "   '" + sym + "'\n"
        output += "\nloop_ \n" + \
                   "   _atom_site_label \n" + \
                   "   _atom_site_occupancy \n" + \
                   "   _atom_site_fract_x \n" + \
                   "   _atom_site_fract_y \n" + \
                   "   _atom_site_fract_z \n" + \
                   "   _atom_site_adp_type \n" + \
                   "   _atom_site_B_iso_or_equiv \n" + \
                   "   _atom_site_type_symbol \n"
        for atom, sites in self.sites_by_atom().items():
            for site in sites:
                output += "   " + site.cif_record() + "\n"
        return output