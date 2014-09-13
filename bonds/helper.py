from __future__ import print_function, division, absolute_import
from .globals import *
from .utils import *
from .site import Site
from .cif import read_cif
from .multi_lattice import MultiLattice as ML

import logging

def cif_compare_s2s(cif1, cif2, coord):
    print('compare "{:s}" to "{:s}":'.format(cif1, cif2))
    
    # Generate a structure with for cif1
    lp, syms, sites = read_cif(cif1)
    latt1 = ML(lp2B(*lp), sites, syms)
    lp, syms, sites = read_cif(cif2)
    grp = latt1.all_sites_by_atom()
    
    mathchingsites = []
    sitesmap = {}
    # rename sites
    logging.debug("Best match of a site in cif2 to one in cif1:")
    for site in sites:
        equivpos = ML.all_in_nb(site, syms)  
        matchingsite = None
        matchingerr = float('inf')
        for pos in equivpos:
            sort_sites2 = sorted(grp[site.atom()], key=lambda s: la.norm(s.pos() - pos))
            best = sort_sites2[0]
            if la.norm(best.pos() - pos) < matchingerr:
                matchingsite = best.label()
                matchingerr = la.norm(best.pos() - pos)
        sitesmap[matchingsite] = site.label()
        logging.debug("  {:s} => {:s}".format(site.label(), matchingsite))
        mathchingsites.append(Site(site.atom(), site.pos(), matchingsite, site.ocpy()))
    
    # Generate a structure with matching site labels for cif2
    latt2 = ML(lp2B(*lp), mathchingsites, syms)
    
    # generate local envs
    envs1 = latt1.all_local_envs(coord)
    envs2 = latt2.all_local_envs(coord)
    envs1 = sorted(envs1, key=lambda env: env._center.label())
    envs2 = sorted(envs2, key=lambda env: env._center.label())
    
    for i in range(len(envs1)):
        env1 = envs1[i]
        env2 = envs2[i]
        grp2 = env2.bonds_by_coord()
        
        label = env1._center.label()        
        msg = "\n\"{:s}\" site at ".format(label)
        msg += "({:>6.3f}, {:>6.3f}, {:>6.3f}) :".format(*env1._centerpos)
        print(msg + "\n" + "=" * len(msg))        
        for atom, bonds in env1.bonds_by_coord().items():
            bonds2 = grp2[atom]
            
                        
            msg = "{:s}-{:s} bonds:".format(env1._center.atom(), atom)
            print(msg)
            print("-"*len(msg))
            for bond in bonds:
                output = "  {:s}-{:s}: ".format(label, bond.label())
                bondlen = la.norm(bond.pos())
                output += "{:>9.6f}".format(bondlen)
                
                # find the matching bond
                idx = -1
                vec = bond.pos()
                diff = float('inf')
                for j, cand in enumerate(bonds2):
                    if cand.label() == bond.label() and la.norm(cand.pos() - vec) < diff:
                        diff = la.norm(cand.pos() - vec)
                        idx = j
                if idx > -1:
                    bond2 = bonds2.pop(idx)
                    bondlen2 = la.norm(bond2.pos())
                    percent = 100 * (bondlen - bondlen2)/bondlen2
                    output += ", {:+6.2f} %".format(percent)
                    output += " ( compare to {:s}-{:s}, {:>9.6f})".format(sitesmap[env2._center.label()], sitesmap[bond2.label()], bondlen2)
                    if abs(percent) >= 10:
                        output += " WARNING: changed more than 10 % !!"
                
                print(output)
    
            
    