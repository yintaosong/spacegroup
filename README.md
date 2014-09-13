# Space Group Project

This is package for multi_lattice analysis based on the notion of sites, bonds, and local environments.
It supports the conversion from and to CIF file formats.

This is a temporary project. It will eventually been merged into the PyStrucTrans project.


## Installation

Clone or download/unzip the project into your working directory.
Then to use it in your local python scripts, start with  `import bonds` 
or `from bonds import *`. 
A useful short hand may be `from bonds import MultiLattice as ML`.

## Example

### Generate multilattice

Let's generate a multilattice with P 2<sub>1</sub> space group, and two occupied *2a* sites, one for Ti, one for Ni. 

Import the package first

    from bonds import *
    from bonds import MultiLattice as ML

First, we need a base matrix. You can generate it by PyStrucTrans. If it is a monoclinic lattice, here is a helper function in `bonds.utils`

    a = 2.91760
    b = 4.04300
    c = 4.64030
    beta = 97.8320
    B = utils.lp2B(a, b, c, beta)

Then we need to define our two sites, using the constructor `Site(atom, pos, label, occupancy)`

    NiSite = Site("Ni", (0.54421, 0.25000, 0.17586), "Ni", 1)
    TiSite = Site("Ti", (0.08316, 0.75000, 0.28152), "Ti", 1)

Finally, we can generate a `MultiLattice` object

    syms = ('x, y, z', '-x, y+1/2, -z')
    latt = ML(B, [NiSite, TiSite], syms)

Now, you can check its sites:

    print(latt.sites())

It will give the two sites you initially assigned,

    [Ni(Ni) @ (0.544210, 0.250000, 0.175860), Ti(Ti) @ (0.083160, 0.750000, 0.281520)]

To check if the symmetry works, print all the sites in the unit cell

    print(latt.all_sites_by_atom())

The output is

    {'Ti': [Ti(Ti) @ (0.083160, 0.750000, 0.281520), Ti(Ti) @ (0.916840, 0.250000, 0.718480)], 'Ni': [Ni(Ni) @ (0.544210, 0.250000, 0.175860), Ni(Ni) @ (0.455790, 0.750000, 0.824140)]}

Clearly, it generates all four atoms in the cell correctly.

### CIF write and read

Now let's export it to a CIF file. The method `to_vesta_cif` comes in handy:

    # call the phase "New structure", and assign the space group No. 4 (P 21) to it
    content = latt.to_vesta_cif('New structure', "P 21", 4)

    # write it to a file
    file = open("NiTi_4atm_P21.cif", "w")
    file.write(content)
    file.close()

Now load the generated file in VESTA to view the multi lattice.

In fact, you can also load a structure from a cif file.

    lp, syms, sites = cif.read_cif("NiTi_4atm_P21.cif")
    B = utils.lp2B(*lp)
    latt2 = ML(B, sites, syms)
    print(latt2.all_sites_by_atom())

Again, we see

    {'Ni': [Ni(Ni) @ (0.544210, 0.250000, 0.175860), Ni(Ni) @ (0.455790, 0.750000, 0.824140)], 'Ti': [Ti(Ti) @ (0.083160, 0.750000, 0.281520), Ti(Ti) @ (0.916840, 0.250000, 0.718480)]}

### Rebase

A powerful tool come with `MultiLattice` is rebase.
Through rebase, you can change the description of the same lattice,
by change the symmetry operations (space group) or the size of the unit cell.

First, let's look at change symmetry. Say, we want a P1 descrition of the NiTi lattice defined above.

    import numpy as np

    L = np.eye(3)                 # not change the unit cell
    O = NiSite.pos()    # move the origin to Ni site
    syms = ['x, y, z']            # P1 symmetry
    lattP1 = latt.rebase(L, O, ['x, y, z'])
    print(lattP1.all_sites_by_atom())

<<<<<<< HEAD
We get again 4 atoms per cell. Notice that because the two Ti sites are considered different, they have different labels now.
=======
We get again 4 atoms per cell.
>>>>>>> 891c13a... init

    {'Ti': [Ti(Ti1) @ (0.372630, 0.000000, 0.542620), Ti(Ti2) @ (0.538950, 0.500000, 0.105660)], 'Ni': [Ni(Ni1) @ (0.000000, 0.000000, 0.000000), Ni(Ni2) @ (0.911580, 0.500000, 0.648280)]}

Export it to a cif file and check in VESTA

    content = lattP1.to_vesta_cif('New structure', "P 1", 1)
    file = open("NiTi_4atm_P1.cif", "w")
    file.write(content)
    file.close()

If we want to enlarge the cell but still keep P 21 symmetry, we can do

    L = np.array([[2, 0, 0], [0, 1, 0], [0, 0, 2]])  # not change the unit cell
    O = [0, 0, 0]                                    # not move the lattice point
    latt_16atm = latt.rebase(L, O)                   # by default using the original symmetry

Again, export to CIF and check in VESTA

    content = lattP1.to_vesta_cif('New structure', "P 21", 4)
    file = open("NiTi_16atm_P21.cif", "w")
    file.write(content)
    file.close()

### Bonds

<<<<<<< HEAD
The reason that the main package is called `bonds` is that it can calculate all the shortes bonds (nearest neighbors) for each site in a multilattice.
=======
The reason that this package is called `bonds`, is that it can calculator all the shortes bonds (nearest neighbors) for each site in a multilattice.
>>>>>>> 891c13a... init

    # define the coordination numbers for sites
    coords = dict()
    coords[('Ti', 'Ti')] = coords[('Ni', 'Ni')] = 6
    coords[('Ti', 'Ni')] = coords[('Ni', 'Ti')] = 8

    envs = latt.all_local_envs(coords)
    for env in envs:
        print(env.tostr())

It outputs

    "Ni" site at ( 0.544,  0.250,  0.176) :
    =======================================
    Ni-Ti bonds:
    ------------
      Ni-Ti:  2.513032
      Ni-Ti:  2.513175
      Ni-Ti:  2.513175
      Ni-Ti:  2.566964
      Ni-Ti:  2.566964
      Ni-Ti:  2.603043
      Ni-Ti:  2.606943
      Ni-Ti:  3.308537

    Ni-Ni bonds:
    ------------
      Ni-Ni:  2.588816
      Ni-Ni:  2.588816
      Ni-Ni:  2.917600
      Ni-Ni:  2.917600
      Ni-Ni:  3.662494
      Ni-Ni:  3.662494

    "Ti" site at ( 0.083,  0.750,  0.282) :
    =======================================
    Ti-Ti bonds:
    ------------
      Ti-Ti:  2.917600
      Ti-Ti:  2.917600
      Ti-Ti:  2.949806
      Ti-Ti:  2.949806
      Ti-Ti:  3.286712
      Ti-Ti:  3.286712

    Ti-Ni bonds:
    ------------
      Ti-Ni:  2.513032
      Ti-Ni:  2.513175
      Ti-Ni:  2.513175
      Ti-Ni:  2.566964
      Ti-Ni:  2.566964
      Ti-Ni:  2.603043
      Ti-Ni:  2.606943
<<<<<<< HEAD
      Ti-Ni:  3.308537
=======
      Ti-Ni:  3.308537
>>>>>>> 891c13a... init
