:py:class:`pychemiq.Molecules`
==================================

Classes
----------

.. py:class:: Molecules(geometry=None, basis=None, multiplicity=None, charge=0, active=None, nfrozen=None)

    Initialize the electronic structure parameters of molecules, including charge, basis-set, atomic coordinates, spin multiplicity, etc. 

   :param str geometry: Enter the element symbol and coordinates of atoms in the molecule. It can be a string type or a string list. For example: geometry="H 0 0 0, H 0 0 0.74" or geometry=["H 0 0 0", "H 0 0 0.74"].
   :param str basis: Set the basis set level required for computation: [MINI / sto-3G / sto-6G / 3-21G / 6-31G / ...]. Starting from version V1.1.0, pyChemiQ supports over 600 commonly used basis sets as well as user-defined basis sets. It supports most basis sets available on the Basis Set Exchange website (excluding a few excessively large ones), including but not limited to STO-nG basis sets, Pople series basis sets, Ahlrichs' def-series basis sets, Dunning series basis sets, and pseudopotential basis sets. In previous versions, supported basis sets included: MINI, STO-3G, STO-6G, 3-21G, and 6-31G. Other basis sets can be directly entered by name. The specific string input rules for basis set names please refer to :doc:`./06API/Configs` . 
   :param int multiplicity: Set the spin multiplicity (M=2S+1), with a default value of 1.
   :param int charge: Set the number of charges in the system. When the number of charges is positive, there is no positive sign, and when it is negative, a negative sign is written. The default value is 0.
   :param int,int active: Sets the active space, specified as two comma-separated integers. Default: not set. For example, \textit{active = 4, 4} means 4 active spatial orbitals with 4 electrons in the active space.
   :param int nfrozen: Sets the number of frozen orbitals. Default: 0. Note: `active` and `nfrozen` cannot be used simultaneously.

   :return: None.


   **Attributes**

   .. py:attribute:: n_atoms

      Obtaining the number of atoms in a molecular system.

   .. py:attribute:: n_electrons

      Obtaining the number of electrons in the molecular system.

   .. py:attribute:: n_orbitals

      Obtaining the total molecular orbital number of the molecular system.

   .. py:attribute:: n_qubits

      Obtaining the total number of quantum bits required for calculation (i.e. number of spin orbitals, 2 * number of molecular orbitals)

   .. py:attribute:: hf_energy

      Obtaining the Hartree-Fock energy(unit: Hartree).

   .. py:attribute:: nuclear_repulsion

      Obtaining the inter nuclear repulsion of the molecular system (unit: Hartree)

   .. py:attribute:: canonical_orbitals

      Obtaining the canonical orbital coefficients (i.e. molecular orbital coefficients) of the molecular system.

   .. py:attribute:: orbital_energies

      Obtaining the energy of each molecular orbital in the system.
      
   .. py:attribute:: overlap_integrals

      Obtaining overlapping integrals of a molecular system.

   .. py:attribute:: one_body_integrals

      Obtaining single electron integrals for a molecular system.

   .. py:attribute:: two_body_integrals

      Obtaining the double electron integral of a molecular system.



   **Methods**

   .. py:method:: get_molecular_hamiltonian()

      Obtaining the Hamiltonian of the initialized molecular system.


---------

**Interface example:**

.. code:: 

    from pychemiq import Molecules
    multiplicity = 1
    charge = 0
    basis =  "sto-3g"
    geom = "H 0 0 0,H 0 0 0.74"
    mol = Molecules(
          geometry = geom,
          basis    = basis,
          multiplicity = multiplicity,
          charge = charge)

Call the following interface to obtain information about the molecular system:

.. code:: 

    print("The number of atoms is", mol.n_atoms)
    print("The number of electrons is", mol.n_electrons)
    print("The number of orbitals is", mol.n_orbitals)
    print("The number of qubits is", mol.n_qubits)
    print("The Hartree-Fock energy is", mol.hf_energy)
    print("The nuclear repulsion is", mol.nuclear_repulsion)


.. parsed-literal::

    The number of atoms is 2
    The number of electrons is 2
    The number of orbitals is 2
    The number of qubits is 4
    The Hartree-Fock energy is -1.1167593072992057
    The nuclear repulsion is 0.7151043390810812


.. code:: 

    print("The canonical orbitals are\n", mol.canonical_orbitals)
    print("The orbital energies are", mol.orbital_energies)
    print("The overlap integrals are\n", mol.overlap_integrals)


.. parsed-literal::

    The canonical orbitals are
     [[-0.54884228  1.21245192]
     [-0.54884228 -1.21245192]]
     
    The orbital energies are [-0.57855386  0.67114349]

    The overlap integrals are
     [[1.         0.65987312]
     [0.65987312 1.        ]]


.. code:: 

    print("The one body integrals are\n", mol.one_body_integrals)
    print("The two body integrals are\n", mol.two_body_integrals)


.. parsed-literal::

    The one body integrals are
     [[-1.25330979e+00  0.00000000e+00]
     [ 4.16333634e-17 -4.75068849e-01]]

    The two body integrals are
     [[[[ 6.74755927e-01 -1.11022302e-16]
       [-8.32667268e-17  6.63711401e-01]]
    
      [[-3.46944695e-17  1.81210462e-01]
       [ 1.81210462e-01  0.00000000e+00]]]
    
    
     [[[-4.85722573e-17  1.81210462e-01]
       [ 1.81210462e-01 -2.22044605e-16]]
    
      [[ 6.63711401e-01 -2.22044605e-16]
       [-1.66533454e-16  6.97651504e-01]]]]

.. code:: 

    print("The molecular hamiltonian is", mol.get_molecular_hamiltonian())


.. parsed-literal::

    The molecular hamiltonian is {
    : 0.715104
    0+ 0 : -1.253310
    1+ 0+ 1 0 : -0.674756
    1+ 0+ 3 2 : -0.181210
    1+ 1 : -1.253310
    2+ 0+ 2 0 : -0.482501
    2+ 1+ 2 1 : -0.663711
    2+ 1+ 3 0 : 0.181210
    2+ 2 : -0.475069
    3+ 0+ 2 1 : 0.181210
    3+ 0+ 3 0 : -0.663711
    3+ 1+ 3 1 : -0.482501
    3+ 2+ 1 0 : -0.181210
    3+ 2+ 3 2 : -0.697652
    3+ 3 : -0.475069
    }
    