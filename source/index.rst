.. pychemiq documentation master file, created by
   sphinx-quickstart on Mon Nov 14 14:20:39 2022.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

pyChemiQ
====================================
pyChemiQ is a Python software library developed by Origin Quantum. Based on the Python interface encapsulated by ChemiQ, it enables the simulation and computation of fermionic systems on real quantum computers or virtual machines. It supports user-defined mappings, ansatzes, and optimizers, facilitating secondary development by users. This package provides a simple, lightweight, and efficient platform for quantum chemical calculations and method development. pyChemiQ can be used to simulate molecules on quantum computers using mean-field and post-mean-field methods to solve quantum chemistry problems. By simplifying the conversion between molecular structure input and quantum circuits, pyChemiQ minimizes the domain expertise required to enter the field, helping  learners to address and research electronic structure problems on quantum computers. 

Currently, pyChemiQ supports inputting molecular structures to obtain the second-quantized Fermionic Hamiltonian. In terms of mapping, pyChemiQ supports the Jordan-Wigner (JW) transformation, Bravyi-Kitaev (BK) transformation, Parity transformation, and Multilayer Segmented Parity (MSP) transformation methods to map the second-quantized Fermionic Hamiltonian operator into a Pauli Hamiltonian operator on a quantum computer. Regarding ansatzes, pyChemiQ also supports different ansatzes for constructing quantum circuits, such as Unitary Coupled Cluster (UCC) ansatze, Hardware-Efficient ansatze, and symmetry-preserved ansatze. For optimizers, pyChemiQ provides the following classical optimizers for variational parameter optimization: NELDER-MEAD, POWELL, COBYLA, L-BFGS-B, SLSQP, and Gradient-Descent. By customizing mapping and ansatz methods, users can also construct or optimize quantum circuits themselves to obtain the optimal ground-state energy solution for electronic structure problems. 

* Welcome to fork our project on GitHub or provide feedback at `pyChemiQ Project <https://github.com/OriginQ/pyChemiQ/>`_.
* To customize pyChemiQ for specific quantum chemistry problems, please contact Mr. Chen at Origin Quantum via (+86)18221003869 or send us an email: `dqa@originqc.com <mailto:dqa@originqc.com>`_.
* To experience all the latest features in our visual teaching software ChemiQ, please download at `official website <https://originqc.com.cn/product/en/Chemiq>`_.
* To cite pyChemiQ or ChemiQ, please use the following format: Wang Q, Liu H Y, Li Q S, et al. Chemiq: A chemistry simulator for quantum computer[J]. arXiv preprint arXiv:2106.10162, 2021. 





.. toctree::
   :maxdepth: 2
   :caption: Installation
   
   01install/install.rst


.. toctree::
   :maxdepth: 2
   :caption: Quick Start
   
   02start/quickstart.rst


.. toctree::
   :maxdepth: 2
   :caption: Basic Tutorial
   
   03basis/hamiltonian.rst
   03basis/mapping.rst
   03basis/ansatz.rst
   03basis/algorithm.rst
   03basis/function.rst

.. toctree::
   :maxdepth: 2
   :caption: Advanced Tutorial
   
   04advanced/fermionpauliop.rst
   04advanced/optimizer.rst
   04advanced/circuit.rst

.. toctree::
   :maxdepth: 1
   :caption: API

   06API/index
   06API/Configs.rst
   06API/Examples.rst


.. toctree::
   :maxdepth: 2
   :caption: Feedback
   
   07FAQ/feedback.rst
   
