Installation Introduction
====================================

We provide pre-compiled Python packages for installation on Linux, macOS, and Windows. Currently, pyChemiQ supports Python versions 3.8, 3.9, and 3.10.

If you have already installed the Python environment and pip tool, pyChemiQ can be installed using the following command:

.. code-block::

   pip install pychemiq


pyChemiQ provides versatile and flexible interfaces for generating Hamiltonians, creating custom ansatz circuits, and optimizing circuit parameters. In addition to performing calculations using pyChemiQ's basic APIs, you can also run simulations directly via configuration files. Advanced features require license authorization and are invoked through the configuration file.

You can reduce quantum circuit depth and execution time by using built-in optimization methods. The interface also offers a range of features such as setting the number of slices, ansatz truncation, and MP2 initial parameter configuration.

If you wish to trial these advanced features of pyChemiQ, please visit the `Origin Quantum Store <https://mall.originqc.com.cn>`_ to purchase a license. The license for ChemiQ and pyChemiQ is universal. If you already have a ChemiQ license, simply enter it into the license parameter in the global settings section of the configuration file.

For detailed parameter configuration, please refer to the Configuration File Parameters section :doc:`../06API/Configs`.