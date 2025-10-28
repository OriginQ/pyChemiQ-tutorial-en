:py:class:`pychemiq.ChemiQ`
=============================

Classes
----------

.. py:class:: ChemiQ()

   The packaged quantum circuit class is a part of the VQE algorithm that is performed on quantum computers/virtual machines. This includes building parameterized quantum circuits to prepare trial wavefunctions and measuring and summing the various sub terms of the Hamiltonian.

   .. py:function:: prepare_vqe(machine_type, mapping_type, n_elec, pauli_size, n_qubits)

   Prepare the VQE algorithm quantum circuit function.

   :param QMachineType machine_type: Enter the type of quantum simulator. Currently, pyChemiQ only supports single threaded CPUs, namely QMachineType.CPU_SINGLE_THREAD. The integration of noisy quantum simulators is still ongoing. Please refer to `pyqpanda.QMachineType <https://pyqpanda-toturial.readthedocs.io/zh/latest/autoapi/pyqpanda/index.html#pyqpanda.QMachineType>`_  for an introduction to this category.

   :param MappingType mapping_type:  Enter the mapping type. Please refer to pychemiq.Transform.Mapping for details.
   :param int n_elec: Enter the number of electrons in the molecular system.
   :param int pauli_size: Enter the number of terms for the Pauli Hamiltonian.
   :param int n_qubits: Enter the total number of qubits required for calculation.

   :return: void



   .. py:function:: get_expectation_value(index, fcalls, task_index, qvec, hamiltonian, para, ansatz, extra_measure)

   Measure the quantum circuit to obtain the expected value of the Hamiltonian.

   :param int index: Enter the system number.
   :param int fcalls: Enter the number of function calls.
   :param int task_index: Task number.
   :param QVec qvec: An array that stores quantum bits. The introduction of this class is detailed in `pyqpanda.QVec <https://pyqpanda-toturial.readthedocs.io/zh/latest/autoapi/pyqpanda/index.html#pyqpanda.QVec>`_ .
   :param Hamiltonian hamiltonian: Enter the Hamiltonian class. In PauliOperator, the Hamiltonian is stored in string format, which is not convenient for subsequent processing. In contrast, the Hamiltonian stores Pauli operators by converting them into a custom Hamiltonian class, making it easier to extract information for each term.
   :param list[float] para: Specify the initial parameters to be optimized.
   :param AbstractAnsatz ansatz: Specify the proposed class. Please refer to pychemiq.Circuit.Ansatz for details.
   :param bool extra_ Measure: Used to distinguish between noisy and non noisy simulations. Boolean value. Set False to non noise simulation.

   :return: The expected value of the Hamiltonian, which is the ground state energy. Double-precision floating-point number.

   .. py:function:: get_loss_func_value(index, para, grad, iters, fcalls, pauli, qvec, ansatz)


   Obtain the value of the loss function.

    :param int index: Enter the system number.
    :param list[float] para: Specify the initial parameters to be optimized.
    :param list[float] gradient: Specify the initial gradient to be optimized.
    :param int iters: Enter the number of iterations for the function.
    :param int fcalls: Number of input function calls.
    :param PauliOperator Pauli: The Pauli Hamiltonian of the specified molecule. Pauli operator class. Please refer to pychemiq.PauliOperator for details.
    :param QVec qvec: An array that stores quantum bits. Please refer to for an introduction to this category.
    :param AbstractAnsatz ansatz: Specify the proposed class. Please refer to pychemiq.Circuit.Ansatz for details.

   :return: The value of the loss function. Dict type.


   .. py:function:: get_qstate()

   :return: Output the final quantum state, i.e., directly print a vector of size :math:`2^N` after optimization. Array of double-precision floating-point numbers.      


   .. py:function:: get_energy_history()

   :return: The energy value after each iteration of the function. Double-precision floating-point array.

---------


**Interface example:**

.. code:: 

      from pychemiq import Molecules,ChemiQ,QMachineType,PauliOperator
      from pychemiq.Transform.Mapping import MappingType
      from pychemiq.Circuit.Ansatz import UCC

      chemiq = ChemiQ()
      machine_type = QMachineType.CPU_SINGLE_THREAD
      mapping_type = MappingType.Jordan_Wigner
      chemiq.prepare_vqe(machine_type,mapping_type,2,1,4)

      ansatz = UCC("UCCD",2,mapping_type,chemiq=chemiq)
      pauli = PauliOperator("Z0 Z1 ",0.1)
      # Expected Value Function and Cost Function
      H = pauli.to_hamiltonian(True)
      result1 = chemiq.getExpectationValue(0,0,0,chemiq.qvec,H,[0],ansatz,False)
      result2 = chemiq.getLossFuncValue(0,[0],[0],0,0,pauli,chemiq.qvec,ansatz)
      energies = chemiq.get_energy_history()

      print(result1)
      print(result2)
      print(energies)

The printed result is:

.. parsed-literal::

      0.1
      ('', 0.1)
      [0.1]