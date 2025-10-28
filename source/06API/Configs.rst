Introduction to Configuration File Parameters
====================================================

In addition to performing calculations using pyChemiQ's basic APIs, you can also run simulations directly via configuration files. Advanced features require license authorization and are invoked through the configuration file. You can reduce quantum circuit depth and execution time by using built-in optimization methods. The interface also offers a range of features such as setting the number of slices, ansatz truncation, and MP2 initial parameter configuration. 

If you wish to trial these advanced features of pyChemiQ, please visit the `Origin Quantum Store <https://mall.originqc.com.cn>`_ to purchase a license. The license for ChemiQ and pyChemiQ is universal. If you already have a ChemiQ license, simply enter it into the license parameter in the global settings section of the configuration file. After setting the configuration file, enter the following command line on the terminal to run the task:

.. code-block::

    from pychemiq import directly_run_config
    directly_run_config("test.chemiq")  # Replace test with the name of your configuration file


The configuration file typically has the extension `.chemiq` (it needs to be placed in the code execution path). More examples of configuration files will be provided in the next section. The settings are primarily divided into six aspects, and the detailed parameter descriptions along with their default values are as follows:

1. General settings
   
Required parameters:

    - task(str): Set the calculation type [energy/PES/MD], that is, single point energy calculation/potential energy surface/molecular dynamics simulation. The default value is energy. 

    - backend(str): Set the backend type for executing calculations [CPU_SINGLE_THREAD]. Currently, only a single-threaded CPU is supported. The development of more backend types is still in progress.

    - license(str): Set the authorization serial number.

Required parameters when the computation type is PES:

    - PES_atoms(int): Specify the atomic indices involved in the potential energy surface (PES) calculation. When scanning bond length, bond angle, or dihedral angle, 2, 3, or 4 atoms are involved, respectively. The atomic index corresponds to the row number of the atom in the molecular coordinate list — the atom in the first row has index 1, and so on. 

    - PES_fast_generation(bool): Specify the method for generating the scanning coordinates of the potential energy curve: when set to T, coordinates are generated automatically after assigning three values in `PES_values`; when set to F, they are generated in a user-defined manner by directly appointing specific bond length or angle in `PES_values`. 

    - PES_values(float): When the `PES_fast_generation` is set to True, three values must be specified: start value, end value, and number of points (int). When the `PES_fast_generation` is set to False, the specific scanning coordinates must be provided. For bond length scans, the unit is angstrom (Å); For bond angle scans, the unit is degrees (°). There is no limit on the number of scanning points if the `PES_fast_generation` is set to False. 

Optional parameters:

    - chem_method(str): Specify which classical computational method to use: [HF/CCSD]. Default is HF.

    - hamiltonian(str): This keyword allows you to define a custom Hamiltonian for computation. Note that when using a custom Hamiltonian, you must set the keyword `nelec`, which specifies the number of electrons.

    - print_out(bool): Set whether to print the scf iteration process. The default is F.

    - print_iters(bool): Set whether to print the optimizer iteration process. The default is F.

    - console_level(int): Set the log level for terminal output: 0 means output enabled, 6 means no output. The default value is 0.

    - logfile_level(int): Set the log level for file output: 0 means output enabled, 6 means no output. The default value is 6.

    - logfile_name(str): Set the log file name; default is empty. For example, if the specified log name is `chemiq.log`, the final output log file will be named `chemiq-YYYYMMDD.log` (where YYYYMMDD is the current date). A log file will only be generated if a log file name is explicitly set.


2. Molecular specification setting
   
Required parameters:

    - mole: Configure the parameters related to the molecular model. The default value is empty. This is followed by sub-parameter settings.

    - geoms(str): Set the molecular coordinates. The first column contains the element symbol (and serial number), followed by the atomic coordinates (x, y, z coordinates) separated by spaces.

    - charge(int): Set the number of charges in the system. When the number of charges is positive, there is no positive sign, and when it is negative, a negative sign is written. The default value is 0.

    - spin(int): Set the spin multiplicity (M=2S+1), with a default value of 1.

    - basis(str): Set the basis set level required for computation: [MINI / sto-3G / sto-6G / 3-21G / 6-31G / ...]. Starting from version V1.1.0, pyChemiQ supports over 600 commonly used basis sets as well as user-defined basis sets. It supports most basis sets available on the Basis Set Exchange website (excluding a few excessively large ones), including but not limited to STO-nG basis sets, Pople series basis sets, Ahlrichs' def-series basis sets, Dunning series basis sets, and pseudopotential basis sets. In previous versions, supported basis sets included: MINI, STO-3G, STO-6G, 3-21G, and 6-31G. Other basis sets can be directly entered by name. The specific string input rules for basis set names are as follows:
  
        - Convert uppercase letters to lowercase. Eg: STO-3G → sto-3g

        - Preserve the plus sign + in diffuse basis sets. Eg: 6-31+G → 6-31+g

        - Retain parentheses `(` and `)` in polarized basis sets; convert asterisk * to underscore _.  Eg: 6-31G(d) → 6-31g(d) ; 6-31G* → 6-31g\_ 

        - Remove spaces, hyphens -, and forward slashes /. Eg: def2-SVP → def2svp

    To use a custom basis set, please place the `.g94` formatted basis set file in the \\ChemiQ\\basis folder under the ChemiQ installation directory (default path: C:\\Users\\YOURUSERNAME\\AppData\\Local\\Programs\\ChemiQ). Please ensure that the basis set name entered in the basis set field exactly matches the filename of the basis set file; otherwise, the program will return an error: `XXX.g94 not exist!`

Optional parameters: 

    - diis(str): Specifies whether to use DIIS (Direct Inversion in the Iterative Subspace) to accelerate SCF convergence. Options: [cdiis / none(empty)]. 

    - diis_n(int): The history length for cdiis — i.e., the number of previous density matrices used to compute the next one. Default: 8. Only effective when `diis = cdiis`.

    - diis_thre(float): The RMSD threshold between consecutive density matrices below which cdiis is activated. Default: 0.1. Only effective when `diis = cdiis`.

    - pauli_group(str): Specifies whether to use Pauli grouping [native / none(empty)]. Default: none. Pauli grouping is a technique that reduces the number of measurements required to estimate the expectation value of the Hamiltonian by analyzing commutation relations among its Pauli terms.

    - pauli_reverse(bool): Sets the qubit ordering to standard or reverse order. Default: T (True), meaning `q3-q2-q1-q0` order.

    - bohr(bool): Specifies whether coordinates are given in bohr units. Default: F (False), meaning angstroms are used.

    - pure(bool): Specifies whether to use spherical harmonic-type or Cartesian-type Gaussian functions. Default: T (True), meaning spherical harmonic-type Gaussians are used.

    - local(bool): Specifies whether to localize HF orbitals. Default: F (False), meaning HF orbitals are not localized.

    - active(int, int): Sets the active space, specified as two comma-separated integers. Default: not set. For example, `active = 4, 4` means 4 active spatial orbitals with 4 electrons in the active space.

    - nfrozen(int): Sets the number of frozen orbitals. Default: 0. Note: `active` and `nfrozen` cannot be used simultaneously.

    - mix_scf(float): This parameter can effectively resolve SCF non-convergence issues using a damping method. Specifically, the density matrix for step n+1, D(n+1), is modified to: :math:`w \times D(n-1) + (1-w) \times D(n)` , where w is the user-specified mixing parameter. Averaging the density matrix smoothens changes between iterations, promoting convergence. Valid range: [0.0, 1.0]. Default: 0.5.

3. Ansatz parameter settings
   
Required parameters:

    - ansatz(str): Specifies the type of quantum circuit ansatz. Options: [UCC/Hardware-efficient/Symmetry-preserved/User-define]. For the first three types, the quantum circuit is automatically generated. For the `User-define` option, you must either provide a quantum circuit in originIR format via the `circuit` parameter, or define it using Pauli operators via the `pauli` parameter. For details on the originIR format, see:  `originIR Format Introduction <http://10.10.10.138/qpanda-3/d7/d65/tutorial_quantum_program_of_content__origin_i_r.html>`_ . 

    - mapping(str): Set mapping method. Options: [JW/P/BK/SP]. These mapping methods are Jordan Wigner Transform, Parity Transform, Bravyi Kitaev Transform, and Segment Parity Transform.

    - excited_level(str): Set the excitation level. Options: [S/D/SD]. Required when ansatz is set to `UCC`.
  
    - circuit(str): Set the quantum circuit using an originIR string. This parameter is required when ansatz is set to `User-define`.

Optional parameters: 

    - restricted(bool): Restrict excitation terms to prepare a superposition of fewer configuration states, thereby shortening the circuit. Default: T (True). Only effective when ansatz is set to `UCC`.

    - cutoff(bool): Truncate the excitation terms in the UCC ansatz based on initial parameters from MP2. Only effective when ansatz is set to `UCC` and init_para_type is set to `MP2`. Default: F (False). 

    - reorder(bool): Arrange qubits in order such that the first half encode spin-up orbitals and the second half encode spin-down orbitals. Setting this parameter to T (True) can reduce the number of qubits required. This option is effective only when the mapping is set to Parity or Bravyi-Kitaev (BK). Default: F (False). 


4. Optimizer parameter settings

Required parameters:

    - Optimizer(str): Set the classic optimizer type. Options: [Nelder-Mead/Powell/Gradient-Descent/COBYLA/L-BFGS-B/SLSQP/GAQPSO].

    - init_Para_Type(str): Set the method for constructing initial parameters. Options: [Zero/Random/input/MP2/CCSD], where `Zero` represents that all initial parameters are set to zero, `Random` represents initial parameters are randomly sampled from the interval [0, 1), `input` represents User-defined initial parameters, and `MP2` represents initial parameters derived from second-order perturbation theory. `CCSD` represents initial parameters obtained from singles and doubles coupled-cluster calculations. `MP2` and `CCSD` are only available when the ansatz is set to `UCCD` or `UCCSD`. Default: Zero.

Optional parameters: 

    - slices(int): Set the number of slices, i.e. the number of quantum circuit repetitions, with a default value of 1.

    - learning_Rate(float): Set the learning rate. The default value is 0.1.

    - iters(int): Set the number of iterations, with a default value of 1000.

    - fcalls(int): Set the number of function calls, with a default value of 1000.

    - xatol(float): Set the variable convergence threshold, with a default value of  :math:`10^4` .

    - fatol(float): Set the expected value convergence threshold, with a default value of :math:`10^4` .

5. Molecular dynamics parameter settings

Required parameters:
   
    - MD: Set the correlation sampling method. The default is 1.

Optional parameters: 

    - axis(str): Set the system to move in a specific direction in the form of a string, in the format of "x y z".

    - save_trajectory(str): Set the name of the saved molecular coordinate file. The default is `traj.csv`.

    - save_topology(str): Set the name of the saved molecular topology file. The default is `topology.txt`.

    - velocity(float): Set the initial velocity of atoms, separated by commas, "0.1 0.2 0.3, -0.1-0.2-0.3", in units of angstom/fs, with default values of all 0.

    - step_size(float): Set the step size, must be greater than 0, in fs, with a default of 0.2.

    - step_number(int): Set the total number of steps, must be greater than 1, with a default of 100.

    - delta_r(float): Set the size of the differential coordinate, must be greater than 0, with a default of 0.001.

6. Real quantum chip computing settings

Required parameters:

    - cloud_api_key (str): API key for the cloud platform. You can query the remaining computing time on the `API key usage query <https://console.originqc.com.cn/en/usage>`_ and copy your APIKEY on `API key copy and refresh <https://console.originqc.com.cn/en/apikey>`_ .
    
    - chip_id (str): specifies which superconducting quantum chip to use, at present the 72 and 102 are both available. Check  `origin quantum cloud platform <https://console.originqc.com.cn/en/services>`_ for detailed real-time information for each superconducting quantum chip. 


Optional parameters: 

    - chip_mode (str): Set the chip task mode: [wait / submit / query / none(empty)]. `wait` indicates to submit the task and wait for the result. It will query the task status every 2 seconds for up to 1 minute(query count can be set using `wait` parameter). If no result is obtained, the backend will stop querying and return the message: "The current task is still in queuing." `submit` means to submit the task only without querying the result. `query` determines to query the task status only — in this case, task ID `chip_task_id` has to be provided. The default mode is submit.

    - chip_task_id (str): Task ID for query. Required only when chip_mode is `query`.

    - cloud_url (str): URL of the cloud platform. Default: https://pyqanda-admin.qpanda.cn.

    - shots (int): refers to the number of sampling times for measurement of quantum circuits on real quantum computers. The higher the sampling times, the smaller the statistical error, but the longer the calculation time. The default sampling times is 1000.


In this section, we provide some examples of calculations using configuration file. This includes calculations of single-point energy, potential energy surface and molecular dynamic simulations. 
In the first example, single-point energy calculation of the hydrogen molecule is performed. The basis set is set as STO-3G, the ansatz is set as UCCSD, the mapping method using BK, and the optimizer is NELDER-MEAD. Initial parameters are set to MP2.

.. code-block::

    general = {
        task    = energy
        backend = CPU_SINGLE_THREAD
        license = XXXXX
    }

    mole = {
        geoms = {
            H 0 0 0
            H 0 0 0.74
        }
        bohr    = F
        charge  = 0
        spin    = 1 
        basis   = sto-3G
        pure    = T 
        local   = F 
    }

    ansatz = UCC {
        excited_level = SD
        restricted    = T
        cutoff        = T
        mapping       = BK
        reorder       = F
    }

    optimizer = NELDER-MEAD {
        learning_rate                 = 0.1 
        init_para_type                = MP2
        slices                        = 1 
        iters                         = 1000 
        fcalls                        = 1000 
        xatol                         = 1e-6 
        fatol                         = 1e-6 
    }


In the second example, we compute the potential energy curve of the hydrogen molecule. Here, we scan 5 points, starting from a bond length of 0.6 angstrom, and ending with 1 angstrom, making step size as 0.1 angstrom. We use the STO-3G basis set, a user-defined ansatz circuit, the Parity mapping, and the SLSQP optimizer. Initial parameters are set to zero.

.. code-block::

    general = {
        task    = PES
        backend = CPU_SINGLE_THREAD
        license = XXXXX
        PES_atoms = 1,2
        PES_fast_generation = T
        PES_values = 0.6,1,5
    }

    mole = {
        geoms = {
            H 0 0 0
            H 0 0 0.54
        }
        charge  = 0
        spin    = 1
        basis   = sto-3G
    }

    ansatz = User-define {
        circuit = {
            QINIT 4
            CREG 4
            CNOT q[1],q[0]
            CNOT q[2],q[1]
            CNOT q[3],q[2]
            H q[1]
            H q[3]
            S q[1]
    }
        mapping       = P
        reorder       = T
    }

    optimizer = SLSQP {
        learning_rate                 = 0.1
        init_para_type                = Zero
        slices                        = 1
        iters                         = 1000
        fcalls                        = 1000
        xatol                         = 1e-6
        fatol                         = 1e-6
    }


In the third example, we calculate the molecular dynamics trajectory of lithium hydride molecules. We select 3-21G as basis-set, the active space uses [4,4], Hardware-efficient type ansate is proposed to be used, mapping uses JW, and optimizer uses L-BFGS-B. The initial parameter is set as random number.

.. code-block::

    general = {
        task    = MD
        backend = CPU_SINGLE_THREAD
        license = XXXXX
    }

    mole = {
        geoms = {
            H 0 0 0.38
            Li 0 0 -1.13
        }
        bohr    = F
        charge  = 0
        spin    = 1 
        basis   = 3-21G
        pure    = T 
        local   = F 
        active = 4,4
    }

    ansatz = Hardware-efficient {
        mapping       = JW
        reorder       = F
    }

    optimizer = L-BFGS-B {
        learning_rate                 = 0.1 
        init_para_type                = Random
        slices                        = 1  
        iters                         = 1000 
        fcalls                        = 1000 
        xatol                         = 1e-6 
        fatol                         = 1e-6 
    }

    MD = 1 {
        velocity           = 0.0
        step_size          = 0.2
        step_number        = 100 
        delta_r            = 0.001
    }

In the fourth example, we use a configuration file to invoke a real quantum chip for calculating the ground-state energy of the hydrogen molecule. The `license` and `cloud_api_key` keywords below must be filled in by the user.

.. code-block::

    general = {
        license = XXXXX
        task = energy {
        chip_mode        = wait
        cloud_url        = https://pyqanda-admin.qpanda.cn
        cloud_api_key    = XXXXX
        shots            = 1000
        chip_id          = 72
        chip_amend       = T
        chip_mapping     = T
        chip_circuit_opt = T
        }
    }

    mole = {
        geoms = {
            H 0 0 0
            H 0 0 0.74
        }
        charge  = 0
        spin    = 1
        basis   = sto-3G
    }

    ansatz = Hardware-efficient {
        mapping       = BK
    }

    optimizer = NELDER-MEAD {
        learning_rate                 = 0.1
        init_para_type                = Random
        slices                        = 1
        iters                         = 1000
        fcalls                        = 1000
        xatol                         = 1e-6
        fatol                         = 1e-6
    }



