Configuration File Example
============================

In the previous section, we provided examples of several common configuration files. In this section, we will present more configuration file examples.


1. When constructing the quantum circuit, we can use Pauli operators to customize an ansatz by using **User-define** keywords with **pauli**.

.. code-block::

    general = {
        task    = energy
        backend = CPU_SINGLE_THREAD
        license = XXXXX
    }

    mole = {
        geoms = {
            H 0.000000 0.000000 0.356115
            H 0.000000 0.000000 -0.356115
        }
        charge  = 0
        spin    = 1
        basis   = sto-3g
    }

    ansatz = User-define{
        pauli = {
            X0 Y2 : -0.125000
            X0 Y2 Z3 : -0.125000
            X0 Z1 Y2 : -0.125000
            X0 Z1 Y2 Z3 : -0.125000
            Y0 X2 : 0.125000
            Y0 X2 Z3 : 0.125000
            Y0 Z1 X2 : 0.125000
            Y0 Z1 X2 Z3 : 0.125000
        }
        mapping  = BK
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

2. We can also constrain the optimizing parameters via linear combination for the Pauli terms in your customized circuit. By bundling parameters, the number of parameters in the classical optimizer is reduced, which facilitates the optimizer in obtaining results more quickly.

.. code-block::

    general = {
        task    = energy
        backend = CPU_SINGLE_THREAD
        license = XXXXX
    }

    mole = {
        geoms = {
            H 0 0 0.38
            Li 0 0 -1.13
        }
        charge  = 0
        spin    = 1 
        basis   = sto-3g
        active  = 4,4
    }

    ansatz = User-define{
        pauli = {
            X0 Y2 : -0.125
            X0 Y2 Z3 : -0.125
            X0 Z1 Y2 : -0.125
            X0 Z1 Y2 Z3 : -0.125
            Y0 X2 : 0.125
            Y0 X2 Z3 : 0.125
            Y0 Z1 X2 : 0.125
            Y0 Z1 X2 Z3 : 0.125
        }
        mapping       = BK
    }
    optimizer = NELDER_MEAD {
        learning_rate                 = 0.1
        init_para_type                = Zero
        # There are 8 parameters x1,x2,...,x8 since there are 8 pauli items
        # only two variables t1,t2 is used for optimizing
        # x1,x2,x3,x4 = t2
        # x5,x6,x7,x8 = t1
        parameter_matrix = {
            0 0 0 0 1 1 1 1
            1 1 1 1 0 0 0 0
        }
        slices = 1
        iters                         = 200
        fcalls                        = 200
        xatol                         = 1e-4
        fatol                         = 1e-4
    }


3. In addition to customizing the quantum circuit ansatz using Pauli operators, we can also employ originIR for this purpose. 

.. code-block::

    general = {
        task    = energy
        backend = CPU_SINGLE_THREAD
        license = XXXXX
    }

    mole = {
        geoms = {
            H 0.000000 0.000000 0.356115
            H 0.000000 0.000000 -0.356115
        }
        charge  = 0
        spin    = 1
        basis   = sto-3g
    }

    ansatz = User-define{
        circuit = {
            H q[0]
            RX q[2],(1.5707963)
            CNOT q[0],q[3]
            CNOT q[1],q[3]
            CNOT q[2],q[3]
            RZ q[3],(1.5707963)
            CNOT q[0],q[3]
            CNOT q[1],q[3]
            CNOT q[2],q[3]
            DAGGER
            H q[0]
            RX q[2],(1.5707963)
            ENDDAGGER
        }
        mapping  = BK
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


4. When using a custom molecular configuration, you can directly specify the molecule for calculation by providing a custom Hamiltonian. In the example below, the custom Hamiltonian is defined in the bracket behind the keywords **hamiltonian**, where the string before the quotation marks denotes the Hamiltonian as a Pauli string, and the number after it is the corresponding coefficient.

.. code-block::

    general = {
        task    = energy
        backend = CPU_SINGLE_THREAD
        license = XXXXX
    }

    mole = {
        hamiltonian = {
            : -0.059663
            X0 Z1 X2 : 0.044918
            X0 Z1 X2 Z3 : 0.044918
            Y0 Z1 Y2 : 0.044918
            Y0 Z1 Y2 Z3 : 0.044918
            Z0 : 0.175755
            Z0 Z1 : 0.175755
            Z0 Z1 Z2 : 0.167143
            Z0 Z1 Z2 Z3 : 0.167143
            Z0 Z2 : 0.122225
            Z0 Z2 Z3 : 0.122225
            Z1 : 0.170014
            Z1 Z2 Z3 : -0.236656
            Z1 Z3 : 0.175702
            Z2 : -0.236656
        }
        # nelec is needed when the hamiltonian is user-defined
        nelec = 2
    }
    ansatz = UCC {
        excited_level = SD
        mapping       = BK
    }
    optimizer = NELDER_MEAD {
        learning_rate                 = 0.1
        init_para_type                = Zero
        slices = 1
        iters                         = 200
        fcalls                        = 200
        xatol                         = 1e-4
        fatol                         = 1e-4
    }

5. When performing a potential energy surface (PES) scan, there are three variables can be scanned: bond length, bond angle and dihedral angle. For example, if one wants to scan bond angle as PES, the index of three atoms needs to be defined in `PES_atoms` and the specific angles have to be provided in `PES_values`. 

.. code-block::

    general = {
        task    = PES
        backend = CPU_SINGLE_THREAD
        license = XXXXX
        PES_atoms = 1,2,3
        PES_values = 30,60,90,120
    }

    mole = {
        geoms = {
            H 0.625276 0.625276 0.625276
            C 0.000000 0.000000 0.000000
            H -0.625276 -0.625276 0.625276
            H -0.625276 0.625276 -0.625276
            H 0.625276 -0.625276 -0.625276
        }
        charge  = 0
        spin    = 1
        basis   = sto-3g
        active  = 4,4
    }
    ansatz = UCC{
        excited_level = SD
        mapping       = BK
    }
    optimizer = SLSQP {
        slices = 1
        learning_rate                 = 0.1
        init_para_type                = MP2
        iters                         = 200
        fcalls                        = 200
        xatol                         = 1e-4
        fatol                         = 1e-4
    }

6. In this example, we utilize the GAQPSO optimizer for variational quantum circuit optimization. This algorithm requires the configuration of several specific parameters:  
`pso_wi`: Maximum inertia weight.  
`pso_we`: Minimum inertia weight.   
`pso_c1`: Cognitive coefficient, the coefficient directing a particle toward its personal best position.   
`pso_c2`: Social coefficient, the coefficient directing a particle toward the global best position.   
`pso_glr`: Learning rate used when updating particle positions with GD or SPSA.   
`pso_deltap`: Perturbation variable for the GD method.   
`pso_thres / pso_thresf`: Threshold for the difference in the objective function value between the current and previous iteration.   
`pso_nearenough`: Distance threshold for a particle to be considered sufficiently close to the global optimum G.   
`pso_cmax`: Maximum number of iterations for the local search.  
`pso_alpha / pso_alphae`: Expansion-contraction coefficient.   

The following parameters are global settings for GAQPSO algorithm:  
`pso_repeatn`: Maximum number of iterations for the algorithm.  
`pso_seed`: Random seed for result reproducibility.    
`pso_n`: Population size.   
`pso_prefix`: Default directory name for storing output.   

.. code-block::

    general = {
        task    = energy
        backend = CPU_SINGLE_THREAD
        license = XXXXX
    }

    mole = {
        geoms = {
            H 0.000000 0.000000 0.356115
            H 0.000000 0.000000 -0.356115
        }
        charge  = 0
        spin    = 1
        basis   = sto-3g
    }

    ansatz = UCC {
        excited_level = SD
        mapping       = BK
    }

    optimizer = GAQPSO {
        learning_rate                 = 0.01
        init_para_type                = MP2

        pso_wi         = 0.8
        pso_we         = 0.8
        pso_c1         = 2.0
        pso_c2         = 3.0
        pso_glr        = 0.01
        pso_deltap     = 1e-4
        pso_thres      = 0.001
        pso_thresf     = 1e-5
        pso_nearenough = 1e-2
        pso_alpha      = 1.0
        pso_alphae     = 0.5
        pso_repeatn    = 1
        pso_seed       = 300
        pso_n          = 20
        pso_cmax       = 2
        pso_prefix     = pso_result

        hamiltonian_simulation_slices = 1
        iters                         = 100
        xatol                         = 1e-8
        fatol                         = 1e-8 
    }
