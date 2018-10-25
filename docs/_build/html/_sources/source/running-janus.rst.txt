Running janus
=================================

Structure of an input file
--------------------------

System
__________________________

Molecular Dynamics
--------------------------

Required keywords
~~~~~~~~~~~~~~~~~

**md_simulation_program**
    - :Description: Specifies what program to use for MD simulation
    - :DataType: String
    - :Possible Values: OpenMM

**start_qmmm**
    - :Description: Specifies at which step to start the QM/MM or adaptive QM/MM approach
    - :DataType: Int
    - :Notes: QM/MM can be started after taking some MD steps so the system can equilibrate

**end_qmmm**
    - :Description: Specifies how many 
    - :DataType: Int
    - :Notes: end_qmmm - start_qmmm = total number of MD steps that will use QM/MM forces
    
Optional keywords
~~~~~~~~~~~~~~~~~

**step_size**
    - :Description: The step size of the MD simulation in femtoseconds
    - :DataType: Int 
    - :Default: 1

**md_ensemble**
    - :Description: Ensemble of MD simulation
    - :DataType: String or List of Strings
    - :Possible Values: NVT, NVE
    - :Default: NVE
    - :Notes: If more than one ensemble is desired (i.e., NVT run before NVE) a list can be created in the 
              order of what is run

**md_steps**
    - :Description: Specifies how many total steps to take for the MD simulation
    - :DataType: Int or List of Ints
    - :Default: end_qmmm - start_qmmm 
    - :Notes: If more than one md_ensemble is desired, the steps for each can be specified in a list where
              each element will correspond to the steps run in each ensemble specified in md_ensemble. The step number specified 
              with start_qmmm will be taken as the step number of the last ensemble specified at which to start QM/MM


Examples
~~~~~~~~


Supported codes
-----------------------
