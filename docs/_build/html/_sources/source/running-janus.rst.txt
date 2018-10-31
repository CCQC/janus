Running janus
=================================

Structure of an input file
--------------------------

System
--------------------------

Required keywords
_________________

----------------------------

**system_info**
    :Description: Files for describing the positions and topology of the initial system
    :DataType: List of strings
    :Values: inpcrd and prmtop files needed for Amber files, gro and top files needed for Gromacs files

Optional keywords
_________________

----------------------------

**system_info_format**
    :Description: Specifies the format of the files given in system_info
    :DataType: String
    :Default: pdb
    :Values: pdb, Amber, Gromacs
    
QMMM
--------------------------

Required keywords
_________________

----------------------------

**run_aqmmm**
    :Description: Specifies whether to run adaptive QM/MM
    :DataType: Bool

**qm_atoms**
    :Description: Specifies what atoms to treat with QM (or other high level treatment)
    :DataType: List
    :Notes: If run_aqmmm is true, the qm_atoms is the same as the qm_center

Optional keywords
_________________

----------------------------

**ll_program**
    :Description: Specifies what program to use for the lower level computation
    :DataType: String
    :Values: OpenMM
    :Default: OpenMM

**hl_program**
    :Description: Specifies what program to use for the high level computations
    :DataType: String
    :Values: OpenMM, Psi4
    :Default: Psi4

**embedding_method**
    :Description: Specifies what QM/MM embedding method to use
    :DataType: String
    :Values: Electrostatic, Mechanical
    :Default: Mechanical

**qmmm_scheme**
    :Description: Specifies what energy scheme to use for computing the QM/MM energy
    :DataType: String
    :Values: subtractive
    :Default: subtractive

**boundary_treatment**
    :Description: Specifies the scheme to use for the treatment of dangling bonds 
    :DataType: String
    :Values: link_atom, RC, RCD
    :Default: link_atom

**link_atom_element**
    :Description: Specifies what atom to use for the link atom
    :DataType: String
    :Values: H
    :Default: H

Examples
_________________

----------------------------

AQMMM
--------------------------

Required keywords
_________________

----------------------------

**qm_center**
    :Description: Specifies what atoms to designate as the center for the high level treatment
    :DataType: List

Optional keywords
_________________

----------------------------

**aqmmm_scheme**
    :Description: Specifies what adaptive QM/MM approach to use
    :DataType: String
    :Values: ONIOM-XS, Hot-Spot, PAP, SAP, DAS
    :Default: ONIOM-XS

**partition_scheme**
    :Description: Specifies how to define the buffer zone atoms
    :DataType: String
    :Values: distance
    :Default: distance

**Rmin**
    :Description: Specifies the radius from the qm center to the inner boundary of the buffer zone in distance partitioning in angstroms
    :DataType: Float
    :Default: 4.0

**Rmin**
    :Description: Specifies the radius from the qm center to the outer boundary of the buffer zone in distance partitioning in angstroms
    :DataType: Float
    :Default: 4.5

**modified_variant**
    :Description: Specifies whether to use the modified variant of either the PAP or SAP schemes
    :DataType: Bool
    :Default: False

Examples
_________________

----------------------------

Molecular Dynamics
--------------------------

High Level 
--------------------------

Low Level
--------------------------

Required keywords
_________________

----------------------------

**md_simulation_program**
    :Description: Specifies what program to use for MD simulation
    :DataType: String
    :Values: OpenMM

**start_qmmm**
    :Description: Specifies at which step to start the QM/MM or adaptive QM/MM approach
    :DataType: Int
    :Notes: QM/MM can be started after taking some MD steps so the system can equilibrate

**end_qmmm**
    :Description: Specifies how many 
    :DataType: Int
    :Notes: end_qmmm - start_qmmm = total number of MD steps that will use QM/MM forces
    
Optional keywords
_________________

----------------------------

**step_size**
    :Description: The step size of the MD simulation in femtoseconds
    :DataType: Int 
    :Default: 1

**md_ensemble**
    :Description: Ensemble of MD simulation
    :DataType: String or List of Strings
    :Values: NVT, NVE
    :Default: NVE
    :Notes: If more than one ensemble is desired (i.e., NVT run before NVE) a list can be created in the 
              order of what is run

**md_steps**
    :Description: Specifies how many total steps to take for the MD simulation
    :DataType: Int or List of Ints
    :Default: end_qmmm 
    :Notes: If more than one md_ensemble is desired, the steps for each can be specified in a list where
              each element will correspond to the steps run in each ensemble specified in md_ensemble. The step number specified 
              with start_qmmm will be taken as the step number of the last ensemble specified at which to start QM/MM

**return_trajectory**
    :Description: Whether to return the trajectory of the MD simulation. Keyword value lists the frame interval to save.
    :DataType: Int 
    :Default: 0 (trajectory not returned)

**return_trajectory_filename**
    :description: name of trajectory file to return
    :datatype: string
    :default: output

**trajectory_format**
    :Description: The format of the trajectory file to return
    :DataType: String
    :Values: NetCDF,
    :Default: NetCDF

**return_system**
    :Description: Whether to return the final position and topology of the system in a pdb file
    :DataType: Bool
    :Default: False

**return_system_filename**
    :description: name of system file to return
    :datatype: string
    :default: final.pdb

**return_info**
    :Description: Whether to return system information such as energy and temperature
    :DataType: List of strings with values to return, will be returned in file "info.dat"
    :Values: potentialEnergy, kineticEnergy, totalEnergy, temperature, density
    :Default: []

**return_info_interval**
    :Description: The frame interval for saving energy, etc. information.
    :DataType: Int
    :Default: 0 (info not returned)



Examples
_________________

----------------------------


Supported codes
-----------------------
Janus only supports Psi4 for quantum mechanics computations and
OpenMM for molecular mechanics and molecular dynamics.
