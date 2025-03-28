import streamlit as st

def main():
    st.title("Compact VASP INCAR Generator")

    # 1) Dictionary of all calculation blocks:
    #    Each key is the short code, and the value is the full text block.
    calc_blocks = {
        "ST": """###Global Parameters
ISTART =  1            (Read existing wavefunction, if there)
ISPIN  =  1            (Non-Spin polarised DFT)
# ICHARG =  11         (Non-self-consistent: GGA/LDA band structures)
LREAL  = .FALSE.       (Projection operators: automatic)
# ENCUT  =  400        (Cut-off energy for plane wave basis set, in eV)
# PREC   =  Accurate   (Precision level: Normal or Accurate, set Accurate when perform structure lattice relaxation calculation)
LWAVE  = .TRUE.        (Write WAVECAR or not)
LCHARG = .TRUE.        (Write CHGCAR or not)
ADDGRID= .TRUE.        (Increase grid, helps GGA convergence)
LASPH  = .TRUE.        (Give more accurate total energies and band structure calculations)
PREC   = Accurate      (Accurate strictly avoids any aliasing or wrap around errors)
# LVTOT  = .TRUE.      (Write total electrostatic potential into LOCPOT or not)
# LVHAR  = .TRUE.      (Write ionic + Hartree electrostatic potential into LOCPOT or not)
# NELECT =             (No. of electrons: charged cells, be careful)
# LPLANE = .TRUE.      (Real space distribution, supercells)
# NWRITE = 2           (Medium-level output)
# KPAR   = 2           (Divides k-grid into separate groups)
# NGXF    = 300        (FFT grid mesh density for nice charge/potential plots)
# NGYF    = 300        (FFT grid mesh density for nice charge/potential plots)
# NGZF    = 300        (FFT grid mesh density for nice charge/potential plots)

Static Calculation
ISMEAR =  0            (gaussian smearing method)
SIGMA  =  0.05         (please check the width of the smearing)
LORBIT =  11           (PAW radii for projected DOS)
NEDOS  =  2001         (DOSCAR points)
NELM   =  60           (Max electronic SCF steps)
EDIFF  =  1E-08        (SCF energy convergence, in eV)
""",

        "MG": """###Collinear Magnetic Calculation
ISPIN      =  2        (Spin polarised DFT)
# MAGMOM   =           (Set this parameters manually)
LASPH      = .TRUE.    (Non-spherical elements, d/f convergence)
GGA_COMPAT = .FALSE.   (Apply spherical cutoff on gradient field)
VOSKOWN    =  1        (Enhances the magnetic moments and the magnetic energies)
LMAXMIX    =  4        (For d elements increase LMAXMIX to 4, f: LMAXMIX = 6)
# AMIX       =  0.2    (Mixing parameter to control SCF convergence)
# BMIX       =  0.0001 (Mixing parameter to control SCF convergence)
# AMIX_MAG   =  0.4    (Mixing parameter to control SCF convergence)
# BMIX_MAG   =  0.0001 (Mixing parameter to control SCF convergence)
""",

        "D3": """###DFT-D3 Correction
IVDW   =  11           (DFT-D3 method of method with no damping)
""",

        "PU": """###DFT+U Calculation
LDAU    = .TRUE.        (Activate DFT+U)
LDAUTYPE=  2            (Dudarev, only U-J matters)
LDAUL   =  2 -1         (Orbitals for each species)
LDAUU   =  2  0         (U for each species)
LDAUJ   =  0  0         (J for each species)
LMAXMIX =  4            (Mixing cut-off, 4-d, 6-f)
""",

        "GW": """###Global Parameters
ISTART =  1            (Read existing wavefunction, if there)
ISPIN  =  1            (Non-Spin polarised DFT)
# ICHARG =  11         (Non-self-consistent: GGA/LDA band structures)
LREAL  = .FALSE.       (Projection operators: automatic)
# ENCUT  =  400        (Cut-off energy for plane wave basis set, in eV)
# PREC   =  Accurate   (Precision level: Normal or Accurate, set Accurate when perform structure lattice relaxation calculation)
LWAVE  = .TRUE.        (Write WAVECAR or not)
LCHARG = .TRUE.        (Write CHGCAR or not)
ADDGRID= .TRUE.        (Increase grid, helps GGA convergence)
LASPH  = .TRUE.        (Give more accurate total energies and band structure calculations)
PREC   = Accurate      (Accurate strictly avoids any aliasing or wrap around errors)
# LVTOT  = .TRUE.      (Write total electrostatic potential into LOCPOT or not)
# LVHAR  = .TRUE.      (Write ionic + Hartree electrostatic potential into LOCPOT or not)
# NELECT =             (No. of electrons: charged cells, be careful)
# LPLANE = .TRUE.      (Real space distribution, supercells)
# NWRITE = 2           (Medium-level output)
# KPAR   = 2           (Divides k-grid into separate groups)
# NGXF    = 300        (FFT grid mesh density for nice charge/potential plots)
# NGYF    = 300        (FFT grid mesh density for nice charge/potential plots)
# NGZF    = 300        (FFT grid mesh density for nice charge/potential plots)

Electronic Relaxation
ISMEAR =  0
SIGMA  =  0.05
EDIFF  =  1E-08
""",

        "DC": """###Global Parameters
ISTART =  1            (Read existing wavefunction, if there)
ISPIN  =  1            (Non-Spin polarised DFT)
# ICHARG =  11         (Non-self-consistent: GGA/LDA band structures)
LREAL  = .FALSE.       (Projection operators: automatic)
# ENCUT  =  400        (Cut-off energy for plane wave basis set, in eV)
# PREC   =  Accurate   (Precision level: Normal or Accurate, set Accurate when perform structure lattice relaxation calculation)
LWAVE  = .TRUE.        (Write WAVECAR or not)
LCHARG = .TRUE.        (Write CHGCAR or not)
ADDGRID= .TRUE.        (Increase grid, helps GGA convergence)
LASPH  = .TRUE.        (Give more accurate total energies and band structure calculations)
PREC   = Accurate      (Accurate strictly avoids any aliasing or wrap around errors)
# LVTOT  = .TRUE.      (Write total electrostatic potential into LOCPOT or not)
# LVHAR  = .TRUE.      (Write ionic + Hartree electrostatic potential into LOCPOT or not)
# NELECT =             (No. of electrons: charged cells, be careful)
# LPLANE = .TRUE.      (Real space distribution, supercells)
# NWRITE = 2           (Medium-level output)
# KPAR   = 2           (Divides k-grid into separate groups)
# NGXF    = 300        (FFT grid mesh density for nice charge/potential plots)
# NGYF    = 300        (FFT grid mesh density for nice charge/potential plots)
# NGZF    = 300        (FFT grid mesh density for nice charge/potential plots)

Elastic constants Calculation
# NGZF    = 300        (FFT grid mesh density for nice charge/potential plots)
IBRION    =  6         (Determine the Hessian matrix)
NFREE     =  2         (How many displacements are used for each direction, 2-4)
ISIF      =  3         (Stress/relaxation: 3-Shape/Ions/V)
NSW       =  1         (Max ionic steps)
POTIM     = 0.015      (Displacement of each ion)
# ENCUT   =  700       (1.3 ~ 1.5 * default cutoff, need to check convergence)
""",

        "BD": """###Bader Charge Analysis
LAECHG     = .TRUE.    (Write core charge into CHGCAR file)
LCHARG     = .TRUE.    (Write CHGCAR file)
""",

        "EC": """###Global Parameters
ISTART =  1            (Read existing wavefunction, if there)
ISPIN  =  1            (Non-Spin polarised DFT)
# ICHARG =  11         (Non-self-consistent: GGA/LDA band structures)
LREAL  = .FALSE.       (Projection operators: automatic)
# ENCUT  =  400        (Cut-off energy for plane wave basis set, in eV)
# PREC   =  Accurate   (Precision level: Normal or Accurate, set Accurate when perform structure lattice relaxation calculation)
LWAVE  = .TRUE.        (Write WAVECAR or not)
LCHARG = .TRUE.        (Write CHGCAR or not)
ADDGRID= .TRUE.        (Increase grid, helps GGA convergence)
LASPH  = .TRUE.        (Give more accurate total energies and band structure calculations)
PREC   = Accurate      (Accurate strictly avoids any aliasing or wrap around errors)
# LVTOT  = .TRUE.      (Write total electrostatic potential into LOCPOT or not)
# LVHAR  = .TRUE.      (Write ionic + Hartree electrostatic potential into LOCPOT or not)
# NELECT =             (No. of electrons: charged cells, be careful)
# LPLANE = .TRUE.      (Real space distribution, supercells)
# NWRITE = 2           (Medium-level output)
# KPAR   = 2           (Divides k-grid into separate groups)
# NGXF    = 300        (FFT grid mesh density for nice charge/potential plots)
# NGYF    = 300        (FFT grid mesh density for nice charge/potential plots)
# NGZF    = 300        (FFT grid mesh density for nice charge/potential plots)

Electronic Relaxation
ISMEAR =  0
SIGMA  =  0.05
EDIFF  =  1E-08

LEPSILON = .TRUE.      (Determined by density functional perturbation theory)
# LCALCEPS = .TRUE.    (OR Determined by self-consistent response of the system to a finite electric field)
LPEAD    = .TRUE.      (Determined derivative of the cell-periodic part of the orbitals using finite differences)
""",

        "PH": """###ISMEAR =  0            (Gaussian smearing)
SIGMA  =  0.01         (Smearing value in eV)
EDIFF  =  1E-08        (SCF energy convergence, in eV)
PREC   =  Accurate     (Precision level)
# ENCUT  =  500        (Cut-off energy for plane wave basis set, in eV)
IALGO  =  38           (Davidson block iteration scheme)
LREAL  = .FALSE.       (Projection operators: false)
LWAVE  = .FALSE.       (Write WAVECAR or not)
LCHARG = .FALSE.       (Write CHGCAR or not)
ADDGRID= .TRUE.        (Increase grid, helps GGA convergence)
IBRION =  6            (Calculate phonon frequencies using finite differences approach)
LPHON_DISPERSION = T   (Calculate phonon dispersion along the q-point path supplied in file QPOINTS)
PHON_NWRITE = -3       (Both phonon eigenvectors and phonon frequencies are written for each q point)
""",

        "NE": """###Global Parameters
ISTART =  1            (Read existing wavefunction, if there)
ISPIN  =  1            (Non-Spin polarised DFT)
# ICHARG =  11         (Non-self-consistent: GGA/LDA band structures)
LREAL  = .FALSE.       (Projection operators: automatic)
# ENCUT  =  400        (Cut-off energy for plane wave basis set, in eV)
# PREC   =  Accurate   (Precision level: Normal or Accurate, set Accurate when perform structure lattice relaxation calculation)
LWAVE  = .TRUE.        (Write WAVECAR or not)
LCHARG = .TRUE.        (Write CHGCAR or not)
ADDGRID= .TRUE.        (Increase grid, helps GGA convergence)
LASPH  = .TRUE.        (Give more accurate total energies and band structure calculations)
PREC   = Accurate      (Accurate strictly avoids any aliasing or wrap around errors)
# LVTOT  = .TRUE.      (Write total electrostatic potential into LOCPOT or not)
# LVHAR  = .TRUE.      (Write ionic + Hartree electrostatic potential into LOCPOT or not)
# NELECT =             (No. of electrons: charged cells, be careful)
# LPLANE = .TRUE.      (Real space distribution, supercells)
# NWRITE = 2           (Medium-level output)
# KPAR   = 2           (Divides k-grid into separate groups)
# NGXF    = 300        (FFT grid mesh density for nice charge/potential plots)
# NGYF    = 300        (FFT grid mesh density for nice charge/potential plots)
# NGZF    = 300        (FFT grid mesh density for nice charge/potential plots)

Budged Elastic Band (NEB)
IMAGES =  5            (no. of images excluding two endpoints, set a NEB run no. of nodes, must be dividable by no. of images. each group of nodes works on one image.   PLEASE CHANG
E IT TO YOUR CHOICE NEBMAKE
.PL)
NSW    =  500          (number of ionic steps)
ISMEAR =  0            (gaussian smearing method)
SIGMA  =  0.05         (please check the width of the smearing)
IBRION =  3            (do MD with a zero time step)
POTIM  =  0            (Zero time step so that VASP does not move the ions)
SPRING =  -5.0         (spring force (eV/A2) between images)
LCLIMB =  .TRUE.       (turn on the climbing image algorithm)
ICHAIN =  0            (Indicates which method to run. NEB (ICHAIN=0) is the default)
IOPT   =  1            (LBFGS = Limited-memory Broyden-Fletcher-Goldfarb-Shanno)
""",

        "FQ": """###Global Parameters
ISTART =  1            (Read existing wavefunction, if there)
ISPIN  =  1            (Non-Spin polarised DFT)
# ICHARG =  11         (Non-self-consistent: GGA/LDA band structures)
LREAL  = .FALSE.       (Projection operators: automatic)
# ENCUT  =  400        (Cut-off energy for plane wave basis set, in eV)
# PREC   =  Accurate   (Precision level: Normal or Accurate, set Accurate when perform structure lattice relaxation calculation)
LWAVE  = .TRUE.        (Write WAVECAR or not)
LCHARG = .TRUE.        (Write CHGCAR or not)
ADDGRID= .TRUE.        (Increase grid, helps GGA convergence)
LASPH  = .TRUE.        (Give more accurate total energies and band structure calculations)
PREC   = Accurate      (Accurate strictly avoids any aliasing or wrap around errors)
# LVTOT  = .TRUE.      (Write total electrostatic potential into LOCPOT or not)
# LVHAR  = .TRUE.      (Write ionic + Hartree electrostatic potential into LOCPOT or not)
# NELECT =             (No. of electrons: charged cells, be careful)
# LPLANE = .TRUE.      (Real space distribution, supercells)
# NWRITE = 2           (Medium-level output)
# KPAR   = 2           (Divides k-grid into separate groups)
# NGXF    = 300        (FFT grid mesh density for nice charge/potential plots)
# NGYF    = 300        (FFT grid mesh density for nice charge/potential plots)
# NGZF    = 300        (FFT grid mesh density for nice charge/potential plots)

Frequence Calculations
NSW    =  1            (number of ionic steps. Make it odd.)
ISMEAR =  0            (gaussian smearing method)
SIGMA  =  0.05         (please check the width of the smearing)
IBRION =  5            (frequence calculation)
POTIM  =  0.02         (displacement step)
NFREE  =  2            (displacement freedom)
""",

        "LR": """###Global Parameters
ISTART =  1            (Read existing wavefunction, if there)
ISPIN  =  1            (Non-Spin polarised DFT)
# ICHARG =  11         (Non-self-consistent: GGA/LDA band structures)
LREAL  = .FALSE.       (Projection operators: automatic)
# ENCUT  =  400        (Cut-off energy for plane wave basis set, in eV)
# PREC   =  Accurate   (Precision level: Normal or Accurate, set Accurate when perform structure lattice relaxation calculation)
LWAVE  = .TRUE.        (Write WAVECAR or not)
LCHARG = .TRUE.        (Write CHGCAR or not)
ADDGRID= .TRUE.        (Increase grid, helps GGA convergence)
LASPH  = .TRUE.        (Give more accurate total energies and band structure calculations)
PREC   = Accurate      (Accurate strictly avoids any aliasing or wrap around errors)
# LVTOT  = .TRUE.      (Write total electrostatic potential into LOCPOT or not)
# LVHAR  = .TRUE.      (Write ionic + Hartree electrostatic potential into LOCPOT or not)
# NELECT =             (No. of electrons: charged cells, be careful)
# LPLANE = .TRUE.      (Real space distribution, supercells)
# NWRITE = 2           (Medium-level output)
# KPAR   = 2           (Divides k-grid into separate groups)
# NGXF    = 300        (FFT grid mesh density for nice charge/potential plots)
# NGYF    = 300        (FFT grid mesh density for nice charge/potential plots)
# NGZF    = 300        (FFT grid mesh density for nice charge/potential plots)

Lattice Relaxation
NSW    =  300          (number of ionic steps)
ISMEAR =  0            (gaussian smearing method )
SIGMA  =  0.05         (please check the width of the smearing)
IBRION =  2            (Algorithm: 0-MD, 1-Quasi-New, 2-CG)
ISIF   =  3            (optimize atomic coordinates and lattice parameters)
EDIFFG = -1.5E-02      (Ionic convergence, eV/A)
""",

        "MT": """###Global Parameters
ISTART =  1            (Read existing wavefunction, if there)
ISPIN  =  1            (Non-Spin polarised DFT)
# ICHARG =  11         (Non-self-consistent: GGA/LDA band structures)
LREAL  = .FALSE.       (Projection operators: automatic)
# ENCUT  =  400        (Cut-off energy for plane wave basis set, in eV)
# PREC   =  Accurate   (Precision level: Normal or Accurate, set Accurate when perform structure lattice relaxation calculation)
LWAVE  = .TRUE.        (Write WAVECAR or not)
LCHARG = .TRUE.        (Write CHGCAR or not)
ADDGRID= .TRUE.        (Increase grid, helps GGA convergence)
LASPH  = .TRUE.        (Give more accurate total energies and band structure calculations)
PREC   = Accurate      (Accurate strictly avoids any aliasing or wrap around errors)
# LVTOT  = .TRUE.      (Write total electrostatic potential into LOCPOT or not)
# LVHAR  = .TRUE.      (Write ionic + Hartree electrostatic potential into LOCPOT or not)
# NELECT =             (No. of electrons: charged cells, be careful)
# LPLANE = .TRUE.      (Real space distribution, supercells)
# NWRITE = 2           (Medium-level output)
# KPAR   = 2           (Divides k-grid into separate groups)
# NGXF    = 300        (FFT grid mesh density for nice charge/potential plots)
# NGYF    = 300        (FFT grid mesh density for nice charge/potential plots)
# NGZF    = 300        (FFT grid mesh density for nice charge/potential plots)

Meta-GGA Calculation
METAGGA = MBJ       (SCAN | RTPSS | MBJ | LIBXC |)
LASPH   = .TRUE.    (include non-spherical contributions in the PAW spheres)
LMAXMIX = 4         (increase for d/f-electrons to obtain fast convergence)
ALGO   =  ALL       (Electronic Minimisation Algorithm, ALGO=58)
TIME   =  0.1       (Timestep for IALGO5X)
""",

        "SR": """###Global Parameters
ISTART =  1            (Read existing wavefunction, if there)
ISPIN  =  1            (Non-Spin polarised DFT)
# ICHARG =  11         (Non-self-consistent: GGA/LDA band structures)
LREAL  = .FALSE.       (Projection operators: automatic)
# ENCUT  =  400        (Cut-off energy for plane wave basis set, in eV)
# PREC   =  Accurate   (Precision level: Normal or Accurate, set Accurate when perform structure lattice relaxation calculation)
LWAVE  = .TRUE.        (Write WAVECAR or not)
LCHARG = .TRUE.        (Write CHGCAR or not)
ADDGRID= .TRUE.        (Increase grid, helps GGA convergence)
LASPH  = .TRUE.        (Give more accurate total energies and band structure calculations)
PREC   = Accurate      (Accurate strictly avoids any aliasing or wrap around errors)
# LVTOT  = .TRUE.      (Write total electrostatic potential into LOCPOT or not)
# LVHAR  = .TRUE.      (Write ionic + Hartree electrostatic potential into LOCPOT or not)
# NELECT =             (No. of electrons: charged cells, be careful)
# LPLANE = .TRUE.      (Real space distribution, supercells)
# NWRITE = 2           (Medium-level output)
# KPAR   = 2           (Divides k-grid into separate groups)
# NGXF    = 300        (FFT grid mesh density for nice charge/potential plots)
# NGYF    = 300        (FFT grid mesh density for nice charge/potential plots)
# NGZF    = 300        (FFT grid mesh density for nice charge/potential plots)

Electronic Relaxation
ISMEAR =  0            (Gaussian smearing, metals:1)
SIGMA  =  0.05         (Smearing value in eV, metals:0.2)
NELM   =  90           (Max electronic SCF steps)
NELMIN =  6            (Min electronic SCF steps)
EDIFF  =  1E-08        (SCF energy convergence, in eV)
# GGA  =  PS           (PBEsol exchange-correlation)

Ionic Relaxation
NSW    =  100          (Max ionic steps)
IBRION =  2            (Algorithm: 0-MD, 1-Quasi-New, 2-CG)
ISIF   =  2            (Stress/relaxation: 2-Ions, 3-Shape/Ions/V, 4-Shape/Ions)
EDIFFG = -2E-02        (Ionic convergence, eV/AA)
# ISYM =  2            (Symmetry: 0=none, 2=GGA, 3=hybrids)
""",

        "SO": """###Spin-Orbit Coupling Calculation
LSORBIT    = .TRUE.    (Activate SOC)
GGA_COMPAT = .FALSE.   (Apply spherical cutoff on gradient field)
VOSKOWN    =  1        (Enhances the magnetic moments and the magnetic energies)
LMAXMIX    =  4        (For d elements increase LMAXMIX to 4, f: LMAXMIX = 6)
ISYM       =  -1       (Switch symmetry off)
# SAXIS    =  0 0 1    (Direction of the magnetic field)
# MAGMOM   =  0 0 3    (Local magnetic moment parallel to SAXIS)
# NBANDS   =           (2 * number of bands of collinear-run)
""",

        "H6": """###HSE06 Calculation
LHFCALC= .TRUE.       (Activate HF)
AEXX   =  0.25        (25% HF exact exchange, adjusted to reproduce experimental band gap)
HFSCREEN= 0.2         (Switch to screened exchange, e.g. HSE06)
ALGO   =  ALL         (Electronic Minimisation Algorithm, ALGO=58)
TIME   =  0.4         (Timestep for IALGO5X)
PRECFOCK= N           (HF FFT grid)
# NKRED    = 2         (Reduce k-grid-even only)
# HFLMAX   = 4         (HF cut-off: 4d, 6f)
# LDIAG    = .TRUE.    (Diagnolise Eigenvalues)
""",

        "MD": """###Global Parameters
ISTART =  1            (Read existing wavefunction, if there)
ISPIN  =  1            (Non-Spin polarised DFT)
# ICHARG =  11         (Non-self-consistent: GGA/LDA band structures)
LREAL  = .FALSE.       (Projection operators: automatic)
# ENCUT  =  400        (Cut-off energy for plane wave basis set, in eV)
# PREC   =  Accurate   (Precision level: Normal or Accurate, set Accurate when perform structure lattice relaxation calculation)
LWAVE  = .TRUE.        (Write WAVECAR or not)
LCHARG = .TRUE.        (Write CHGCAR or not)
ADDGRID= .TRUE.        (Increase grid, helps GGA convergence)
LASPH  = .TRUE.        (Give more accurate total energies and band structure calculations)
PREC   = Accurate      (Accurate strictly avoids any aliasing or wrap around errors)
# LVTOT  = .TRUE.      (Write total electrostatic potential into LOCPOT or not)
# LVHAR  = .TRUE.      (Write ionic + Hartree electrostatic potential into LOCPOT or not)
# NELECT =             (No. of electrons: charged cells, be careful)
# LPLANE = .TRUE.      (Real space distribution, supercells)
# NWRITE = 2           (Medium-level output)
# KPAR   = 2           (Divides k-grid into separate groups)
# NGXF    = 300        (FFT grid mesh density for nice charge/potential plots)
# NGYF    = 300        (FFT grid mesh density for nice charge/potential plots)
# NGZF    = 300        (FFT grid mesh density for nice charge/potential plots)

Electronic Relaxation
ISMEAR =  0
SIGMA  =  0.05
EDIFF  =  1E-08

Molecular Dynamics
IBRION =  0            (Activate MD)
NSW    =  100          (Max ionic steps)
EDIFFG = -1E-02        (Ionic convergence, eV/A)
POTIM  =  1            (Timestep in fs)
SMASS  =  0            (MD Algorithm: -3=microcanonical, 0=canonical)
# TEBEG  =     100     (Start temperature K)
# TEEND  =     100     (Final temperature K)
# MDALGO =  1          (Andersen Thermostat)
# ISYM   =  0          (Switch symmetry off)
NWRITE =  0            (For long MD-runs use NWRITE=0 or 1)
""",

        "BS": """###Global Parameters
ISTART =  1            (Read existing wavefunction, if there)
ISPIN  =  1            (Non-Spin polarised DFT)
# ICHARG =  11         (Non-self-consistent: GGA/LDA band structures)
LREAL  = .FALSE.       (Projection operators: automatic)
# ENCUT  =  400        (Cut-off energy for plane wave basis set, in eV)
# PREC   =  Accurate   (Precision level: Normal or Accurate, set Accurate when perform structure lattice relaxation calculation)
LWAVE  = .TRUE.        (Write WAVECAR or not)
LCHARG = .TRUE.        (Write CHGCAR or not)
ADDGRID= .TRUE.        (Increase grid, helps GGA convergence)
LASPH  = .TRUE.        (Give more accurate total energies and band structure calculations)
PREC   = Accurate      (Accurate strictly avoids any aliasing or wrap around errors)
# LVTOT  = .TRUE.      (Write total electrostatic potential into LOCPOT or not)
# LVHAR  = .TRUE.      (Write ionic + Hartree electrostatic potential into LOCPOT or not)
# NELECT =             (No. of electrons: charged cells, be careful)
# LPLANE = .TRUE.      (Real space distribution, supercells)
# NWRITE = 2           (Medium-level output)
# KPAR   = 2           (Divides k-grid into separate groups)
# NGXF    = 300        (FFT grid mesh density for nice charge/potential plots)
# NGYF    = 300        (FFT grid mesh density for nice charge/potential plots)
# NGZF    = 300        (FFT grid mesh density for nice charge/potential plots)

Electronic Relaxation
ISMEAR =  0
SIGMA  =  0.05
EDIFF  =  1E-08
""",

        "EL": """###Electron Localization Function
ISTART =  1            (Read existing wavefunction, if there)
LELF   = .TRUE.        (Activate ELF)
""",

        "OP": """###Optical properties
ALGO     =  Exact
NBANDS   =             (Set this parameters manually)
LOPTICS  = .TRUE.
CSHIFT   =  0.100
NEDOS    =  2000
ISMEAR   =  0
SIGMA    =  0.01
EDIFF    =  1.E-8
# LPEAD  = .TRUE.
""",

        "PC": """###Decomposed Charge Density
ISTART =      1        (Job: 0-new  1-cont  2-samecut)
ICHARG =      1        (Read charge: 1-file 2-atom 10-const)
LPARD  = .TRUE.        (Activate decomposed charge density)
LSEPB  = .TRUE.        (Separately write PARCHG.nb by every band)
LSEPK  = .TRUE.        (Separately write PARCHG.nk by every kpoint)

Method I: Partial Charge for the specified BANDS and KPOINTS
IBAND  = 20 21 22 23   (Set this parameters manually)
KPUSE  = 1 2 3 4       (Set this parameters manually)

# Method II: Partial Charge in the energy rang of [-10.3 -5.1]
# EINT = -10.3 -5.1

# Method III: Partial Charge in the energy rang of [EF-1 -EF]
# NBMOD=-3
# EINT = -1

*********** Notes *************
(1) Copy IBZKPT as KPOINTS for static calculation,
(2) Band structure calculation.
""",

        "PY": """###ISMEAR =  0            (Gaussian smearing)
SIGMA  =  0.05         (Smearing value in eV)
# VASP + phonopy with force constants using density-functional-perturbation
IBRION =  8            (Determines the Hessian matrix using DFPT)
# VASP + phonopy with force constants using finite displacement method
# IBRION =  -1         (The ions are not moved)
EDIFF  =  1E-08        (SCF energy convergence, in eV)
PREC   =  Accurate     (Precision level)
# ENCUT  =  500        (Cut-off energy for plane wave basis set, in eV)
IALGO  =  38           (Davidson block iteration scheme)
LREAL  = .FALSE.       (Projection operators: false)
LWAVE  = .FALSE.       (Write WAVECAR or not)
LCHARG = .FALSE.       (Write CHGCAR or not)
ADDGRID= .TRUE.        (Increase grid, helps GGA convergence)
""",

        "DM": """###Global Parameters
ISTART =  1            (Read existing wavefunction, if there)
ISPIN  =  1            (Non-Spin polarised DFT)
# ICHARG =  11         (Non-self-consistent: GGA/LDA band structures)
LREAL  = .FALSE.       (Projection operators: automatic)
# ENCUT  =  400        (Cut-off energy for plane wave basis set, in eV)
# PREC   =  Accurate   (Precision level: Normal or Accurate, set Accurate when perform structure lattice relaxation calculation)
LWAVE  = .TRUE.        (Write WAVECAR or not)
LCHARG = .TRUE.        (Write CHGCAR or not)
ADDGRID= .TRUE.        (Increase grid, helps GGA convergence)
LASPH  = .TRUE.        (Give more accurate total energies and band structure calculations)
PREC   = Accurate      (Accurate strictly avoids any aliasing or wrap around errors)
# LVTOT  = .TRUE.      (Write total electrostatic potential into LOCPOT or not)
# LVHAR  = .TRUE.      (Write ionic + Hartree electrostatic potential into LOCPOT or not)
# NELECT =             (No. of electrons: charged cells, be careful)
# LPLANE = .TRUE.      (Real space distribution, supercells)
# NWRITE = 2           (Medium-level output)
# KPAR   = 2           (Divides k-grid into separate groups)
# NGXF    = 300        (FFT grid mesh density for nice charge/potential plots)
# NGYF    = 300        (FFT grid mesh density for nice charge/potential plots)
# NGZF    = 300        (FFT grid mesh density for nice charge/potential plots)

The dimer method
NSW    =  500          (number of ionic steps)
ISMEAR =  0            (gaussian smearing method )
SIGMA  =  0.05         (please check the width of the smearing)
IBRION =  3            (do MD with a zero time step)
POTIM  =  0            (Zero time step so that VASP does not move the ions)
ICHAIN =  2            (Use the dimer method required for the latest code)
DdR    =  0.005        (The dimer separation)
DRotMax  =  1          (Maximum rotation steps per translation step)
DFNMin =  0.01         (Magnitude of the rotational force below which the dimer is not rotated)
DFNMax =  1.0          (Magnitude of the rotational force below which dimer rotation stops)
IOPT   =  2            (CG = Conjugate Gradient)
""",

        "PZ": """###Global Parameters
ISTART =  1            (Read existing wavefunction, if there)
ISPIN  =  1            (Non-Spin polarised DFT)
# ICHARG =  11         (Non-self-consistent: GGA/LDA band structures)
LREAL  = .FALSE.       (Projection operators: automatic)
# ENCUT  =  400        (Cut-off energy for plane wave basis set, in eV)
# PREC   =  Accurate   (Precision level: Normal or Accurate, set Accurate when perform structure lattice relaxation calculation)
LWAVE  = .TRUE.        (Write WAVECAR or not)
LCHARG = .TRUE.        (Write CHGCAR or not)
ADDGRID= .TRUE.        (Increase grid, helps GGA convergence)
LASPH  = .TRUE.        (Give more accurate total energies and band structure calculations)
PREC   = Accurate      (Accurate strictly avoids any aliasing or wrap around errors)
# LVTOT  = .TRUE.      (Write total electrostatic potential into LOCPOT or not)
# LVHAR  = .TRUE.      (Write ionic + Hartree electrostatic potential into LOCPOT or not)
# NELECT =             (No. of electrons: charged cells, be careful)
# LPLANE = .TRUE.      (Real space distribution, supercells)
# NWRITE = 2           (Medium-level output)
# KPAR   = 2           (Divides k-grid into separate groups)
# NGXF    = 300        (FFT grid mesh density for nice charge/potential plots)
# NGYF    = 300        (FFT grid mesh density for nice charge/potential plots)
# NGZF    = 300        (FFT grid mesh density for nice charge/potential plots)

Piezoelectric Calculation
IBRION   = 8       (second derivatives of total energy with respect to the position using DFPT)
LEPSILON = .TRUE.  (Determines the static ion-clamped dielectric matrix using DFPT)
LPEAD    = .TRUE.  (Derivative of the cell-periodic part of orbitals)
"""
    }

    calc_labels = {
            "ST": "ST) Static-Calculation",
            "SR": "SR) Standard Relaxation",
            "MG": "MG) Magnetic Properties",
            "SO": "SO) Spin-Orbit Coupling",
            "D3": "D3) DFT-D3 no-damping Correction",
            "H6": "H6) HSE06 Calculation",
            "PU": "PU) DFT+U Calculation",
            "MD": "MD) Molecular Dynamics",
            "GW": "GW) GW0 Calculation",
            "BS": "BS) BSE Calculation",
            "DC": "DC) Elastic Constant",
            "EL": "EL) ELF Calculation",
            "BD": "BD) Bader Charge Analysis",
            "OP": "OP) Optical Properties",
            "EC": "EC) Static Dielectric Constant",
            "PC": "PC) Decomposed Charge Density",
            "PH": "PH) Phonon-Calculation",
            "PY": "PY) Phonon with Phonopy",
            "NE": "NE) Nudged Elastic Band (NEB)",
            "DM": "DM) The Dimer Method",
            "FQ": "FQ) Frequence Calculation",
            "LR": "LR) Lattice Relaxation",
            "MT": "MT) Meta-GGA Calculation",
            "PZ": "PZ) Piezoelectric Calculation"
        }

    # ------------------------------------------------------------------
    # 3) The order to display them in the sidebar
    # ------------------------------------------------------------------
    calc_order = [
        "ST", "SR", "MG", "SO", "D3", "H6", "PU", "MD",
        "GW", "BS", "DC", "EL", "BD", "OP", "EC", "PC",
        "PH", "PY", "NE", "DM", "FQ", "LR", "MT", "PZ"
    ]
    st.markdown("Select one or more calculation blocks **on the left sidebar** to combine their parameters.")
    # ------------------------------------------------------------------
    # 4) Checkboxes in sidebar: each key is code, label is from calc_labels
    # ------------------------------------------------------------------
    st.sidebar.header("Calculation Types")
    selected_codes = []
    for code in calc_order:
        # Show checkbox only if code is in our dictionary
        if code in calc_blocks:
            label = calc_labels[code]
            if st.sidebar.checkbox(label, value=False):
                selected_codes.append(code)
    # ------------------------------------------------------------------
    # 5) Combine the selected blocks into final_incar text
    # ------------------------------------------------------------------
    final_incar = ""
    for code in selected_codes:
        block_text = calc_blocks[code].strip()
        final_incar += block_text + "\n\n"
    final_incar = final_incar.strip()
    # ------------------------------------------------------------------
    # 6) Display and download
    # ------------------------------------------------------------------
    if final_incar:
        st.subheader("Generated INCAR Content")
        st.code(final_incar, language="text")
        st.download_button(
            label="Download INCAR",
            data=final_incar,
            file_name="INCAR",
            mime="text/plain"
        )
    else:
        st.info("No blocks selected yet. Pick one or more from the left sidebar.")

if __name__ == "__main__":
    main()