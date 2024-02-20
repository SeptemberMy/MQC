# Overview

**MQC** is based on **Fortran90**. To maintain the readability,code is not parallel in principle,and should be consistent with formulas as much as possible. As a result,the efficiency may be low,but will remain within an acceptable range.

**MQC** can  run on **Linux** and **Windows** platforms and used some linear algebra libraries,like **BLAS/LAPACK** ,so you need to configure these environments before using it.
# Package Features
- **Intgral**
    - *Overlap integrals **(S)***
    - *Kinetic integrals **(T)***
    - *Nuclear attraction integrals **(NAI)***
    - *Electron-electron repulsion integrals **(ERI)***
- **Single Point Energy**
  - *Restricted-Hartree-Fock **(RHF)***
  - *Møller–Plesset perturbation theory **(MP2,MP2.5,MP3)*** 
  - *Coupled cluster singles and doubles **(CCSD)*** 
  - *CCSD with perturbative triples correction **CCSD(T)***
- **Wavefunction Analysis and Property**
  - *Mulliken charges*
  - *Löwdin charges*
  - *dipole moment*
