<p align="center">
    <img src="Doc/source/picture/logo.svg">
</p>

# My Quantum Chemistry(version:demo)
## Overview
**My Quantum Chemistry(MQC)** is an open-source package based on   Ab-initio Quantum Chemistry(mainly) methods.This package aims to provide a simple and convenient platform for quickly implementing new theoretical  methods, which also can be a teaching tool (toy model) for beginners to learn Quantum Chemistry.

**MQC** is based on **Fortran90**. To maintain the readability,code is not parallel in principle,and should be consistent with formulas as much as possible. As a result,the efficiency may be low,but will remain within an acceptable range.


## Installation
**MQC** can  run on Linux and Windows platforms and used some linear algebra libraries,like [**BLAS**](https://www.netlib.org/blas/) [**LAPACK**](https://www.netlib.org/lapack/) ,so you need to configure these environments before using it.

The package provides dynamic linked library (```libblas.dll``` and  ```liblapack.dll```) of  **LAPACK 3.11**  compiled by **gfortan 8.1.0**(MinGW) for Windows and static library(```librefblas.a``` and ```libreflapack.a```) compiled by **gfortan 11.3.0** for Linux in the ```lapack``` folder.

### There are two ways to use MQC
#### 1.Compile from source code
You need fortran compiler,both **gfortran** and **ifort** are acceptable.
Using **Makefile** can easily compile binary files,but you should make sure that the LAPACK or MKL was added to environment variables.
You can also specify the library directory by yourself.

To install the package,you need to first enter the directory where the source code is located and execute the following command in the terminal:

**for Windows**
```bash
make
```
**for Linux**
```bash
make mqclinux
```
then you can get the binary file ```mqc```.
You can use the following command to clean the auxiliary files:

**for Windows**
```bash
make clean
```
**for Linux**
```bash
make cleanlinux
```
Done!
#### 2.Directly use the precompiled binary file
You can get the precompiled binary file  ```mqc``` from github release.
## Quickstart
The running mode of mqc is ```mqc``` +input file path,default to the same path as **mqc** binary file.

For example,  ```mqc .\INPUT.inp``` is the same as ```mqc```.


You need to prepare your input file ```INPUT.inp``` following this type firstly.

The structure of input file  has the following forms :

```
INPUT title
! System
# Block name 1
  keyword11 option11
  keyword12 option12
  ...
END  
# Block name 2
  keyword21 option21
  ...
END  
...
# MOLECULE  
  Charge   Spin-multiplicity
  Unit name
  Num of atoms
  Z(atomic number which must be interger)         Coordinates( X       Y       Z)
END
...
Other blocks
...
> END
```

All the keywords are either **all uppercase** or **all lowercase** .

Consider that you want to calculate the electronic energy of a $H_2O$ molecule. 
```c++
Hello water !
! WINDOWS
# TASK
  ENERGY  HF
END
# BASIS
  STO-3G
END
# MOLECULE 
  0 1
  BOHR
  3
  8   0.0000000000      0.0000000000      0.0000000000
  1   0.0000000000      0.0000000000      2.0579117598
  1   1.9923609397      0.0000000000     -0.5152602865
END
>END
```
Keep the input file```INPUT.inp```  and the binary file  ```mqc``` in the same directory,then run ```mqc```.
Done!

For more specific content, please refer to the **document**.
## Package Features
- **Integral**
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
#### 
## Document
**It is easy to generate document of MQC.**
The **Doc** folder is a [*sphnix* ](https://www.sphinx-doc.org/en/master/)project.
You can use the following command to generate **html** or **$latex$**  project of document.
```bash
make html
```
```bash
make latex
```
You can use ```make cleandoc``` to clean up the document files.
See [how to generate document](./Doc/README.md).
## Reference
[1] [*Programming Projects* from Crawford Group](https://github.com/CrawfordGroup/ProgrammingProjects) 
[2] Q. Sun. J. Comput. Chem. 2015, 36, 1664–1671. DOI: 10.1002/jcc.23981
