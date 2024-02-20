# Installation

## Prerequisites
To compile MQC, please make sure that the following prerequisites are present:

- Fortran compiler.  You can use [Intel® Fortran Compiler](https://www.intel.com/content/www/us/en/developer/tools/oneapi/fortran-compiler.html) or [GFortran](https://gcc.gnu.org/fortran/).
- [BLAS](http://www.netlib.org/blas/).
- [LAPACK](http://www.netlib.org/lapack/)
- [(optional)libcint](https://github.com/sunqm/libcint)
  
The package provides dynamic linked library (```libblas.dll``` and  ```liblapack.dll```) of  **LAPACK 3.11**  compiled by **gfortan 8.1.0**(MinGW) for Windows and static library(```libcint.a```,```librefblas.a``` and ```libreflapack.a```) compiled by **gfortan 11.3.0** for Linux in the ```lib``` folder.

## There are two ways to use MQC
### 1.Compile from source code
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
### 2.Directly use the precompiled binary file
You can get the precompiled binary file  ```mqc``` from github release.