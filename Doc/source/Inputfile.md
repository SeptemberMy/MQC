# Input files 
## Structure
In a nutshell,  **MQC**  can be controlled by keywords in input file.
The input file  should follow the following rules :
- All the keywords are either **all uppercase** or **all lowercase**.Mixing uppercase and lowercase is not accepted.
- Input file must start writing from the first line,end with `>END`.
The content after `>END` will not be read and can be used as comments.
- The first line of input file is considered a `title` and does not work.
- The control related to the system is usually written immediately before the title, starting with `!`.
- Keywords related to computation are given in the form of modules(blocks),starting with `#`,end with `END`.
- Although spaces generally have little impact, it is still recommended not to add too many spaces.
- Except for modules TASK and MOLECULE, all other modules are generally in the form of `[keyword] [values]`.

The structure of input file has the following forms :

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
Modules(Blocks) and keywords can be found in the next section.
## Modules(Blocks) 
Currently, there are modules available in **MQC**.

| Blocks name | outline | optional |
| :---: | :-----------: | :--------: | 
|TASK|Provide basic calculation information|No|
|MOLECULE|Provide basic molecular structure information|No|
|BASIS|Provide basis set information|No|
|INT| Control the calculation of general molecular integrals|Yes|
|SCF| Control the SCF procedure|Yes|
|POSTHF| Control the post-Hatree-fock procedure|Yes|
|OUTPUT| Control program output content|Yes|
|BASPATH| Provide the basis set file path|Yes|
## Key words
### TASK block
Keywords can be written in the same line without order.
`[task] [method]`
|Key words of task |*Describe* | 
| :---: | :--------: | 
| energy |*Calculate **Single Point Energy***|

|Key words of method |*Describe*| 
| :---: | :--------: | 
|hf|*Hartree-Fock method*|
|mp2|*Second-order Møller–Plesset perturbation theory*|
|mp3|*Third-order Møller–Plesset perturbation theory*|
|ccsd|*Coupled cluster singles and doubles*|
|ccsd(t)|*CCSD with perturbative triples correction*|
### MOLECULE block
MOLECULE block should follow the given forms.

**NOTE**:
-  The value of  *Charge* and  *Spin-multiplicity* should be interger in the same line.
-  *Unit name* is `bohr` or `angstrom`.
-  *The number of atoms* should be same as the number of following line except `END` .
-  Atomic number `Z` must be interger.
### BASIS block
For convenience, **MQC** did not define its own basis set file.

Basis set file of **MQC** using **Gaussian type** can be easily obtained from [Basis set exchange](https://www.basissetexchange.org/).
And it can be easily customized.

**NOTE**:
 - If you want to use user-defined basis set,please follow the Gasussian type and add  ****  line above the first element.
 - The key word is the name of basis set file which can specify path, see **BASPATH block**.
eq:
```c++
****
H     0
S    3   1.00
      0.3425250914D+01       0.1543289673D+00
      0.6239137298D+00       0.5353281423D+00
      0.1688554040D+00       0.4446345422D+00
****
He     0
S    3   1.00
      0.6362421394D+01       0.1543289673D+00
      0.1158922999D+01       0.5353281423D+00
      0.3136497915D+00       0.4446345422D+00
****
...

```
The keywords for available basis sets in MQC :
#### Pople's family basis sets
 |keyword | availability |
 | :---: | :--------: | 
 |sto-2g |  H-Xe| 
 |sto-3g |  H-Xe| 
 |sto-4g |  H-Xe|
 |sto-5g |  H-Xe|
 |sto-6g |  H-Xe|
 |3-21g  |  H-Xe|
 |4-31g  |  H-He,B-Ne,P-Cl|
 |5-21g  |  H-Be|
 |6-21g  |  H-Ar|
 |6-31g  |  H-Kr|
 |6-31g(2df,p)  |  H-Ar| 
 |6-311g |  H-Ca,Ga-Kr,I| 
 |6-311gss| H-Ca,Ga-Kr,I| 
 |6-311g(2df,2pd)  |  H-Ne,K-Ca| 
 |6-311++g  |  H-Ca| 
#### Dunning's family basis sets(Correlation-consistent basis sets)
 |keyword | availability |
 | :---: | :--------: | 
 |cc-pvdz | H-Ar,Ca-Kr| 
 |cc-pvtz | H-Ar,Ca-Kr| 
 |cc-pvqz | H-Ar,Ca-Kr|
 |cc-pv5z | H-Ar,Ca-Kr|
 |cc-pv6z | H-He,Be-Ne,Al-Ar|
 |cc-pv8z | H,Ne|
 |cc-pv9z | Ne|
 |aug-cc-pvdz | H-Ar,Sc-Kr| 
 |aug-cc-pvtz | H-Ar,Sc-Kr| 
 |aug-cc-pvqz | H-Ar,Sc-Kr| 
 #### Karlsruhe’s family basis sets(The def/def2 basis sets)
 |keyword | availability |
 | :---: | :--------: | 
 |def2-sv(p) | H-Rn| 
 |def2-svp | H-Rn| 
 |def2-svpd | H-La,Hf-Rn|
 |def2-tzvp | H-Rn| 
 |def2-tzvpp | H-Rn|
 |def2-tzvpd | H-La,Hf-Rn|
 |def2-tzvppd | H-La,Hf-Rn|
 |def2-qzvp | H-Rn| 
 |def2-qzvpp | H-Rn| 
 |def2-qzvpd | H-La,Hf-Rn| 
 |def2-qzvppd | H-La,Hf-Rn|
### INT block
 |keyword | default|**options**|
 | :---: | :--------: | ---|
 |METHODS | 1 |**1**:Gauss-Hermite numerical integration <br>**2**:Analytical algorithm|
 |METHODT | 1|**1**:Direct method <br>**2**:Integral by parts |
 |METHODV | 2|**1**:Analytical algorithm <br>**2**:Dupuis-Rys-King Method <br>**3**:Obara-Saika Recursive Method|
 |METHODERI |4|**1**:Dupuis-Rys-King method <br>**2**:Obara-Saika Recursive method <br>**3**:Ohata method<br>**4**:Dupuis-Rys-King method(fast)|
### SCF block
|keyword | default|*Describe* and **options**|
 | :---: | :--------: | ---|
| damp | 0 |*Set damp factor* <br> $D_{n+1}=damp*D_{n}+(1-damp)*D_{n+1}$ <br>**range: 0~1.0**|
| maxcycle | 200 |*The maximum number of SCF cycles* <br> **N is a positive integer.**|
| convere | 12 |*The SCF energy convergence threshold value* <br> **Set threshold to $10^{-N}$** |
| converrms | 12 |*The SCF RMS Density Matrix change convergence threshold value* <br> **Set threshold to $10^{-N}$**  |
| diis | 1 |*Pulay's Direct Inversion in the Iterative Subspace (DIIS) extrapolation method* <br> **0:off,1:on** |
| ndiis | 7 | *The size of DIIS space* <br>**N is a positive integer.**|
| guess | 1 |*The initial guess for the SCF* <br>**1:core Hamiltonian,2:Generalized Wolfsberg−Helmholz (GWH)** |
### POSTHF block
|keyword | default|**options**|
 | :---: | :--------: | ---|
| damp | 0 |*Set damp factor* <br> $D_{n+1}=damp*D_{n}+(1-damp)*D_{n+1}$ <br>**range: 0~1.0** |
| maxcycle | 200 |*The maximum number of CC Iterative process* <br> **N is a positive integer.**|
| convere | 12 |*The CC Iterative process energy convergence threshold value* <br> **Set threshold to $10^{-N}$** |
| converrms | 12 |*The CC Iterative process RMS Density Matrix change convergence threshold value* <br> **Set threshold to $10^{-N}$**|
| diis | 1 |DIIS for CC Iterative process <br> **0:off,1:on**|
| ndiis | 7 | *The size of DIIS space* <br>**N is a positive integer.**|
### OUTPUT block
The document number used in MQC is shown in the table below.
 |File id |name| content | type |
| :---:| :---:|:---: |:---: | 
| 33  | INPUT.inp |input file| formatted |
| 99  | OUTPUT.log |output file| formatted |
| 55  |~ |basis set file | formatted |
| 56  |~|aux basis set file | formatted |
| 66  |Structure.xyz|Storing the structure(xyz form)  | formatted |
| 70  |2eint.bin| two-electron intgral | unformatted |
### BASPATH block
  This block is used to specify the file paths of the basis set and auxiliary basis set, and can be used to customize the basis set.
  |keyword | availability |
 | :---: | :--------: | 
 |basfpath | Absolute path or relative path| 




