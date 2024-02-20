# Quick Example
* How to run MQC

The running mode of mqc is ```mqc + input file path + output file path ```.
The default path is the same as **mqc** binary file.

For example,  ```mqc INPUT.inp``` is the same as ```mqc```.

**Note that if you want to specify the output file path,you must specify the input file path at the same time.**

**The default path of basis set files is ```./basis```,which can be customized by input file.**

For example,  ```mqc INPUT.inp OUTPUT.log```.The file names can be customized.

You need to prepare your input file ```INPUT.inp``` following this type firstly.

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