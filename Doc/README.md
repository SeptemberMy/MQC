# Document
## How to generate document
**It is easy to generate document of MQC.**
The **Doc** folder is a *sphnix* project.

Several types document can be automatically generated through ```Make.bat``` (Windows)and```Makefile```(Linux).

Just use 
```bash
make html
```
or 
```bash
make latex
```
then you can find html or $Latex$ files in **build** folder.

You also can clean the **build** folder use
```bash
make clean
```
If you want to get PDF,Please compile tex files using texlive(xelatex).

For more information, refer to the the [documentation of *sphnix*](https://www.sphinx-doc.org/en/master/).




