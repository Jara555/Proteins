### Classes
The following folder contains the main classes used for protein folding.

### AminoAcid
This class contains all properties and methods of a amino acid. AminoAcid objects have a type: H(ydrophobic), P(olair) or (C)ysteine and x, y, z coordinates. A list of AminoAcids forms a protein.

### Protein
This class contains all properties and methods of a protein. A protein is implemented as being a list of AminoAcid objects. The methods in the Protein class are related to its properties, for example folding by changing the coordinates or calculating its stability.

### Algorithm
This class contains all properties and methods which are shared by different algorithms. The class is a super class of the algorithm subclasses, which can be found in the algorithms folder. This class manages for example the writing of .log and .csv result files and the printing of progress in the terminal.  

