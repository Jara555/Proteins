# Protein folding
Proteins are strings of aminoacids. (Dis)fuctioning of proteins is known to be determined by the dimensional folding structure. In a simplified model, the aminoacids can be represented by H (hydrophobic) and P (polar). In this model the stability of a folding structure is measured by the amount of H-bonds present in the structure. An H-bond is defined as two Hs that are adjecent to each other in the structure, while not in the aminoacid-sequence. Each H-bonds has a score of -1. The lower the score of the total structure, the more stable the protein is and the better the folding structure is (see figure).

![best random folding pattern of protein1](https://github.com/Jara555/Proteins/blob/master/Figures/random_protein1.png)

This project tries to find an algorithm that finds the best possible folding structure for a benchmark of proteins:
* HHPHHHPH
* HHPHHHPHPHHHPH 
* HPHPPHHPHPPHPHHPPHPH 
* PPPHHPPHHPPPPPHHHHHHHPPHHPPPPHHPPHPP 
* HHPHPHPHPHHHHPHPPPHPPPHPPPPHPPPHPPPHPHHHHPHPHPHPHH 

## Getting Started

### Prerequisites
This code is written in [Python3.6.3](https://www.python.org/downloads/). 
Furthermore the code uses the following packages: matpotlib. This is easy to install using the following instruction:

```
pip3 install matplotlib
```

### Structure

This code used three two classes: Protein and AminoAcid. These can together with their corresponding methods be found in the map new_classes. All benchmark proteins can be found in the map data and the results are safed as a .csv file in the results map.

### Testing

To run the code with the standardconfiguration (e.g. protein HHPHHHPH) use the following instruction:

```
python application.py
```

## Authors

* Amber Brands 
* Jara Linders
* Marrit Leenstra

## Acknowledgments

* Minor programmeren van de UvA
* Special thanks to our TA Bram van den Heuvel

