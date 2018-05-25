# Protein folding
Proteins are strings of aminoacids. (Dis)fuctioning of proteins is known to be determined by the dimensional folding structure. In a simplified model, the aminoacids can be represented by H (hydrophobic) and P (polar). In this model the stability of a folding structure is measured by the amount of H-bonds present in the structure. An H-bond is defined as two Hs that are adjecent to each other in the structure, while not in the aminoacid-sequence. Each H-bonds has a score of -1. A more advanced protein can also have aminoacid represented by a C (cysteine). These kan form C-bonds with a stability score of -5 and C-H bonds with a stability score of -1. The lower the stability score of the total structure, the more stable the protein is and the better the folding structure is (see figure).

![best random folding pattern of protein1](https://github.com/Jara555/Proteins/blob/master/doc/randomizer_protein7.png)

This program implements algorithms, all with the aim to find the best possible folding structure for a benchmark of proteins as listed below. Proteins can be fold in 2D as well as 3D. Furthermore, besides these benchmark proteins, the program can be run for any custum created protein.
1) HHPHHHPH
2) HHPHHHPHPHHHPH 
3) HPHPPHHPHPPHPHHPPHPH 
4) PPPHHPPHHPPPPPHHHHHHHPPHHPPPPHHPPHPP 
5) HHPHPHPHPHHHHPHPPPHPPPHPPPPHPPPHPPPHPHHHHPHPHPHPHH 
6) PPCHHPPCHPPPPCHHHHCHHPPHHPPPPHHPPHPP
7) CPPCHPPCHPPCPPHHHHHHCCPCHPPCPCHPPHPC
8) HCPHPCPHPCHCHPHPPPHPPPHPPPPHPCPHPPPHPHHHCCHCHCHCHH
9) HCPHPHPHCHHHHPCCPPHPPPHPPPPCPPPHPPPHPHHHHCHPHPHPHH

## Getting Started

### Prerequisites
This code is written in [Python3.6.3](https://www.python.org/downloads/). 
Furthermore some packages are used. They are easy to install using the following instruction:

```
pip install -r requirements.txt
```

### Structure

This code uses three classes: AminoAcid, Protein and Algorithm. These can together with their corresponding methods be found in the map classes. The map algorithms contains al algorithm subclasses of the superclass Algorithm. All benchmark proteins can be found in the map data and the results will be saved as a .csv or .log files in the results map. The main script for running the program can be found in proteins.py.

### Algorithms
Proteins are folded according to the following folding algorithms:
- Randomizer: this algorithm chooses a random folding pattern and checks stability for a number of iterations. 
- Hill Climber: this algorithm starts with folding pattern extracted from the randomizer algorithm. Subsequently, it changes the direction from a random choosen amino acids and checks stability for a number of iterations.
- Simulated Annealing: this algorithm allows overlap and a degradation of stability for a maximal amount of times in order to escape local minima/maxima. Temperature is cooling down every run, as is the allowed amount of overlap and degradation.
- Depth First: this algorithm checks all folding patterns possible by a depth first search.
- Branch 'n Bound: this algorithm starts with folding pattern extracted from the randomizer algorithm. Subsequently, it checks folding patterns by a depth first search and excludes specific patterns based on the number of Hbonds already found.

### Running

To run the program with the standardconfigurations use the following instructions:

```
python proteins.py -a <algorithm> -p <protein> -d <dimensions> -i <iterations> -c <csv>

```

        :argument -a <algorithm> : First letters of algorithm names
                                        R = Randomizer
                                        HC = HillClimber
                                        SA = SimulatedAnnealing
                                        DF = DepthFirst
                                        BB = BranchNBound
        :argument -p <protein> : Protein to be fold
                                 Can be an integer reflecting 1 of the benchmark proteins 
                                 Or a custom string reflecting a list of amino acids
                                        int options = 1/2/3/ ... 9
                                        str options  = H/P/C
        :argument -d <dimensions> : dimensions to be fold in
                                        2 = 2D
                                        3 = 3D (default)
        :argument -i <iterations> : Maximal iterations to be run 
                                    Required for R and HC, optional for DF and BB
                                        default = None 
        :argument -c <csv> : Write solutions to .csv file
                                        ON = write
                                        OFF = not write (default)

### ExperimentAlgorithms

After each run a .log file is created containing the running info of the algorithm for the specific protein. To analyze the log results the program experiment.py can be used, found in the map /experiment. The experiment program will analyze the running info of all algorithms for which .log files exist.  

To run the experiment with the standardconfigurations use the following instructions:

```
python experimentAlgorithms.py -p <protein> -d <dimensions>

```
        :argument -p <protein> : Protein to be analyzed
                                 Integer reflecting 1 of the proteins in the /data folder
        :argument -d <dimensions> : dimension to be analyzed
                                        2 = 2D
                                        3 = 3D
                             
                                        
### ExperimetProteins

To run the experiment with the standardconfigurations use the following instructions:

```
python experimentProteins.py -p <protein> -d <dimensions>

```
        :argument -a <algorithm> : First letters of algorithm names
                R = Randomizer
                HC = HillClimber
                SA = SimulatedAnnealing
                DF = DepthFirst
                BB = BranchNBound (default)
        :argument -d <dimensions> : dimensions to be fold in
                2 = 2D
                3 = 3D (default)
        :argument -i <iterations> : Maximal iterations to be run (required for R and HC, optional for DF and BB)
        :argument -n <number> : Number of proteins to be created (default: 100)
        :argument -l <length> : Lenght of the protein to be crated (default: 14)
        :argument -f <fixed h-number> : Fixed number of H's in the protein

## Authors

* Amber Brands 
* Jara Linders
* Marrit Leenstra

## Acknowledgments

* Minor programmeren van de UvA
* Special thanks to our TA Bram van den Heuvel

