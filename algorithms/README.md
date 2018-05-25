### Algorithms
The following folder contains the algorithms used for folding proteins. The proteins can be folded according to the following algorithms:

### Randomizer
1. A random folding pattern is generated.
2. The algorithm computes the stability score and compares it with the maximal stability score.
3. If the stability score is lower than the max stability, save pattern and change max stability, else go to (1)

### Hill Climber 
1. Starts with a folding pattern extracted from the randomizer algorithm. 
2. It selects a random amino acid and alters its direction in all possible directions (overlap is not allowed).
3. If the stability score is lower or eqeual than the max stability, save pattern and change max stability, else go to (1)

### Hill Climber with Degradation
1. Starts with a folding pattern extracted from the randomizer algorithm.
2. It selects a random amino acid and alters its direction randomly in all poissible directions. Overlap is allowed depending on how many iterations already have passed.
3. If the stability score is lower or eqeual than the max stability, save pattern and change max stability, else go to (1)

### Depth First 
This algorithm checks all folding patterns possible by a depth first search.

### Branch 'n Bound 
This algorithm starts with folding pattern extracted from the randomizer algorithm. Subsequently, it checks folding patterns by a depth first search and excludes specific patterns based on the number of Hbonds already found.
