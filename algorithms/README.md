### Algorithms
The following folder contains the algorithms used for folding proteins. The proteins can be folded according to the following algorithms:

### Randomizer
1. A random folding pattern is generated.
2. The algorithm computes the stability score and compares it with the maximal stability score.
3. If the stability score is lower than the max stability, save pattern and change max stability, else go to (1)

### Hill Climber 
1. The algorithm starts with a folding pattern extracted from the randomizer algorithm. 
2. It selects a random amino acid and alters it direction randomly.
3. If the stability score is lower than the max stability, save pattern and change max stability, else go to (1)

### Simulated Annealing 
This algorithm allows overlap and a degradation of stability for a maximal amount of times in order to escape local minima/maxima. Temperature is cooling down every run, as is the allowed amount of overlap and degradation.

### Depth First 
This algorithm checks all folding patterns possible by a depth first search.

### Branch 'n Bound 
This algorithm starts with folding pattern extracted from the randomizer algorithm. Subsequently, it checks folding patterns by a depth first search and excludes specific patterns based on the number of Hbonds already found.
