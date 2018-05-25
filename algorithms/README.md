### Algorithms
The following folder contains the algorithms used for protein folding. All algorithms share properties and methods of the Algorithm super class. The algorithms have the following approaches:

### Randomizer
1. A random folding pattern is generated.
2. The algorithm computes the stability score and compares it with the maximal stability score.
3. If the stability score is lower than the max stability, save pattern and change max stability. Go back to (1).

### Hill Climber 
1. Starts with a folding pattern extracted from the randomizer algorithm. 
2. Selects an amino acid and alters its orientation in all possible directions (overlap is not allowed).
3. If the stability score is better or equal than the max stability, save pattern and change max stability. Go back to (2).

### Simulated Annealing
1. Starts with a folding pattern extracted from the randomizer algorithm.
2. Selects an amino acid and alters its orientation in all possible directions.
4. A probability is calculated for the acceptance of overlap and degradations in stability score. This calculation is based on a temperature which is cooling down during every iteration.
5. If the stability score is better or equal than the max stability, save pattern and change max stability. Go back to (2).

### Depth First 
1. This algorithm checks all folding patterns possible by a depth first search.
2. Pruning/exclusing overlapping patterns.

### Branch 'n Bound 
1. Starts with a folding pattern extracted from the randomizer algorithm.
2. Checks folding patterns by a depth first search.
4. Prunes overlapping patterns.
3. Prunes when more than 3 aminoacids are placed in a straight line.
4. Prunes when potential stability is lower than best stability already found.
