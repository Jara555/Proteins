## 
The following folder can contain .csv and .log files. The structure of the two will be as follows

 ### log files
 Title: P <protein number> - <dimensions> D - <algorithm>.log
 The files contain the summarized data of an algorithm:
 
 ```
 _________________________________________________



SimulatedAnnealing

_________________________________________________



Protein 3: HPHPPHHPHPPHPHHPPHPH



Best Stability Found: -7.0



Best Folding Pattern Found:

"['0', '+Y', '+Z', '+Z', '-Y', '-Z', '+X', '+Z', '-Y', '-Z', '-Z', '+Y', '-Z', '-X', '-Z', '+Y', '-Z', '-Y', '-Y', '+Z']"



Started with stability of:  -3.0

Total Iterations:  99960

Total Overlap:     41420

Found in Run:      0

Elapsed Time:      22.6794

```

### csv files Algorithm Experiment
The files contain more detailed information in 3 columns inclusing:
- iteration number 
- stability 
- folding pattern per iteration

```
 2,-4.0,"['0', '+Y', '+X', '+X', '+Y', '+Y', '+X', '+X', '+X', '+Y', '+X', '-Y', '-Y', '-Y', '-X', '-X', '-X', '+Y', '+X', '+X']"
 ```
 ### csv files Protein Experiment
 Title: experimentProteins -l <length> -d <dimensions> -f <number of fixed H's>.csv
 
 The files contain 5 columns:
  - Number of H's
  - Max cluster length
  - Number of clusters
  - Upperbound stability 
  - Stability found with Branch 'n Bound algorithm
 
```
7,2,5,-6,-5.0
```
