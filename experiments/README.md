### Experiments

During this project two experiments were carried out. The first analysis the results of the different algorithms, whereas the second analyses the effects of different protein properties on the stability. The results are saved in respectively the /AlgorithmResults folder and the /ProteinResults folder. Furthermore the ProteinRandomizer class can be found in this folder, since it is only used in the Proteins experiment. Lastly the statespace.xls document shows the summerized results of all algorithms.

### AlgorithmResults
Each protein has been runned with different algorithms. The performance of the algorithms has been analysed with the Algorithm experiment. The results are visualised in graphs in this folder.

Title structure: Protein \<number\> - \<dimensions\> D.png

Example:

![algorithm results of protein 2](https://raw.githubusercontent.com/Jara555/Proteins/master/experiments/AlgorithmResults/Protein2-3D.png)
 
### ProteinResults
In the Protein experiment the effects of 3 different protein properties on the stability have been analysed: 
- the number of H's in the protein 
- the lenght of the longest cluster in the protein
- the number of clusters in the protein

Title structure: \<Property\> -l \<lenght\> -d \<dimensions\> -f \<fixed number of H's\>

Examples:
![HCount example](https://raw.githubusercontent.com/Jara555/Proteins/master/experiments/ProteinResults/HCount-l14-d2.png)
![Cluster example](https://raw.githubusercontent.com/Jara555/Proteins/master/experiments/ProteinResults/Cluster-l14-d2-f7.png)
