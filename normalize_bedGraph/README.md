
## normalize_bedGraph

This script allows a user to normalize a bedGraph file by any desired value. This normalization factor is appended to the new file name for documentation purposes. 

You can run the *normalize_bedGraph* script as follows:
```
perl normalize_bedGraph [factor] [multiply|divide] [bedGraph file 1] [bedGraph file 2] ... [bedGraph file N]
```

Here is a specific example:
 
```
perl normalize_bedGraph 0.85867 multiply Dm_S2_PRO-seq_F.bedGraph Dm_S2_PRO-seq_R.bedGraph
```
This generates two files named Dm_S2_PRO-seq_F_normM0.86.bedGraph and Dm_S2_PRO-seq_R_normM0.86.bedGraph.
