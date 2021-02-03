
## normalize_bedGraph

You can run the *normalize_bedGraph* script as follows:
```
perl normalize_bedGraph [factor] [multiply|divide] [bedGraph file 1] [bedGraph file 2] ... [bedGraph file N]
```

Here is a specific example:
 
```
perl normalize_bedGraph 1.3222567 multiply Dm_S2_Start-seq_Mock-tr_A_5prRNA_forward.bedGraph Dm_S2_Start-seq_Mock-tr_A_5prRNA_reverse.bedGraph
```
This generates two files named Dm_S2_Start-seq_Mock-tr_A_5prRNA_forward_normM1.32.bedGraph and Dm_S2_Start-seq_Mock-tr_A_5prRNA_reverse_normM1.32.bedGraph.
