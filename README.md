# Ligand-Receptor-Pipeline 
First put your Seurat obj in the InputObj folder

Update the config csv with :

1.	name of Seurat clusters obj

3.	vector of signal receiving cluster separated by space

5.	vector of signal sending cluster separated by space

7.	Output Plot Path

9.	Project Name

the Graphs obj will be crated in the outputObj folder

Now you can run the FLOW analysis with:

pyton3 main.py

In order to run the TF analysis first add the cluster annotation to the Annot.csv file.

Now run -> python3 tf_network.py




