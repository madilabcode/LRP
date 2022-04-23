# Ligand-Receptor-Pipeline 
First put your Seurat obj in the InputObj folder

Update the config csv with :

  1. name of seurat clusterd obj 

  2. vector of signal reciving cluster seprated by space 

  3. vector of signal sending  cluster seprated by space

  4. Output Plot Path

  6. Project Name  

the Graphs obj will be crated in the outputObj folder 

Now you can run the FLOW anaylsis with:

pyton3 main.py 

In order to run the TF anaylsis first add the cluster annotation to the Annot.csv file.

Now run -> python3 tf_network.py 




