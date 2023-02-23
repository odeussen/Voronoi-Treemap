# Voronoi-Treemap
Processing code to create a Voronoi Treemap from a CSV-file. Output is to a PDF-File and an SVG-File. 
In contrast to other code here the lines between cells are curves. Other code typically creates power 
diagrams, which are much faster to compute, but have straight lines.

Please install processing from processing.org for running and editing the program. After installation 
open the editor and select "voronoi_treemap.pde" as the main file. 

If you don't want to install processing, running versions for mac, windows, linux are in the zip files

An example of an input CSV-file is data/flare-2.csv, a standard file for visualizations of hierachical data.
Important: the file has to start with "name,size" to indicate the pairs in each line
Here the first lines of the CSV-file: 
```
name,size
flare,
flare.analytics,
flare.analytics.cluster,
flare.analytics.cluster.AgglomerativeCluster,3938
flare.analytics.cluster.CommunityStructure,3812
flare.analytics.cluster.HierarchicalCluster,6714
flare.analytics.cluster.MergeEdge,743
flare.analytics.graph,
...
```
The file descrbes a tree structure with leaves having a size. Sizes are relative in the Voronoi Treemap, 
so the sizes of all leaves belonging to an internal node (in case of describing a file system, a directory) 
are summed up for the directory. The whole tree is shown as partition of a circle.

The output of the program for data/flare-2.csv are in data/flare-2.pdf and data/flare-2.svg 

![data/flare-2.svg](/data/flare-2.svg)
