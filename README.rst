
Clustering protein-protein interaction(PPI) network using Louvain community detection algorithm
===============================================================================================


This repository contain python scripts used -
     
*  to detect statistically significant communities from PPI network.
*  to perform centrality analysis


Requirements
============

* python 3.0
* networkX 2.4
* qstest 1.1.0
* pandas 1.0.4

Usage
=====

For significant community detection
-----------------------------------

::

     python python/find__significant_module.py -n example/example_network.txt -g example/example_input.txt -o output_dir

For details::

     python python/find__significant_module.py -h

For centrality analysis
-----------------------
::

     python python/perform_centrality_analysis.py -n example/example_network.txt -g example/example_input.txt -o output_dir

For details::

     python python/perform_centrality_analysis.py -h

References
==========

.. [#] Blondel, V. D., Guillaume, J. L., Lambiotte, R. & Lefebvre, E. Fast unfolding of communities in large networks. J. Stat. Mech. Theory Exp. 2008, (2008).
.. [#] Kojaku, S. & Masuda, N. A generalised significance test for individual communities in networks. Sci. Rep. 8, 1–10 (2018).
