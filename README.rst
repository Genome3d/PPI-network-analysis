Louvain clustering and centrality analysis on protein-protein interaction(PPI) network 
======================================================================================

This repository contain python scripts
     
*  to detect statistically significant communities from PPI network.
*  to perform centrality analysis


Requirements
============

* python 3.0
* networkx 2.4
* qstest 1.1.0
* pandas 1.0.4

Usage
=====

For significant community detection
***********************************

.. code-block:: bash

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
