Read clustering
===============
You can also cluster reads by modification status. 
Note, this currently cluster entire reference sequences, 
so it's reasonable to run it only with individual transcritps.
Definitely, avoid clustering reads from entire chromosomes. 

.. code-block:: bash

   ~/src/modPhred/src/mod_cluster.py --minfreq 0.01 -i mod.gz -e pdf


This will produce plots similar to this one

.. image:: cluster.png
   :align: center
   :width: 100%
	   
