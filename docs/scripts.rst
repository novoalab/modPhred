Scripts
=======

mods_from_bams
--------------

If you happen to have BAM files with modifications already encoded in them
and just want to process them with modPhred
all you need to do is:
- make new directory ie ``mods_from_bams/CC.raw_1x``

- copy ``modPhred.pkl`` there (it stores info about mods that are present in BAMs)

.. code-block:: bash

   rsync -av modPhred/raw_1x.m6A/modPhred.pkl mods_from_bams/CC.raw_1x

- prepare your BAM files in this output dir - BAM files have to be in ``mods_from_bams/CC.raw_1x/minimap2``

- and run the script
  
.. code-block:: bash

   ~/src/modPhred/src/mods_from_bams.py -o mods_from_bams/CC.raw_1x \
-i mods_from_bams/CC.raw_1x/minimap2/*.bam -f ~/cluster/rna_mods/ref/curlcake.fa

