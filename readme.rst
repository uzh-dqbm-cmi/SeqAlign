SeqAlign: Sequence Alignment
============================

Sequence alignment done right -- Modules implementing sequence alignment algorithms correctly.. :). See this `paper by Flouri et al. <http://www.biorxiv.org/content/biorxiv/early/2015/11/12/031500.full.pdf>`__
for further discussion about the correctness and completeness of sequence alignment algorithms.

Moreover, the implemented alignment algorithms not only operate on sequences composed of character elements generated from a given alphabet, 
but also on generic sequence elements (i.e. nodes) where each element/node has some feature vector representation. 
In fact, elements could have arbitrary representation, but what is important is to define/implement a similarity measure function that is able to compute elements similarity in the sequences being aligned.

-------------------------------------------

TODO
-----

	#. Provide examples for non character elements in sequence alignment
	#. Provide examples for character sequence alignment (for now check the tests module under ``src`` folder)

-------------------------------------

Non temporal Aligner
--------------------


Supported alignment methods
++++++++++++++++++++++++++++

	#. ``global`` alignment (i.e. Needlman-Wunch or more precisely a variant of David Sankoff algorithm). 
	   See  Sankoff, D. (1972). "Matching sequences under deletion-insertion constraints". Proceedings of the National Academy of Sciences of the United States of America. 69 (1): 4–6
	#. ``local`` alignment (i.e. Smith-Waterman algorithm)
	#. ``semi-global`` alignment (shorter to longer sequences)
	#. ``end-gap-free`` alignment (i.e. detecting overlap)
 
Supported alignment penalties
++++++++++++++++++++++++++++++

The current implementation allows for specifying:

	- linear gap penalty
	- affine gap penalty
	
-------------------------------------


Temporal Aligner
-----------------

Supported alignment methods
++++++++++++++++++++++++++++

	#. ``global`` alignment (i.e. Needlman-Wunch or more precisely a variant of David Sankoff algorithm). 
	   See  Sankoff, D. (1972). "Matching sequences under deletion-insertion constraints". Proceedings of the National Academy of Sciences of the United States of America. 69 (1): 4–6
	#. ``local`` alignment (i.e. Smith-Waterman algorithm)
	#. ``semi-global`` alignment (shorter to longer sequences)
	#. ``end-gap-free`` alignment (i.e. detecting overlap)
 
 
Supported alignment penalties
+++++++++++++++++++++++++++++

The current implementation allows for specifying:

	- linear gap penalty

-------------------------------------

``NB``: All implemented algorithms return alignment score with all the possible alignment paths having the optimal score.

