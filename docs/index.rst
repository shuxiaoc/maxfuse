MaxFuse Documentation
==================================

``MaxFuse`` (``Ma``\ tching ``x``\ cross modalities via ``Fu``\ zzy ``s``\ moothed ``e``\ mbeddings) 
is a python package for integrating single-cell datasets from different modalities with no overlapping features and/or under low signal-to-noise ratio regimes. 
For most single-cell cross modality integration methods, the feasibility of cross-modal integration relies on the existence of highly correlated, a priori 'linked' features. 
When such linked features are few or uninformative, a scenario that we call 'weak linkage', existing methods fail. 
We developed MaxFuse, a cross-modal data integration method that, through iterative co-embedding, data smoothing, and cell matching, leverages all information in each modality to obtain high-quality integration. 
A prototypical example of weak linkage is the integration of spatial proteomic data with single-cell sequencing data. 
For details, please refer to `the manuscript <https://www.biorxiv.org/content/10.1101/2023.01.12.523851v2>`__.


***************
Getting started
***************

The ``MaxFuse`` package can also be installed via pip:

.. code-block:: bash
    :linenos:

    conda create -n maxfuse python=3.8
    conda activate maxfuse
    python -m pip install maxfuse

.. note::
    To avoid potential dependency conflicts, we recommend
    installing within Python virtual environment such as conda.

Now you are all set! Please proceed to `tutorials <tutorials.rst>`__
for a list of examples.

Note in cases when integrating single cell data across protein and RNA modalities, 
many times the nomenclature of features are different (e.g., mRNA ITGAM could be named as CD11b-1 when used as antibody). 
We gathered a `.csv <https://github.com/shuxiaoc/maxfuse/blob/main/docs/protein_gene_conversion.csv>`__ file that covers many of such naming conversions and used during the MaxFuse process. 
Of course, this is not a complete conversion, and users should manually add in new naming conversions if they were not included in this .csv file.


***************
Code archive
***************
The analysis presented in `the manuscript <https://www.biorxiv.org/content/10.1101/2023.01.12.523851v2>`__ was also 
deposited in `this <https://github.com/shuxiaoc/maxfuse>`__ GitHub repository, under `this <https://github.com/shuxiaoc/maxfuse/tree/main/Archive>`__ folder. 
Note in the manuscript we used a development version of MaxFuse with slightly different grammar and can also be found there. 
If you require additional information on the analysis/data, please contact Zongming Ma (zongming@wharton.upenn.edu).


***************
License
***************
MaxFuse is under the `Academic Software License Agreement <https://github.com/shuxiaoc/maxfuse/blob/main/LICENSE>`__, please use accordingly.



.. toctree::
   :maxdepth: 2
   :caption: Contents

   tutorials
   api

******************
Indices and tables
******************

* :ref:`genindex`
* :ref:`modindex`
