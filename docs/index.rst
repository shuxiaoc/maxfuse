MaxFuse Documentation
==================================

``MaxFuse`` (``Ma``\ tching ``x``\ cross modalities via ``Fu``\ zzy ``s``\ moothed ``e``\ mbeddings)
is an algorithmic pipeline for integrating multi-model single-cell data.

MaxFuse leverages all features available in each modality of interest
to strengthen the weak linkage across them via a sequence of
coembedding-smoothing-matching iterations,
and it can provide high-quality final matching and integration results
under low signal-to-noise-ratio regimes.

For more details about MaxFuse,
please check out `our paper <https://google.com/>`__.
See `this paper <https://google.com/>`__ for an application of MaxFuse
on human intestine data from the `HuBMAP <https://hubmapconsortium.org/>`__
consortium.

***************
Getting started
***************

The ``MaxFuse`` package can be installed via conda:

.. code-block:: bash
    :linenos:

    conda install -c conda-forge maxfuse

Or, it can also be installed via pip:

.. code-block:: bash
    :linenos:

    pip install maxfuse

.. note::
    To avoid potential dependency conflicts, we recommend
    installing within Python virtual environment such as conda.

Now you are all set! Please proceed to `tutorials <tutorials.rst>`__
for a list of examples.


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
