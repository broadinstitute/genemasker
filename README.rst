genemasker
**********

This application produces masked group/gene variant inclusion files for rare variant analysis.

Installation
************

This application has been tested using the following version of Python and required modules.
   - python=3.9.21
   - numpy=2.0.2
   - pandas=2.2.3
   - scikit-learn=1.6.1
   - scipy=1.13.1
   - psutil=6.1.1
   - dask=2024.8.0

It can be installed in a conda environment via the provided yml file as below.

.. code-block:: bash

   conda env create -f env/environment.yml
   source activate genemasker

With your environment activated, you are ready to install genemasker.

.. code-block:: bash

   python setup.py install

Another option is to run genemasker via Docker (see the included Dockerfile).

Getting Started
***************

.. code-block:: bash

   source activate genemasker

Verify that genemasker is functional using the following command to display help.

.. code-block:: bash

   genemasker --help
