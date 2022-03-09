.. UniOnco documentation master file, created by
   sphinx-quickstart on Fri Mar  4 01:43:10 2022.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

UniOnco: Unified Cancer Ontology
================================
Team: Jackson Michuda, Abby Pandya, Vinodh Rajapakse, Joseph Wakim



Quickstart
----------

1. **Download Source Code:**

.. parsed-literal::

   $ git clone https://github.com/jmichuda/bmi-214-final-project.git

2. **Install Poetry:**

.. parsed-literal::

   $ curl -sSL https://raw.githubusercontent.com/python-poetry/poetry/master/get-poetry.py | python -

3. **Install Environment:**

.. parsed-literal::

   $ cd bmi-214-final-project
   $ poetry install

4. **Generate Ontology by Running Build Pipeline:**

.. parsed-literal::

   $ poetry run python src/ontology.py


.. Tip:: See **Demos** for usage of API tools and query functionality


.. toctree::
   :hidden:
   :maxdepth: 2
   :caption: Introduction

   Abstract <abstract>

   Background <background>


.. toctree::
   :hidden:
   :maxdepth: 2
   :caption: Demos

   CIViC API Usage <CIViC_API_Usage.ipynb>

   cBioPortal API Usage <cBioPortal_API_Usage.ipynb>

   Query Ontology <query_therapy_regimen.ipynb>


.. toctree::
   :hidden:
   :maxdepth: 4
   :caption: Additional Resources

   Modules <modules>

   References <references>



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
