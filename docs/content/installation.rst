Installation
=============

Remember -- deepTools are available for **command line usage** as well as for
**integration into Galaxy servers** !

.. contents:: 
    :local:

Command line installation using ``conda``
-----------------------------------------

The recommended way to install deepTools (including its requirements) is via `miniconda <https://docs.conda.io/projects/miniconda/en/latest/>`_ or `anaconda <https://www.anaconda.com/>`_.

.. code:: bash

    $ conda install -c bioconda deeptools

Command line installation using ``pip``
---------------------------------------

deepTools can also be installed using `pip <https://pip.pypa.io/en/stable/>`_.
You can either install the latest release from `pypi <https://pypi.org/>`_:

.. code:: bash

	$ pip install deeptools

or a specific version with:

.. code:: bash

	$ pip install deeptools==3.5.3

In case you would like to install an unreleased or development version, deepTools can also be installed from the repository:

.. code:: bash

	$ git clone https://github.com/deeptools/deepTools.git
	$ cd deepTools
	$ pip install .

Galaxy installation
--------------------

deepTools can be easily integrated into a local `Galaxy <https://galaxyproject.org>`_.
All wrappers and dependencies are available in the `Galaxy Tool
Shed <https://toolshed.g2.bx.psu.edu/>`_.

Installation via Galaxy API (recommended)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

First generate an `API Key <https://wiki.galaxyproject.org/Admin/API#Generate_the_Admin_Account_API_Key>`_
for your admin user and run the the installation script:
::

	$ python ./scripts/api/install_tool_shed_repositories.py \
		--api YOUR_API_KEY -l http://localhost/ \
		--url https://toolshed.g2.bx.psu.edu/ \
		-o bgruening -r <revision> --name suite_deeptools \
		--tool-deps --repository-deps --panel-section-name deepTools

The ``-r`` argument specifies the version of deepTools. You can get the
latest revision number from the test tool shed or with the following
command:
::

	$ hg identify https://toolshed.g2.bx.psu.edu/repos/bgruening/suite_deeptools

You can watch the installation status under: Top Panel --> Admin --> Manage
installed tool shed repositories

Installation via web browser
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

-  go to the `admin page <http://localhost:8080/admin>`_
-  select *Search and browse tool sheds*
-  Galaxy tool shed --> Sequence Analysis --> deeptools
-  install deeptools
