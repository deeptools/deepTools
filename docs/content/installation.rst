Installation
=============

Remember -- deepTools are available for **command line usage** as well as for
**integration into Galaxy servers**!

.. contents:: 
    :local:

Requirements
-------------

* Python 2.7 or Python 3.x
* numpy >= 1.8.0
* scipy >= 0.17.0
* py2bit >= 0.1.0
* pyBigWig >= 0.2.1
* pysam >= 0.8
* matplotlib >= 1.4.0

The fastest way to obtain **Python 2.7 or Python 3.x together with numpy and scipy** is
via the `Anaconda Scientific Python
Distribution <https://store.continuum.io/cshop/anaconda/>`_.
Just download the version that's suitable for your operating system and
follow the directions for its installation. All of the requirements for deepTools can be installed in Anaconda with:

.. code:: bash

    $ conda install -c bioconda deeptools

Command line installation using ``pip``
-----------------------------------------

Install deepTools using the following command:
::

	$ pip install deeptools

All python requirements should be automatically installed.

If you need to specify a specific path for the installation of the tools, make use of `pip install`'s numerous options:

.. code:: bash

    $ pip install --install-option="--prefix=/MyPath/Tools/deepTools2.0" git+https://github.com/deeptools/deepTools.git


Command line installation without ``pip``
-------------------------------------------

You are highly recommended to use `pip` rather than these more complicated steps.

1. Install the requirements listed above in the "requirements" section. This is done automatically by `pip`.

2. Download source code
::

	$ git clone https://github.com/deeptools/deepTools.git

or if you want a particular release, choose one from https://github.com/deeptools/deepTools/releases:
::

	$ wget https://github.com/deeptools/deepTools/archive/1.5.12.tar.gz
	$ tar -xzvf

3. install the source code (if you don't have root permission, you can set
a specific folder using the ``--prefix`` option)
::

	$ python setup.py install --prefix /User/Tools/deepTools2.0

Galaxy installation
--------------------

deepTools can be easily integrated into a local `Galaxy <http://galaxyproject.org>`_.
All wrappers and dependencies are available in the `Galaxy Tool
Shed <http://toolshed.g2.bx.psu.edu/view/bgruening/deeptools>`_.

Installation via Galaxy API (recommended)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

First generate an `API Key <http://wiki.galaxyproject.org/Admin/API#Generate_the_Admin_Account_API_Key>`_
for your admin user and run the the installation script:
::

	$ python ./scripts/api/install_tool_shed_repositories.py \
		--api YOUR_API_KEY -l http://localhost/ \
		--url http://toolshed.g2.bx.psu.edu/ \
		-o bgruening -r <revision> --name suite_deeptools \
		--tool-deps --repository-deps --panel-section-name deepTools

The ``-r`` argument specifies the version of deepTools. You can get the
latest revision number from the test tool shed or with the following
command:
::

	$ hg identify http://toolshed.g2.bx.psu.edu/repos/bgruening/suite_deeptools

You can watch the installation status under: Top Panel --> Admin --> Manage
installed tool shed repositories

Installation via web browser
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

-  go to the `admin page <http://localhost:8080/admin>`_
-  select *Search and browse tool sheds*
-  Galaxy tool shed --> Sequence Analysis --> deeptools
-  install deeptools

Installation with Docker
^^^^^^^^^^^^^^^^^^^^^^^^

The deepTools Galaxy instance is also available as a docker container, for those wishing to use the Galaxy framework but who also prefer a virtualized solution. This container is quite simple to install:
::

    $ sudo docker pull quay.io/bgruening/galaxy-deeptools

To start and otherwise modify this container, please see the instructions on `the docker-galaxy-stable github repository <https://github.com/bgruening/docker-galaxy-stable>`__. Note that you must use `bgruening/galaxy-deeptools` in place of `bgruening/galaxy-stable` in the examples, as the deepTools Galaxy container is built on top of the galaxy-stable container.

.. tip:: For support or questions please make a post on `Biostars <http://biostars.org>`__. For feature requests or bug reports please open an issue `on github <http://github.com/deeptools/deeptools>`__.
