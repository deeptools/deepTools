Installation
=============

Remember -- deepTools are available for **command line usage** as well as for
**integration into Galaxy servers**!

.. contents:: 
    :local:

Requirements
-------------

* Python 2.7
* numpy, scipy, bx-python, and pyBigWig
* pysam >= 0.8

The fastet way to obtain **Python 2.7 together with numpy and scipy** is
via the `Anaconda Scientific Python
Distribution <https://store.continuum.io/cshop/anaconda/>`_.
Just download the version that's suitable for your operating system and
follow the directions for its installation. All of the requirements for deepTools can be installed in Anaconda with:

    $ conda install -c https://conda.anaconda.org/bioconda pysam pyBigWig bx-python

Command line installation using ``pip``
-----------------------------------------

Install deepTools using the following command:
::

	$ pip install git+https://github.com/fidelram/deepTools.git

All python requirements are automatically installed.


Command line installation without ``pip``
-------------------------------------------

1. Download source code
::

	$ git clone https://github.com/fidelram/deepTools.git

or if you want a particular release, choose one from https://github.com/fidelram/deepTools/releases:
::

	$ wget https://github.com/fidelram/deepTools/archive/1.5.12.tar.gz
	$ tar -xzvf

2. The config file will tell you what deepTools expects to be installed properly:
::

	$ cat deepTools/deeptools/config/deeptools.cfg
	
	[external_tools]
	sort: sort
	
	[general]
	# if set to max/2 (no quotes around)
	# half the available processors will
	# be used
	default_proc_number: max/2
	test_root: ../deeptools/test/

	# temporary dir:
	# deepTools bamCoverage, bamCompare and correctGCbias
	# write files to a temporary dir before merging them
	# and creating a final file. This can be speed up
	# by writting to /dev/shm but for this a large
	# physical memory of the server is required. If
	# this is the case in your system, uncomment
	# the following line. Otherwise, setting the
	# variable to 'default', deepTools will use the
	# temporary file configured in the system.
	# Any other path that wants to be used for temporary
	# files can by given as well (ie, /tmp)
	#tmp_dir: /dev/shm
	tmp_dir: default

3. install the source code (if you don't have root permission, you can set
a specific folder using the ``--prefix`` option)
::

	$ python setup.py install --prefix /Users/frd2007/Tools/deepTools

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
		--api YOUR_API_KEY -l http://localhost:8080 \
		--url http://toolshed.g2.bx.psu.edu/ \
		-o bgruening -r <revision> --name deeptools \
		--tool-deps --repository-deps --panel-section-name deepTools

The ``-r`` argument specifies the version of deepTools. You can get the
latest revsion number from the test tool shed or with the following
command:
::

	$ hg identify http://toolshed.g2.bx.psu.edu/view/bgruening/deeptools

You can watch the installation status under: Top Panel --> Admin --> Manage
installed tool shed repositories

Installation via web browser
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

-  go to the `admin page <http://localhost:8080/admin>`_
-  select *Search and browse tool sheds*
-  Galaxy tool shed --> Sequence Analysis --> deeptools
-  install deeptools

remember: for support, questions, or feature requests contact:
deeptools@googlegroups.com
