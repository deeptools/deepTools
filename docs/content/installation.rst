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
* samtools (preferably in your PATH)
* 1 UCSC tool: bedGraphToBigWig (preferably in your PATH) 

The fastet way to obtain **Python 2.7 together with numpy and scipy** is
via the `Anaconda Scientific Python
Distribution <https://store.continuum.io/cshop/anaconda/>`_.
Just download the version that's suitable for your operating system and
follow the directions for its installation.

Installing samtools
^^^^^^^^^^^^^^^^^^^^

1. Download source code and unzip it
::

	$ wget http://sourceforge.net/projects/samtools/files/samtools/1.2/samtools-1.2.tar.bz2/download -O samtools-1.2.tar.bz2
	$ bunzip2 samtools-1.2.tar.bz2 
	$ tar -xvf samtools-1.2.tar
	$ cd samtools-1.2

2. Compile
::

	$ make

3. Check whether the tool is running
::

	$ ./samtools

4. Add the location where it was installed to your PATH variable
::

	$ export PATH=$PATH:/Users/frd2007/Tools/samtools-1.2

Installing UCSC tools (aka Kent's tools)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

1. Figure out which operating system version you have
::

	$ uname -a

2. Download the compiled binaries from the corresponding folder (shown here for the Linux server,
check the website for the correct folder for your operating system) and make them executable.
::

	$ mkdir UCSCtools
	$ cd UCSCtools
	$ wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bedGraphToBigWig
	$ chmod +x bedGraphToBigWig

4. Check whether the tools work (you should get a usage message)
::

	$ ./bedGraphToBigWig

5. Add them to the PATH variable
::

	$ export PATH=$PATH:/Users/frd2007/Tools/UCSCtools


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
	samtools: samtools
	bedgraph_to_bigwig: bedGraphToBigWig
	
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

As you can see, deepTools expects samtools to be available via the command ``samtools``,
as well as the UCSC tool ``bedGraphToBigWig``.
You can either specify the path where you installed the tools in the *.cfg file
or add the tools to your PATH in your .bashrc file
(see above for details of the installation of these tools).

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
