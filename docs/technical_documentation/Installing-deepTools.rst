`**WIKI-START** <Home>`__ > `**deepTools technical
documentation** <Technical-documentation>`__ > `deepTools
installation <installing-deepTools>`__

Dependencies
------------

The following software is required for proper installation and usage of
deepTools:

-  Python 2.7
-  numpy
-  scipy
-  pysam (>= 0.7.7)

| 
| Installation
| ---------------

-  `Installing deepTools within a local *Galaxy instance* <#galaxy>`__
-  `Installing deepTools with a *Galaxy instance* using
   *docker* <#dockergalaxy>`__
-  `Installing deepTools from source (Linux) <#linux>`__
-  `Installing deepTools on a Mac <#mac>`__
-  `Troubleshooting <#trouble>`__

--------------

Galaxy Installation
^^^^^^^^^^^^^^^^^^^

| deepTools can be easily integrated into the
`Galaxy <http://galaxyproject.org>`__ framework. All wrappers and
dependencies are
| available in the `Galaxy Tool
Shed <http://toolshed.g2.bx.psu.edu/view/bgruening/deeptools>`__.

Installation via Galaxy API (recommended)
'''''''''''''''''''''''''''''''''''''''''

| First, generate an `API
Key <http://wiki.galaxyproject.org/Admin/API#Generate_the_Admin_Account_API_Key>`__
for your admin
| user and run the the installation script:

::

    python ./scripts/api/install_tool_shed_repositories.py --api YOUR_API_KEY -l http://localhost:8080 --url http://toolshed.g2.bx.psu.edu/ -o bgruening -r <revision> --name deeptools --tool-deps --repository-deps --panel-section-name deepTools

The -r argument specifies the version of deepTools. You can get the
latest revision number from the test tool shed or via the command line
with the following command:

::

    hg identify http://toolshed.g2.bx.psu.edu/view/bgruening/deeptools

You can watch the installation status within Galaxy: Top Panel → Admin →
Manage installed tool shed repositories

Installation via web browser
''''''''''''''''''''''''''''

-  go to the `admin page <http://localhost:8080/admin>`__
-  select *Search and browse tool sheds*
-  Galaxy tool shed → Sequence Analysis → deeptools
-  install deeptools

Using deepTools and Galaxy from a docker Image
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

At first you need to install `docker <http://docker.io/>`__. Please
follow the instruction on
https://www.docker.io/gettingstarted/#h_installation

After the successful installation, all what you need to do is:

``docker run -d -p 8080:80 bgruening/galaxy-deeptools``

I will shortly explain the meaning of all the parameters. For a more
detailed describtion please consult the `docker
manual <http://docs.docker.io/>`__, it's really worth reading.

Let's start: ``docker run`` will run the Image/Container for you. In
case you do not have the Container stored locally, docker will download
it for you. ``-p 8080:80`` will make the port 80 (inside of the
container) available on port 8080 on your host. Inside the container a
Apache Webserver is running on port 80 and that port can be bound to a
local port on your host computer. With this parameter you can access
your Galaxy instance via ``http://localhost:8080`` immediately after
executing the command above. ``bgruening/galaxy-deeptools`` is the
Image/Container name, that directs docker to the correct path in the
`docker index <https://index.docker.io/u/bgruening/galaxy-stable/>`__.
``-d`` will start the docker container in daemon mode.

For more information please refer to the `project
site <https://github.com/bgruening/docker-recipes/tree/master/galaxy-deeptools>`__.

Installation from source (Linux, command line)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The easiest way to install deepTools is by **downloading the source file
and using python pip** or easy\_install tools:

Requirements: Python 2.7, numpy, scipy installed

Commands:

::

      $ cd ~
      $ export PYTHONPATH=$PYTHONPATH:~/lib/python2.7/site-packages
      $ export PATH=$PATH:~/bin:~/.local/bin

If pip is not already available, install it with:

::

      $ easy_install --prefix=~ pip

Install deepTools and dependencies with pip:

::

      $ pip install --user deeptools

Done!

**Another option is to clone the repository:**

::

    $ git clone https://github.com/fidelram/deepTools

Then go to the deepTools directory, edit the ``deepTools.cfg`` file with
your favorite editor (e.g. vim) and run the install script:

::

    $ cd deepTools/
    $ vim config/deepTools.cfg
    $ python setup.py install

| By default, the script will install python library and executable
| codes globally, which means you need to be root or administrator of
| the machine to complete the installation. If you need to
| provide a non-standard install prefix, or any other non-standard
| options, you can provide many command line options to the install
| script.

::

    $ python setup.py --help

For example, to install deepTools under a specific location use:

::

    $ python setup.py install --prefix <target directory>

Installation on a MAC
~~~~~~~~~~~~~~~~~~~~~

| Although the installation of deepTools itself is quite simple,
| the installation of the required modules SciPy and NumPy might demand
| a bit of extra work.

| The easiest way to install them is together with the
| `Anaconda Scientific Python
Distribution <https://store.continuum.io/cshop/anaconda/>`__. After
installation, open a terminal ("Applications" --> "Terminal") and type:

::

     $ pip install deeptools

Done!

| If you prefer to install the dependencies individually, follow
| these steps (Python 2.7 must already be installed):

Download the packages and install them using dmg images:

-  http://sourceforge.net/projects/numpy/files/NumPy/
-  http://sourceforge.net/projects/scipy/files/scipy/

Then install deepTools via the terminal ("Applications" --> "Terminal"):

::

     $ cd ~
     $ export PYTHONPATH=$PYTHONPATH:~/lib/python2.7/site-packages
     $ export PATH=$PATH:~/bin:~/.local/bin:~/Library/Python/2.7/bin

If pip is not already available, install with:

::

     $ easy_install --prefix=~ pip

Install deepTools and dependencies with pip:

::

     $ pip install --user deeptools

Troubleshooting
'''''''''''''''

| The easy\_install command is provided by the python package
setuptools.
| You can download the package from
`https://pypi.python.org/pypi/setuptools <https://pypi.python.org/pypi/setuptools>`__

::

     $ wget https://bitbucket.org/pypa/setuptools/raw/bootstrap/ez_setup.py -O - | python
     

or the user-specific way:

::

     $ wget https://bitbucket.org/pypa/setuptools/raw/bootstrap/ez_setup.py
     $ python ez_setup.py --user

| Numpy/Scipy Installation:
| see
`http://www.scipy.org/install.html <http://www.scipy.org/install.html>`__

for support, questions, or feature requests contact:
deeptools@googlegroups.com
