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

    $ conda install -c conda-forge -c bioconda deeptools

Note that for ARM architecture (e.g. M1 on OSX) you could go via the pip installation (see below), or install via the osx-64 env:

.. code:: bash

    $ CONDA_SUBDIR=osx-64 conda create -c conda-forge -c bioconda -n deeptools deeptools


Command line installation using ``pip``
---------------------------------------

deepTools can also be installed using `pip <https://pip.pypa.io/en/stable/>`_.
You can either install the latest release from `pypi <https://pypi.org/>`_:

.. code:: bash

	$ pip install deeptools

or a specific version with:

.. code:: bash

	$ pip install deeptools==3.5.6

In case you would like to install an unreleased or development version, deepTools can also be installed from the repository:

.. code:: bash

	$ git clone https://github.com/deeptools/deepTools.git
	$ cd deepTools
	$ pip install .

Command line installation using ``uv``
--------------------------------------

`uv <https://docs.astral.sh/uv/>`_ can install deepTools as a standalone command line
application, fully isolated from your other Python environments:

.. code:: bash

    $ uv tool install deeptools

The ``deeptools`` commands (``bamCoverage``, ``computeMatrix``, ...) are then available on
your ``PATH`` without activating any environment. Upgrade or remove with
``uv tool upgrade deeptools`` / ``uv tool uninstall deeptools``.

To run deepTools once without installing it, use:

.. code:: bash

    $ uvx --from deeptools bamCoverage --help

Inside a project you can instead add it as a dependency with ``uv add deeptools`` and run
the tools via ``uv run bamCoverage ...``.

Command line installation using ``pipx``
----------------------------------------

`pipx <https://pipx.pypa.io/>`_ installs deepTools into its own isolated environment while
exposing the command line tools globally:

.. code:: bash

    $ pipx install deeptools

Run a single tool without a persistent install using:

.. code:: bash

    $ pipx run --spec deeptools bamCoverage --help

Upgrade with ``pipx upgrade deeptools`` and remove with ``pipx uninstall deeptools``.

Command line installation using ``pixi``
----------------------------------------

`pixi <https://pixi.sh/>`_ can install deepTools globally from the ``bioconda`` channel,
which ships prebuilt binaries (no compiler needed):

.. code:: bash

    $ pixi global install -c conda-forge -c bioconda deeptools

To hack on deepTools from a checkout, the repository ships a pixi workspace that pins the
build toolchain (Rust, ``maturin``, ``libclang``, ``htslib``, ``pkg-config``, ``openssl`` --
see ``pyproject.toml``'s ``[tool.pixi.dependencies]``). This builds the Rust extension from
source and installs deepTools in editable mode:

.. code:: bash

    $ git clone https://github.com/deeptools/deepTools.git
    $ cd deepTools
    $ pixi install
    $ pixi run build     # maturin develop --release
    $ pixi run test      # pytest pydeeptools/deeptools/test/

Building from source
--------------------

Installs from `pypi <https://pypi.org/>`_ (``pip``, ``uv``, ``pipx``) use a prebuilt wheel
when one is available for your platform, so **no compiler is needed**. deepTools is only
built from source when no matching wheel exists (an unsupported platform, or a git
checkout / source distribution). In that case its Rust extension has to be compiled, which
needs a build toolchain: Rust (see ``Cargo.toml``'s ``rust-version`` for the minimum
supported version), ``maturin``, ``libclang``, ``pkg-config``, ``perl`` (to build the
vendored OpenSSL) and the usual HTSlib system libraries (``zlib``, ``bzip2``, ``liblzma``,
``libcurl``).

There are two equally supported ways to provide that toolchain -- pick whichever suits
your setup:

1. **Self-contained, via the pixi workspace (recommended for reproducible builds).**
   ``pixi install`` pins Rust, ``maturin``, ``libclang``, ``htslib``, ``pkg-config`` and
   ``openssl`` from conda-forge -- see ``pyproject.toml``'s ``[tool.pixi.dependencies]``
   for the exact pins -- so nothing has to be present on the host, and it works identically
   on Linux and macOS. The workspace also links conda-forge's OpenSSL
   (``OPENSSL_NO_VENDOR``), so no ``perl`` toolchain is needed on this route. Use the
   ``pixi install`` / ``pixi run build`` steps shown above.

2. **System toolchain, via ``pip`` / ``uv`` / ``pipx``.** Install Rust yourself (for
   example with `rustup <https://rustup.rs/>`_ or Homebrew) plus the system libraries
   listed above, then build with any of:

   .. code:: bash

       $ pip install .          # from a checkout
       $ uv tool install deeptools
       $ pipx install deeptools

   This is the same path deepTools' CI uses, so it is well exercised across platforms.

Both routes produce the same package; they differ only in where the build toolchain comes
from.

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
