Installation
============

For installing _PyForFluids_ you just need to:

.. code-block:: sh

        pip install pyforfluids

Make sure to check the requirements first!

Requirements
------------

Be sure to install `numpy` and a fortran compiler previously, since both are
needed for the compilation of `Fortran` code.

NumPy
~~~~~

.. code-block:: sh

        pip install numpy

Fortran Compiler
~~~~~~~~~~~~~~~~

**Linux**

- **Debian-based** (Debian, Ubuntu, Mint,...)

.. code-block:: sh

        sudo apt install gfortran

- **Arch-based** (Arch, Manjaro, Garuda, ...)

.. code-block:: sh

   sudo pacman -S gfortran

**Windows**

We recommended using the Windows Subsystem for Linux and following the Linux
instructions.

`WSL <https://www.windowscentral.com/install-windows-subsystem-linux-windows-10>`_

If WSL ain't being used, the native Windows wheels will be download instead,
so no need to worry!

**MacOS**

.. code-block:: sh

        brew install gfortran


.. toctree::
   :maxdepth: 4

   pyforfluids
   pyforfluids.models
