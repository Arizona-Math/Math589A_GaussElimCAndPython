Gaussian Elimination and LU-decomposition in C and Python
=========================================================

C code
------

These are implementations of Gaussian Elimination and LU-decomposition in C.
The extra feature (beyond the typical textbook version) is the ability to
perform the algorithms 'in place', so that no new memory is allocated on the C side.


Python code
-----------

The module file `gauss_solve.py` implement functions which accomplish their
tasks by calling C via the `ctypes` module.

Makefile
--------

The project is managed by `make`, available on all UNIX systems. The
file `Makefile` contains all recipes to build the dynamic library
(shared object file) `libgauss.so` which is callable from
Python. Also, a basic driver `main.c` is included, to allow testing C
code from C itself (without Python).
