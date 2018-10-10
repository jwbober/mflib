This is supposed to be a C++ library and associated programs for computing
spaces of classical modular forms, with functionality over the complex number
and modulo primes.

It's probably not very user friendly at the moment, and has so far only been
used by the author.

Some executables built by default
=================================

`newforms_acb`
--------------

The main useful program that exists right now is `newforms_acb`, with command
line

`newforms_acb level weight chi ncoeffs outpath prec targetaccuracy [nthreads] [verbose]`

When run, this creats the file `outpath/level/weight/level.weight.chi.mfdb`
with coefficients of that level, weight, and character. All computations are
with `prec` working precision, and `ncoeffs` will be computed and output to an
absolute accuracy at most 2^`targetaccuracy`, though if the working precision
was not high enough they may be output with lower accuracy (and as a special
exception they may be output with infinite accuracy when we can determine that
they are integers).

The program prints a lot to standard output. If the strings `error` or `ohno`
don't appear anywhere, then everything should have gone well.

Some parts of the computation are multithreaded, and will use as many threads
as specified on the command line (default 1), though if many threads are used
the bottleneck may be in single-threaded parts of the computation. Given
infinite memory, for computing many spaces of forms it is best to just run one
single threaded process for each, but sometimes the computation can take a
lot of memory and then it is better to run fewer processes and allow more
threads. To compute just a single space, it should be the case that more
threads are always better.

`print-mfdb`
------------

Prints contents of a `.mfdb` file produced by the above.

`newform-dimension`
-------------------

`newform-dimension N k` prints the dimensions of each `S_k^new(N, chi)`, while
`newform-dimension N1 N2 k1 k2` prints the same for each N1 <= N <= N2 and k1
<= k <= k2.

`hecke-polynomials-from-mfdb`
-----------------------------

Given an input mfdb file, this program computes the decomposition of each
character space into Galois orbits over the rational numbers, and puts this
information into a `polydb` output file. The program expects that if there
is data about one character in a Galois orbit then there is data about all
the characters in that Galois orbit, except that only one of each pair of
complex conjugate characters should be present, and it should be the one
with a lower index. So the output from the program `newforms_acb` needs to
be joined together into one database file before being input to the program.

The usage looks like

`hecke-polynomials-from-mfdb level weight mfdb polydb [overwrite] [nthreads]`

`level` and `weight` are integers (the input mfdb file could in principle contain
coefficients for many different levels and weights) and `mfdb` and `polydb` are
filenames. If the output file already exists, then setting `overwrite` to a
nonzero integer will cause _everything_ in the output file to be erased. (That's
probably not the best way of doing things, but it is the way things are now.)

Only some parts of the program are multithreaded, and in particular it will
only use as many threads as there are characters in the Galois orbit.

The program print a ridulous amount to standard output right now. If all
goes well it should terminate with the message `main function exiting normally.`.

In the misc/ directory there are a few programs for reading the output of this
program.

libmodform
==========

See the various source files that build executables for usage. The main way to
create a space of cusp forms is with the `get_cuspforms_acb` or
`get_cuspforms_modp` functions. (A space of cusp forms needs to know about
spaces of lower levels, so there is a clunky caching feature to make sure that
each space is only created once.)

Some of the useful functions in the objects created are `newspace_basis`, `trace`,
`newforms` (only available over the complex numbers) and `hecke_matrix` (only
available mod p).

How to Build
============

Install dependencies, type `make` and cross your fingers. Maybe modified some
hardcoded paths in the Makefile. `make link` to put symlinks in `~/bin`,
`~/include` and `~/lib`. Type `make` again in the misc directory to try to
build a bunch of examples.

Dependencies
============

flint, arb, gmp or similar, NTL, sqlite3. Also, a class number table and a
factorization table. (See below.)

Not every library is needed for every feature. Maybe a program that only
works mod p doesn't need to link to arb, for example. NTL is currently
only used in the `hecke_polynomials` programs (for its excellent polynomial
factorization). sqlite3 is only needed for reading and writing the default
file format.

Factorization table
-------------------

This is not strictly necessary for using the library, but it is expected
to be present for most of the executables, which expect to find a file
`~/include/int-factorization-table`. At the (12n)th byte of this file, there
should be a 12 byte struct

```
struct int_factor_t {
    int p;
    int e;
    int f; // f == p^e
};
```

where p is the largest prime dividing n, e is the power to which it divides n,
and f is p^e. A file of this sort can be produced by
`misc/build-factor-table.cc`. A largest possible table is 25769803764 bytes in
size. (Maybe everything really should have been unsigned here.)

Class number table
------------------

We expect to find a (binary) file ~/include/imaginary-quadratic-classnumbers
which, at the (4n)th byte, should have a integer equal to the class number of
the imaginary quadratic order of discriminant -n. (Whenever such a order
exists). In other words, we mmap this file to an array `int * classnumbers` and
expect that `classnumbers[d]` is h(-d).

A table of this sort can be built by the code in
`misc/build-classnumber-table.cc` (which requires pari). Also, note that the
executable built by that source file can be run multiple times in parallel, and
the results concatenated together. (d < 1000000000 has so far been sufficient
for my purposes, and only takes 4 billion bytes to store.)
