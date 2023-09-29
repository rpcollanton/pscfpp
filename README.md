
# PSCF - Polymer Self-Consistent Field Theory (C++/CUDA)

PSCF is a package of software for field-theoretic predictions of 
equilibrium properties of inhomogeneous polymer liquids, with a focus 
on microstructures formed by block polymers. The package provides an 
extensive set of tools for self-consistent field theory (SCFT) 
calculations.  The current version is the first to also provides tools
for field-theoretic Monte-Carlo (FTMC) simulations, which are based on 
a partial-saddle point approximation to the formulation of the 
partitition function of the material as a functional integral.
The acronym PSCF stands for "Polymer Self-Consistent Field", which 
was the original purpose of the package. The version described here 
is written primarily in C++, with GPU accelerated code in CUDA.

## History

This C++/CUDA version of PSCF is intended to supersede an older Fortran 
SCFT program of the same name. The Fortran PSCF program is maintained in 
a separate github.com repository at https://github.com/dmorse/pscf.
This current C++/CUDA version provides almost all of the capabilities of
the legacy Fortran program, and several important new capabilities, as
discussed below.

Differences between the current C++/CUDA version of PSCF and the legacy
Fortran version include:

   - The new version is an extensible package of several different 
     programs designed for use with different geometries, boundary 
     conditions, algorithms or hardware, implemented within a common
     framework.

   - The new version enables simulations of mixtures containing both
     linear and acyclic branched block polymers, while the Fortran PSCF 
     code allowed only linear polymers.

   - The new version enables use of graphics processing units (GPUs) to
     dramatically accelerate some programs.

   - Starting with this release, the C++/CUDA version provides tools
     for stochastic field-theoretic Monte-Carlo (FTMC) simulations.

## Programs

PSCF currently contains code for three programs:

   - **pscf_fd** : The program pscf_fd is is designed to performed SCFT
     calculations on one-dimensional problems in Cartesian, cylindrical 
     or spherical coordinates. A simple finite difference method is used
     to solve the underlying partial differential equation. It is useful 
     for treating problems involving flat and curved interfaces, as well 
     as cylindrical or spherical copolymer micelles.  The suffix "fd" 
     stands for "finite difference".
     
   - **pscf_pc** : The pscf_pc rogram is designed to treat structures
     that are periodic in 1, 2 or 3 dimensions, using a pseudo-spectral
     algorithm and conventional CPU hardware. This program provides 
     capabilities for SCFT calculations analogous to those of the older 
     PSCF Fortran program, as well as tools for FTMC simulations. The
     The suffix "pc" stands for "periodic CPU".  
     
   - **pscf_pg** : The pscf_pg program is a GPU-accelerated variant
     of pscf_pc that provides tools for SCFT and FTMC simulations. It 
     is based on algorithms similar to those used in pscf_pc, with 
     almost identical features but provides higher performance for 
     large systems. The suffix "pg" stands for "periodic GPU". 

## Features

All PSCF programs are designed to treat an incompressible mixture 
containing any number of block polymers, homopolymers and small molecule 
(point-like) solvent molecular species. All polymeric species are 
treated using the standard Gaussian model of each polymer block as a 
continuous random walk.

Features that are common to all three PSCF programs that apply to both
SCFT and FTMC calculations include:

  - Ability to treat mixtures of any number of block polymers, 
    homopolymers, and solvent species

  - Arbitrarily complex acyclic branched block polymers, in addition to 
    linear block polymers

  - Canonical, grand-canonical or mixed statistical ensembles - users 
    may specify either a volume fraction or a chemical potential for 
    each molecular species

  - Thorough user and developer documentation provided as a web manual

  - Open source code written in object oriented C++

Features for SCFT calculations that are common to all three PSCF programs 
include:

  - Efficient Anderson-mixing iteration algorithms

  - Efficient SCFT calculations for sequences of parameter choices along 
    a path in parameter space ("sweeps"), using extrapolation to 
    construct initial guesses

  - Examples of converged SCFT solutions and input scripts

  - Python tools for data analysis

Features specific to the pscf_pc and pscf_pg programs for periodic 
structures include:

  - Ability to perform both SCFT and FTMC calculations

  - Accurate pseudo-spectral solution of the modified diffusion equation

  - Ordered phases with 1, 2 or 3 dimensional periodicity

  - All possible 2D and 3D lattice systems (i.e., orthorhombic, monoclinic,
    etc.)

  - Automatic optimization of unit parameters in SCFT so as to minimize 
    free energy density

  - Imposition of any user-selected space-group symmetry in SCFT

  - Built-in database of symmetry operations for all 230 3D space groups 
    and 17 2D plane groups for use in SCFT

  - Availability of a companion package 
    [Polymer Visual](<https://github.com/kdorfmanUMN/polymer_visual/>)
    for visualization of periodic structures

Features specific to SCFT calculations performed with the pscf_pc CPU 
program are:

  - Inhomogeneous density constraints (a "mask")

  - External fields

  - Thin polymer films

The mask and external field features are used in pscf_pc to implement 
simulations of thin films. A mask is used to constrain a polymer 
material to a slit within a periodic supercell, and localized external 
fields are used to represent selective interactions with the top and 
bottom surfaces.

## Getting the source code

The current PSCF source code is maintained in the github repository

   <https://github.com/dmorse/pscfpp>.

The source code may be obtained by using a git version control system
client to clone the repository. To do so on a machine with a git client,
enter the command:
```
git clone --recursive https://github.com/dmorse/pscfpp.git
```
Note the use of the --recursive option to the git clone command. This
option is necessary to clone several git submodules that are maintained 
in separate repositories. This command will create a new directory 
called pscfpp/ that contains all of the source code and associated
documentation, including all required git submodules.

We do *not* recommend that users obtain the source code by simply
downloading a zip or tar file of a tagged release from the PSCF github 
repository, because this file will not contain source code for the other 
git repositories that are automatically downloaded and installed as 
submodules by the above git command. It is possible to install the
required submodules after the fact, but requires some additional 
steps.

## Documentation

PSCF is distributed with source files for an html web manual. A copy of 
the manual for the most recent numbered release is available online at

  <https://dmorse.github.io/pscfpp-man>

After cloning the source code, users can also use the 
[doxygen](<https://www.doxygen.nl/index.html>) documentation generation 
program to generate a local copy of the web manual.  To do this, the 
doxygen application must be installed on your computer, and the directory 
containing the doxygen executable must be in your unix command search 
PATH. To generate documentation:

   - Change directory (cd) to the pscfpp/ root directory

   - Enter "make html"

This should create many html files in the docs/html subdirectory of the 
pscfpp/ root directory.  To begin reading the resulting documentation, 
point a web browser at the file pscfpp/docs/html/index.html .  This is 
the main page of the web manual.

## Dependencies

PSCF has been developed on both linux and and Mac OS X operating systems,
and is designed to run on these or other unix-like systems.

The pscf_fd and pscf_pc CPU programs depend on the following external 
libraries:

  - [GSL](<https://www.gnu.org/software/gsl/>) - GNU Scientific Library

  - [FFTW](<https://www.fftw.org/>) Fast Fourier Transform library

The pscf_fd one-dimensional program pscf_fd requires only GSL, while
pscf_pc program relies on both GSL and FFTW.

The GPU-accelerated pscf_pg program can only be compiled and run on a 
computer with an appropriate NVIDIA GPU and an NVIDIA CUDA development kit.

Procedures for installing these dependencies are different for different 
operating system environments and different package managers. Instructions 
for some common environments are given in the web manual.

## Compiling

The PSCF source code is written using C++ as the primary language, with
CUDA used in pscf_pg for GPU acceleration. PSCF is only provided in 
source file format, and so must be compiled from source.

The build system used to compile and install PSCF relies on standard unix 
utilities that are available in any unix-like command-line environment. 
To compile and run this or other unix software on Mac OS X, the user must 
first install the Mac OS X unix command line tools. 

Complete directions for compiling and installing PSCF are provided in
section 2 of the [web manual](<https://dmorse.github.io/pscfpp-man>).
A brief summary is given below of instructions for steps that must be 
taken after cloning the git repository and installing all of the 
required dependencies:

   - Add the pscfpp/bin directory to your linux command search PATH
     environment variable. This is where executables will be installed.

   - Add the pscfpp/lib/python directory to your PYTHONPATH environment 
     variable. This provides access to python scripts used by the build 
     system.

   - Run the pscfpp/setup script by entering "./setup" from the pscfpp/ 
     root directory, with or without an optional filename argument. 
     This script installs default versions of several files that are 
     required by the build system. The setup script usually only needs 
     to be run once, before compiling for the first time.  In a linux 
     environment, it is usually sufficient to run the setup script 
     without a filename argument.  See Section 2.6 of the 
     [web manual](<https://dmorse.github.io/pscfpp-man>) for more 
     information about how to invoke the setup script for Apple OS X.

   - To compile and install only CPU-based programs in the package on
     a computer that does not have an appropriate NVIDA GPU, enter
     "make all" from the pscfpp/ root directory immediately after 
     running the setup script.

   - To compile all PSCF programs on a machine that has an appropriate
     NVIDA GPU and a CUDA development kit, first enter "./configure -c1" 
     from the pscfpp/ root directory to enable compilation of CUDA code, 
     then use the -a option of the configure script to set the correct 
     GPU architecture for your system, and enter "make all" to compile. 

The "make all" command installs all resulting executable files in the 
pscfpp/bin directory.

## Command line usage

PSCF is a package containing three different SCFT programs (pscf_fd,
pscf_pc, and pscf_pg) that are designed for different geometries, using 
different algorithms or different computer hardware.  To perform a
calculation, Each of these programs must read a parameter file and a 
command file. The parameter file, which is processed first, is a fixed
format file that contains parameters required to describe the physical 
system of interest and initialize the program before performing any 
computation. The command file is a script that contains a sequence of 
commands that are read and executed in the order that they appear.  
Contents and syntax for PSCF parameter and commands files are discussed 
in Sec. 3 of the web manual. 

The command file usually contains the names of several input and output 
field data files as arguments to commands that read or write these files. 
Specifically, a command file to perform either a standard SCFT calculation 
or an FTMC simulation usually contains a command to read a specific file
that contains an initial guess for the monomer chemical potential fields.
A command file for an SCFT calculation normally also contains commands 
that write final chemical potential and monomer concentration fields to 
specific files. 

The usual command line syntax for invoking the pscf_fd program for 
one-dimensional SCFT programs is
```
pscf_fd -e -p param -c command
```
where "param" denotes the name of a parameter file, and "command" denotes 
the name of a command file.  The -p and -c options are required.

The "-e" command line option in the above example is not required, but 
causes the program to "echo" the parameter file to standard output as 
this file is being processed. Use of this option makes it easier to 
debug any failure arising from a syntax error in this file. 

The pscf_pc and pscf_pg programs require a value for the spatial 
dimension of the structure as an additional command line parameter, 
which specifies the number 1, 2, or 3 of coordinates along which the 
structure is periodic.  This must be provided as an integer argument of 
the -d option.  The usual syntax for invoking pscf_pc for a three 
dimensionally periodic structure (e.g., a network phase or an 
arrangement of spheres), while also using the -e  option, is thus
```
pscf_pc -d 3 -e -p param -c command
```
The syntax for a simulation of one-dimensional lamellar phase would 
instead use 1 as the argument of the -d option. The syntax for invoking 
pscf_pg is the same as that for pccf_pc, except for the executable 
program name.

Programs that are invoked as shown above would both write log output 
that is produced during execution to the user's screen (i.e., to standard 
output).  This log output may instead be redirected to a file by using 
the unix ">" standard output redirection operator. For example, a 
command such as 
```
pscf_pc -d 3 -e -p param -c command > log
```
would redirect log output that is written during simulation of a 
3-dimensionally periodic structure to a file named "log" in the current 
working directory. Standard output is normally redirected to a file when
a program is run in a queue on a shared computer cluster.

## Examples

The quickest way to learn the syntax of the PSCF parameter and command
files is by examining examples.  The directory pscfpp/examples contains 
examples of input files for many different types of SCFT calculations.  
Top level subdirectories of the examples/ directory contain examples 
for different PSCF programs (pscf_fd, pscf_pc, or pscf_pg). 

Subdirectory examples/fd contains examples of SCFT calculations for the 
1D finite-difference program pscf_fd. Top level subdirectories of the
directory examples/fd contain examples for planar, cylindrical and 
spherical geometries, as indicated by the subdirectory names. One or
more example is given for each geometry.

Subdirectory examples/pc contains examples for the pscf_pc CPU program. 
Top level subdirectories of examples/pc contain solutions for a particular 
type of physical system. For example, the directory examples/pc/diblock 
contains examples for a diblock copolymer melt.  Subdirectories of 
examples/pc/diblock contain examples for lamellar, hexagonal, and BCC 
structures, among others.

Subdirectory examples/pg contains examples for the pscf_pg program for 
periodic structures. Subdirectories are organized in a manner similar 
to that used for the examples/pc directory tree.

We refer to a directory that contains all of the input files for a single 
example or a set of very closely related examples as an example directory.  
Each such example directory contains at least one parameter file (usually 
named "param"), at least one command file (usually named "command"), and 
an input chemical potential field (w field) file. Example directories 
that contain input files for two or more examples usually contain a single
file that contains an initial guess for the chemical potential field, but 
may contain different parameter or command files for use in different 
examples.

Example directories that contain files for a single example usually
contain a shell script named "run" that can be executed to run the 
example using the supplied input files.  The simplest way to run such an 
example is thus to change directory (cd) to the relevant example directory 
and enter "./run" from within that directory. (Note the use of the 
prefix "./", which tells the operating system to look for the run script 
in the current working directory).  Users may also inspect the text of 
such a run script to see the exact command used to run the example. 
Example directories that contain input files for two or more closely 
related examples may contain several such scripts with names that are 
variants of "run", each of which can be execute to run a specific 
example. 

Each example directory also usually contains a script named "clean" that 
can be executed to remove all output files that are created by running 
the example.

## License

The C++/CUDA version of PSCF is free, open source software. It is
distributed under the terms of the GNU General Public License (GPL)
as published by the Free Software Foundation, either version 3 of the
License or (at your option) any later version.  PSCF is distributed
without any warranty, without even the implied warranty of merchantability
or fitness for a particular purpose.  See the LICENSE file or the
<a href=http://www.gnu.org/licenses/>gnu web page</a> for details.

## Support

Development of PSCF is currently supported by the National Science 
Foundation program for Cyberinfrastructure for Sustained Scientific 
Development (CSSI) under Grant No. 2103627.

## Contributors

- David Morse
- Guo Kang Cheong
- Anshul Chawla
- Ryan Collanton
- Ben Magruder
- Kexin Chen
- Ying Zheng

