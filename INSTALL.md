# Installing Hdef


## Prerequisites


### Required Prerequisites

- GNU Fortran (gfortran) and C (gcc) compilers
- CMake (3.12+)
- Make


### Recommended Prerequisites

- Git (for easily downloading Hdef and keeping it updated)
- LAPACK (Hdef will have fewer features if built without LAPACK)
- Generic Mapping Tools ([GMT](https://www.generic-mapping-tools.org)) (for plotting Hdef outputs)

I recommend using a package manager to install the prerequisites listed above. On macOS, I use [Homebrew](https://brew.sh) and on Debian-based Linux systems I use [APT](https://wiki.debian.org/Apt). For example, to install the required program `cmake`:

On macOS with Homebrew:

```
brew install cmake
```

On Ubuntu with APT:

```
apt install cmake
```


## Installation Instructions


### Download Hdef from [GitHub](https://github.com/mherman09/Hdef)

In your terminal, navigate to the directory where you keep external programs. For example, I have a Research/ subdirectory in my home directory. I recommend using `git` to download Hdef:

```
git clone https://github.com/mherman09/hdef ./hdef
```

If you download Hdef as a zip file instead of using `git`, make sure you unzip it before moving on to the next steps.



### Configure Hdef

From the directory where you downloaded Hdef, create a build/ directory where you will configure and compile Hdef:

```
cd hdef
mkdir build
cd build
```

On **Linux** operating systems with the GNU C and Fortran compilers, type:

```
cmake ..
```

On **macOS**, you will need to specify the GNU C compiler that was installed with your package manager. For example, if Homebrew installed GCC version 15, you would type:

```
cmake --DCMAKE_C_COMPILER=gcc-15 ..
```


### Build Hdef

Compile the Hdef programs, test functions, and scripts in your build/ directory with:

```
cmake --build .
```


### Test Hdef (optional, but recommended)

Hdef includes numerous unit tests and test scripts that can be run after building the programs. To test *all the Hdef subroutines and executables*, navigate to the build/ directory and type:

```
cmake --build . -t test
```

To run *only the executable test scripts*, type:

```
cmake --build . -t test-exec
```

If the codes are compiled and installed correctly, the tests *should* complete with no errors. However, the tests are still a bit of a work in progress. So please report testing issues to me if you run into them!



### Install Hdef (optional)

You only need to follow these instructions if you wish to have the Hdef executables, scripts, libraries, and tests installed in a particular location on your computer such as /opt/hdef or /usr/local/hdef. This can be useful if multiple users will be working with the code.

For example, to install Hdef in the directory /opt/hdef, type:

```
cmake --install . --prefix /opt/hdef
```

Note that you may need to run this command with `sudo` before `cmake`.



### Set PATH Variable

Finally, you should add the directory where you installed the Hdef programs to your PATH variable so your computer can find the programs. For example, if your shell is `bash` and you installed Hdef in `/opt/hdef`, add the following to the `.bashrc` file in your home directory:

```
export PATH=$PATH:/opt/hdef/bin
```

Close your terminal window and reopen it, and type:

```
o92util
```

If everything is set up properly, you should see the usage statement for the program `o92util`, which looks something like this:

```
 Usage: o92util ...options...

 Input fault options
 -ffm FFMFILE         Fault file in USGS .param format
 -fsp FSPFILE         Fault file in SRCMOD FSP format
 -mag MAGFILE         Fault file in "psmeca -Sa" format (...mag)
 -flt FLTFILE         Fault file with slip and dimensions (...slip wid len)
 -tns TNSFILE         Tensile source file (IN DEVELOPMENT)
 -fn|-pt              Treat faults as finite rectangular (default) or point
 -empirical OPT       Empirical scaling relation
 -thr THR             Set low slip to zero

 Input target/receiver options
 -sta STAFILE         Station/receiver locations
 -auto h DEPTH N      Generate horizontal location grid
 -auto v AZ N         Generate vertical location grid (through centroid)
 -auto:thr DISP       Displacement threshold for auto grids
 -trg TRGFILE         Target/receiver geometry

 Input half-space options
 -haf HAFSPCFILE      Elastic half-space properties

 Output options
 -disp DSPFILE        Displacement (E N Z)
 -strain STNFILE      Strain matrix (EE NN ZZ EN EZ NZ)
 -rotation ROTFILE    Rotation matrix (EE NN ZZ EN EZ NZ)
 -stress STSFILE      Stress matrix (EE NN ZZ EN EZ NZ)
 -estress ESTSFILE    Effective (maximum) shear stress
 -normal NORFILE      Normal traction on target faults (requires -trg)
 -shear SHRFILE       Shear traction on target faults (requires -trg)
 -coul COULFILE       Coulomb stress on target faults (requires -trg)
 --keep-all-lines     Keep comment (#), segment (>), and blank lines

 Miscellaneous options
 -parallel [NTHREADS] Calculate deformation in parallel
 -geo|-xy             Use geographic (default) or cartesian coordinates
 -az                  Displacement vector outputs (AZ HMAG Z)
 -gmt FILE            GMT psxy -SJ file (lon lat str len proj_wid)
 -prog                Turn on progress indicator
 -v LVL               Turn on verbose mode
 -debug [ROUTINE]     Turn on debugging

 See o92util man page for details
```


## Questions, Comments, or Concerns?

Contact me at matthew.w.herman@gmail.com.
