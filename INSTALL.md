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

To install Hdef in the directory /opt/hdef, type:

```
cmake --install . --prefix /opt/hdef
```

Note that you may need to run this command with `sudo`.




## Questions, Comments, or Concerns?

Contact me at matthew.w.herman@gmail.com.
