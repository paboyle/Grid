# Using Xcode with Grid on Mac OS

This guide explains how to use Xcode as an IDE on Mac OS.

# Initial setup

For first time setup of the Xcode and Grid build environment on Mac OS, you will need to do the following, in this order:

1. Install Xcode and the Xcode command-line utilities
2. Set Grid environment variables
3. Install and build Open MPI ***optional***
4. Install and build Grid pre-requisites
5. Install and Build Grid

Apple's [Xcode website][Xcode] is the go-to reference for 1, and the definitive reference for 4 and 5 is the [Grid Documentation][GridDoc].

[Xcode]: https://developer.apple.com/xcode/
[GridDoc]: https://github.com/paboyle/Grid/blob/develop/documentation/Grid.pdf

The following sections explain these steps in more detail

## 1. Install Xcode and the Xcode command-line utilities

See Apple's [Xcode website][Xcode] for instructions on installing Xcode.

Once Xcode is installed, install the Xcode command-line utilities using:

    xcode-select --install

*NB: the screenshots from this guide were generated from Xcode 10.1.*

## 2. Set Grid environment variables

To make sure we can share Xcode projects via git and have them work without requiring modification, we will define Grid environment variables. To make sure these environment variables will be available to the Xcode build system, issue the following command:

    defaults write com.apple.dt.Xcode UseSanitizedBuildSystemEnvironment -bool NO

These are the environment variables we will define for Grid:

Variable | Typical Value | Use
--- | --- | ---
`Grid` | `/Users/user_id/src/Grid` | Path to grid source
`GridPre` | `/Users/user_id/bin` | Path to install directory containing grid pre-requisites built from source
`GridPkg` | **MacPorts**=`/opt/local`, **Homebrew**=`/usr/local` | Path to package manager install directory
`GridDebug` | `mpidebug` | Grid environment for debug configurations
`GridRelease` | `mpirelease` | Grid environment for release configurations

Choose either of the following ways to do this:

### Method 1 -- Apple Script

* Start *Script Editor* (cmd-space, *script editor*)
* Click on *New Document*. Paste the following into the new script, editing the paths appropriately (just replace `user_id` with your *user_id* if you are unsure):

```apple script
do shell script "launchctl setenv Grid $HOME/src/Grid
launchctl setenv GridDebug mpidebug
launchctl setenv GridRelease mpirelease
launchctl setenv GridPre $HOME/bin
launchctl setenv GridPkg /opt/local"
```

* Save the script inside `~/Applications` and give it the name `GridEnv.app`.
* Open `System Preferences`, `Users & Groups`
* Click on `Login Items`
* Click the plus sign to add a new login item
* Select the `~/Applications` folder and select `GridEnv.app`

Log out and in again.

### Method 2 -- `environment.plist`

Make the file `environment.plist` in `~/Library/LaunchAgents` with the following contents, editing the paths appropriately (just replace `user_id` with your *user_id* if you are unsure):

```html
<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE plist PUBLIC "-//Apple//DTD PLIST 1.0//EN" "http://www.apple.com/DTDs/PropertyList-1.0.dtd">
<plist version="1.0">
<dict>
<key>Label</key>
<string>Grid.startup</string>
<key>ProgramArguments</key>
<array>
<string>sh</string>
<string>-c</string>
<string>launchctl setenv Grid $HOME/src/Grid
launchctl setenv GridDebug mpidebug
launchctl setenv GridRelease mpirelease
launchctl setenv GridPre $HOME/bin
launchctl setenv GridPkg /opt/local</string>

</array>
<key>RunAtLoad</key>
<true/>
</dict>
</plist>
```

These environment variables will be set each time your machine boots, so either reboot, or load them manually with:

    launchctl load  ~/Library/LaunchAgents/environment.plist
    launchctl start ~/Library/LaunchAgents/environment.plist

NB: if they are already loaded, you will need to unload them first, with:

    launchctl stop   ~/Library/LaunchAgents/environment.plist
    launchctl unload ~/Library/LaunchAgents/environment.plist

## 3. Install and build Open MPI -- ***optional***

Download the latest version of [Open MPI][OMPI] version 3.1 (I used 3.1.3).

NB: Grid does not have any dependencies on fortran, however many standard scientific packages do, so you may as well download GNU fortran (e.g. MacPorts ``gfortran`` package) and build Open MPI like so:

[OMPI]: https://www.open-mpi.org/software/ompi/v3.1/

    ../configure CC=clang CXX=clang++ F77=gfortran FC=gfortran CXXFLAGS=-g --prefix=$GridPre/openmpi-3.1.3
    make -j 4 all install

(If you don't want to bother with fortran bindings, just don't include the F77 and FC flags)

## 4. Install and build Grid pre-requisites

To simplify the installation of **Grid pre-requisites**, you can use your favourite package manager, e.g.:

### 1. [MacPorts][MacPorts]

[MacPorts]: https://www.macports.org "MacPorts package manager"

Install [MacPorts][MacPorts] if you haven't done so already, and then install packages with:

    sudo port install <portname>

These are the `portname`s for mandatory Grid libraries:

* git
* gmp
* mpfr

and these are the `portname`s for optional Grid libraries:

* fftw-3
* hdf5
* lapack
* doxygen
* OpenBLAS

***Please update this list with any packages I've missed! ... and double-check whether OpenBLAS is really for Grid***

### 2. [Homebrew][Homebrew]

[Homebrew]: https://brew.sh "Homebrew package manager"

Install [Homebrew][Homebrew] if you haven't done so already, and then install packages with:

    sudo brew install <packagename>

The same packages are available as from MacPorts.

### Install LIME ***optional***

There isn't currently a port for [C-LIME][C-LIME], so download the source and then build it:

[C-LIME]: https://usqcd-software.github.io/c-lime/ "C-language API for Lattice QCD Interchange Message Encapsulation / Large Internet Message Encapsulation"

../configure --prefix=$GridPre/lime-1.3.2 CC=clang
make -j 4
make install

## 5. Install and Build Grid

### Install Grid

[Grid]: https://github.com/paboyle/Grid

Start by cloning [Grid (from GitHub)][Grid] ***into the directory you specified in***  `$Grid`. Bear in mind that git will create the `Grid` subdirectory to place Grid in, so for example if `$Grid` is set to `~/src/Grid` then install Grid with:

    cd ~/src

followed by either:

    git clone git@github.com:paboyle/Grid.git

or

    git clone https://github.com/paboyle/Grid.git

depending on whether you are using https or ssh.

### Build Grid

The Xcode build system supports ***debug*** and ***release*** configurations for each project (more configurations can be defined). We will create separate Grid build directories for each configuration, using the Grid **Autoconf** build system to make each configuration. NB: it is **not** necessary to run `make install` on them once they are built (IDE features such as *jump to definition* will work better of you don't).

Below are shown the ``configure`` script invocations below for debug and release configurations, with and without MPI. You are unlikely to need all four combinations, so as a minimum, just build ``build_mpidebug`` (or ``build_debug`` if you don't want MPI). You probably want to build the corresponding release configuration (``build_mpirelease`` or ``build_release``), but this is optional.

For each configuration, run the `configure` script (shown below) ***first***, then make the configuration with:

    make -j 4

NB: you **do not** need to run ``make install`` for these to work with Xcode.

### 1. build_mpidebug

This is the build for every day developing with Xcode. It uses the Xcode clang c++ compiler, with MPI and debug symbols. This is the version of Grid Xcode links to for **debug** configurations.

    ../configure --with-hdf5=$GridPkg --with-gmp=$GridPkg --with-mpfr=$GridPkg --with-fftw=$GridPkg --with-lime=$GridPre/lime-1.3.2 --enable-simd=GEN --enable-precision=double CXX=clang++ --prefix=$GridPre/GridMPIDebug --enable-comms=mpi-auto MPICXX=$GridPre/openmpi-3.1.3/bin/mpicxx CXXFLAGS=-g --enable-doxygen-doc

### 2. build_mpirelease

Much the same as `build_mpidebug`, except without debug symbols for **release** configurations (it can be handy to use `#ifdef DEBUG` while developing, and it's useful to be able to compile and test the alternate case).

    ../configure --with-hdf5=$GridPkg --with-gmp=$GridPkg --with-mpfr=$GridPkg --with-fftw=$GridPkg --with-lime=$GridPre/lime-1.3.2 --enable-simd=GEN --enable-precision=double CXX=clang++ --prefix=$GridPre/GridMPIRelease --enable-comms=mpi-auto MPICXX=$GridPre/openmpi-3.1.3/bin/mpicxx CXXFLAGS=-g --enable-doxygen-doc

### 3. build_debug

Debug configuration without MPI:

    ../configure --with-hdf5=$GridPkg --with-gmp=$GridPkg --with-mpfr=$GridPkg --with-fftw=$GridPkg --with-lime=$GridPre/lime-1.3.2 --enable-simd=GEN --enable-precision=double CXX=clang++ --prefix=$GridPre/GridDebug --enable-comms=none CXXFLAGS=-g --enable-doxygen-doc

### 4. build_release

Release configuration without MPI:

    ../configure --with-hdf5=$GridPkg --with-gmp=$GridPkg --with-mpfr=$GridPkg --with-fftw=$GridPkg --with-lime=$GridPre/lime-1.3.2 --enable-simd=GEN --enable-precision=double CXX=clang++ --prefix=$GridPre/GridRelease --enable-comms=none

# Make a new application using Grid

NB: Instead of following the instructions in this section, you can clone `HelloGrid` from the [University of Edinburgh GitLab site][HelloGrid].

[HelloGrid]: https://git.ecdf.ed.ac.uk/s1786208/HelloGrid

## Make a new application

To make a hello world application for Grid:

* Start Xcode
* Click 'Create a new project'
* Click ‘macOS’, then in the ‘Application’ section choose ‘Command Line Tool’, then click ‘Next’
* Choose options for your new project:
  * Product Name: HelloGrid
  * Team: None
  * Organisation Name: sopa
  * Organisation Identifier: uk.ac.ed.ph
  * Language: C++
  * ... then click ‘Next’
* Choose a location for your project, e.g. `$HOME/src`. NB: The project and all it’s files will be created inside `$HOME/src/HelloGrid`. If you are using Git, you can put the new project under Git source control immediately, if you like. Now click ‘Create’.

## Configure your new application to use Grid

Click the project name (`HelloGrid`) in the project navigator pane on the left (command-1 if it's not visible), then click the project name (`HelloGrid`) under `PROJECT` in the second pane. Click the `Build Settings` tab on the right, then under that click `All` and `Combined`. You should see:

![Project settings](GridXcFig1.png)

We now need to make changes to two sections (these are listed in alphabetical order), bearing in mind that if you are not using MPI (or you gave your build directories different names) replace `build_mpidebug` and `build_mpirelease` with the directory names you used.

### 1. Search Paths

#### HEADER_SEARCH_PATHS

Obtain a list of header locations required by Grid by running the following from your Grid build directory:

    ./grid-config --cxxflags

Output should look similar to:

    -I$GridPre/openmpi-3.1.3/include -I$GridPkg/include -I$GridPre/lime-1.3.2/include -I$GridPkg/include -I$GridPkg/include -I$GridPkg/include -O3 -g -std=c++11

The header locations follow the `-I` switches. You can ignore the other switches, and you can ignore duplicate entries, which just mean that your package manager has installed multiple packages in the same location.

*Note: `grid-config` will output absolute paths. Make sure to replace absolute paths with environment variables (such as `$GridPre`) in your settings, so that the project will work unmodified for other collaborators downloading the same project from git.*

Set the **Debug** HEADER_SEARCH_PATHS to:

    $Grid/build_$GridDebug/Grid
    $Grid

followed by (***the order is important***) the locations reported by `grid-config --cxxflags`, ignoring duplicates, e.g.:

    $GridPre/openmpi-3.1.3/include
    $GridPkg/include
    $GridPre/lime-1.3.2/include

**Note: the easiest way to set this value is to put it all on one line, space separated, and edit the text to the right of `Debug`**, i.e.:

    $Grid/build_$GridDebug/Grid $Grid $GridPre/openmpi-3.1.3/include $GridPkg/include $GridPre/lime-1.3.2/include

Similarly, set the **Release** HEADER_SEARCH_PATHS to exactly the same settings, replacing `debug` in the first entry with `release`, e.g.:

    $Grid/build_$GridRelease/Grid $Grid $GridPre/openmpi-3.1.3/include $GridPkg/include $GridPre/lime-1.3.2/include

#### LIBRARY_SEARCH_PATHS

Obtain a list of library locations required by Grid by running the following from your Grid build directory:

    ./grid-config --ldflags

Output should look similar to:

    -L$GridPre/openmpi-3.1.3/lib -L$GridPkg/lib -L$GridPre/lime-1.3.2/lib -L$GridPkg/lib -L$GridPkg/lib -L$GridPkg/lib

Set the **Debug** LIBRARY_SEARCH_PATHS to:

    $Grid/build_$GridDebug/Grid
    $Grid/build_$GridDebug/Hadrons

followed by the locations reported by `grid-config --ldflags`, ignoring duplicates (again, the order is important), e.g.:

    $GridPre/openmpi-3.1.3/lib
    $GridPkg/lib
    $GridPre/lime-1.3.2/lib

Again, this can be done all on one line:

    $Grid/build_$GridDebug/Grid $Grid/build_$GridDebug/Hadrons $GridPre/openmpi-3.1.3/lib $GridPkg/lib $GridPre/lime-1.3.2/lib

Similarly, set the **Release** LIBRARY_SEARCH_PATHS to exactly the same settings, replacing `debug` in the first two entries with `release`, e.g.:

    $Grid/build_$GridRelease/Grid $Grid/build_$GridRelease/Hadrons $GridPre/openmpi-3.1.3/lib $GridPkg/lib $GridPre/lime-1.3.2/lib

### 2. Linking

#### OTHER_LDFLAGS

The easiest way to link to all required libraries is to obtain a list of all libraries required by Grid by running the following from your Grid build directory:

    ./grid-config --libs

and pasting the output ***with `-lGrid -lHadrons ` prepended*** (including the `-l` switches) directly into `OTHER_LDFLAGS`, e.g.:

    -lGrid -lHadrons -lmpi -lhdf5_cpp -lz -lcrypto -llime -lfftw3f -lfftw3 -lmpfr -lgmp -lstdc++ -lm -lz -lhdf5

NB: The library names should be the same for your **debug** and **release** configurations, so you can set this once for both configurations.

## Edit your source code

A hello world for grid is:

```c++
#include <Grid/Grid.h>
using namespace Grid;

int main(int argc, char * argv[]) {
  Grid_init(&argc,&argv);
  std::cout << GridLogMessage << "Hello Grid" << std::endl;
  Grid_finalize();
  return 0;
}
```

## Create a `.gitignore` file for Xcode

You can create an up-to-date .gitignore file to ignore all the Xcode temporary build files using [gitignore.io][GIO].

[GIO]: https://www.gitignore.io/api/xcode

NB: If you let Xcode add your project to git when you created it, you probably want to remove your personal scheme selection from git:

    git rm --cached HelloGrid.xcodeproj/xcuserdata/$USER.xcuserdatad/xcschemes/xcschememanagement.plist

## Run your program under the Xcode debugger

First, specify command-line arguments. From the menu, select `Product`, then `Scheme`, then `Edit Scheme`. Select `Run` on the left, then select the `Arguments` tab on the right. Add the following to `Arguments passed on Launch`:

    --grid 4.4.4.8

If your program will be manipulating files, it's a good idea to specify the working directory on the `Options` tab under `Use Custom Working Directory` (by default, Xcode launches the program inside the Xcode build folder). 

Then click `Close`.

Let's set a breakpoint by clicking on:

    Grid_finalize();

then from the menu selecting `Debug`, then `Breakpoints`, then `Add Breakpoint at Current Line`.

Now click on the `Play` button (the right pointing triangle just to the right of the maximise button) to run your program under the debugger. (You may see dialog boxes the first couple of times asking whether to allow MPI to receive network requests - say yes to these.)

The debug output pane opens at the bottom of Xcode, with output on the right (ending with `Hello Grid`) and local variables on the left i.e.:

![Running under the debugger](GridXcFig2.png)

See the Xcode documentation to learn about the debugger. When you're done, press `ctl-cmd-Y` to let the program run to completion. 

## Running multiple MPI processes

### Use Xcode to launch `mpirun`

You can tell Xcode to use mpirun to launch multiple copies of a target executable, however if you do this the debugger will attach to mpirun - not your target process. This can be useful - so long as you do not need to debug.

To do this, edit the Scheme again (cmd-<), click on the `info` tab, then under Executable select `Other...` and enter the full path to your `mpirun`, e.g.:

    $GridPre/openmpi-3.1.3/bin

NB: if `mpirun` is a link, don't be surprised if Xcode saves the name of the linked-to executable.

Click on the `Arguments` tab, then under `Arguments Passed on Launch` add:

    -np 2 $(TARGET_BUILD_DIR)/$(TARGETNAME)

and make sure this is the first argument - i.e. drag it to the top of the list of arguments. NB: You probably want to specify more arguments for the MPI run, but this should get you started.

Then click `Close`.

### Use Xcode to debug multple instances of your target

From the `Debug` menu, select `Attach to Process by PID or Name ...`.  In the `PID or Process Name` field, enter the name of your target. Then click `Attach`.

From a terminal session, locate and run your executable using mpirun (*the mangled name of the project build products will not be exactly the same as this example*):

    $GridPre/openmpi-3.1.3/bin/mpirun -np 2 ~/Library/Developer/Xcode/DerivedData/HelloGrid-fiyyuveptaqelbbvllomcgjyvghr/Build/Products/Debug/HelloGrid --grid 4.4.4.8

The Xcode debugger will attach to the first process. From the `Debug` menu in Xcode, select `Attach to Process`, and other running instances of your application will appear at the top of the list. Attach to as many instances as you wish to debug.

You are now debugging multiple MPI instances, and the Xcode debugger should look similar to this:   

![Debugging multiple MPI instances under the Xcode debugger](GridXcFig3.png)
