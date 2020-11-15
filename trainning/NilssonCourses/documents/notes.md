## A quick look at the source code of applications

Use the pre-defined alias app to go to the applications directory: $FOAM_APP

Allwmake is used to compile all the applications.
• solvers contains the source code of the solvers.
• utilities contains the source code of the utilities.
• test contains source code for testing specific features of OpenFOAM

## Source code and binary file directory organization, browsing, name conventions, and compilation of installation

### directory organization 
• etc: The OpenFOAM environment and global settings, needed for everything.
• tutorials: Usage examples. Assumed to be already known.
• applications: Source code of solvers and utilities.
• src: Source code of libraries, used by the applications.
• wmake: Command and instructions for compilation of applications and libraries.
Followed by a look at Allwmake and the Make directories.
• build: Intermediate compilation files, for applications and libraries.
• platforms: Final binaries for applications and libraries.
• bin: Executable bash scripts.
• doc: Doxygen source files and coding style instructions
• modules: Additional modules (not in slides)

### change compilation configuration at different level
Have a look in that file and see that the prefs.sh file must be found by foamEtcFile. You
can list the appropriate paths by:
foamEtcFile -list

Those are searched in order, and you can therefore do the modifications at different levels (for
only you, at a site, or for all users).

copy the file etc/config.sh/example/prefs.sh to your wanted etc level, and set the
environment variables in that file

### doc directory 
The doc directory contains source code for the Doxygen documentation of OpenFOAM. It needs
to be compiled to generate all the html pages. Go to  Doxygen folder and read README to compile the Doxygen User Guide in local .


## User directory organization, and compilation as a user

### User directory organization

• applications will be used for our own developed solvers and utilities
• src will be used for our own developed libraries
• run will be used for our cases, including running the original tutorials (assumed to be
known)

### Check PATH or LD_LIBRARY_PATH

For applications we discussed that the environment variable $PATH is used to find the executable files. The environment variable LD_LIBRARY_PATH is similarly used to find the libraries when dynamically linking:

You can check which libraries the solver is linking to by:
ldd `which icoFoam`


##  Find solver and utility tutorials in the source code and learn how to use them
• Some functionality that was earlier available as utilities for post processing have been
reorganized using functionObjects (discussed more later).
• We will have a look at how to extract data for plotting and visualization in the coming
slides.
• You can list the available options by typing
postProcess -list
Use the -help flag for more information, as usual.
• You can figure out how to do other kinds of postProcessing by looking at info and examples in $WM_PROJECT_DIR/etc/caseDicts/postProcessing/
- Gnuplot at http://www.gnuplot.info/
- postProcess dict: /etc/caseDicts/postProcessing/
- function objects code: $FOAM_SRC/functionObjects

## Debug 
 There are many debugging approaches, and we will discuss three alternatives here:
• Info statements in the code
• Built-in DebugSwitch option in OpenFOAM 
  - each class thus has its own DebugSwitch.
  - DebugSwitches set to zero will produce no debug information
  - Different levels of debugging can be chosen by setting a DebugSwitch to 1, 2, 3 ...
  - It is recommended to make a copy of that file to a specific location and make any modifications there. This file will override the original file

```cmd
mkdir -p $HOME/.OpenFOAM/$WM_PROJECT_VERSION
cp $WM_PROJECT_DIR/etc/controlDict $HOME/.OpenFOAM/$WM_PROJECT_VERSION
```

• Debugging using the Gnu debugger, GDB
    - In $WM_PROJECT_DIR/etc/bashrc you find an environment variable WM_COMPILE_OPTION
that can be set to Debug. That is what you need to do if you want to compile using the debug
flag, or use the Debug version. 
    - copy the file etc/config.sh/example/prefs.sh to your wanted etc level, and set the
Debug option in that file
    - If you only want to debug your own development, you can compile only that in debug mode, for example:
     Copy the icoFoam solver, re-name the executable name and location, and modify the first
    line in Make/options to:
    `EXE_INC = -O0 -fdefault-inline -ggdb3 -DFULLDEBUG \`
    With this you can debug only this application.