YGSL
====

YGSL is a Yorick plug-in to bring some special functions of the GSL
(GNU Scientific Library) into Yorick.


Installation
------------

In short, building and installing the plug-in can be as quick as:
````{.sh}
cd $BUILD_DIR
$SRC_DIR/configure
make
make install
````
where `$BUILD_DIR` is the build directory (at your convenience) and
`$SRC_DIR` is the source directory of the plug-in code.  The build and
source directories can be the same in which case, call `./configure` to
configure for building.

If the plug-in has been properly installed, it is sufficient to use any
function of YGSL to automatically load the plug-in.  You may force the
loading of the plug-in by something like:
````{.sh}
#include "gsl.i"
````
or
````{.sh}
require, "gsl.i";
````
in your code.

More detailled installation explanations are given below.

1. You must have Yorick and the GSL (GNU Scientific Library) installed
   on your machine.  (See the *"Links"* section below.)

2. Unpack the software code somewhere or clone the Git repository.

3. Configure for compilation.  There are two possibilities:

   For an **in-place build**, go to the source directory, say `$SRC_DIR`, of
   the plug-in code and run the configuration script:
   ````{.sh}
   cd $SRC_DIR
   ./configure
   ````
   To see the configuration options, call:
   ````{.sh}
   ./configure --help
   ````

   To compile in a **different build directory**, say `$BUILD_DIR`, create the
   build directory, go to the build directory and run the configuration
   script:
   ````{.sh}
   mkdir -p $BUILD_DIR
   cd $BUILD_DIR
   $SRC_DIR/configure
   ````
   where `$SRC_DIR` is the path to the source directory of the plug-in code.
   To see the configuration options, call:
   ````{.sh}
   $SRC_DIR/configure --help
   ````

4. Compile the code:
   ````{.sh}
   make
   ````

4. Install the plug-in in Yorick directories:
   ````{.sh}
   make install
   ````


License
-------

YGSL is open source sofware released under the CeCILL-C license
<http://www.cecill.info/index.en.html>.


History
-------

YGSL was a component of Yeti (a group of Yorick plugins), it is now a
standalone plug-in. You can find more informations about Yeti at
<http://www-obs.univ-lyon1.fr/labo/perso/eric.thiebaut/yeti.html>.


Links
-----

 * Yorick: <http://github.com/LLNL/yorick/>;
 * GSL (GNU Scientific Library): <http://www.gnu.org/software/gsl/>;
