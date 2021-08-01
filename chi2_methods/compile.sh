#!/bin/sh

# the name of the code and exectutalbe
code=chi2.cxx
exec=chi2.exe

rm -f $exec
# Get compiler and flags used by root, needed a small "hack" to work at feynman (lib64)
cxx="$(root-config --cxx)"
flagsNlibs="-L/usr/local/lib64 $(root-config --cflags --glibs)"

# print some info to the screen
echo "Will compile $code using the $cxx compiler"
echo "Executable to be produced: $exec"

# Now compile!
$cxx $flagsNlibs -o $exec $code && echo "Now run by ./$exec" || echo "Compilation failed"

####
