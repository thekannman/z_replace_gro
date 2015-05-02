This program allows the replacement of molecules in the .gro files used by the GROMACS simulation pacakge. 

THIS PROGRAM HAS NOT BEEN THOROUGHLY TESTED!!!

Current limitations of the program include:
* Only monatomic molecules can be added, but polyatomic molecules can be removed.

The following libraries are required:
* The Boost program_options library.
* The Armadillo matrix library.
* The xdrfile library for reading/writing GROMACS .xtc and .trr files.
