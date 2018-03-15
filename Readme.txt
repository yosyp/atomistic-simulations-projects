MSE6270_MD code, C++ version
Details about the code are available in the online documentation (you should get a link from the professor).

To make an executable, go to src, type "module load intel", "make", and "cd ..". 
As a result you will get two executables "mse627" and "CG", MD code and crystal generator.
After that copy "mse627" to hw3 directory, where you can run it.

For people using Microsoft Visual Studio, open "MSE6270_MD.vcxproj" project. When compile do not forget 
to switch to "Release". The executable will be in Release directory. You can copy it to a target directory 
(hw3) to run. Alternatively, you can put input files to the project directory to run your code from 
Visual Studio.

The code also includes some useful files: "Convert_alloy_tab.cpp" and "main_direct.cpp". 
The first one can be used to build a code for conversion of *.alloy EAM tables to *tab tables.
The second file shows how to use the code as a library, without building MD system based on 
"md.rc" and "md.input" files. It provides mode flexibility but requires more advanced programming skills.
You can compile these files by adding corresponding modifications to the Makefile.

Warning: The code may require understanding of object oriented concept and some experience in C++/Java. 
For beginners it is recommended to use Fortran version of MD code, which is simpler. 