INSTRUCTIONS TO INSTALL THE NURBS TOOLBOX MEX-FILES IN MATLAB

- Copy all the files (the files, not the folder) in the folder of the NURBS toolbox (your_folder/nurbs/inst/)
- Start Matlab, and go to the same folder mentioned above.
- Run the script compile.m. This will compile the files, and add the folder to your path.

After compiling, check that Matlab is using the compiled functions by typing the command "which bspeval", for instance.

It is necessary to have a MEX-compiler to compile the MEX-functions. If you don't have a MEX-compiler, you will get some error like
 ??? Undefined function or method 'mex' for input arguments of type 'char'.
