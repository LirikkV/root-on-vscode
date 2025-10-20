How to run macros in VS code with debugger?
we need Makefile and Macro.c

1) Into Makefile we need DEBUGFLAG = -ggdb -g3 //this will compile Macro.C with debug information
2) Compilated file after "make -j4" needs to be already exist to recompile it 
    so for first time you need to compile it by hands with "make -j4"