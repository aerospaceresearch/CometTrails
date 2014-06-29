This directory should contain the contents of the downloaded cspice archive. Necessary are the following subdirectories:

/include/
/lib/

To re-built cspice with statically linked libraries (and be able to use the CMake-option LinkStatic), do the following:

- go to /src/cspice/mkproduct.bat
- replace both occurences of "set cl= /c" with "set cl= /c /MT".
- go to /src/csupport/mkproduct.bat
- replace one occurence of "set cl= /c" with "set cl= /c /MT".
- open the VS2013 x86 Native Tools Command Prompt
- run makeall.bat