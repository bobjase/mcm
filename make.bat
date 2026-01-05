@echo off
echo Compiling...
del mcm.exe
"C:\temp_mingw\mingw64\bin\g++.exe" -pthread -DNDEBUG -O3 -fomit-frame-pointer -march=x86-64 -mtune=generic -std=c++17 -Wno-deprecated-declarations -o mcm Archive.cpp Huffman.cpp MCM.cpp Memory.cpp Util.cpp Compressor.cpp File.cpp LZ.cpp Tests.cpp
if errorlevel 1 (
  echo Compilation failed.
) else (
  echo mcm.exe created.
  .\mcm.exe --phase1-profile testFiles/large.txt
)