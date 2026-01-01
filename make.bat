@echo off
echo Compiling...
C:\msys64\ucrt64\bin\g++ -pthread -DNDEBUG -O3 -fomit-frame-pointer -march=x86-64 -mtune=generic -flto -std=c++14 -Wno-deprecated-declarations -o mcm Archive.cpp Huffman.cpp MCM.cpp Memory.cpp Util.cpp Compressor.cpp File.cpp LZ.cpp Tests.cpp
if errorlevel 1 (
  echo Compilation failed.
) else (
  echo mcm.exe created.
)