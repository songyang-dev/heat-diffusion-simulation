@echo off
if "%1"=="release" (
    cmake .. -G"NMake Makefiles" -DCMAKE_BUILD_TYPE=Release
) else (
    cmake .. -G"NMake Makefiles" -DCMAKE_BUILD_TYPE=Debug
)
nmake