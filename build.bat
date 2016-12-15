@echo off

SET OPTS=/W4 /wd4310 /wd4100 /wd4201 /wd4505 /wd4996 /wd4127 /wd4510 /wd4512 /wd4610 /wd4390 /WX
SET OPTS=%OPTS% /GR- /EHa- /nologo /FC
SET DEBUG=/Zi
SET OPTIM=/O2

pushd ..\build
REM cl %OPTS% %DEBUG% ..\code\wfc.c /Fewfc
cl %OPTS% %OPTIM% ..\code\wfc.c /Fewfc
popd
