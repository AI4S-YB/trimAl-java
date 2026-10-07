@echo off
setlocal
if defined JAVA_HOME (
  "%JAVA_HOME%\bin\java.exe" -cp "%~dp0target\trimal.jar" trimal.ReadAl %*
) else (
  java -XX:TieredStopAtLevel=1 -cp "%~dp0target\trimal.jar" trimal.ReadAl %*
)
exit /b %errorlevel%
