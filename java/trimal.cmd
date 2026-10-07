@echo off
setlocal
if defined JAVA_HOME (
  "%JAVA_HOME%\bin\java.exe" -cp "%~dp0target\trimal.jar" trimal.TrimAl %*
) else (
  java -XX:TieredStopAtLevel=1 -cp "%~dp0target\trimal.jar" trimal.TrimAl %*
)
exit /b %errorlevel%
