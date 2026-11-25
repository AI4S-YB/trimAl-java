@echo off
rem Launcher of the trimAl v1.5.1 Java port (GPL-3.0-or-later, see LICENSE and NOTICE.md).
setlocal
if defined JAVA_HOME (
  "%JAVA_HOME%\bin\java.exe" -cp "%~dp0target\trimal.jar" trimal.StatAl %*
) else (
  java -XX:TieredStopAtLevel=1 -cp "%~dp0target\trimal.jar" trimal.StatAl %*
)
exit /b %errorlevel%
