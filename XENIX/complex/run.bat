@echo off
echo. | date > runtime
echo. | time >> runtime
call %1
echo. >> runtime
echo. | date >> runtime
echo. | time >> runtime
