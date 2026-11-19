@echo off
chcp 65001 >nul
cd /d "%~dp0"
echo [1/4] EIS fitting (several minutes)...
python step1_extract.py || goto err
echo [2/4] path test...
python step2_paths.py || goto err
echo [3/4] ML-1...
python step3_ml1.py || goto err
echo [4/4] plot...
python step4_plot.py || goto err
echo Done. See out\results.png
pause
exit /b 0
:err
echo Error occurred. Copy the message above and send it to Claude.
pause
