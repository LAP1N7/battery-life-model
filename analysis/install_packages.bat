@echo off
chcp 65001 >nul
cd /d "%~dp0"
rem PyCharm 프로젝트에 지정한 Python(3.14)에 패키지를 한 번에 설치합니다.
set PY=C:\Users\qkrdm\AppData\Local\Python\pythoncore-3.14-64\python.exe
if not exist "%PY%" set PY=py

echo [1/3] pip 업데이트...
"%PY%" -m pip install --upgrade pip

echo [2/3] 필수 패키지 설치 (numpy, pandas, scipy, scikit-learn, matplotlib)...
"%PY%" -m pip install -r requirements.txt || goto err

echo [3/3] 선택 패키지 설치 (shap: ML 해석용, 실패해도 분석은 진행 가능)...
"%PY%" -m pip install shap || echo   shap 설치 실패 - 지금 단계에서는 필요 없으니 넘어가도 됩니다.

echo.
echo 설치 확인:
"%PY%" -c "import numpy, pandas, scipy, sklearn, matplotlib; print('필수 패키지 OK')" || goto err
echo 완료. PyCharm으로 돌아가서 실행하세요.
pause
exit /b 0
:err
echo 오류가 났습니다. 위 메시지를 복사해서 Claude에게 보내 주세요.
pause
