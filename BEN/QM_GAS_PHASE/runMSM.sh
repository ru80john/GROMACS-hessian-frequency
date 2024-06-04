#!/bin/bash

# 獲取腳本所在目錄的絕對路徑
script_dir="$(cd "$(dirname "$0")" && pwd)"
python_script='../../../Python_Modified_Seminario_Method/modified_Seminario_method.py'

#vibrational scaling vector
v_scaling_vec='0.957'

python $python_script $script_dir/ $script_dir/ $v_scaling_vec
