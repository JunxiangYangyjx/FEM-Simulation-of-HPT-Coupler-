import json
import numpy as np
import os
import random
import pickle
import shutil
from utils_MAXWELL import *
from utils_HFSS import *
import re
import pandas as pd
from scipy.signal import argrelextrema
import numpy as np

def searching_boundary():
    w1_bound = [1.0 , 4.0]  # mm
    k1_bound = [0.85, 1.10]  
    k2_bound = [0.85, 1.10]  
    space_bound = [1, 5]  # mm
    n1_bound = [5, 15]
    n2_bound = [10, 25]
    
    return w1_bound, k1_bound, k2_bound, space_bound, n1_bound, n2_bound

def convert_to_config(w1, k1, k2, space, n1, n2, save_dir, index=""):
    config_file = "config_template.json" # a template config file
    with open(config_file, "r+") as f:
        kwargs = json.load(f)
    kwargs["w1"] = w1
    kwargs["space"] = space
    kwargs["wouts"] = [w1 / (k1 ** i) for i in range(0, n1)]
    kwargs["wins"] = [w1/ (k2 ** i) for i in range(0, n2)]
    kwargs["k1"] = k1
    kwargs["k2"] = k2
    kwargs["n1"] = n1
    kwargs["n2"] = n2
    kwargs["res_path"] = save_dir
    config_file = f'{save_dir}/config{index}.json'
    with open(config_file, 'w+') as f:
        json.dump(kwargs, f)
    return config_file

def get_random():
    w1_bound, k1_bound, k2_bound, space_bound, n1_bound, n2_bound = searching_boundary()
    w1 = random.uniform(*w1_bound)
    k1 = random.uniform(*k1_bound)
    k2 = random.uniform(*k2_bound)
    space = random.uniform(*space_bound)
    n1 = random.sample(list(range(n1_bound[0], n1_bound[1]+1)), 1)[0]
    n2 = random.sample(list(range(n2_bound[0], n2_bound[1]+1)), 1)[0]
    return w1, k1, k2, space, n1, n2

def check_rad(config_file):
    Flag = True
    with open(config_file, "r+") as f:
        kwargs = json.load(f)
    wins = kwargs["wins"]
    wouts = kwargs["wouts"] 
    radins1 = [kwargs["rout_out"] - wouts[0]]    
    radouts1 = [kwargs["rout_out"]]                     

    for i in range(1, len(wouts)):
        current_w = wouts[i]          
        previous_w = wouts[i - 1]     

        radin = radins1[-1] - kwargs["space"] - current_w
        radout = radouts1[-1] - kwargs["space"] - previous_w

        radins1.append(radin)
        radouts1.append(radout)
            
    radins2 = [kwargs["rout_in"] - wins[0]]   
    radouts2 = [kwargs["rout_in"]]  

    for i in range(1, len(wins)):
        current_w = wins[i]         
        previous_w = wins[i - 1]    

        radin = radins2[-1] - kwargs["space"] - current_w
        radout = radouts2[-1] - kwargs["space"] - previous_w

        radins2.append(radin)
        radouts2.append(radout)

    if radins1[-1] <= 90 :   # 外圈内径小于90mm（内圈外径）不合理
        Flag = False
        
    if radins2[-1] <= 10: # rin小于0不合理，给一个设定范围作为保护
        Flag = False
    S_out = 0 
    S_in = 0
    for i in range(kwargs["n1"]):
        S_out = S_out + radouts1[i]**2 - radins1[i]**2

    for i in range(kwargs["n2"]):
        S_in = radouts2[i]**2 - radins2[i]**2 + S_in
    
    if np.abs(S_out-S_in) > 400:
        Flag = False
    return Flag

def latin(num=200, restore=False):

    check_point_file = "./checkpoint.pickle"

    if restore:
        assert os.path.exists(check_point_file), \
            "the checkpoint does not exist, please rerun."

    bounds = [*searching_boundary()]
    is_cat = [False]*4 + [True]*2

    d = 1.0 / num
    D = 6
    save_dir = "./latin"

    if not os.path.exists(save_dir):
        os.makedirs(save_dir)

    if not restore:
        params = [[] for _ in range(num)]

        for i in range(D):
            if not is_cat[i]:
                temp = []
                for j in range(num):
                    tmp = np.random.uniform(low=j*d, high=(j+1)*d)
                    tmp = tmp*(bounds[i][1]-bounds[i][0]) + bounds[i][0]
                    temp.append(tmp)
            else:
                tmp = list(range(bounds[i][0], bounds[i][1]+1))
                np.random.shuffle(tmp)
                temp = [tmp[j % len(tmp)] for j in range(num)]

            np.random.shuffle(temp)

            for j, item in enumerate(temp):
                params[j].append(item)

        j_start = 0

    else:
        with open(check_point_file, "rb") as f:
            j_start, params = pickle.load(f)

    for j in range(j_start, len(params)):

        with open(check_point_file, "wb+") as f:
            pickle.dump([j, params], f)

        w1, k1, k2, space, n1, n2 = params[j]

        _j = 0
        while True:

            if _j != 0:
                w1 = random.uniform(*bounds[0])
                k1 = random.uniform(*bounds[1])
                k2 = random.uniform(*bounds[2])
                space = random.uniform(*bounds[3])
                n1 = random.randint(bounds[4][0], bounds[4][1])
                n2 = random.randint(bounds[5][0], bounds[5][1])

            config_file = convert_to_config(
                w1, k1, k2, space, n1, n2,
                save_dir,
                index=j
            )

            if not check_rad(config_file):
                print(f"Invalid configuration with radius: {config_file}")
                _j += 1
                continue

            result = evaluate_simulation(
                w1, k1, k2, space, n1, n2,
                save_dir=save_dir,
                index=j
            )

            del_cache(config_file)

            if not result["success"]:
                print(f"Simulation failed at index {j}")
                _j += 1
                continue

            print(f"Sample {j} finished.")
            break


def del_cache(file):
    with open(file, "r+") as f:
        kwargs = json.load(f)
    dir_cache = kwargs['project_path']+kwargs['project_id']+".aedtresults"
    if os.path.exists(dir_cache):
        shutil.rmtree(dir_cache)
        
# 把文件转移到特定文件夹中
def move_file(files, save_dir):
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)
    for file in files:
        with open(file, "r+") as f:
            kwargs = json.load(f)
        shutil.move(kwargs['project_path']+kwargs['project_id']+".aedt", save_dir+"/"+kwargs['project_id']+".aedt")


def evaluate_simulation(w1, k1, k2, space, n1, n2, save_dir=".", index=0, desktop_instance=None):
    desktop_created_inside = False
    try:
        if desktop_instance is None:
            desktop_instance = Desktop(
                specified_version="2023.1",
                non_graphical=False,
                close_on_exit=False,
                student_version=False
            )
            desktop_created_inside = True

        config_file = convert_to_config(w1, k1, k2, space, n1, n2, save_dir, index=index)

        _, success_maxwell = run_maxwell(config_file, desktop_instance=desktop_instance, index=index)

        if not success_maxwell:
            print("Maxwell failed.")
            return {
                "success": False,
                "index": index
            }

        _, success_hfss = run_hfss(config_file, desktop_instance=desktop_instance, index=index)

        if not success_hfss:
            print("HFSS failed.")
            return {
                "success": False,
                "index": index
            }

        return {
            "success": True,
            "index": index
        }

    except Exception as e:
        print(f"Evaluation error: {e}")
        return {
            "success": False,
            "index": index
        }

    finally:
        if desktop_created_inside and desktop_instance:
            desktop_instance.release_desktop(close_projects=True)
            


import os
import re
import json
import math
import numpy as np
from scipy.signal import argrelextrema


def aggregate_results(folder_path, output_name="all_results.json"):

    if not os.path.exists(folder_path):
        raise ValueError("Folder does not exist.")

    files = os.listdir(folder_path)
    indices = set()

    for file in files:
        match = re.search(r'index(\d+)', file)
        if match:
            indices.add(int(match.group(1)))

        match_cfg = re.search(r'config(\d+)', file)
        if match_cfg:
            indices.add(int(match_cfg.group(1)))

    indices = sorted(indices)
    results = {}

    for idx in indices:

        entry = {}

        # =====================================================
        # 1️⃣ CONFIG
        # =====================================================
        config_path = os.path.join(folder_path, f"config{idx}.json")

        if os.path.exists(config_path):
            with open(config_path, "r") as f:
                config_data = json.load(f)

            entry["geometry"] = {
                "w1": config_data.get("w1"),
                "space": config_data.get("space"),
                "k1": config_data.get("k1"),
                "k2": config_data.get("k2"),
                "n1": config_data.get("n1"),
                "n2": config_data.get("n2"),
                "rout_out": config_data.get("rout_out"),
                "rout_in": config_data.get("rout_in")
            }

        # =====================================================
        # 2️⃣ L MATRIX
        # =====================================================
        L_path = os.path.join(folder_path, f"index{idx}_L.csv")

        if os.path.exists(L_path):
            data = np.loadtxt(L_path, delimiter=",", skiprows=1)
            row = data if data.ndim == 1 else data[0]

            entry["L_matrix"] = {
                "L_inner": float(row[1]),
                "L_outer": float(row[2]),
                "M_samePCB": float(row[3]),
                "M_inner_inner": float(row[4]),
                "M_cross": float(row[5]),
                "M_outer_outer": float(row[6])
            }

        # =====================================================
        # 3️⃣ SINGLE PCB
        # =====================================================
        single_path = os.path.join(folder_path, f"index{idx}-singlePCB.csv")

        if os.path.exists(single_path):

            data = np.loadtxt(single_path, delimiter=",", skiprows=1)
            freq = data[:, 0]      # MHz
            dB20 = data[:, 1]

            min_indices = argrelextrema(dB20, np.less)[0]

            if len(min_indices) > 0:
                series_index = min_indices[0]
            else:
                series_index = np.argmin(dB20)

            entry["singlePCB"] = {
                "series_freq_MHz": float(freq[series_index]),
                "series_dB20": float(dB20[series_index])
            }

        # =====================================================
        # 4️⃣ DOUBLE PCB（扁平结构）
        # =====================================================
        double_path = os.path.join(folder_path, f"index{idx}-doublePCB.csv")

        if os.path.exists(double_path):

            data = np.loadtxt(double_path, delimiter=",", skiprows=1)
            freq = data[:, 0]
            dB20 = data[:, 1]

            # ---- 串联谐振 ----
            min_indices = argrelextrema(dB20, np.less)[0]
            if len(min_indices) > 0:
                series_index = min_indices[0]
            else:
                series_index = np.argmin(dB20)

            f_series = freq[series_index]

            # ---- 并联谐振（只在串联之后找）----
            mask = freq > f_series

            if np.any(mask):
                freq_after = freq[mask]
                dB_after = dB20[mask]

                max_indices = argrelextrema(dB_after, np.greater)[0]

                if len(max_indices) > 0:
                    # 找到局部极大值
                    parallel_index = np.where(mask)[0][max_indices[0]]
                else:
                    # 没找到极大值 → 设为扫描上限
                    parallel_index = len(freq) - 1
            else:
                # 极端情况（串联已经是最大频率）
                parallel_index = len(freq) - 1

            entry["doublePCB"] = {
                "series_freq_MHz": float(freq[series_index]),
                "series_dB20": float(dB20[series_index]),
                "parallel_freq_MHz": float(freq[parallel_index]),
                "parallel_dB20": float(dB20[parallel_index])
            }

        # =====================================================
        # 5️⃣ Calculate（整合进来）
        # =====================================================
        if all(k in entry for k in ["geometry", "L_matrix", "singlePCB"]):

            geom = entry["geometry"]
            Lmat = entry["L_matrix"]
            single = entry["singlePCB"]

            # ---- R_total from dB20 ----
            dB20 = single["series_dB20"]
            R_total = 10**(dB20 / 20)

            # ---- coil lengths ----
            w = geom["w1"]
            space = geom["space"]
            pitch = w + space

            n1 = geom["n1"]
            rout_out = geom["rout_out"]
            l_out = 2*math.pi*(n1*rout_out - (n1-1)*n1/2*pitch)

            n2 = geom["n2"]
            rout_in = geom["rout_in"]
            l_in = 2*math.pi*(n2*rout_in - (n2-1)*n2/2*pitch)

            k_ratio = l_out / l_in

            R_outer = k_ratio/(1+k_ratio) * R_total
            R_inner = 1/(1+k_ratio) * R_total

            # ---- Cs calculation ----
            f_s = single["series_freq_MHz"] * 1e6   # MHz → Hz

            L_outer = Lmat["L_outer"] * 1e-6
            L_inner = Lmat["L_inner"] * 1e-6
            M = Lmat["M_samePCB"] * 1e-6

            L_eq = L_outer + L_inner + 2*M
            omega = 2*math.pi*f_s
            Cs = 1/(omega**2 * L_eq)

            entry["Calculate"] = {
                "R_outer_est": R_outer,
                "R_inner_est": R_inner,
                "R_total_est": R_total,
                "Cs_est_pF": Cs * 1e12
            }

        results[f"index{idx}"] = entry

    output_path = os.path.join(folder_path, output_name)

    with open(output_path, "w") as f:
        json.dump(results, f, indent=4)

    print(f"Aggregation complete → {output_path}")

    return results