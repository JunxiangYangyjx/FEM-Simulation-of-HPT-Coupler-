# Import all dependent modules
import numpy as np
import os
import re
import json
import random
import pickle
from optimizer import get_random, convert_to_config, del_cache, obj_func, check_rad
from utils_MAXWELL import HPT, run
from sklearn.manifold import TSNE


# Function: get_x
# Description: get the geometric parameters
# Arguments: kwargs
# Returns: various geometric parameters
def get_x(kwargs):
    return (kwargs['w1'], kwargs['k1'], kwargs['k2'],
            kwargs['n1'], kwargs['n2'], kwargs['space'])


# Function: get_y
# Description: get the objective values
# Arguments: results
# Returns: various objective values
def get_y(results):
    return (results['SRFout1'], results['SRFin1'], results['R1'], results['R2'], results['obj'])


# Function: load_all
# Description: load all results from the folder
# Arguments: path
# Returns: a dictionary all_results summarizing all simulations
# 得到一个包含所有仿真配置和结果的字典all_results
def load_all(path):
    all_results = {}
    files = os.listdir(path)    # 获取指定路径下的所有文件名列表    
    f = re.compile(r"config(\d+).json") # 使用正则表达式匹配文件名，寻找形如 "configX.json" 的文件
    indices = [int(f.match(file)[1]) for file in files if f.match(file)]    # 遍历文件列表，提取所有匹配的文件的数字索引，indices是一个包含所有找到的配置文件索引的列表
    for idx in indices:
        file_full = os.path.join(path, f"config{idx}.json") # 构造配置文件的完整路径
        files_full = [os.path.join(path, file) for file in [f"index{idx}-reZout.csv",     # 构造当前索引对应的结果文件列表
                                                            f"index{idx}-reZin.csv",
                                                            f"index{idx}-dB20Z.csv",]]
        with open(file_full, "r") as f: # 打开当前索引对应的配置文件，读取其中的内容并将其加载到 kwargs 字典中。
            kwargs = json.load(f)
        cir_pcb = HPT(kwargs, None)   # 创建一个 IPTCoil 对象，用于处理当前配置的仿真结果。
        cir_pcb.report_files.extend(files_full) # 将结果文件的完整路径添加到 cir_pcb 对象的 report_files 属性中，以便在后续的解析中使用这些文件。
        cir_pcb.parse_results()                 # 解析结果文件中的数据，并将结果存储在 cir_pcb.parsed_results 中。
        # obj_value = obj_func(cir_pcb.parsed_results)
        # all_results[idx] = {"kwargs": kwargs, "results": cir_pcb.parsed_results, "obj": obj_value}    # 将当前索引对应的配置和解析结果存储在 all_results 字典中
        all_results[idx] = {"kwargs": kwargs, "results": cir_pcb.parsed_results}
        del cir_pcb # 清除变量
    return all_results


# Function: t_sne_dist
# Description: convert high-dimensional data to 2-dimensional space using t-SNE (reduce sparsity)
# Arguments: x
# Returns: x_transformed
# 目前的 t-SNE 使用默认参数，有时可能无法获得最佳结果。可以考虑在函数中添加更多参数，以控制 t-SNE 的行为，例如 perplexity、learning_rate 和 n_iter 等。
# 用t-sne给x降维
def t_sne_dist(x):
    n_samples = x.shape[0]
    perplexity_value = min(30, n_samples - 1)
    tsne_sim_x = TSNE(n_components=2, perplexity = perplexity_value)   # 用TSNE讲指定数据降维到二维空间
    x_transformed = tsne_sim_x.fit_transform(x) # 对数据 x 执行 t-SNE 降维，并返回降维后的数据
    return x_transformed


# Function: dict_search
# Description: using dictionary-based search for an improved Q with the required self-resonant frequency
# Arguments: aug_num, target_SRF, prev_result_path, 
#                num, threshold, weight,
#                check_point_file, restore
# Returns: results
# 用字典检索方法获得一个接近SRF的高Q值
def dict_search(aug_num, prev_result_path, 
                num=20, threshold=0.1, weight=10.0,
                check_point_file = "./checkpoint-dict.pickle", restore=False):
    if restore:
        assert os.path.exists(check_point_file), "the checkpoint does not exist, please rerun"
    
    save_dir = "./latin/dict_search"
    if not os.path.exists(save_dir):    # 创建文件夹
        os.makedirs(save_dir)
        
    all_results = load_all(prev_result_path)    # 加载之前已经保存的仿真结果
    selected_data = [(get_x(all_results[key]["kwargs"]),                            
                  get_y(all_results[key]["results"])) for key in all_results]       # selected_data 是一个列表，列表中的每个元素是一个元组，包含两个部分：仿真参数的 x 值和仿真结果的 y 值
    selected_data = sorted(selected_data, key=lambda x: x[1][4], reverse=True)      # obj最大的排在前面
    dict_search_data = np.concatenate([np.concatenate(item, axis=0)[None]           # 把get_x 和get_y整合成的一个元祖selected_data再转化为二维数组(K,8)，其中有几个仿真结果就有几行，每一行都包括八个元素（输入和结果）
                                       for item in selected_data], axis=0)
    dict_path = "dict_search_contents.csv"
    if not os.path.exists(dict_path):
        np.savetxt(dict_path, dict_search_data, delimiter=",")
    # 这两行代码将 dict_search_data 按列划分为两部分
    _x = dict_search_data[:, :6]    # 前6行是参数，相当于输入
    _y = dict_search_data[:, 6:]    # 后面的是lable
    results = {key:{"obj": None, "all": None} for key in range(aug_num)}
    if restore:
        with open(check_point_file, "rb") as f:
            j_start = pickle.load(f)[0]
    else:
        j_start = 0 # 没有中断的话默认是0
    
    for j in range(j_start, aug_num):
        with open(check_point_file, "wb+") as f:
            pickle.dump([j], f)
        
        _j = 0
        while True:
            #################################################
            # get a random param set based on dictionary search
            #################################################
            while True:
                new_x = np.zeros((0, 6))    # 创建一个形状为 (0, 6) 的 NumPy 数组
                all_configs = []
                for _ in range(num):
                    w1, k1, k2, space, n1, n2 = get_random()
                    all_configs.append((w1, k1, k2, space, n1, n2))
                    file = convert_to_config(w1, k1, k2, space, n1, n2, ".", index="_tmp") # config_tmp.json
                    with open(file, "r") as f:
                        kwargs = json.load(f)
                        _new_x = np.array([get_x(kwargs)])  # 使用 get_x(kwargs) 函数从配置中提取特征向量，并将其转换为 NumPy 数组（形状为 (1, 6)）。
                        new_x = np.concatenate([new_x, _new_x], axis=0) # 使用 np.concatenate 将新的特征向量 _new_x 添加到 new_x 数组中，以构建本次循环中生成的所有新特征向量的集合。
                x = np.concatenate([_x, new_x], axis=0) # 将之前已有的特征数据 _x 与本次生成的特征 new_x 进行拼接，更新为新的 x 数组。
                new_y = np.zeros((num, _y[-1].shape[0]))    # 创建一个num行，5列的全0数组
                y = np.concatenate([_y, new_y], axis=0)
                idx = np.abs(y[:, 0]-y[:, 0]) <= threshold # 计算 y 数组中第一个特征SRF1与第二个特征SRF2之间的绝对差值，并判断哪些样本满足条件：差值小于等于指定的阈值 threshold = 0.1。
                inputs_sim_transformed = t_sne_dist(x)
                dist = np.sqrt(((inputs_sim_transformed[-num:][:, None]-inputs_sim_transformed[idx])**2).sum(axis=-1))  # 计算新生成的样本与满足条件的样本在降维空间中的欧氏距离。
                if np.any(dist < threshold):    # 如果找到一个满足条件的样本，使用 break 退出 while True 循环
                    break
            
            min_dist = dist.min(axis=1) # 计算最小距离
            min_dist[min_dist == 0] = 1e-10  # 避免除以零的错误
            prob = (weight/min_dist).cumsum()/(weight/min_dist).sum()   # 用一个常数 weight 除以每个样本的最小距离。距离越小，结果越大，意味着更靠近其他样本的样本将具有更高的选择概率。并计算累计和以及归一化
            rand_idx = np.where(prob>=random.random())[0][0] # roulette method    找到累积概率 prob 中第一个大于等于随机数的位置。这个位置对应的是选中的样本索引。
            # rand_idx = np.argmin(min_dist) # deterministic
            #################################################
            
            w1, k1, k2, space, n1, n2 = all_configs[rand_idx]
            config_file = convert_to_config(w1, k1, k2, space, n1, n2, save_dir, index=j)
            if not check_rad(config_file):
                print(f"Invalid configuration with radius: {config_file}")
                _j += 1
                continue
            else:
                print(config_file, j)
            cir_pcb, if_success = run(config_file, index=j)
            del_cache(config_file) # EXTENSION
            if not if_success:
                _j += 1
                print(cir_pcb.error_log)
                continue
            # obj_value = obj_func(cir_pcb.parsed_results) # EXTENTION
            # results[j]["obj"] = obj_value
            # results[j]["all"] = cir_pcb.parsed_results
            results[j] = cir_pcb.parsed_results
            break
        dict_search_data = update_dict(dict_path, cir_pcb.kwargs, cir_pcb.parsed_results)
    with open(f"{save_dir}/all_results.json", "w") as f:
        json.dump(results, f)
    return results


# Function: update_dict
# Description: update the stored dictionary after each FEM simulation
# Arguments: dict_path, kwargs, results
# Returns: dict_search_data
# 更新字典
def update_dict(dict_path, kwargs, results):
    dict_search_data = np.loadtxt(dict_path, delimiter=",")                 # 加载现有数据
    new_xy = np.array([get_x(kwargs)+get_y(results)])                       # 将组合后的新样本转化为二维 NumPy 数组
    dict_search_data = np.concatenate((dict_search_data, new_xy), axis=0)   # 将新的样本 new_xy 追加到现有数据 dict_search_data 的末尾，沿着行的方向进行拼接。
    np.savetxt(dict_path, dict_search_data, delimiter=",")                  # 保存更新后的数据
    return dict_search_data
