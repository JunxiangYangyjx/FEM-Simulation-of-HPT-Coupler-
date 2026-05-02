from utils_MAXWELL import *
from utils_HFSS import *
from optimizer import *
# from dict_search import dict_search

if __name__ == "__main__":  # 确保这段代码只在直接运行脚本时执行，而不是在被导入时执行
    
    w1, k1, k2, space, n1, n2, save_dir, index = 2, 1.11, 1.05, 2.5, 9, 20, ".", 0
    # w1, k1, k2, space, n1, n2, save_dir, index = 2.41, 0.91, 1.01, 1.14, 6, 19, ".", 0
    result = evaluate_simulation(w1, k1, k2, space, n1, n2, save_dir=save_dir, index=index)
    # print(result)
    
    # print(cir_pcb.parse_results())
    
    # Stage 1
    # num = 100
    # latin(num=num, restore=False)
    # latin(num=num, restore=True)
    # aggregate_results(folder_path="./temp")
    # Stage 2
    # aug_num = 100
    # prev_result_path = "latin"
    # dict_search(aug_num, prev_result_path, num=20, threshold=0.1, weight=10.0, restore=False)
    # dict_search(aug_num, prev_result_path, num=20, threshold=0.1, weight=10.0, restore=True)
    