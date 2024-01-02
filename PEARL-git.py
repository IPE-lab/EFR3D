import numpy as np

from tool_func import *
import numpy as np

from utils import *
from plane_detection import *
# matplotlib.use("Agg")  # Force matplotlib to not use any Xwindows backend.
import gco  # noqa: E402
import matplotlib.pyplot as plt 
from scipy.spatial import Delaunay
import mpl_toolkits.mplot3d as a3


#生成初始标签
def initial_model_gen(points, min_ratio=0.00, threshold=0.01, iterations=1000):
    planes_fun, planes_idx = DetectMultiPlanes(points, min_ratio=min_ratio, threshold=threshold, iterations=iterations, idx=True)
    ini_inliers = points[list(set(chain.from_iterable(planes_idx)))]
    ini_outliers = points[list(set(range(len(points))) - set(chain.from_iterable(planes_idx)))]
    # print(len(ini_inliers), len(points), len(ini_outliers))
    # print(datacost.size, points.size, datacost)
    ini_labels = []
    for i_ in range(len(planes_idx)):
        print(len(planes_idx[i_]))
        for j_ in range(len(planes_idx[i_])):
            ini_labels.append(i_)
        
    return planes_fun, planes_idx, ini_inliers, ini_outliers, ini_labels
    
#能量最小化
def energy_min(planes_fun, inliers, outliers,ini_labels, ksi=5, algorithm='expansion', n_inter=-1):
    # def energy_opt(datcost, smooth, edges, edge_weights,algorithm='expansion', n_inter=-1):
    
    initialize = 0
    extend_points = np.ones((len(inliers), 4))
    #nx4矩阵 xi yi zi 1
    extend_points[:, 0:3] = inliers
    #4xm矩阵 ai bi ci di.T
    planes_para = np.array(planes_fun).T
    # print(extend_points)
    # print(planes_para)
    #计算当前datacost
    datacost = abs(np.dot(extend_points, planes_para))
    
    #三角剖分
    tri = Delaunay(inliers)
    #获得四面体顶点 每行 按从小到大排列
    vertices = np.sort(tri.simplices)
    #权重6列 代表四面6个边 01 02 03 12 13 23
    weights = np.ones((len(vertices), 6))
    # w=exp(-||p-q||^2/ksi^2)
    w0_1 = inliers[vertices[:, 0]] - inliers[vertices[:, 1]]
    weights[:, 0] = np.exp(- np.linalg.norm(w0_1, axis=1) ** 2 / ksi**2)
    w0_2 = inliers[vertices[:, 0]] - inliers[vertices[:, 2]]
    weights[:, 1] = np.exp(- np.linalg.norm(w0_2, axis=1) ** 2 / ksi**2)
    w0_3 = inliers[vertices[:, 0]] - inliers[vertices[:, 3]]
    weights[:, 2] = np.exp(- np.linalg.norm(w0_3, axis=1) ** 2 / ksi**2)
    w1_2 = inliers[vertices[:, 1]] - inliers[vertices[:, 2]]
    weights[:, 3] = np.exp(- np.linalg.norm(w1_2, axis=1) ** 2 / ksi**2)
    w1_3 = inliers[vertices[:, 1]] - inliers[vertices[:, 3]]
    weights[:, 4] = np.exp(- np.linalg.norm(w1_3, axis=1) ** 2 / ksi**2)
    w2_3 = inliers[vertices[:, 2]] - inliers[vertices[:, 3]]
    weights[:, 5] = np.exp(- np.linalg.norm(w2_3, axis=1) ** 2 / ksi**2)
    
    #边的左顶点 要求做顶点的索引小于右边的
    edge0 = np.ones((len(vertices), 6))
    edge0[:, 0] = vertices[:, 0]
    edge0[:, 1] = vertices[:, 0]
    edge0[:, 2] = vertices[:, 0]
    edge0[:, 3] = vertices[:, 1]
    edge0[:, 4] = vertices[:, 1]
    edge0[:, 5] = vertices[:, 2]
    #边的左顶点
    edge1 = np.ones((len(vertices), 6))
    edge1[:, 0] = vertices[:, 1]
    edge1[:, 1] = vertices[:, 2]
    edge1[:, 2] = vertices[:, 3]
    edge1[:, 3] = vertices[:, 2]
    edge1[:, 4] = vertices[:, 3]
    edge1[:, 5] = vertices[:, 3]
    
    # 不同类为1 同类为0 nlabelxnlabel
    smooth = smooth = (1 - np.eye(len(planes_fun))).astype(np.float64)
    # print(datacost, ini_labels)
    # print(vertices, edge0.T.flatten(), edge1.T.flatten())
    # print(weights.T.flatten())
    
    #设定datacost neighbour smoothcost 
    # 防止内存溢出 缩放矩阵
    max_arr = max(np.abs(datacost).max(), np.abs(weights).max() * smooth.max())
    down_weight_factor = max_arr + 1e-10
    print('max_arr', max_arr)
    gc = gco.GCO()
    
    return labels

#融合平面
def merge(planes):
      
      plane, point_idx = map(list, zip(*planes))
      # print(res1, res2)

      dbscan = DBSCAN(eps=0.001, min_samples=1)
      dbscan.fit_predict(plane)
      label_pred = dbscan.labels_
      
      plane = np.array(plane)
      point_idx = np.array(point_idx, dtype=object)
      
      planes_merged = []
      # 看看有多少平面离的比较近 聚类
      for i in list(set(label_pred)):
            # print(i, list(set(label_pred)))
            # 对于每一簇平面 选取第一个为最终平面
            plane_i = plane[label_pred == i][0]
            # 将该簇平面附近的点融合为一个列表
            point_idx_i = point_idx[label_pred == i]
            point_idx_i = list(set(chain.from_iterable(list(point_idx_i))))
            planes_merged.append([list(plane_i), point_idx_i])
      return planes_merged
    
#outlier recheck    
def ocheck(planes_fun, inliers, outliers,ini_labels):
    jj = 0
#迭代直至能量不减少
def energy_opt(planes_fun, inliers, outliers,ini_labels):
    tt = 0

