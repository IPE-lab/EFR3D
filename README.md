# **An Energy-based 3D Fracture Reconstruction method (EFR3D)**
## Running screenshots show
- **results**
  - <img src="img/results.png" width="400" />
***
## Paper Support
- Original information: Energy-Based Fracture Network Reconstruction of Shale Gas Reservoir
- Recruitment Journal: Asia Pacific Unconventional Resources Symposium held in Brisbane, Australia on 14 – 15 November 2023. 
- Original DOI: https://doi.org/10.2118/217311-MS
***
## Description of the project
Microseismic monitoring is a commonly used technique in characterizing hydraulic fractures. However, fracture network reconstruction remains challenging due to the heterogeneity and complex stress field of shale gas reservoirs. An Energy-based 3D Fracture Reconstruction method (EFR3D) is proposed to derive the complex fracture network from microseismic data in a shale gas reservoir. The EFR3D method mainly combines the Propose Expand and Re-Estimate Labels algorithm (PEARL), the Density-Based Spatial Clustering of Applications with Noise algorithm (DBSCAN), and the Alpha-shape algorithm to detect the fracture orientation and fracture shape. This method formulates fracture orientation detection as an energy minimization task to improve reconstruction accuracy. The effectiveness of the proposed method was evaluated by performing a verification procedure against different fracture numbers, fracture orientations, and fracture scales using the Monte Carlo simulation. The results show that the proposed method has good adaptability and high accuracy in various fracture configurations. Furthermore, a field application of six horizontal wells located in the southern Sichuan basin, southwest China, was conducted to illustrate the robustness and practicability of the proposed method. The proposed method can serve as a practical and reliable approach to characterize hydraulic fractures.
***
## Functions of the project
PEARL git. py is used for mesh reconstruction and contains 5 functions, of which functions 3 and 5 are not publicly disclosed.
1. Preliminary identification of crack surface and generation of initial labels
    Initial_ Model_ Gen (points, min-ratio=0.00, threshold=0.01, iterations=1000):
2. Minimize the energy of the current label
    Energy_ Min (planes'fun, inliers, outliers, ini'labels, ksi=5, algorithm='expansion ', n-inter=-1):
3. Check if the current plane contains enough interior points
    Ocheck (planes'fun, inliers, outliers, ini'labels):
4. Fusion of insufficient inner points on the plane and re marking as outer points
    Merge (planes):
5. Re identify external points, calculate energy minimization based on the new and previous labels, and iterate until the energy does not decrease:
    Energy_ Opt (planes'fun, inliers, outliers, ini'labels):
***
## The operating environment of the project
-	Python == 3.7.0
-	pandas == 1.1.5
-	numpy == 1.21.6
-	matplotlib == 3.5.1
-	scipy == 1.1.2
-	gco == 1.1.2
-	plane_detection == 1.0.0
***
## How to use the project
#### 1、Read the coordinates of microseismic events.

#### 2、Set coordinates as points variables.

#### 3、Call functions 1 to 5 in sequence.
