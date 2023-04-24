import os
import numpy as np
from pathlib import Path
from numpy import sin, cos
from scipy.spatial.transform import Rotation
import cv2
from collections import OrderedDict
def read_qinpose(qinpose_path):
    lines = [l.strip() for l in  qinpose_path.read_text().splitlines()]
    num_image = int(lines[0])
    poses = OrderedDict()
    for li, l in enumerate(lines[1:]):
        toks = l.split(' ')
        if len(toks)!=14:
            continue
        filename = toks[0]
        f = float(toks[1])
        cx = float(toks[2])
        cy = float(toks[3])
        pxsz_x = float(toks[4])
        pxsz_y = float(toks[5])
        width = int(toks[6])
        height = int(toks[7])
        X, Y, Z= float(toks[8]), float(toks[9]),float(toks[10])
        omega,phi,kappa = float(toks[11]),float(toks[12]),float(toks[13])
        id = os.path.splitext(filename)[0]

        
        K = np.array([f / pxsz_x, 0, float(width-1) / 2.0 + cx / pxsz_x, 0, f / pxsz_y, float(height-1) / 2.0 - cy / pxsz_y, 0, 0, 1]).reshape(3,3)
        C = np.array([X,Y,Z])
        R = np.array([cos(kappa) * cos(phi), cos(omega) * sin(kappa) + cos(kappa) * sin(omega) * sin(phi), sin(kappa) * sin(omega) - cos(kappa) * cos(omega) * sin(phi), cos(phi) * sin(kappa),
              sin(kappa) * sin(omega) * sin(phi) - cos(kappa) * cos(omega), -cos(kappa) * sin(omega) - cos(omega) * sin(kappa) * sin(phi), -sin(phi), cos(phi) * sin(omega), -cos(omega) * cos(phi)]).reshape(3,3)
        poses[id]=dict(filename=filename, idx=li,
                      f=f,cx=cx,cy=cy,pxsz_x=pxsz_x, pxsz_y=pxsz_y,
                      width=width, height=height, X=X, Y=Y,Z=Z,
                      omega=omega,phi=phi, kappa=kappa, K=K, C=C, R=R
                    )
    return poses

def read_correspondence_pix(correspondence_pix_path: str):
    with open(str(correspondence_pix_path),'rb') as f:
        data = f.read()
        numpts = np.frombuffer(data, dtype=np.uint32, count=1,offset=0)[0]
        decdata = np.frombuffer(data,dtype=np.float32,count=numpts*7,offset=4).reshape(numpts,7)
        srcUV = decdata[:,:2]
        tgtUV = decdata[:,2:4]
        srcXYZ = decdata[:,4:]
        return srcUV, tgtUV, srcXYZ


work_dir = Path(r"M:\MoAve2\data_for_moAve_2\ver2")
qinpose_path = work_dir / "undistorted_Img" / "Ori_QIN.qinv2"
res_dir = work_dir / "results"
neighbor_path = res_dir / "Point_clouds"/ "point_clouds" / "Neighborhood_file.bin"
tmp_dir = res_dir / "Point_clouds"/ "temp_folder_cluster"

assert qinpose_path.exists(), 'qinpose_path does not exist'
assert neighbor_path.exists(), 'neighbor_path does not exist'

def read_neighbor(neighbor_path):
    with open(str(neighbor_path),'rb') as f:
        data = f.read()
        numimg = np.frombuffer(data, dtype=np.uint32, count=1,offset=0)[0]
        maxpair = np.frombuffer(data, dtype=np.uint32, count=1,offset=4)[0]
        neighborbin = np.frombuffer(data,dtype=np.uint32,count=numimg*(maxpair+1),offset=8).reshape(numimg,maxpair+1)
        return numimg, maxpair, neighborbin
  
numimg, maxpair, neighborbin = read_neighbor(neighbor_path)
poses = read_qinpose(qinpose_path)
imgnames  = list(poses.keys())

for srcIdx in range(numimg):
    for tgtIdx in neighborbin[srcIdx][1:]:
        if tgtIdx <=0: continue
        srcName = imgnames[srcIdx]
        tgtName = imgnames[tgtIdx]
        print(f'Pair: {srcName} {tgtName}')
        bin_path = tmp_dir / srcName / f"{srcName}_{tgtName}_correspondence_pix.bin"
        assert bin_path.exists(), f'{bin_path} does not exist'

        srcK = poses[srcName]['K']
        srcR = poses[srcName]['R']
        srcC = poses[srcName]['C']
        
        tgtK = poses[tgtName]['K']
        tgtR = poses[tgtName]['R']
        tgtC = poses[tgtName]['C']
        
        srcUV, tgtUV, srcXYZ = read_correspondence_pix(bin_path)
        
        reproject_tgtUV = (srcXYZ + srcC - tgtC) @ tgtR.T @ tgtK.T
        
        reproject_tgtUV = reproject_tgtUV[:,:2] / reproject_tgtUV[:,2:]
        
        success, rvec, tvec = cv2.solvePnP(objectPoints=srcXYZ.astype(np.float64),
                                           imagePoints=srcUV.astype(np.float64),
                                           cameraMatrix=srcK.astype(np.float64),
                                           distCoeffs=None)
        assert success, 'solvePnP failed' 
        est_srcR, _ = cv2.Rodrigues(rvec)
        est_srcC = -est_srcR.T @ tvec
        
        success, rvec, tvec = cv2.solvePnP(objectPoints=srcXYZ.astype(np.float64),
                                           imagePoints=tgtUV.astype(np.float64),
                                           cameraMatrix=tgtK.astype(np.float64),
                                           distCoeffs=None)
        assert success, 'solvePnP failed' 
        est_tgtR, _ = cv2.Rodrigues(rvec)
        est_rel_C = -est_tgtR.T @ tvec
        
        rel_R = srcR @ tgtR
        est_rel_R = est_srcR @ est_tgtR
        
        rel_euler = Rotation.from_matrix(rel_R).as_euler('xyz', degrees=True)
        est_rel_euler = Rotation.from_matrix(est_rel_R).as_euler('xyz', degrees=True)
        
        # print(rel_euler, est_rel_euler)
        print('Err. Euler (deg): ', rel_euler - est_rel_euler)
        
        # d_rvec = cv2.Rodrigues(tgtR @ est_tgtR.T)[0]
        # d_pos = (tgtC - srcC).ravel() - est_tgtC_srcC.ravel()
        
        # print("rotation error: ", d_rvec.ravel())
        # print("translation error: ", d_pos.ravel())
        # print(tgtR, est_tgtR)
        # print(tgtC.ravel() - srcC.ravel(), est_tgtC_srcC.ravel())
       
    #     break
    # break




