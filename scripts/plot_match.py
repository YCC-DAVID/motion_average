import os
import os.path
import numpy as np
import matplotlib.pyplot as plt
import rasterio
from rasterio.windows import Window
from pathlib import Path
import imageio
from rasterio.enums import Resampling
np.set_printoptions(suppress=True)

def read_correspondence_pix(correspondence_pix_path):
  with open(correspondence_pix_path,'rb') as f:
    data = f.read()
    numpts = np.frombuffer(data, dtype=np.uint32, count=1,offset=0)[0]
    decdata = np.frombuffer(data,dtype=np.float32,count=numpts*7,offset=4).reshape(numpts,7)
    srcUV = decdata[:,:2]
    tgtUV = decdata[:,2:4]
    srcXYZ = decdata[:,4:]
    return srcUV, tgtUV, srcXYZ
  
img_folder = Path(r'M:\MoAve\Dublin\Area2\images_undist')
qinpose_path = Path(img_folder,'pose.qinv2')
assert qinpose_path.exists(), f'{qinpose_path} not found'

def read_qinpose(qinpose_path):
  from numpy import sin, cos
  lines = [l.strip() for l in  qinpose_path.read_text().splitlines()]
  num_image = int(lines[0])
  poses = dict()
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

poses = read_qinpose(qinpose_path)
print(poses)


tmp_folder = Path("M:\\MoAve\\Dublin\\Area2_test\\res\\Point_clouds\\temp_folder_cluster")
srcname = "125872_id1269c1_20150326125028"
tgtname = "125874_id1271c1_20150326125030"

srcC = poses[srcname]['C']
srcR = poses[srcname]['R']
srcK = poses[srcname]['K']

tgtC = poses[tgtname]['C']
tgtR = poses[tgtname]['R']
tgtK = poses[tgtname]['K']

srcpath = list(img_folder.glob(f'{srcname}.*'))[0]
tgtpath = list(img_folder.glob(f'{tgtname}.*'))[0]
assert Path(srcpath).exists(), f'{srcpath} not found.'
assert Path(tgtpath).exists(), f'{tgtpath} not found.'
corrpix_path = Path(tmp_folder, srcname, f'{srcname}_{tgtname}_correspondence_pix.bin')
assert corrpix_path.exists(), f'{corrpix_path} not found.'
srcgrids_path = [ Path(tmp_folder, srcname, f'{srcname}_{tgtname}_temgrid{i}.tif') for i in 'XYZ']
tgtgrids_path = [ Path(tmp_folder, tgtname, f'{tgtname}_{srcname}_temgrid{i}.tif') for i in 'XYZ']
for i in srcgrids_path+tgtgrids_path:
  assert i.exists(), f'{i} not found.'


srcUV, tgtUV, srcBinXYZ = read_correspondence_pix(corrpix_path)
color = np.random.rand(len(srcUV))

if 0:
  srcimg = imageio.imread(srcpath)
  tgtimg = imageio.imread(tgtpath)
  _s = np.s_[::100]
  fig1, ax1 = plt.subplots(1,1)
  ax1.imshow(srcimg)
  ax1.scatter(srcUV[_s,0], srcUV[_s,1], c=color[_s])
  fig2, ax2 = plt.subplots(1,1)
  ax2.imshow(tgtimg)
  ax2.scatter(tgtUV[_s,0], tgtUV[_s,1], c=color[_s])
  plt.show()

if 0:
  idx = np.argmin(np.linalg.norm(srcUV - np.array([[4300,625]]), axis=1))
  print(idx)
  print(srcUV[idx,:])

if 1:
  idx = 0
  srcX, srcY = srcUV[idx,:]
  tgtX, tgtY = tgtUV[idx,:]
  srcbinxyz = srcBinXYZ[idx,:]
  
  srcXYZ = []
  tgtXYZ = []
  for gridpath in srcgrids_path:
    with rasterio.open(gridpath,'r') as tp:
      result = tp.read(1, window = Window(col_off=srcX, row_off=srcY,width=1,height=1), resampling=Resampling.bilinear)
      srcXYZ.append(result.ravel()[0])
      
  for gridpath in tgtgrids_path:
    with rasterio.open(gridpath,'r') as tp:
      result = tp.read(1, window = Window(col_off=tgtX, row_off=tgtY,width=1,height=1))
      tgtXYZ.append(result.ravel()[0])
  
  print('----')
  print(srcbinxyz, srcXYZ, tgtXYZ )
  print(srcC+srcXYZ)
  print(tgtC+tgtXYZ)
  
  _aa = srcK@srcR@np.array(srcXYZ)
  _aa /= _aa[2]
  print(f'Check srcUV {_aa[0]:.3f}, {_aa[1]:.3f}')
  print(f'            {srcX:.3f}, {srcY:.3f}')
  
  _aa = tgtK@tgtR@np.array(tgtXYZ)
  _aa /= _aa[2]
  print(f'Check tgtUV {_aa[0]:.3f}, {_aa[1]:.3f}')
  print(f'            {tgtX:.3f}, {tgtY:.3f}')
  
  srcWld = srcC+srcXYZ
  tgtWld = tgtC+tgtXYZ
  
  src_in_tgt = tgtR @ (srcWld - tgtC)
  tgt_in_tgt = tgtR @ (tgtWld - tgtC)
  print('----------------')
  print('Reproj Depth in TGT:', src_in_tgt[2])
  print('Depth in TGT:', tgt_in_tgt[2])
  print('----------------')
  print('Length src2tgt', np.linalg.norm(src_in_tgt))
  print('Length ray', np.linalg.norm(tgt_in_tgt))
  print('----------------')
  print('src_in_tgt', src_in_tgt)
  print('tgt_in_tgt', tgt_in_tgt)
  print('----------------')
  print('src2tgt px', tgtK@src_in_tgt/src_in_tgt[2])
  print('tgt px', tgtK@tgt_in_tgt/tgt_in_tgt[2])
  print('match tgt', tgtX, tgtY)
      
