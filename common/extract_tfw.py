import os
import os.path as osp
import numpy as np
import gdal
import argparse

np.set_printoptions(suppress=True)

def Xfm2Numpy(xfm):
    return np.array([[xfm[1],xfm[2],xfm[0]],
                     [xfm[4],xfm[5],xfm[3]],
                     [0,0,1]],dtype=np.float64)
    
def Numpy2Xfm(mat):
    return (mat[0,2],mat[0,0],mat[0,1], mat[1,2], mat[1,0], mat[1,1])

def str_tfw(mat):
    return f'{mat[0,0]}\n{mat[1,0]}\n{mat[0,1]}\n{mat[1,1]}\n{mat[0,2]}\n{mat[1,2]}\n'
Geo2Tfw = np.array([[1, 0, 0.5],
                     [0, 1, 0.5],
                     [0, 0, 1]],dtype=np.float64)

Tfw2Geo = np.array([[1, 0, -0.5],
                     [0, 1, -0.5],
                     [0, 0, 1]],dtype=np.float64)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--xml",type=str, nargs='?')
    parser.add_argument("--glob", type=str, nargs='?')
    args = parser.parse_args()

    tiflist = []
    if args.xml is not None:
        import re
        content = open(args.xml,'r').read()
        _xmllist = re.findall('<Image_Name>(.*?)</Image_Name>',content)
        tiflist += _xmllist

    if args.glob is not None:
        from glob import glob
        _globlist = list(glob(args.glob,recursive=True))
        tiflist += _globlist
    
    for tifname in tiflist:
        ds = gdal.Open(tifname, gdal.GA_ReadOnly)
        tfwname = tifname[:-4]+'.tfw'
        geoAff = Xfm2Numpy(ds.GetGeoTransform())
        tfwAff =  geoAff @ Geo2Tfw
        open(tfwname,'w').write(str_tfw(tfwAff))
        
    print('Done')
    





    
