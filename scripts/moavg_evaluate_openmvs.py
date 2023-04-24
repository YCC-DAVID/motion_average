from pyrsp import rsp_utils
from pathlib import Path
import numpy as np
import rasterio
import matplotlib.pyplot as plt

eval_dir = Path(r'M:\MoAve\Dublin\Area1\res\colmap')

cell_size = 0.1
before_las_dir = eval_dir / 'before_las'
after_las_dir = eval_dir / 'after_las'
gt_before_las_dir = eval_dir / 'gt_before_las'
gt_after_las_dir = eval_dir / 'gt_after_las'

before_path = before_las_dir / 'dsm.tif'
after_path = after_las_dir / 'dsm.tif'
gt_before_path = gt_before_las_dir / 'dsm.tif'
gt_after_path = gt_after_las_dir / 'dsm.tif'


if not before_path.exists():
  rsp_utils.rsp_las2dsm(str(before_las_dir),str(before_las_dir), cell_size=cell_size)

if not after_path.exists():
  rsp_utils.rsp_las2dsm(str(after_las_dir),str(after_las_dir), cell_size=cell_size)
  
if not gt_before_path.exists():
  rsp_utils.rsp_las2dsm(str(gt_before_las_dir),str(gt_before_las_dir), cell_size=cell_size)
  
if not gt_after_path.exists():
  rsp_utils.rsp_las2dsm(str(gt_after_las_dir),str(gt_after_las_dir), cell_size=cell_size)
  
#######################
do_whitten = True

# %%
gtbefore = rasterio.open(gt_before_path).read(1)
gtafter = rasterio.open(gt_after_path).read(1)
after = rasterio.open(after_path).read(1)
before = rasterio.open(before_path).read(1)
print('Shapes:', gtbefore.shape, gtafter.shape, before.shape, after.shape)
######### Filter NAN INF etc
finite_threashold = 1000.0
second_threshold = 5.0
third_threshold = 0.5

validmask = np.logical_and.reduce([
    np.logical_and(before>-finite_threashold, before<finite_threashold),
    np.logical_and(after>-finite_threashold, after<finite_threashold),
    np.logical_and(gtbefore>-finite_threashold, gtbefore<finite_threashold),
    np.logical_and(gtafter>-finite_threashold, gtafter<finite_threashold)])

print('  valid {}/{}'.format(validmask.sum(), np.prod(validmask.shape)))
err_bef_before = (gtbefore - before)[validmask]
err_after_before = (gtbefore - after)[validmask]
del gtbefore
err_bef_after = (gtafter - before)[validmask]
err_after_after = (gtafter - after)[validmask]
del gtafter
del after
del before

validmask2 = np.logical_and.reduce([np.abs(err_bef_before)<second_threshold,
                                    np.abs(err_after_before)<second_threshold,
                                    np.abs(err_bef_after)<second_threshold,
                                    np.abs(err_after_after)<second_threshold,])

validmask3 = np.logical_and.reduce([np.abs(err_bef_before)<third_threshold,
                                    np.abs(err_after_before)<third_threshold,
                                    np.abs(err_bef_after)<third_threshold,
                                    np.abs(err_after_after)<third_threshold,])

print('  valid2 {}/{}'.format(validmask2.sum(), len(validmask2)))
print('  valid3 {}/{}'.format(validmask3.sum(), len(validmask3)))

if do_whitten:
    err_bef_before -= np.mean(err_bef_before[validmask2])
    err_after_before -= np.mean(err_after_before[validmask2])

    err_bef_after -= np.mean(err_bef_after[validmask2])
    err_after_after -= np.mean(err_after_after[validmask2])
# %%
fig, ax = plt.subplots()
ax.hist(err_bef_before, bins=64, range=(-3,3), color='red', ls='dashed', lw=3, alpha=0.5, label='Before Adjustment')
ax.hist(err_after_before, bins=64, range=(-3,3), color='green', ls='dashed', lw=3,  alpha=0.5, label='After Adjustment')
ax.set_ylabel('Frequency')
ax.set_xlabel('Vertical Error')
ax.set_title("Compare with LiDAR Register Before")
ax.legend()

fig, ax = plt.subplots()
ax.hist(err_bef_after, bins=64, range=(-3,3), color='red', ls='dashed', lw=3, alpha=0.5, label='Before Adjustment')
ax.hist(err_after_after, bins=64, range=(-3,3), color='green', ls='dashed', lw=3,  alpha=0.5, label='After Adjustment')
ax.set_ylabel('Frequency')
ax.set_xlabel('Vertical Error')
ax.set_title("Compare with LiDAR Register After")
ax.legend()
plt.show()

######### Filter Outliers

### Compare with Before:
def compute_metrics(err_bef, err_after, threshold):    
    mean_before = np.mean(err_bef[validmask3])
    mae_before = np.mean(np.abs(err_bef[validmask3]))
    rmse_before = np.sqrt(np.mean(np.power(err_bef[validmask3],2)))
    print('  Before: Mean {}, MAE {}, RMSE {}'.format(mean_before, mae_before, rmse_before))

    mean_after = np.mean(err_after[validmask3])
    mae_after = np.mean(np.abs(err_after[validmask3]))
    rmse_after = np.sqrt(np.mean(np.power(err_after[validmask3],2)))

    print('  After: Mean {}, MAE {}, RMSE {}'.format(mean_after, mae_after, rmse_after))

    print(f'Threashold  : {third_threshold} m')
    print(f'MAE  Changed: {(mae_before-mae_after)/mae_before*100:+.3f}%')
    print(f'RMSE Changed: {(rmse_before-rmse_after)/rmse_before*100:+.3f}%')
 
print("Compute Metrics with GT Reigstered with Before")
compute_metrics(err_bef_before,err_after_before,second_threshold)
print("Compute Metrics with GT Reigstered with After")
compute_metrics(err_bef_after,err_after_after,second_threshold)
  
  


