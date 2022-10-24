import enum
from mimetypes import common_types
from pickle import FALSE
import numpy as np
import os.path as osp
import matplotlib.pyplot as plt
from matplotlib import cm
from imageio import imread
from glob import glob
from plot_trans3d_graph import Node, Edge, TranslateGraph

import sys
import argparse
np.set_printoptions(suppress=True)

def multivariate_gaussian(x, mean, cov):
  dis_x = x-mean
  num = np.exp((-0.5*(dis_x @ np.linalg.inv(cov))*dis_x).sum(axis=1))
  denom = np.sqrt((2*np.pi)**3*np.linalg.det(cov))
  return num/denom

def parseArgs():
  parser = argparse.ArgumentParser(description="Filter correspondence offset")
  parser.add_argument("input_graph_pri", type=str, help="input prior graph")
  parser.add_argument("input_graph_post", type=str, help="input posterior graph")
  parser.add_argument("--imgdir", type=str, help="img dir")
  
  args = parser.parse_args()
  return args

if __name__ == '__main__':
  args = parseArgs()
  print(args)

  input_graph_pri = args.input_graph_pri
  input_graph_post = args.input_graph_post
  imgdir = args.imgdir
  bPlotImage = imgdir is not None

  # pri_graph = TranslateGraph()
  post_graph = TranslateGraph()
  # pri_graph.read(input_graph_pri)
  post_graph.read(input_graph_post)
  rootdir = post_graph.basepath
  tempdir = osp.join(rootdir, 'Point_clouds','temp_folder_cluster')

  name_to_nid = dict()
  nid_to_name = dict()
  for nid,n in enumerate(post_graph.nodes):
    nodename = osp.splitext(n.name)[0]
    name_to_nid[nodename] = nid
    nid_to_name[nid] = nodename

  for ei, e in enumerate(post_graph.edges):
    srcId = e.source
    tgtId = e.target
    srcName = nid_to_name[srcId]
    tgtName = nid_to_name[tgtId]
    
    print(f'Edge {ei}: {srcId},{tgtId} {srcName}, {tgtName}')

    correspondence_pix_path = osp.join(tempdir, srcName, f'{srcName}_{tgtName}_correspondence_pix.bin')
    out_weight_filtered_path = osp.join(tempdir, srcName, f'{srcName}_{tgtName}_correspondence_pix.weight.bin')
    out_posterior_filtered_path = osp.join(tempdir, srcName, f'{srcName}_{tgtName}_correspondence_pix.graphopt.bin')
    correspondence_offset_path = osp.join(tempdir, srcName, f'{srcName}_{tgtName}_correspondence_offset.bin')

    with open(correspondence_pix_path,'rb') as f:
      data = f.read()
      numpts = np.frombuffer(data, dtype=np.uint32, count=1,offset=0)[0]
      decdata = np.frombuffer(data,dtype=np.float32,count=numpts*4,offset=4).reshape(numpts,4)
      srcUV = decdata[:,:2]
      tgtUV = decdata[:,2:]
      
    with open(correspondence_offset_path,'rb') as f:
      data = f.read()
      numpts = np.frombuffer(data, dtype=np.uint32, count=1,offset=0)[0]
      decdata = np.frombuffer(data,dtype=np.float32,count=numpts*4,offset=4).reshape(numpts,4)
      offset = decdata[:,:3]
      weights = decdata[:,3]
      validmask = weights>0
      assert(not np.isinf(offset[validmask,:]).any())
    # Approach #1: Inlier by comparing MedianFilterGrid vs TempGrid
    if not osp.exists(out_weight_filtered_path):
      w_hist, w_edges = np.histogram(weights[validmask],bins=512)
      w_cdf = np.cumsum(w_hist) / w_hist.sum()
      w_threshold = w_edges[np.searchsorted(w_cdf,0.3)]
      mask_by_weight = weights > w_threshold
      print(f'masking_by_weight {mask_by_weight.sum()}/{len(weights)}')

      with open(out_weight_filtered_path,'wb') as f:
        numpts = np.array([mask_by_weight.sum()],dtype=np.uint32)
        encdata = np.hstack([srcUV[mask_by_weight,:], tgtUV[mask_by_weight,:]]).astype(dtype=np.float32)
        f.write(numpts.tobytes())
        f.write(encdata.tobytes())
        
    # Approach #2:
    if not osp.exists(out_posterior_filtered_path):
      # Get statistic from measurement
      observed_mean = np.average(offset[validmask,:], axis=0, weights=weights[validmask])
      observed_cov = np.cov(offset[validmask,:], rowvar=False, aweights = weights[validmask])
      # print('obs. mean', observed_mean)
      # print('obs. cov', observed_cov)

      post_xfm = -post_graph.nodes[srcId].xfm + post_graph.nodes[tgtId].xfm
      post_cov = np.eye(3)*0.01
      # print('edge_xfm', e.xfm)
      # print('post_xfm:', post_xfm)
    
      p_prior = np.zeros_like(weights)
      p_post = np.zeros_like(weights)
      
      p_prior[validmask] = multivariate_gaussian(offset[validmask,:], observed_mean, observed_cov)
      p_post[validmask] = multivariate_gaussian(offset[validmask,:], post_xfm, observed_cov)
      
      # p_prior[validmask] = np.log(np.linalg.norm(offset[validmask,:] - observed_mean, axis=1))
      # p_post[validmask] = np.log(np.linalg.norm(offset[validmask,:] - post_xfm, axis=1))
      
      p_hist, p_edges = np.histogram(p_post[validmask],bins=512)
      p_cdf = np.cumsum(p_hist) / p_hist.sum()
      p_threshold = p_edges[np.searchsorted(p_cdf,0.3)]
      mask_by_posterior = p_post > p_threshold
      print(f'masking_by_posterior {mask_by_posterior.sum()}/{len(weights)}')

      with open(out_posterior_filtered_path,'wb') as f:
        numpts = np.array([mask_by_posterior.sum()],dtype=np.uint32)
        encdata = np.hstack([srcUV[mask_by_posterior,:], tgtUV[mask_by_posterior,:]]).astype(dtype=np.float32)
        f.write(numpts.tobytes())
        f.write(encdata.tobytes())

      if bPlotImage:
        srcImgPath = glob(osp.join(imgdir,f'{srcName}.*'))[0]
        tgtImgPath = glob(osp.join(imgdir,f'{tgtName}.*'))[0]
        srcImg = imread(srcImgPath)
        tgtImg = imread(tgtImgPath)
        
        splr = np.s_[::5]
        show_prior = True
        if show_prior:      
          prob = p_prior
          label = 'Probability of observed distribution'
        else:
          prob = p_post
          label = 'Probability of graph-optimized distribution'
        
        fig, ax = plt.subplots(1,2)
        ax[0].imshow(srcImg)
        pcm = ax[0].scatter(srcUV[validmask,0][splr], srcUV[validmask,1][splr],c=p_prior[validmask][splr], s=10, alpha=0.6)
        ax[0].set_title("SrcImage")
        ax[1].imshow(tgtImg)
        pcm = ax[1].scatter(tgtUV[validmask,0][splr], tgtUV[validmask,1][splr],c=p_prior[validmask][splr], s=10, alpha=0.6)
        ax[1].set_title("TgtImage")
        fig.colorbar(pcm, ax = ax, location='bottom')
        fig.suptitle(f'Observed \n({srcName}, {tgtName})')
        
        fig, ax = plt.subplots(1,2)
        ax[0].imshow(srcImg)
        pcm = ax[0].scatter(srcUV[validmask,0][splr], srcUV[validmask,1][splr],c=p_post[validmask][splr], s=10, alpha=0.6)
        ax[0].set_title("SrcImage")
        ax[1].imshow(tgtImg)
        pcm = ax[1].scatter(tgtUV[validmask,0][splr], tgtUV[validmask,1][splr],c=p_post[validmask][splr], s=10, alpha=0.6)
        ax[1].set_title("TgtImage")
        fig.colorbar(pcm, ax = ax, location='bottom')
        fig.suptitle(f'Posterior \n({srcName}, {tgtName})')
        
      if bPlotImage:
        fig, ax = plt.subplots(2,1)
        pri_bins, pri_edges, _ = ax[0].hist(p_prior[validmask], bins=256)
        pri_cdf = np.concatenate([[0],np.cumsum(pri_bins)]) / pri_bins.sum()
        ax0_twinx=ax[0].twinx()
        ax0_twinx.plot(pri_edges, pri_cdf, c='r')
        ax[0].set_xlabel("p(x)")
        ax[0].set_ylabel("Count")
        ax[0].set_title("Observed Distribution")
        ax0_twinx.set_ylabel("F(p(x))")
        
        post_bins, post_edges, _ = ax[1].hist(p_post[validmask], bins=256)
        post_cdf =  np.concatenate([[0],np.cumsum(post_bins)]) / post_bins.sum()
        ax1_twinx = ax[1].twinx()
        ax1_twinx.plot(post_edges, post_cdf, c='r')
        ax[1].set_title("Graph-Optimized Distribution")
        ax[1].set_xlabel("p(x)")
        ax[1].set_ylabel("Count")
        ax1_twinx.set_ylabel("F(p(x))")
        
        fig.suptitle(f'Histogram of Probability\n({srcName}[{srcId}], {tgtName}[{tgtId}])')
        
        for perc in [0.1, 0.5, 0.9]:
          perc_idx = np.searchsorted(pri_cdf, perc)
          ax0_twinx.axvline(pri_edges[perc_idx], color='g')
          ax0_twinx.text(pri_edges[perc_idx], 0.7, f' {perc*100:.1f}%\n @{perc_idx}')
        for perc in [0.1, 0.5, 0.9]:
          perc_idx = np.searchsorted(post_cdf, perc)
          ax1_twinx.axvline(post_edges[perc_idx], color='g')
          ax1_twinx.text(post_edges[perc_idx], 0.7, f' {perc*100:.1f}%\n @{perc_idx}')
        
        plt.show()
      bPlotImage = False
