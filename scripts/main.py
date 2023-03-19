import os
import os.path as osp
from subprocess import Popen, PIPE, STDOUT, DEVNULL

IMG_FOLDER= r'M:\MoAve\NUS\images_undist'
RES_FOLDER= r'M:\MoAve\NUS\res'
GCPFILE = r'M:\MoAve\NUS\GCP_undist.json'
MINRAYS=3
SIGMA=1
GRAPH_NAME = r'graph_msp_gcp.txt'
SOL_GRAPH_NAME = GRAPH_NAME[:-4]+'_SOL' + GRAPH_NAME + GRAPH_NAME[-4:]

EXE_CREATE_GRAPH = r'C:\ssdev\method\motionaverage\x64\Release\create_translate_graph_MSP.exe'
EXE_SOLVE_GRAPH = r'C:\ssdev\method\motionaverage\x64\Release\solve_translate_graph.exe'
EXE_APPLY_GRAPH = r'C:\ssdev\method\motionaverage\x64\Release\apply_translate_graph_MSP.exe'

input_qinpose = osp.join(IMG_FOLDER, 'pose.qinv2')
cmd_create_graph = [EXE_CREATE_GRAPH, RES_FOLDER,input_qinpose,  GRAPH_NAME, '-v', '0', '--minrays' ,str(MINRAYS), '--sigma', str(SIGMA)]
if GCPFILE:
  cmd_create_graph += ['--gcpfile' ,GCPFILE]

print(' '.join(cmd_create_graph))

p = Popen(cmd_create_graph,stdout=STDOUT, stdin=DEVNULL)
p.wait()

cmd_solve_graph = [EXE_SOLVE_GRAPH, osp.join(RES_FOLDER, GRAPH_NAME), osp.join(RES_FOLDER, SOL_GRAPH_NAME)]
print(' '.join(cmd_solve_graph))

p = Popen(cmd_solve_graph,stdout=STDOUT, stdin=DEVNULL)
p.wait()

output_qinpose = osp.join(IMG_FOLDER, 'pose_sol.qinv2')
cmd_apply_graph = [EXE_APPLY_GRAPH,  '--inqinv2', input_qinpose, '--outqinv2', output_qinpose]
print(' '.join(cmd_apply_graph))

p = Popen(cmd_apply_graph,stdout=STDOUT, stdin=DEVNULL)
p.wait()
