import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.collections import LineCollection
import numpy as np
from glob import glob
import re
import struct
import os

gtpath = r'N:\dump_graph\gt.bin'
itpath = list(glob(r'N:\dump_graph\iter*.bin'))
itpath = sorted(itpath, key=lambda x: int(re.search('iter(.*)\.bin', x).group(1)))

gtpathtimestamp = os.path.getmtime(gtpath)
itpath = filter(lambda p: os.path.getmtime(p)>gtpathtimestamp, itpath)


def read_nodes(f):
    num_nodes = struct.unpack('N', f.read(8))[0]
    nodes = np.asarray(struct.unpack(f'{2*num_nodes}d', f.read(num_nodes*2*8))).reshape((num_nodes, -1))
    return nodes

def read_edges(f):
    num_edges = struct.unpack('N', f.read(8))[0]
    edges = []
    for ni in range(num_edges):
        edges.append(dict())
        num_nei = struct.unpack('N', f.read(8))[0]
        for ei in range(num_nei):
            oid = struct.unpack('i', f.read(4))[0]
            edata = struct.unpack('2d', f.read(16))
            edges[ni][oid] = np.array(edata)
    return edges

with open(gtpath,'rb') as f:
    gtnodes = read_nodes(f)
    gtedges = read_edges(f)

snaps=[]
for snap in itpath:
    with open(snap,'rb') as f:
        snaps.append(read_nodes(f))
num_frames = len(snaps)
#################3

def draw_nodes(nodes, with_anno=True, ann_fmt='{}', **kwargs):
    ax.scatter(nodes[:,0], nodes[:,1], **kwargs)
    if not with_anno:
        return
    for i, pos in enumerate(nodes):
        ax.annotate(ann_fmt.format(i), pos)


def update(cur_frame, init=False):
    if not init:
        _xlim = ax.get_xlim()
        _ylim = ax.get_ylim()
    ax.clear()
    ax.set_title(f'Iter {cur_frame}')
    draw_nodes(gtnodes, True, ann_fmt='{}|', c='r')
    draw_nodes(snaps[cur_frame], True, c='g')

    lines = []
    cur_nodes = snaps[cur_frame]
    for ni, pos in enumerate(cur_nodes):
        for oni, e in gtedges[ni].items():
            opos = cur_nodes[oni]
            opos_e = opos + e
            lines.append([(opos[0], opos[1]), (opos_e[0], opos_e[1])])
        # ax.annotate('_{}'.format(ni), opos_e)
    lc = LineCollection(lines, lw=1)
    plt.gca().add_collection(lc)

    if not init:
        ax.set_xlim(_xlim)
        ax.set_ylim(_ylim)
    fig.canvas.draw()

def on_press(event):
    global cur_frame
    if event.key == 'left':
        cur_frame = np.clip(cur_frame-1, 0, num_frames-1)
    elif event.key == 'right':
        cur_frame = np.clip(cur_frame+1, 0, num_frames-1)
    else:
        return

    update(cur_frame, init=False)

mpl.rcParams['keymap.back'].remove('left')
mpl.rcParams['keymap.forward'].remove('right')
fig, ax = plt.subplots()
fig.canvas.mpl_connect('key_press_event', on_press)
cur_frame = 0
update(1, init=True)
plt.show()


