import sys
from pyvis.network import Network
import numpy as np

class Node:
    def __init__(self):
        self.name: str
        self.path: str
        self.xfm = np.zeros(3)


class Edge:
    def __init__(self):
        self.name: str
        self.source: int
        self.target: int
        self.xfm = np.zeros(3)
        self.cov = np.eye(3)

class TranslateGraph:
    def __init__(self):
        self.basepath: str
        self.baseID: int
        self.baseGeo = np.zeros(3)
        self.nodes = []
        self.edges = []

    def read(self, path):
        f = open(path)
        lines = [l.strip() for l in f.readlines()]
        f.close()

        ptr = 0
        if lines[ptr] != 'TranslateGraph 3':
            print(f"Error, unknown line type {lines[ptr]}")
        ptr += 1
        self.basepath = lines[ptr]
        ptr += 1
        self.baseID = int(lines[ptr])
        ptr += 1
        self.baseGeo = np.array([float(t) for t in lines[ptr].split(' ')])
        ptr += 1

        num_nodes = int(lines[ptr])
        ptr += 1
        for ni in range(num_nodes):
            n = Node()
            n.name = lines[ptr]
            ptr += 1
            n.path = lines[ptr]
            ptr += 1
            n.xfm = np.array([float(t) for t in lines[ptr].split(',')])
            ptr += 1
            self.nodes.append(n)

        num_edges = int(lines[ptr])
        ptr += 1
        for ei in range(num_edges):
            e = Edge()
            e.name = lines[ptr]
            ptr += 1
            e.target, e.source = [int(t) for t in lines[ptr].split(',')]
            ptr += 1
            e.xfm = np.array([float(t) for t in lines[ptr].split(',')])
            ptr += 1
            e.cov = np.array([float(t) for t in lines[ptr].split(',')]).reshape(3, 3)
            ptr += 1

            self.edges.append(e)

if __name__ == "__main__":
    input_graph1 = r'N:\citymapper_70\res\graph.txt'
    input_graph2 = r'N:\citymapper_70\res\graph_direct.txt'

    g1 = TranslateGraph()
    g1.read(input_graph1)

    g2 = TranslateGraph()
    g2.read(input_graph2)

    with np.printoptions(precision=5, suppress=True):
        for n1, n2 in zip(g1.nodes, g2.nodes):
            offset = n1.xfm - n2.xfm
            dist = np.linalg.norm(offset)
            print(dist,": ",offset)

    sys.exit(0)

    for ei, e in enumerate(g1.edges):
        initxfm = g1.nodes[e.target].xfm - g1.nodes[e.source].xfm
        edgexfm = e.xfm
        diffxfm = edgexfm - initxfm
        confidence = 1/e.cov[0 , 0]
        print(f'src {g1.nodes[e.source].name[:-4]} tgt {g1.nodes[e.target].name[:-4]} conf {confidence}')
        # print(f'Init: {initxfm}')
        # print(f'Edge: {edgexfm}')
        print(f'Diff: {diffxfm}')

    sys.exit(0)

    baseXfm = tg.nodes[0].xfm

    net = Network()
    for ni, n in enumerate(tg.nodes):

        font_color='blue'
        if n.name == '029_065_id3188c81615_124544_Backward.tif':
            font_color = 'red'
        label = n.name[:7]+n.name[-12:-4]
        net.add_node(ni,  label, value=1, x=n.xfm[0]-baseXfm[0], y=n.xfm[1]-baseXfm[1], color=font_color)

    for ei, e in enumerate(tg.edges):
        net.add_edge(e.source, e.target)

    # net.toggle_physics(True)
    net.show_buttons(filter_=['nodes', 'edges'])
    net.show('mygraph.html')