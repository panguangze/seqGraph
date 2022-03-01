import networkx as nx
from algo import *
edges = [(0, 0), (0, 1), (0, 2), (1, 2), (2, 0), (2, 1), (2, 2)]
G = nx.DiGraph(edges)
newg = {}
for line in nx.generate_adjlist(G):
    line_lst = line.split(" ")
    newg[line_lst[0]]=line_lst[1:]

print(list(simple_cycles(newg)))
