import argparse, os
from algo import *
from graphg import solve
from fastg2 import MyGraph
import networkx as nx

from utils import *


def parse_user_input():
    parser = argparse.ArgumentParser(
        description=
        'Extract cycle nodes and edges'
    )
    parser.add_argument('-g', '--graph',
                        help='(spades 3.50+) assembly graph FASTG file to process; recommended for spades 3.5: before_rr.fastg, for spades 3.6+:assembly_graph.fastg',
                        required=True, type=str
                        )
    parser.add_argument('-o', '--output_fastg_file',
                        help='Output file',
                        required=True, type=str
                        )
    return parser.parse_args()


def nodes_to_idxes(nodes):
    res = {}
    for idx, node in enumerate(nodes):
        if node == "EDGE_656670_length_61_cov_396.500000'":
            print("yes")
        res[node] = idx
    return res


if __name__ == '__main__':
    args = parse_user_input()
    fastg = args.graph
    fp = open(fastg, 'r')
    output_fastg_file = args.output_fastg_file

    G = get_fastg_digraph(fastg)
    print("finished read")
    G.remove_nodes_from(list(nx.isolates(G)))
    print("finished remove")
    print("nodes count: ", len(G.nodes))
    print("edges count: ", len(G.edges))
    # print(G.edges)
    node_idxes_map = nodes_to_idxes(G.nodes)
    # print(node_idxes_map)
    # Edges = [ [ 0, 1 ], [ 0, 2 ],
    #           [ 1, 3 ], [ 2, 1 ],
    #           [ 2, 5 ], [ 3, 0 ],
    #           [ 4, 5 ] ];
    Edges = [[0,1],[1,2],[2,3],[3,1],[2,4]]
    # solve(5,5,Edges)
    # solve(len(G.nodes), len(G.edges), list(G.edges), node_idxes_map)
    g = MyGraph(len(G.nodes))
    for link in G.edges:
        g.addEdge(node_idxes_map[link[0]], node_idxes_map[link[1]])
    g.printGraph()
    g.nonCyclicNode()
    # comps = nx.strongly_connected_component_subgraphs(G)
    print("finished saperating")
    # newg = {}
    # for line in nx.generate_adjlist(G):
    #     line_lst = line.split(" ")
    #     newg[line_lst[0]]=line_lst[1:]
    # print(newg)
    # print(list(simple_cycles(newg)))
    # print("all subgraph: ", length(comps))
    # cycles = nx.simple_cycles(G)
    # edges = [(0, 0), (0, 1), (0, 2), (1, 2), (2, 0), (2, 1), (2, 2)]
    # G = nx.DiGraph(edges)
    # minimum_cycle_basis
    # print(list(nx.recursive_simple_cycles(G)))
    # subg = nx.simple_cycles(G)
    # print(list(nx.recursive_simple_cycles(G)))
    # visited = set()
    # for idx, s in enumerate(subg):
    #     print("searching in graph", idx)
    #     for l in s:
    #         print(l)
    #         visited.add(l)
    # print("finished nodes set")
    # print("all nodes:", len(visited))
    # for c in comps:
    # subG = nx.simple_cycles(c)
    # print(nx.simple_cycles(c))
    # for s in subG:
    #     print(len(s))
    # try:
    #
    #     cyc = nx.recursive_simple_cycles(c)
    #     # print(cyc)
    # except:
    #     continue
    # print(list(subG))
    # print(list(nx.simple_cycles(G))))
    print("finished cycling")
    fp.close()
