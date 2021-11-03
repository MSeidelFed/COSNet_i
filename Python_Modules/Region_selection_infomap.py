#!/usr/bin/env python3

"""
RegionSelectionInfomap.py:

python3 RegionSelectionInfomap.py 

"""
## IMPORTS

from sys import argv

import networkx as nx

from infomap import Infomap

import matplotlib.pyplot as plt

## Data

if __name__ == "__main__":

    if len(argv) == 3:

        File_name = argv[1]

        lines = []
        with open(File_name, 'r') as f:
            lines = f.readlines()

        ## Initialize graph

        G = nx.Graph()

        ## fill empty graph with weighted nodes

        for i in lines:

            link_i = i.split()
            G.add_edge(link_i[0], link_i[1], weigt = int(link_i[2]))
        
        ## Arg 2

        Print_graph_PNG = argv[2]

        if Print_graph_PNG == True:

            ## plotting the graph

            pos = nx.spring_layout(G, seed=7)  # positions for all nodes - seed for reproducibility

            ### nodes
            nx.draw_networkx_nodes(G, pos, node_size=700)

            ### edges
            nx.draw_networkx_edges(G, pos, width=6)
            nx.draw_networkx_edges(
                G, pos, width=6, alpha=0.5, edge_color="b", style="dashed"
            )

            ### labels
            nx.draw_networkx_labels(G, pos, font_size=20, font_family="sans-serif")

            ax = plt.gca()
            ax.margins(0.08)
            plt.axis("off")
            plt.tight_layout()
            plt.show()
            plt.savefig(fname = "test.png")

        ## finding communites in the created graph
        ### initializing a new infomap object

        im = Infomap(silent=True)

        ### filling the object with the networkx graph
    
        mapping = im.add_networkx_graph(G)

        ### Infomap analysis

        im.run()
        print(f"Found {im.num_top_modules} modules with codelength: {im.codelength}")

        for node in im.nodes:
            print(node.node_id, node.module_id, node.flow, mapping[node.node_id])

    else:
            print('################################################################################################\n')
            print('    Usage: python3 RegionSelectionInfomap.py [Input_file] [Print_graph_PNG]  \n')
            print('################################################################################################\n')


