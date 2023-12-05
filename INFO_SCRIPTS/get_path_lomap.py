import lomap
import networkx as nx

db_mol = lomap.DBMolecules("tmp/", output=True)

strict, loose = db_mol.build_matrices()

nx_graph = db_mol.build_graph()

# get node values for the parent molecules
for key, value in nx.get_node_attributes(nx_graph, 'fname_comp').items():
     if 'parent_1.sdf' == value:
         node1 = key
     if 'parent_2.sdf' == value:
         node2 = key

# get the shortest path from netowrk
sp = nx.shortest_path(nx_graph, source=node1, target=node2, weight='weight')

for ea in pathGraph.edges():
    print(ea, nx_graph.edges[ea[0], ea[1]])
    print(nx_graph.nodes[ea[0]]['fname_comp'])
    print(nx_graph.nodes[ea[1]]['fname_comp'])



