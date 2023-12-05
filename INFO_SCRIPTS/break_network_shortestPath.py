import networkx as nx
import matplotlib.pyplot as plt
import math

G = nx.Graph()

G.add_nodes_from([1,2,3,4,5])

G.add_edge(1,2)
G.add_edge(1,3)
G.add_edge(1,4)
G.add_edge(2,3)
G.add_edge(3,4)
G.add_edge(3,5)
G.add_edge(2,4)
G.add_edge(4,5)
G.add_edge(2,5)
G.add_edge(1,5)

#print(G.edges())

G[1][3]['weight'] = 0.7
G[5][3]['weight'] = 0.2
G[1][2]['weight'] = 0.8
G[4][2]['weight'] = 0.9
G[3][2]['weight'] = 0.6
G[5][2]['weight'] = 0.05
G[4][5]['weight'] = 0.5
G[1][4]['weight'] = 0.6
G[3][4]['weight'] = 0.1
G[1][5]['weight'] = 0.01

#print(G.get_edge_data(2,3))

# copy the original graph
inverted = G.copy()

for u, v, data in inverted.edges(data=True):
    #print(u, v, data)
    try:
       #data['weight'] = 1 / data['weight']         # one way
       data['weight'] = -math.log(data['weight'])  # second way
    #except ZeroDivisionError:
    except ValueError:
        print("Warning ! You are dividing by zero ")
        print("or something wrong in the logarithm")

source = 1
target = 5

print(f'Shortest path between {source} and {target}')
#print(nx.shortest_path(inverted, source=source, target=target, weight='weight'))
#print(nx.shortest_path(G, source=source, target=target, weight='weight', method='dijkstra'))
print(nx.shortest_path(inverted, source=source, target=target, weight='weight', method='dijkstra'))
print('Shortest path length: ')
print(nx.shortest_path_length(inverted, source=source, target=target, weight='weight', method='dijkstra'))

##----------------------------
## plotting ------------------

pos = nx.spring_layout(G, seed=7)  # positions for all nodes - seed for reproducibility

## nodes
nx.draw_networkx_nodes(G, pos, node_size=700)

## edges
nx.draw_networkx_edges(G, pos, edgelist=G.edges(), width=6)

## node labels
nx.draw_networkx_labels(G, pos, font_size=20, font_family="sans-serif")

## edge weight labels
edge_labels = nx.get_edge_attributes(G, "weight")
nx.draw_networkx_edge_labels(G, pos, edge_labels)

ax = plt.gca()
ax.margins(0.08)
plt.axis("off")
plt.tight_layout()
plt.show()

quit()

