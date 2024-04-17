import networkx as nx
import matplotlib.pyplot as plt

def parse_graph_description(description):
    """
    Parse a graph description and return a directed graph.
    
    Parameters:
        description (str): Multi-line string where each line describes a node and its edges.
    
    Returns:
        nx.DiGraph: A directed graph constructed based on the description.
    """
    graph = nx.DiGraph()
    lines = description.strip().split('\n')
    for line in lines:
        parts = line.split(':')
        node = parts[0].split()[1]  # Obtener el nodo después de 'Node'
        if len(parts) > 1 and parts[1].strip():  # Comprobar si hay aristas
            edges = parts[1].strip().split()
            for edge in edges:
                graph.add_edge(node, edge)
        else:
            graph.add_node(node)  # Añadir el nodo sin aristas
    return graph

def visualize_graph(graph):
    """
    Visualize the given directed graph using matplotlib.
    
    Parameters:
        graph (nx.DiGraph): The graph to be visualized.
    """
    pos = nx.spring_layout(graph)  # positions for all nodes
    nx.draw(graph, pos, with_labels=True, node_color='skyblue', node_size=2000, edge_color='k', linewidths=1, font_size=15, arrows=True)
    plt.show()


# Example description
description = """
Node TA has edges to: AG 
Node GT has edges to: TA 
Node AG has edges to: GT GT 
"""

# Parse the description and visualize the graph
graph = parse_graph_description(description)
visualize_graph(graph)
