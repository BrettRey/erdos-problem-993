
import json
import networkx as nx
from graph6 import parse_graph6

def analyze_structure(g6_str):
    if isinstance(g6_str, str):
        g6_str = g6_str.encode('ascii')
    
    n, adj = parse_graph6(g6_str)
    
    # Convert to NetworkX for easier analysis
    G = nx.Graph()
    G.add_nodes_from(range(n))
    for u, neighbors in enumerate(adj):
        for v in neighbors:
            if u < v:
                G.add_edge(u, v)
                
    degrees = [d for n, d in G.degree()]
    max_deg = max(degrees)
    leaves = degrees.count(1)
    diameter = nx.diameter(G)
    
    # Check if it's a spider (one vertex with degree > 2)
    high_degree_nodes = [v for v, d in G.degree() if d > 2]
    is_spider = len(high_degree_nodes) == 1
    
    structure_type = "Unknown"
    details = ""
    
    if is_spider:
        center = high_degree_nodes[0]
        structure_type = "Spider"
        # Calculate leg lengths
        legs = []
        for neighbor in G.neighbors(center):
            # perform BFS/DFS from neighbor away from center to count length
            curr = neighbor
            prev = center
            length = 1
            while G.degree(curr) > 1:
                # Find next node
                for nbr in G.neighbors(curr):
                    if nbr != prev:
                        prev = curr
                        curr = nbr
                        length += 1
                        break
            legs.append(length)
        legs.sort(reverse=True)
        details = f"Legs: {legs}"
        
    elif len(high_degree_nodes) == 0:
        structure_type = "Path"
    else:
        structure_type = f"Tree (Branch points: {len(high_degree_nodes)})"
        # Maybe check if it's a caterpillar (path spine)
        # Remove leaves, check if path
        H = G.copy()
        while True:
            leaves_nodes = [node for node in H.nodes() if H.degree(node) == 1]
            if not leaves_nodes: break
            if len(H) <= 2: break # It's a path
            # But wait, caterpillar check is: removing leaves once leaves a path.
            H.remove_nodes_from(leaves_nodes)
            break
            
        # Check if H is a path
        if nx.is_connected(H):
            h_degrees = [d for n, d in H.degree()]
            if all(d <= 2 for d in h_degrees):
                structure_type = "Caterpillar"
                details = f"Spine length: {len(H)}"
    
    print(f"Graph: {g6_str.decode('ascii')}")
    print(f"  N={n}, MaxDeg={max_deg}, Leaves={leaves}, Diameter={diameter}")
    print(f"  Type: {structure_type}")
    if details:
        print(f"  {details}")
    print("-" * 40)

def main():
    # LC Failures from analysis_n26.json
    lc_failures = [
        "Y???????????_?O?C??_?A??C??C??A???_??C???O?[?_?F`???^???",
        "Y???????????_?O?C??_?A??C??C??A???_??C?C?O?K@_?F@???|???"
    ]
    
    print("Analyzing LC Failure Structures:")
    for g6 in lc_failures:
        analyze_structure(g6)

if __name__ == "__main__":
    main()
