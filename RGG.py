def random_geometric_graph(n, radius, dim=2, pos=None):
    G=nx.Graph()
    G.name="Random Geometric Graph"
    G.add_nodes_from(range(n)) 
    if pos is None:
        # random positions
        for n in G:
            G.node[n]['pos']=[random.random() for i in range(0,dim)]
    else:
        nx.set_node_attributes(G,'pos',pos)
    # connect nodes within "radius" of each other
    # n^2 algorithm, could use a k-d tree implementation
    nodes = G.nodes(data=True)
    while nodes:
        u,du = nodes.pop()
        pu = du['pos']
        for v,dv in nodes:
            pv = dv['pos']
            d = sum(((a-b)**2 for a,b in zip(pu,pv)))
            if d <= radius**2:
                G.add_edge(u,v)
return G