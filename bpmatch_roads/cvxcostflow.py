
import itertools

import numpy as np

from mygraph import mygraph
from dijkstra import Dijkstra




""" Utility Algorithms """

def Residuals( graph, capacity, flow=None ) :
    rcapacity = {}
    if flow is None : flow = {}
    
    for e in graph.edges() :
        u = capacity.get( e, np.Inf )
        x = flow.get( e, 0. )
        assert u >= 0. and x >= 0. and x <= u
        
        rcapacity[(e,-1)] = x
        rcapacity[(e,+1)] = u - x
        
    return rcapacity


def AugmentCapacity( rcapacity, e, x ) :
    # query residual capacity in both directions
    frwd = (e,+1)
    bkwd = (e,-1)
    
    u_frwd = rcapacity.setdefault( frwd, np.Inf )
    u_bkwd = rcapacity.setdefault( bkwd, np.Inf )
    
    assert x <= u_frwd
    assert -x <= u_bkwd
    
    rcapacity[frwd] -= x
    rcapacity[bkwd] += x
    
    
def MaintainResidualGraphEdge( rgraph, e, rcapacity, graph ) :
    i,j = graph.endpoints(e)
    iter = zip( [ (e,+1), (e,-1) ], [ (i,j), (j,i) ] )
    
    for label, (ii,jj) in iter :
        if rgraph.has_edge(label) : rgraph.remove_edge(label)
        u = rcapacity[label]
        if u > 0. : rgraph.add_edge( label, ii, jj )
        
        
        



""" Convex Cost Flow Algorithm """

def MinCostCovexFlow( network, capacity, supply, cost, epsilon=1 ) :
    """
    graph is a mygraph (above)
    capacity is a dictionary from E -> real capacities
    supply is a dictionary from E -> real supplies
    cost is a dictionary from E -> lambda functions of convex cost edge costs
    """
    flow = {}
    potential = {}
    excess = {}
    reducedcost = {}
    
    for i in network.nodes() :
        pot[i] = 0.
    for e in network.edges() :
        flow[e] = 0.
        
    U = sum([ abs(b) for b in supply.values() ])
    delta,Delta = (1,1)
    while Delta <= U : delta,Delta = Delta,2*Delta
    
    while True :
        
        if Delta <= 1 : break
    
    
    












if __name__ == '__main__' :
    import networkx as nx
    g = mygraph()
    g.add_edge( 'a', 0, 1 )
    g.add_edge( 'b', 0, 2 )
    g.add_edge( 'c', 2, 3 )
    g.add_edge( 'd', 3, 0 )
    
    c = { 'a' : 10., 'b' : 5., 'c' : 1., 'd' : .5 }
    f = { 'a' : 3., 'b' : 1.37 }
    
    if False :
        g = mygraph()
        c = {}
        edgegen = itertools.count()
        
        nxg = nx.DiGraph()
        
        X = [ ( i, np.random.rand(2) ) for i in range(10) ]
        for (i,x), (j,y) in itertools.product( X, X ) :
            cost = np.linalg.norm( y - x )
            if cost < .3 :
                e = edgegen.next()
                g.add_edge( e, i, j )
                c[e] = cost
                
                nxg.add_edge( i, j, weight=cost )
                
        s = 0
        d, pred = Dijkstra( g, c, s )
        nxd = nx.single_source_dijkstra_path_length( nxg, s )
        
        
    cap = {}
    for e in g.edges() : cap[e] = 10.
    rcap = Residuals( g, cap )
        
        
    
    
    
    
    
    
    
    