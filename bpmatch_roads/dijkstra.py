
import numpy as np

from mygraph import mygraph
from priodict import *


def Dijkstra( graph, cost, s ) :
    d = { s : 0. }      # only in here if they are seen... duh!!
    pred = { s : None }
    
    OPEN = priorityDictionary()
    OPEN[s] = 0.
    
    while len( OPEN ) > 0 :
        i = OPEN.smallest()
        del OPEN[i]
        
        for e in graph.V[i] :
            _, j = graph.endpoints(e)
            new_d = d[i] + cost.get( e, np.Inf )
            
            if j not in d : OPEN[j] = new_d
            if new_d < d.get( j, np.Inf ) :
                d[j] = new_d
                pred[j] = i
                
    return d, pred
                
                




""" unit test """

if __name__ == '__main__' :
    import itertools
    
    import numpy as np
    import networkx as nx
    
    g = mygraph()
    g.add_edge( 'a', 0, 1 )
    g.add_edge( 'b', 0, 2 )
    g.add_edge( 'c', 2, 3 )
    g.add_edge( 'd', 3, 0 )
    
    c = { 'a' : 10., 'b' : 5., 'c' : 1., 'd' : .5 }
    f = { 'a' : 3., 'b' : 1.37 }


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
    
    
    
    
    
    
    
    