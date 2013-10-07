
import math
import numpy as np

from mygraph import mygraph
from dijkstra import Dijkstra




""" Utility Algorithms """

def Residual( graph, capacity, flow=None ) :
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
    
    
def ResidualGraph( rcapacity, network, Delta=None ) :
    rgraph = mygraph()
    for e in network.edges() :
        ResidualGraphEdge( rgraph, e, rcapacity, network, Delta )
    return rgraph
    
    
def ResidualGraphEdge( rgraph, e, rcapacity, network, Delta=None ) :
    if Delta is None : Delta = 0.
    
    i,j = network.endpoints(e)
    iter = zip( [ (e,+1), (e,-1) ], [ (i,j), (j,i) ] )
    
    for label, (ii,jj) in iter :
        if rgraph.has_edge(label) : rgraph.remove_edge(label)
        u = rcapacity[label]
        if u > Delta : rgraph.add_edge( label, ii, jj )
        
        
        


def Excess( flow, graph, supply ) :
    excess = {}
    for i in graph.nodes() :
        excess[i] = supply[i]
        
        for e in graph.W[i] :   # edges in
            excess[i] += flow.get(e, 0. )
        for e in graph.V[i] :
            excess[i] -= flow.get(e, 0. )
            
    return excess


def ReducedCost( cost, potential, graph ) :
    rcost = {}
    for e in graph.edges() :
        i,j = graph.endpoints(e)
        rcost[e] = cost[e] + potential.get(j,0.) - potential.get(i,0.)
        
    return rcost


def LinearizeCost( lincost, flow, Delta, network, cost, edge=None ) :
    if edge is None :
        iter = network.edges()
    else :
        iter = [ edge ]
    
    for e in iter :
        x = flow.get(e, 0. )
        cc = cost[e]
        for dir in [ +1, -1 ] :
            lincost[(e,dir)] = cc( x + dir * Delta ) / Delta
                



""" Convex Cost Flow Algorithm """



def MinCostConvexFlow( network, capacity, supply, cost, epsilon=None ) :
    """
    graph is a mygraph (above)
    capacity is a dictionary from E -> real capacities
    supply is a dictionary from E -> real supplies
    cost is a dictionary from E -> lambda functions of convex cost edge costs
    """
    if epsilon is None : epsilon = 1
    
    flow = { e : 0. for e in network.edges() }
    potential = { i : 0. for i in network.nodes() }
    
    U = sum([ abs(b) for b in supply.values() ])
    temp = math.floor( math.log(U,2) )
    Delta = 2.**temp
    
    print 'total supply: %f' % U
    print 'Delta: %d' % Delta
    
    # initialize of intermediate data
    residual = Residual( network, capacity, flow )
    lincost = {}
    
    while Delta >= epsilon :
        print 'new phase: Delta=%f' % Delta
        
        # Delta has changed, so we need to [re-] linearize the costs
        # also will need to update relevant edge, after every augmentation
        LinearizeCost( lincost, flow, Delta, network, cost )    # all edges
        rgraph = ResidualGraph( residual, network )
        redcost = ReducedCost( lincost, potential, rgraph )
        
        # for every arc (i,j) in the residual network G(x)
        for resedge in rgraph.edges() :
            e,dir = resedge
            rcap = residual[resedge]
            # if residual capacity >= Delta and reduced cost < 0
            if rcap >= Delta and redcost[resedge] < 0. :
                pattern = 'residual capacity %f >= %f on residual edge %s with reduced cost %f'
                strdata = ( rcap, Delta, resedge, redcost[resedge] )
                print pattern % strdata
                
                # send rij flow along arc
                flow[e] += dir * rcap
                AugmentCapacity( residual, e, dir * rcap )
                LinearizeCost( lincost, flow, Delta, network, cost, e )     # just modify one edge
                
        # while there are imbalanced nodes
        while True :
            excess = Excess( flow, network, supply )
            print 'excess: %s' % repr(excess)
            
            SS = [ i for i,ex in excess.iteritems() if ex >= Delta ]
            TT = [ i for i,ex in excess.iteritems() if ex <= -Delta ]
            print 'surplus nodes: %s' % repr( SS )
            print 'deficit nodes: %s' % repr( TT )
            if len( SS ) <= 0 : break
            
            s = SS[0] ; t = TT[0]
            print 'augmenting %s to %s' % ( repr(s), repr(t) )
            
            Delta_rgraph = ResidualGraph( residual, network, Delta )
            print Delta_rgraph
            dist, _ = Dijkstra( Delta_rgraph, redcost, s )
            print 'Dijkstra shortest paths: %s' % repr( dist )
            
            # update the potentials and reduced costs
            for i in network.nodes() : potential[i] -= dist.get(i, np.Inf )
            redcost = ReducedCost( lincost, potential, rgraph )
            
            print dist
            print redcost
            
            # find shortest path w.r.t. reduced costs
            PATH = [] ; curr = t
            while curr is not s :
                UP = [ e for e in Delta_rgraph.W[curr] if redcost[e] <= 0. ]
                edge = UP[0]
                i,j = Delta_rgraph.endpoints(edge)
                
                PATH.insert( 0, edge )
                curr = i
                
            # augment Delta flow along the path P
            for edge in PATH :
                e,dir = edge
                flow[e] += dir * Delta
                AugmentCapacity( residual, e, dir * Delta )
                LinearizeCost( lincost, flow, Delta, network, cost, e )
                
            # update stuff; but not really... 
            # flow already done; SS, TT, and Delta_rgraph at next loop begin
            
        # end the iteration
        if Delta <= 1 : break
        Delta = Delta / 2
    
    return flow






if __name__ == '__main__' :
    import itertools
    import networkx as nx
    
    """
    convert linear instances on non-multi graphs to networkx format
    for comparison against nx.min_cost_flow() algorithm
    """
    def mincostflow_nx( network, capacity, supply, weight ) :
        digraph = nx.DiGraph()
        for i in network.nodes() :
            digraph.add_node( i, demand=-supply.get(i, 0. ) )
            
        for e in network.edges() :
            i,j = network.endpoints(e)
            digraph.add_edge( i, j, capacity=capacity.get(e, np.Inf ), weight=weight.get(e, 1. ) )
        return digraph
    
    
    
    
    
    g = mygraph()
    g.add_edge( 'a', 0, 1 )
    g.add_edge( 'b', 1, 2 )
    g.add_edge( 'c', 2, 3 )
    g.add_edge( 'd', 3, 0 )
    
    u = { e : 10. for e in g.edges() }
    supply = { 0 : 1., 1 : 2., 2 : -3., 3 : 0. }
    c = { 'a' : 10., 'b' : 5., 'c' : 1., 'd' : .5 }
    cf = {}
    
    for e in c : cf[e] = lambda x : c[e] * x
    
    flow = MinCostConvexFlow( g, u, supply, cf, epsilon=.01 )
    
    digraph = mincostflow_nx( g, u, supply, c )
    compare = nx.min_cost_flow( digraph )
    
    
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
    rcap = Residual( g, cap )
        
    rg = ResidualGraph( rcap, g )
    
    
    
    
    
    
    
    