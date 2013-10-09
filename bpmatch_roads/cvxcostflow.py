
import math
import numpy as np

from mygraph import mygraph
from dijkstra import Dijkstra
from toposort import toposort




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
    
    
def ResidualGraph( rgraph, rcapacity, network, Delta=None, edge=None ) :
    if Delta is None : Delta = 0.
    if edge is None :
        iter = network.edges()
    else :
        iter = [ edge ]
        
    for e in iter :
        i,j = network.endpoints(e)
        iter = zip( [ (e,+1), (e,-1) ], [ (i,j), (j,i) ] )
        
        for label, (ii,jj) in iter :
            if rgraph.has_edge(label) : rgraph.remove_edge(label)
            u = rcapacity[label]
            if u > Delta : rgraph.add_edge( label, ii, jj )
            
        
        
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






def ReducedCost( cost, potential, graph ) :
    rcost = {}
    for e in graph.edges() :
        i,j = graph.endpoints(e)
        # if node is not in potential, it is assumed to have -inf potential
        rcost[e] = cost[e] + potential.get(j,0.) - potential.get(i,0.)
        #rcost[e] = cost[e] + potential[j] - potential[i]
        
    return rcost


def LinearizeCost( lincost, flow, Delta, network, cost, edge=None ) :
    if edge is None :
        iter = network.edges()
    else :
        iter = [ edge ]
    
    for e in iter :
        x = flow.get(e, 0. )
        cc = cost.get( e, line(0.) )
        for dir in [ +1, -1 ] :
            lincost[(e,dir)] = ( cc( x + dir * Delta ) - cc(x) ) / Delta
            
            
            
            

def Excess( flow, graph, supply ) :
    excess = {}
    for i in graph.nodes() :
        excess[i] = supply.get(i, 0. )
        
        for e in graph.W[i] :   # edges in
            excess[i] += flow.get(e, 0. )
        for e in graph.V[i] :
            excess[i] -= flow.get(e, 0. )
            
    return excess

            
            
"""
an alternative potential update routine by myself, which
(perhaps) more gracefully handles [negatively] infinite potentials
at previous stages (owing to dis-connectivity).
Only updates potentials at the nodes in upstream.
The costs should be the *NON*-reduced costs!

Missing nodes are un-reachable, and therefore have -inf potential implicitly.
"""
def UpdatePotential( potential, upstream, graph, cost ) :
    # construct the ancestor graph
    succ_graph = mygraph()
    for child, umbilical in upstream.iteritems() :
        if umbilical is None : continue
        parent,_ = graph.endpoints(umbilical)
        succ_graph.add_edge( umbilical, parent, child )
    
    # obtain a topological ordering of the nodes
    order = toposort( succ_graph )
    
    # update potentials
    for child in order :
        e = upstream[child]
        # if it's the top node, set potential to 0.
        if e is None :
            potential[child] = 0.
        # otherwise, ensure it's parent potential minus (non-reduced cost of shortest path)
        else :
            parent,_ = graph.endpoints(e)
            potential[child] = potential[parent]
            potential[child] -= cost[e]
            


""" Convex Cost Flow Algorithm """
    
def MinCostConvexFlow( network, capacity, supply, cost, epsilon=None ) :
    """
    see "Fragile" version;
    this wrapper adds robustness:
    pre-processes to ensure strong connectivity of *any* Delta-residual graph;
    edge weights should be prohibitively expensive
    """
    # pre-process
    REGULAR = 0
    AUGMENTING = 1
    
    network_aug = mygraph()
    capacity_rename = {}
    cost_aug = {}
    
    for e in network.edges() :
        i,j = network.endpoints(e)
        newedge = (REGULAR,e)
        
        network_aug.add_edge( newedge, i, j )
        if e in capacity : capacity_rename[ newedge ] = capacity[e]
        if e in cost : cost_aug[ newedge ] = cost[e]
        
    U = sum([ abs(b) for b in supply.values() ])
    # since costs are convex, flow cannot have cost greater than M
    M = max([ c(U) for c in cost.values() ])
    prohibit = line(M)

    # add a doubly-linked chain
    NODES = network.nodes()
    edgegen = itertools.count()
    for i,j in zip( NODES[:-1], NODES[1:] ) :
        frwd = (AUGMENTING, edgegen.next() )
        network_aug.add_edge( frwd, i, j )
        cost_aug[frwd] = prohibit
        
        bkwd = (AUGMENTING, edgegen.next() )
        network_aug.add_edge( bkwd, j, i )
        cost_aug[bkwd] = prohibit
    
    # run the "fragile" version
    flow = FragileMCCF( network_aug, capacity_rename, supply, cost_aug, epsilon )
    
    # prepare output --- perhaps do some feasibility checking in the future
    res = { e : x for (type,e), x in flow.iteritems() if type == REGULAR }
    return res
    #
    
    
    
def FragileMCCF( network, capacity, supply, cost, epsilon=None ) :
    """
    network is a mygraph (above)
    capacity is a dictionary from E -> real capacities
    supply is a dictionary from E -> real supplies
    cost is a dictionary from E -> lambda functions of convex cost edge costs
    
    assumes every Delta-residual graph is strongly connected,
    i.e., there exists a path with inf capacity b/w any two nodes;
    """
    if epsilon is None : epsilon = 1
    
    U = sum([ abs(b) for b in supply.values() ])
    temp = math.floor( math.log(U,2) )
    Delta = 2.**temp
    
    print 'total supply: %f' % U
    print 'Delta: %d' % Delta
    
    """ ALGORITHM """
    flow = { e : 0. for e in network.edges() }
    potential = { i : 0. for i in network.nodes() }
    
    # initialize algorithm data
    residual = Residual( network, capacity, flow )
    rgraph = mygraph()
    Delta_rgraph = mygraph()
    
    lincost = {}
    
    while Delta >= epsilon :
        print '\nnew phase: Delta=%f' % Delta
        
        # Delta is fresh, so we need to [re-] linearize the costs
        # also will need to update relevant edge, after every augmentation
        LinearizeCost( lincost, flow, Delta, network, cost )    # all edges
        print 'local costs: %s' % repr( lincost )
        
        print 'residual cap: %s' % repr( residual )
        ResidualGraph( rgraph, residual, network )
        print 'residual graph: %s' % rgraph.edges()
        
        redcost = ReducedCost( lincost, potential, rgraph )
        print 'reduced cost, phase begin: %s' % repr( redcost )
        
        # for every arc (i,j) in the residual network G(x)
        for resedge in rgraph.edges() :
            e,dir = resedge
            # if residual capacity >= Delta and reduced cost < 0
            # modification WHILE: augment by Delta, not by residual capacity
            # b/c, we must allow re-linearization between each augmentation
            while True :
                rcap = residual[resedge]
                if not ( rcap >= Delta and redcost[resedge] < 0 ) : break
                
                pattern = 'residual capacity %f >= %f on residual edge %s with reduced cost %f'
                strdata = ( rcap, Delta, resedge, redcost[resedge] )
                print pattern % strdata
                
                # send rij flow along arc --- modification: Delta flow
                flow[e] += dir * Delta
                AugmentCapacity( residual, e, dir * rcap )
                ResidualGraph( rgraph, residual, network, edge=e )
                LinearizeCost( lincost, flow, Delta, network, cost, e )     # just modify one edge
                
                redcost = ReducedCost( lincost, potential, rgraph )
                print 'reduced cost, correction: %s' % repr( redcost )
                
        # while there are imbalanced nodes
        while True :
            print 'flow: %s' % repr( flow )
            excess = Excess( flow, network, supply )
            print 'excess: %s' % repr(excess)
            
            SS = [ i for i,ex in excess.iteritems() if ex >= Delta ]
            TT = [ i for i,ex in excess.iteritems() if ex <= -Delta ]
            print 'surplus nodes: %s' % repr( SS )
            print 'deficit nodes: %s' % repr( TT )
            if len( SS ) <= 0 : break
            
            s = SS[0] ; t = TT[0]
            print 'augmenting %s to %s' % ( repr(s), repr(t) )
            
            print 'residual: %s' % repr( residual )
            ResidualGraph( Delta_rgraph, residual, network, Delta )
            print 'Delta residual graph: %s' % repr( Delta_rgraph.edges() )
            
            print 'potential: %s' % repr( potential )
            
            redcost = ReducedCost( lincost, potential, rgraph )
            print 'reduced cost, shortest paths: %s' % repr( redcost )
            
            dist, upstream = Dijkstra( Delta_rgraph, redcost, s )
            print 'Dijkstra shortest paths: %s' % repr( dist )
            print 'Dijkstra upstreams: %s' % repr( upstream )
            
            # update the potentials and reduced costs; by connectivity, should touch *every* node
            for i in Delta_rgraph.nodes() : potential[i] -= dist[i]
            #UpdatePotential( potential, upstream, Delta_rgraph, lincost )
            print 'new potentials: %s' % repr( potential )
            
            redcost = ReducedCost( lincost, potential, rgraph )
            print 'new reduced costs: %s' % repr( redcost )
            
            # find shortest path w.r.t. reduced costs
            PATH = [] ; j = t
            while j is not s :
                e = upstream[j]
                i,_ = Delta_rgraph.endpoints(e)
                PATH.insert( 0, e )
                j = i
            print 'using path: %s' % repr( PATH )
            
            # augment Delta flow along the path P
            for edge in PATH :
                e,dir = edge
                flow[e] += dir * Delta
                AugmentCapacity( residual, e, dir * Delta )
                
                ResidualGraph( Delta_rgraph, residual, network, edge=e )
                LinearizeCost( lincost, flow, Delta, network, cost, e )
                
            # update stuff; but not really... 
            # flow already done; SS, TT, and Delta_rgraph at next loop begin
            
        # end the iteration
        if Delta <= epsilon : break
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
    
    if False :
        g.add_edge( 'a', 0, 1 )
        g.add_edge( 'b', 1, 2 )
        g.add_edge( 'c', 2, 3 )
        g.add_edge( 'd', 3, 0 )
        
        u = { e : 10. for e in g.edges() }
        supply = { 0 : 1., 1 : 2., 2 : -3., 3 : 0. }
        c = { 'a' : 10., 'b' : 5., 'c' : 1., 'd' : .5 }
    else :
        u = {}
        c = {}
        s = {}
        
        s[0] = 10.
        
        g.add_edge( 'a', 0, 1 )
        c['a'] = 1.
        u['a'] = 1.35
        
        #g.add_edge( 'aprime', 0, 1 )
        #c['aprime'] = 1000.
        
        g.add_edge( 'b', 0, 2 )
        c['b'] = 10.
        
        g.add_edge( 'c', 1, 3 )
        g.add_edge( 'd', 2, 3 )
        
        s[3] = -10.
        supply = s
    
    cf = {}
    class line :
        def __init__(self, slope ) :
            self.m = slope
        def __call__(self, x ) :
            return self.m * x
        
    for e in c :
        cf[e] = line( c[e] )
    
    flow = MinCostConvexFlow( g, u, supply, cf, epsilon=.001 )
    
    digraph = mincostflow_nx( g, u, supply, c )
    compare = nx.min_cost_flow( digraph )
    
    
    
    
    