

import itertools

import numpy as np
import bintrees

import networkx as nx

""" my dependencies """
import setiptah.roadgeometry.roadmap_basic as ROAD
import setiptah.roadgeometry.astar_basic as ASTAR

from nxflow.capscaling import SOLVER


# to construct the optimization problem


""" utility """

def get_road_data( road, roadnet ) :
    for _,__,key, data in roadnet.edges_iter( keys=True, data=True ) :
        if key == road : return data

        
        
        


""" ALGORITHM HIGH LEVEL """

def ROADSBIPARTITEMATCH( P, Q, roadnet ) :
    MATCH = []
    
    segment_dict = SEGMENTS( P, Q, roadnet )
    surplus_dict = dict()
    objective_dict = dict()
    
    for road, segment in segment_dict.iteritems() :
        match = PREMATCH( segment )
        MATCH.extend( match )
        
        surplus_dict[road] = SURPLUS( segment )
        
        roadlen = get_road_data( road, roadnet ).get( 'length', 1 )
        measure = MEASURE( segment, roadlen )
        objective_dict[road] = OBJECTIVE( measure )
        #objective_dict[road] = objective
        
    assist = SOLVER( roadnet, surplus_dict, objective_dict )
    
    topograph = TOPOGRAPH( segment_dict, assist, roadnet )
    
    match = TRAVERSE( topograph )
    MATCH.extend( match )
    
    return MATCH


def WRITEOBJECTIVES( P, Q, roadnet ) :
    MATCH = []
    
    segment_dict = SEGMENTS( P, Q, roadnet )
    surplus_dict = dict()
    objective_dict = dict()
    
    for road, segment in segment_dict.iteritems() :
        match = PREMATCH( segment )
        MATCH.extend( match )
        
        surplus_dict[road] = SURPLUS( segment )
        
        roadlen = get_road_data( road, roadnet ).get( 'length', 1 )
        measure = MEASURE( segment, roadlen )
        objective_dict[road] = OBJECTIVE( measure )
        #objective_dict[road] = objective
        
    return objective_dict


""" Algorithm Utilities """

def ensure_road( road, data ) :
    curr = data.setdefault( road )
    if curr is None : data[road] = bintrees.RBTree()
    return data[road]

class TwoQueues() :
    def __init__(self) :
        self.P = []
        self.Q = []
        
    def __repr__(self) :
        return '<P:%s,Q:%s>' % ( repr(self.P), repr(self.Q) )

def ensure_key( key, tree ) :
    curr = tree.set_default( key )
    if curr is None : tree[key] = TwoQueues()
    return tree[key]

class LineData :
    def __init__(self, m,b) :
        self.slope = m
        self.offset = b
        
    def __call__(self, x ) :
        return self.slope * x + self.offset
    
    def __repr__(self) :
        return '<%f z + %f>' % ( self.slope, self.offset )

class terminal :    # simple node type for TRAVERSE
    def __init__(self, q ) :
        self.q = q



""" ALGORITHM SUB-ROUTINES """

def SEGMENTS( P, Q, roadnet ) :
    """
    returns:
    a dictionary whose keys are coordinates and whose values are local (P,Q) index queues 
    """
    segments = dict()
    for _,__,road in roadnet.edges_iter( keys=True ) :
        ensure_road( road, segments )   # these, and only these, roads are allowed
        
    for i, p in enumerate( P ) :
        r,y = p
        tree = segments[r]      # crash by design if r not in segments
        queues = ensure_key( y, tree )
        queues.P.append( i )
        
    for j, q in enumerate( Q ) :
        r,y = q
        tree = segments[r]
        queues = ensure_key( y, tree )
        queues.Q.append( j )
        
    return segments
    
    
def PREMATCH( segment ) :
    match = []
    for y, q in segment.iter_items() :
        annih = min( len( q.P ), len( q.Q ) )
        for k in range( annih ) :
            i = queues.P.pop(0)
            j = queues.Q.pop(0)
            match.append( (i,j) )
            
    return match


def SURPLUS( segment ) :
    deltas = [ len( q.P ) - len( q.Q ) for y,q in segment.iter_items() ]
    return sum( deltas )


def MEASURE( segment, length, rbound=None ) :
    if rbound is not None :
        lbound = length
    else :
        lbound = 0.
        rbound = length
        
    # bintree instead of dict so that it is enumerated in sorted order
    measure = bintrees.RBTree()
    
    posts = [ lbound ] + [ y for y,q in segment.iter_items() ] + [ rbound ]
    intervals = zip( posts[:-1], posts[1:] )
    
    deltas = [0] + [ len(q.P)-len(q.Q) for y,q in segment.iter_items() ]
    F = np.cumsum( deltas )
    
    for (a,b), f in zip( intervals, F ) :
        measure.setdefault( f, 0. )
        measure[f] += b - a
        
    return measure


def INTERVALS( segment ) :      # very similar routine, used to build the walk graph
    res = dict()
    
    posts = [ '-' ] + [ y for y,q in segment.iter_items() ] + [ '+' ]
    intervals = zip( posts[:-1], posts[1:] )
    
    deltas = [0] + [ len(q.P)-len(q.Q) for y,q in segment.iter_items() ]
    F = np.cumsum( deltas )
    
    for I, f in zip( intervals, F ) :
        res.setdefault( f, [] )
        res[f].append( I )
        
    return res


def EDGES( segment ) :      # very similar routine, used to build the walk graph
    edges = dict()
    
    posts = [ '-' ] + [ q for y,q in segment.iter_items() ] + [ '+' ]
    posts = [ terminal(q) for q in posts ]
    intervals = zip( posts[:-1], posts[1:] )
    
    deltas = [0] + [ len(q.P)-len(q.Q) for y,q in segment.iter_items() ]
    F = np.cumsum( deltas )
    
    for I, f in zip( intervals, F ) :
        edges.setdefault( f, [] )
        edges[f].append( I )
        
    return edges


def OBJECTIVE( measure ) :
    # prepare constants kappa and alpha
    ff = np.array( [ -np.inf ] + [ f for f in measure ] )
    PREALPHA = np.array( [ 0. ] + [ w for f,w in measure.items() ] )
    PREKAPPA = np.array( [ 0. ] + [ f*w for f,w in measure.items() ] )
    
    ALPHAM = np.cumsum( PREALPHA )
    ALPHAP = ALPHAM[-1] - ALPHAM
    ALPHA = ALPHAP - ALPHAM
    
    KAPPAM = np.cumsum( PREKAPPA )
    KAPPAP = KAPPAM[-1] - KAPPAM
    KAPPA = KAPPAP - KAPPAM
    
    Czminus = bintrees.RBTree()
    for f, kappa, alpha in zip( ff, KAPPA, ALPHA ) :
        Czminus.insert( f, LineData( alpha, kappa ) )
        
    return Czminus










def TOPOGRAPH( segment_dict, assist, roadnet ) :
    topograph = nx.DiGraph()
    
    special = dict()
    for u in roadnet.nodes_iter() :
        data = TwoQueues()
        node = terminal( data )
        special[u] = node
        
    for u,v, road, data in roadnet.edges_iter( keys=True, data=True ) :
        segment = segment_dict[road]
        z = assist[road]
        
        edges = EDGES( segment )
        for f, intervals in edges.iteritems() :
            h = f + z
            for (ll,rr) in intervals :
                if ll.q == '-' : ll = special[u]
                if rr.q == '+' : rr = special[v]
                
                if h > 0 :
                    topograph.add_edge( ll, rr, weight = h )
                if h < 0 :
                    topograph.add_edge( rr, ll, weight = -h )
                    
    return topograph
    
    
def TRAVERSE( topograph ) :
    match = []
    nodes_ord = nx.topological_sort( topograph )
    
    LISTS = dict()
    for u in nodes_ord : LISTS[u] = []
    
    for u in nodes_ord :
        L = LISTS[u]
        queue = u.q
        
        if len( queue.P ) > 0 : L.extend( queue.P )
        
        for k in range( len( queue.Q ) ) :
            i = L.pop(0)
            match.append( (i,queue.Q[k] ) )
        
        for _,v, data in topograph.out_edges_iter( u, data=True ) :
            w = data.get('weight')
            pre, post = L[:w], L[w:]
            LISTS[v].extend( pre )
            L = post
            
    return match







""" MATCH TESTING Utilities """

def MATCHCOSTS( match, P, Q, roadnet ) :
    costs = []
    for i, j in match :
        p = ROAD.RoadAddress( *P[i] )
        q = ROAD.RoadAddress( *Q[j] )
        d = ROAD.distance( roadnet, p, q, 'length' )
        costs.append( d )
    return costs

def ROADMATCHCOST( match, P, Q, roadnet ) :
    costs = MATCHCOSTS( match, P, Q, roadnet )
    return sum( costs )





class costWrapper :
    def __init__(self, lines ) :
        self.lines = lines
        
    def __call__(self, z ) :
        """ this is an O(log n) query function, can be reduced to O(1) by random access with saturation """
        _, line = self.lines.floor_item( -z )
        return line( z )



def drawCBounds( ZZ, Ctree, ax=None ) :
    if ax is None :
        plt.figure()
        ax = plt.gca()
        
    cost = costWrapper( Ctree )
    Cref = np.array([ cost( z ) for z in ZZ ]) 
    ax.plot( ZZ, Cref, c='b', linewidth=3, alpha=.25 )
    for f, (kappa, alpha) in Ctree.iter_items() :
        C = kappa + alpha * ZZ
        ax.plot( ZZ, C, c='k', linestyle='--' )
    #ax.set_aspect( 'equal' )
    return ax







""" convenient sampling utility for the unit test below, might as well be a package-export, though """

class WeightedSet :
    def __init__(self, weight_dict ) :
        """
        keys are targets, values are weights; needn't sum to 1
        doesn't check for repeats
        """
        targets = weight_dict.keys()
        weights = weight_dict.values()
        scores = np.cumsum( np.array( weights ) )
        
        self._hiscore = scores[-1]
        self._tree = bintrees.RBTree()
        for target, score in zip( targets, scores ) :
            self._tree.insert( score, target )
            
    def sample(self) :
        z = self._hiscore * np.random.rand()
        _, res = self._tree.ceiling_item( z )
        return res
    
    
class UniformDist :
    def __init__(self, roadnet=None, length=None ) :
        if roadnet is not None :
            self.set_roadnet( roadnet, length )
        
    def set_roadnet(self, roadnet, length=None ) :
        if length is None : length = 'length'
        
        weight_dict = dict()
        for _,__, road, data in roadnet.edges_iter( keys=True, data=True ) :
            weight_dict[road] = data.get( length, 1 )
            
        self.roadnet = roadnet
        self.road_sampler = WeightedSet( weight_dict )
        
    def sample(self) :
        road = self.road_sampler.sample()
        L = get_road_data( road, self.roadnet ).get( 'length', 1 )
        y = L * np.random.rand()
        return (road,y)
    
    
    
    
if __name__ == '__main__' :
    import matplotlib.pyplot as plt
    plt.close('all')

    roadnet = nx.MultiDiGraph()
    roadnet.add_edge( 0,1, 'N', length=1. )
    roadnet.add_edge( 1,2, 'E', length=1. )
    roadnet.add_edge( 2,3, 'S', length=1. )
    roadnet.add_edge( 3,0, 'W', length=1. )
    
    sampler = UniformDist( roadnet )
    
    NUMPOINT = 50
    ZZ = np.linspace(-NUMPOINT,NUMPOINT,1000)
    #
    PP = [ sampler.sample() for i in xrange(NUMPOINT) ]
    QQ = [ sampler.sample() for i in xrange(NUMPOINT) ]
    
    z = np.arange(-NUMPOINT,NUMPOINT,.25)
    objective_dict = WRITEOBJECTIVES( PP, QQ, roadnet )
    for road, Cz in objective_dict.iteritems() :
        cost = costWrapper( Cz )
        C = [ cost(zz) for zz in z ]
        plt.figure()
        plt.plot( z, C, '--', marker='x' )
        
        
    
    match = ROADSBIPARTITEMATCH( PP, QQ, roadnet )
    costs = MATCHCOSTS( match, PP, QQ, roadnet )
    cost = ROADMATCHCOST( match, PP, QQ, roadnet )
    print match
    print costs
    print cost
    
    
    # compare to optimal matching
    if True and NUMPOINT <= 50 :
        PR = [ ROAD.RoadAddress( road, y ) for road, y in PP ]
        QR = [ ROAD.RoadAddress( road, y ) for road, y in QQ ]
        
        class pointnode() :
            def __init__(self, point, idx ) :
                self.point = point
                self.idx = idx
                
        RED = [ pointnode(p,i) for i,p in enumerate(PR) ]
        BLU = [ pointnode(q,j) for j,q in enumerate(QR) ]
        graph = nx.Graph()
        match_mat = np.zeros((NUMPOINT,NUMPOINT))
        for (i,r), (j,b) in itertools.product( enumerate(RED), enumerate(BLU) ) :
            w = ROAD.distance( roadnet, r.point, b.point, 'length' )
            graph.add_edge( r, b, weight=-w )
            match_mat[i,j] = w
        match_dict = nx.max_weight_matching( graph, True )
        
        match_brute = [ (r.idx,match_dict[r].idx) for r in RED ]      # match pruning
#        matchstats = [ ( r.point, b.point, ROAD.distance( roadnet, r.point, b.point, 'length' ) )
#                      for r,b in match ]
        costs_brute = MATCHCOSTS( match_brute, PP, QQ, roadnet )
        cost_brute = ROADMATCHCOST( match_brute, PP, QQ, roadnet )
        print match_brute
        print costs_brute
        print cost_brute
        #print 'optimal match has cost: %f' % matchcost

        
        
    if False :       # validate CTREES
        zmin = -NUMPOINT/2
        zmax = NUMPOINT/2
        ZZZ = range( zmin, zmax )
        #
        for road, C in CTREES.iteritems() :
            ax = drawCBounds( ZZ, CTREES[road] )
            #plt.plot( ZZ, Cz, linestyle='--' )
            
            width = get_road_data( road, roadnet ).get( 'length', 1 )
            Cmatch = [ MATCHCOST( P[road], Q[road], width, z ) for z in ZZZ ]
            plt.scatter( ZZZ, Cmatch, marker='x' )
            
            
    if False :
            # validate C's of the different roads
            for road in CTREES :
                ax = drawCBounds( ZZ, CTREES[road] )
                ax.set_title('road=%s' % road )
                zr = assist[road].value
                cost = costWrapper( CTREES[road] )
                Cr = cost( zr )
                #ax.axvline( assist[road].value )
                ax.scatter( [ zr ], [ Cr ], marker='o' )
                ax.scatter( [ zr ], [ cost[road].value ], marker='x' )
                
            dC = dict()
            for road in assist :
                zr = assist[road].value
                f, (kappa,alpha) = CTREES[road].floor_item( -zr )
                dC[road] = alpha
                
            # compute a matching and verify cost
            matchZ = dict()
            for road, var in assist.iteritems() :
                matchZ[road] = int( round( var.value ) )
            the_match = ROADMATCH( PP, QQ, matchZ, roadnet )
            the_match_cost = ROADMATCHCOST( the_match, roadnet )
            print 'constructed matching has cost: %f' % the_match_cost
                
                

            
            
