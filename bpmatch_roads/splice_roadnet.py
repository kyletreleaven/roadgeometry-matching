
""" migrate this to vehicle routing at some point """

import itertools

import numpy as np
import matplotlib.pyplot as plt
plt.close('all')

import networkx as nx
import bintrees

import roadmap_basic as ROAD
import astar_basic as ASTAR
import bm_roadnet




""" utility """

def get_road_data( road, roadnet ) :
    for _,__,key, data in roadnet.edges_iter( keys=True, data=True ) :
        if key == road : return data




""" ALGORITHM BEGINS """
class demand :
    def __init__(self, pick, delv ) :
        self.pick = pick
        self.delv = delv

def ROADSSPLICE( demands, roadnet ) :
    QQ = [ dem.pick for dem in demands ]
    PP = [ dem.delv for dem in demands ]
    
    match = bm_roadnet.ROADSBIPARTITEMATCH( PP, QQ, roadnet )
    #cycles = CYCLEFACTOR( match )
    tour = DOUBLETOUR( roadnet )
    order = ORDERPOINTS( QQ, tour )
    
    res = []
    M = dict( match )
    
    #print order
    #print M
    
    while len( order ) > 0 :
        i = order[0]
        #print 'start at %d', i
        while i in M :
            res.append( i )
            #print i
            
            j = M[i]
            del M[i]
            order.remove( i )
            i = j
            
    return res



def CYCLEFACTOR( match ) :
    res = []
    
    M = dict(match)
    while len( M ) > 0 :
        curr = []
        i = M.iterkeys().next()
        while i in M :
            curr.append( i )
            j = M[i]
            del M[i]
            i = j
            
        res.append( curr )
    
    return res


class traversal :
    def __init__(self, road, forward ) :
        self.road = road
        self.forward = forward
        
    def __repr__(self) :
        dir = 'forward'
        if not self.forward : dir = 'reverse'
        return '(%s,%s)' % ( self.road, dir )

def DOUBLETOUR( roadnet ) :
    eulerian = nx.MultiDiGraph()
    for u,v, road in roadnet.edges_iter( keys=True ) :
        eulerian.add_edge( u,v, label=traversal( road, True ) )
        eulerian.add_edge( v,u, label=traversal( road, False ) )
        
    tour = []
    walk = [ u for u in nx.eulerian_circuit( eulerian ) ]
    for u, v in walk :
        edges = eulerian.edge[u][v]
        key = edges.keys()[0]       # get key of the some (e.g., first) edge from u, v
        data = eulerian.get_edge_data( u, v, key )
        tour.append( data.get('label') )
        eulerian.remove_edge( u, v, key )
        
    return tour
        
        
def ARRANGEMENT( points ) :
    def ensure_road( road, data ) :
        curr = data.setdefault( road )
        if curr is None : data[road] = bintrees.RBTree()
        return data[road]

    def ensure_key( key, tree ) :
        curr = tree.set_default( key )
        if curr is None : tree[key] = []
        return tree[key]

    TREES = dict()
    for idx, (road,y) in enumerate( points ) :
        tree = ensure_road( road, TREES )
        Q = ensure_key( y, tree )
        Q.append( idx )
        
    return TREES

def quantity( arr ) :
    num = 0
    for k, tree in arr.iteritems() :
        Z = [ len(q) for y,q in tree.iter_items() ]
        num += sum( Z )
    return num


def COLLECT( arrangement, tour ) :
    order = []
    
    arr = arrangement.copy()
    #print quantity( arr ) 
    
    for t in tour :
        if t.road not in arr : continue
        
        tree = arr[ t.road ]
        if t.forward :
            items = tree.iter_items()
        else :
            items = tree.iter_items( reverse=True )
        items = [ i for i in items ]        # because it *wasn't* an iterator
            
        for key, Q in items :
            order.extend( Q )
            tree.remove(key)
            
    return order


def ORDERPOINTS( points, tour ) :
    arr = ARRANGEMENT( points )
    order = COLLECT( arr, tour )
    return order








""" MATCH TESTING Utilities """

def CYCLECOST( cycle, demands, roadnet ) :
    addr = lambda p : ROAD.RoadAddress( *p )
    
    carry = [ ROAD.distance( roadnet, addr(dem.pick), addr(dem.delv) ) for dem in demands ]
    
    edges = zip( cycle, cycle[1:] + cycle[:1] )
    
    print edges
    
    edges = [ ( demands[i], demands[j] ) for i,j in edges ] 
    empty = [ ROAD.distance( roadnet, addr(dem1.delv), addr(dem2.pick) ) for dem1,dem2 in edges ]
    return sum( carry ) + sum( empty )


def MATCHCOST( match, demands, roadnet ) :
    cycles = CYCLEFACTOR( match )
    costs = [ CYCLECOST( cycle, demands, roadnet ) for cycle in cycles ]
    return sum( costs )






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
            
        self.road_sampler = WeightedSet( weight_dict )
        
    def sample(self) :
        road = self.road_sampler.sample()
        L = get_road_data( road, roadnet ).get( 'length', 1 )
        y = L * np.random.rand()
        return (road,y)
    
    



if __name__ == '__main__' :
    
    roadnet = nx.MultiDiGraph()
    roadnet.add_edge( 0,1, 'N', length=1. )
    roadnet.add_edge( 1,2, 'E', length=1. )
    roadnet.add_edge( 2,3, 'S', length=1. )
    roadnet.add_edge( 3,0, 'W', length=1. )
    
    
    tour = DOUBLETOUR( roadnet )
    print tour
    
    sampler = UniformDist( roadnet )
    
    NUMPOINT = 100
    ZZ = np.linspace(-NUMPOINT,NUMPOINT,1000)
    #
    PP = [ sampler.sample() for i in xrange(NUMPOINT) ]
    QQ = [ sampler.sample() for i in xrange(NUMPOINT) ]
    DEMANDS = [ demand(p,q) for p,q in zip( PP, QQ ) ]
    
    tour = ROADSSPLICE( DEMANDS, roadnet )
    cost = CYCLECOST( tour, DEMANDS, roadnet )
    
    
    #order = ORDERPOINTS( PP, tour )
    #point_tour = [ PP[i] for i in order ]
    #print point_tour
    
    
    
    
    #match = ROADSBIPARTITEMATCH( PP, QQ, roadnet )
    #costs = MATCHCOSTS( match, PP, QQ, roadnet )
    #cost = ROADMATCHCOST( match, PP, QQ, roadnet )
    
    
            