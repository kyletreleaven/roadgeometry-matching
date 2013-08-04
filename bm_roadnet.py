
import itertools

import numpy as np
import matplotlib.pyplot as plt
plt.close('all')

import networkx as nx
import bintrees

import roadmap_basic as ROAD

# to construct the optimization problem
import cvxpy




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


""" matching """
def MATCH( P, Q, length, z=0 ) :
    PP = sorted(P)
    QQ = sorted(Q)
    if z > 0 : PP = [ 0. for i in xrange( z ) ] + PP
    if z < 0 : QQ = [ 0. for i in xrange( -z ) ] + QQ
    
    surp = len(PP) - len(QQ)
    if surp > 0 : QQ = QQ + [ length for i in xrange( surp ) ]
    if surp < 0 : PP = PP + [ length for i in xrange( -surp ) ]
    
    return zip(PP,QQ)

def MATCHCOST( P, Q, length, z=0 ) :
    return sum([ np.abs(q-p) for p,q in MATCH( P, Q, length, z) ])


def obtain_segment_objective( P, Q, length ) :
    # prepare the points index
    points = bintrees.RBTree()
    points.insert( 0., 0 )
    points.insert( length, 0 )
    
    # insert the problem data
    for p in P : points.insert( p, 1 )
    for q in Q : points.insert( q, -1 )
    
    # prepare F and assign intervals
    posts = [ y for y, sign in points.iter_items() ]
    signs = [ sign for y, sign in points.iter_items() ]
    pairs = zip( posts[:-1], posts[1:] )
    F = np.cumsum( signs[:-1] )
    #
    measure = bintrees.RBTree()
    for (a,b), f in zip( pairs, F ) :
        val = measure.setdefault( f, 0. )
        measure[f] = val + ( b - a )
        
    # prepare constants kappa and alpha
    #temp = sorted( measure.iteritems() )
    ff = np.array( [ -np.inf ] + [ f for f,w in measure.iter_items() ] )
    PREALPHA = np.array( [ 0. ] + [ w for f,w in measure.iter_items() ] )
    PREKAPPA = np.array( [ 0. ] + [ f*w for f,w in measure.iter_items() ] )
    
    ALPHAM = np.cumsum( PREALPHA )
    ALPHAP = ALPHAM[-1] - ALPHAM
    ALPHA = ALPHAP - ALPHAM
    
    KAPPAM = np.cumsum( PREKAPPA )
    KAPPAP = KAPPAM[-1] - KAPPAM
    KAPPA = KAPPAP - KAPPAM
    
    Czminus = bintrees.RBTree()
    for f, kappa, alpha in zip( ff, KAPPA, ALPHA ) :
        Czminus.insert( f, ( kappa, alpha ) )
        
    return Czminus




def evalC( z, Ctree ) :
    """ this is the O(log n) query function """
    _, ( kappa, alpha ) = Ctree.floor_item( -z )
    return kappa + alpha * z

def drawCBounds( ZZ, Ctree, ax=None ) :
    if ax is None :
        plt.figure()
        ax = plt.gca()
        
    Cref = np.array([ evalC( z, Ctree ) for z in ZZ ]) 
    ax.plot( ZZ, Cref, c='b', linewidth=3, alpha=.25 )
    for f, (kappa, alpha) in Ctree.iter_items() :
        C = kappa + alpha * ZZ
        ax.plot( ZZ, C, c='k', linestyle='--' )
    #ax.set_aspect( 'equal' )
    return ax





if __name__ == '__main__' :
    roadnet = nx.MultiDiGraph()
    roadnet.add_edge( 0,1, 'N', length=1. )
    roadnet.add_edge( 1,2, 'E', length=1. )
    roadnet.add_edge( 2,3, 'S', length=1. )
    roadnet.add_edge( 3,0, 'W', length=1. )
    
    def get_road_data( road, roadnet ) :
        for _,__,key, data in roadnet.edges_iter( keys=True, data=True ) :
            if key == road : return data
            
    weight_dict = dict()
    for _,__, road, data in roadnet.edges_iter( keys=True, data=True ) :
        weight_dict[road] = data.get( 'length', 1 )
        
    road_sampler = WeightedSet( weight_dict )
    def sample() :
        road = road_sampler.sample()
        L = get_road_data( road, roadnet ).get( 'length', 1 )
        y = L * np.random.rand()
        return (road,y)
    
    NUMPOINT = 10
    ZZ = np.linspace(-NUMPOINT,NUMPOINT,1000)
    def shatter( points, roadnet ) :
        P = dict()
        for _,__,road in roadnet.edges_iter( keys=True ) : P[road] = []
        for road,y in points : P[road].append( y )
        return P
    #
    PP = [ sample() for i in xrange(NUMPOINT) ]
    QQ = [ sample() for i in xrange(NUMPOINT) ]
    P = shatter( PP, roadnet )
    Q = shatter( QQ, roadnet )
    
    # compute the level intervals of C(z;road)
    CTREES = dict()
    for road in P :
        L = get_road_data( road, roadnet ).get( 'length', 1 )
        CTREES[road] = obtain_segment_objective( P[road], Q[road], L )
    
    """ construct the problem """
    surplus = dict()
    assist = dict()
    cost = dict()
    
    for _,__,road in roadnet.edges_iter( keys=True ) :
        surplus[road] = len( P[road] ) - len( Q[road] )
        assist[road] = cvxpy.variable( name='z_%s' % road )
        cost[road] = cvxpy.variable( name='c_%s' % road )
        
    objfunc = sum( cost.values() )
    OBJECTIVE = cvxpy.minimize( objfunc )
    
    CONSTRAINTS = []
    
    # the flow conservation constraints
    for u in roadnet.nodes_iter() :
        INFLOWS = []
        for _,__,road in roadnet.in_edges( u, keys=True ) :
            INFLOWS.append( assist[road] + surplus[road] )
            
        OUTFLOWS = []
        for _,__,road in roadnet.out_edges( u, keys=True ) :
            OUTFLOWS.append( assist[road] )
            
        #conserve_u = cvxpy.eq( sum(OUTFLOWS), sum(INFLOWS) )
        #
        error_u = sum(OUTFLOWS) - sum(INFLOWS)
        conserve_u = cvxpy.leq( cvxpy.abs( error_u ), .0001 )
        
        CONSTRAINTS.append( conserve_u )
        
    # the cost-form constraints
    for road in cost :
        for f, (kappa,alpha) in CTREES[road].iter_items() : 
            # is this plus or minus alpha?
            LB = cvxpy.geq( cost[road], kappa + alpha * assist[road] )
            CONSTRAINTS.append( LB )
    
    prog = cvxpy.program( OBJECTIVE, CONSTRAINTS )
    
    def showprog( prog ) :
        print 'minimize %s' % str( objfunc )
        print 'subject to:'
        for constr in CONSTRAINTS :
            print str( constr )
    showprog( prog )
    #print dict([ (road,help.value) for road, help in assist ])
    
    
    if True :      # solve the LP
        prog.solve()
        print 'cost "predicted": %f' % prog.objective.value
        
        print 'solution:'
        print [ a.value for a in assist.values() ]
        if False :
            violations = [ ( str(c), c.left.value, c.right.value, c.left.value >= c.right.value )
                          for c in CONSTRAINTS
                          if not c.left.value >= c.right.value ]
            for v in violations : print v
        
        for road in CTREES :
            ax = drawCBounds( ZZ, CTREES[road] )
            ax.set_title('road=%s' % road )
            zr = assist[road].value
            Cr = evalC( zr, CTREES[road] )
            #ax.axvline( assist[road].value )
            ax.scatter( [ zr ], [ Cr ], marker='o' )
            ax.scatter( [ zr ], [ cost[road].value ], marker='x' )
            
        dC = dict()
        for road in assist :
            zr = assist[road].value
            f, (kappa,alpha) = CTREES[road].floor_item( -zr )
            dC[road] = alpha
        
    if True :      # compare to optimal matching
        PR = [ ROAD.RoadAddress( road, y ) for road, y in PP ]
        QR = [ ROAD.RoadAddress( road, y ) for road, y in QQ ]
        
        class pointnode() :
            def __init__(self, point ) :
                self.point = point
        RED = [ pointnode(p) for p in PR ]
        BLU = [ pointnode(q) for q in QR ]
        graph = nx.Graph()
        match_mat = np.zeros((NUMPOINT,NUMPOINT))
        for (i,r), (j,b) in itertools.product( enumerate(RED), enumerate(BLU) ) :
            w = ROAD.distance( roadnet, r.point, b.point, 'length' )
            graph.add_edge( r, b, weight=-w )
            match_mat[i,j] = w
        match_dict = nx.max_weight_matching( graph, True )
        
        match = [ (r,match_dict[r]) for r in RED ]      # match pruning
        matchstats = [ ( r.point, b.point, ROAD.distance( roadnet, r.point, b.point, 'length' ) )
                      for r,b in match ]
        matchcost = sum([ c for r,b,c in matchstats ])
        print 'optimal match has cost: %f' % matchcost

        
        
    if False :       # validate CTREES
        zmin = -NUMPOINT/2
        zmax = NUMPOINT/2
        ZZZ = range( zmin, zmax )
        #
        for road, C in CTREES.iteritems() :
            ax = drawCBounds( ZZ, CTREES[road] )
            #Cz = np.array([ evalC( z, C ) for z in ZZ ])
            #plt.plot( ZZ, Cz, linestyle='--' )
            
            width = get_road_data( road, roadnet ).get( 'length', 1 )
            Cmatch = [ MATCHCOST( P[road], Q[road], width, z ) for z in ZZZ ]
            plt.scatter( ZZZ, Cmatch, marker='x' )
            
            
            
            
            
            
            