
import numpy as np
import matplotlib.pyplot as plt
plt.close('all')

import networkx as nx
import bintrees

# to construct the optimization problem
import cvxpy


def MATCH( P, Q, z=0 ) :
    PP = sorted([ (p,i) for i,p in enumerate(P) ])
    QQ = sorted([ (q,j) for j,q in enumerate(Q) ])
    
    I = [ i for p,i in PP ]
    J = [ j for q,j in QQ ]
    I = I[-z:] + I[:-z]
    return zip( I, J )

def MATCHCOST( P, Q, M, circum ) :
    cost = 0.
    for i,j in M :
        c1 = abs( Q[j] - P[i] )
        c2 = circum - c1
        cost += min( c1, c2 )
    return cost




def obtain_segment_objective( P, Q, length ) :
    # Step 1. Construct F
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








if __name__ == '__main__' :
    # sampling utility
    WIDTH = 1.
    def sample() :
        return WIDTH * np.random.rand()
    
    # obtain a BM problem instance by sampling
    NUMPOINT = 10
    P = [ sample() for i in xrange(NUMPOINT) ]
    Q = [ sample() for i in xrange(NUMPOINT) ]
    
    Czminus = obtain_segment_objective( P, Q, WIDTH )
    def evalC( z ) :
        _, ( kappa, alpha ) = Czminus.floor_item( -z )
        return kappa + alpha * z
    
    if True :
        zmin = -NUMPOINT/2
        zmax = NUMPOINT/2
        zz = np.linspace( zmin, zmax, 1000 )
        C = np.array([ evalC(z) for z in zz ])
        
        zzz = range( zmin, zmax+1 )
        matches = [ MATCH( P, Q, z ) for z in zzz ]
        Cz = [ MATCHCOST( P, Q, match, WIDTH ) for match in matches ]
        
        plt.figure()
        plt.plot(zz,C)
        plt.scatter(zzz,Cz, marker='x')
        
    
    # augment Ctree with the offset as well
    
    
    
    
    