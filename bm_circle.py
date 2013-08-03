
import numpy as np
import matplotlib.pyplot as plt
plt.close('all')

import networkx as nx
import bintrees

# sampling utility
WIDTH = 1.
def sample() :
    return WIDTH * np.random.rand()

# obtain a BM problem instance by sampling
NUMPOINT = 10
P = [ sample() for i in xrange(NUMPOINT) ]
Q = [ sample() for i in xrange(NUMPOINT) ]

def MATCH( P, Q, z=0 ) :
    PP = sorted([ (p,i) for i,p in enumerate(P) ])
    QQ = sorted([ (q,j) for j,q in enumerate(Q) ])
    
    I = [ i for p,i in PP ]
    J = [ j for q,j in QQ ]
    I = I[-z:] + I[:-z]
    return zip( I, J )

def MATCHCOST( P, Q, M, circumf=WIDTH ) :
    cost = 0.
    for i,j in M :
        c1 = abs( Q[j] - P[i] )
        c2 = circumf - c1
        cost += min( c1, c2 )
    return cost



""" start algorithm """
points = bintrees.RBTree()
# insert the endpoints
points.insert( 0., 0 )
points.insert( WIDTH, 0 )

# insert the problem data
for p in P : points.insert( p, 1 )
for q in Q : points.insert( q, -1 )

# prepare F and assign intervals
posts = [ y for y, sign in points.iter_items() ]
signs = [ sign for y, sign in points.iter_items() ]
#
pairs = zip( posts[:-1], posts[1:] )
F = np.cumsum( signs[:-1] )
#
yintervals = bintrees.RBTree()
for (a,b), f in zip( pairs, F ) :
    vals = yintervals.setdefault( f, [] )
    vals.append( (a,b) )

# prepare to plot F(y)
Ftree = bintrees.RBTree()
for lbd, f in zip( posts[:-1], F ) :
    Ftree.insert( lbd, f )
def evalF( y ) :
    _, f = Ftree.floor_item( y )
    return f


# simple plot of F
y = np.linspace(0,WIDTH,1000)
F = np.array([ evalF(yy) for yy in y ])

def drawpoints( points, ax=None, **kwargs ) :
    if ax is None :
        plt.figure()
        ax = plt.gca()
        
    radius = WIDTH / np.pi / 2
    to_rad = lambda y : 2 * np.pi * y / WIDTH
    theta = to_rad( np.array( points ) )
    X = radius * np.cos( theta )
    Y = radius * np.sin( theta )
    ax.scatter( X, Y, **kwargs )
    return ax

plt.figure()
plt.subplot(2,1,1)
ax = plt.gca()
drawpoints( P, ax, marker='x' )
drawpoints( Q, ax, marker='o' )
ax.set_aspect('equal')
ax.set_title( 'Two Sets of Points Around a Circle of Circumference %f' % WIDTH )
#redX = radius * np.array( np.cos( to_rad)) 
#Y = np.zeros(len(P))
#plt.scatter( P, Y, marker='x' )
#plt.scatter( Q, Y, marker='o' )
#plt.xlim(0.,WIDTH)
plt.subplot(2,1,2)
plt.plot(y,F)



# prepare the C(z) value lookup
# Step 1. get support of each value f
measure = dict()
for f, supp in yintervals.iter_items() :
    measure[f] = sum([ b-a for a,b in supp ])    # do the flip here

# prepare constants kappa and alpha
temp = sorted( measure.iteritems() )
ff = np.array( [ -np.inf ] + [ f for f,w in temp ] )
PREALPHA = np.array( [ 0. ] + [ w for f,w in temp ] )
PREKAPPA = np.array( [ 0. ] + [ f*w for f,w in temp ] )

ALPHAM = np.cumsum( PREALPHA )
ALPHAP = ALPHAM[-1] - ALPHAM
ALPHA = ALPHAP - ALPHAM

KAPPAM = np.cumsum( PREKAPPA )
KAPPAP = KAPPAM[-1] - KAPPAM
KAPPA = KAPPAP - KAPPAM

Czminus = bintrees.RBTree()
for f, kappa, alpha in zip( ff, KAPPA, ALPHA ) :
    Czminus.insert( f, ( kappa, alpha ) )
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
    Cz = [ MATCHCOST( P, Q, match ) for match in matches ]
    
    plt.figure()
    plt.plot(zz,C)
    plt.scatter(zzz,Cz, marker='x')
    

# augment Ctree with the offset as well




