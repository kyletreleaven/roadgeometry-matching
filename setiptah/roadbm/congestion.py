
import numpy as np
import networkx as nx

import bintrees

import setiptah.roadbm as roadbm
import setiptah.roadgeometry.roadmap_basic as ROAD

import scipy.signal as sig


if __name__== '__main__' :
    
    import matplotlib.pyplot as plt
    plt.close('all')

    # parameters
    roadmap = nx.MultiDiGraph()
    roadmap.add_edge(0,0,'A', length=1. )       # circle

    """ congestion function """
    #rho = lambda x : np.power( x, 2. )          # square law, why not!
    rho = lambda x : abs(x)                      # linear congestion?

    N = 10


    # instance    
    sampler = roadbm.UniformDist(roadmap)
    
    S = [ sampler.sample() for i in xrange(N) ]
    T = [ sampler.sample() for i in xrange(N) ]
    
    # algorithm
    # stolen form ROADSBIPARTITEMATCH    
    
    MATCH = []
    
    segment_dict = roadbm.SEGMENTS( S, T, roadmap )
    measure_dict = {}
    
    for road, segment in segment_dict.iteritems() :
        match = roadbm.PREMATCH( segment )
        MATCH.extend( match )
        
        #surplus_dict[road] = SURPLUS( segment )
        
        roadlen = ROAD.get_road_data( road, roadmap ).get( 'length', 1 )
        measure = roadbm.MEASURE( segment, roadlen )
        measure_dict[road] = measure
        #objective_dict[road] = OBJECTIVE( measure )
        #objective_dict[road] = objective
    """ measure_dict now contains sequence of W_n """

    
    
    
    """ prepare convolution """
    a = measure.min_key()
    b = measure.max_key()
    
    # serialize
    if True :
        W = [ w for k, w in measure.items() ]
        W.reverse()
        
    else :
        W = np.zeros(b-a + 1 )
        for k, w in measure.items() :
            W[-k] = w   # reverse here, and account for the circular shift!
    
    """
        need C(z) for z \in [-N,N] ??? is this true?
        so, W_{n-z} = 0 for all n < a - N, and for all n > b + N
    """
    # relevant range
    A, B = a - N, b + N
    L = B-A + 1
    
    alpha = {}
    for n in xrange(A,B+1) :
        alpha[n] = abs(n) * rho(n)

    # prepare sequence form 
    if False :
        # w/ circular shift?
        alpha_seq = np.zeros(L)
        for n in xrange(A,B+1) : alpha_seq[n] = alpha[n]
        
        #for n in xrange(B+1) : alpha_seq[n] = alpha[n]
        #for n in xrange(A,0) :
            #alpha_seq[ B+1 + n-A ] = alpha[n]
            
    else :
        alpha_seq = [ alpha[k] for k in xrange(A,B+1) ]
        


    
            
    
        
        
    # long-hand, for verification
    C_manual = {}
    for z in xrange(-N,N+1) :
        C_manual[z] = 0.
        for k, w in measure.items() :
            n = k + z
            C_manual[z] += alpha[n] * w
            
    # efficient, FFT method
    C_fft = sig.fftconvolve( W, alpha_seq )
    
    plt.figure()
    plt.subplot(3,1,1)
    plt.scatter( range(len(W)), W )
    plt.subplot(3,1,2)
    plt.scatter( range(len(alpha_seq)), alpha_seq )
    plt.subplot(3,1,3)
    plt.scatter( range(len(C_fft)), C_fft )
    
    
    #C_right = C_fft[:B+1]
    #C_left = C_fft[A:]
    #C_fix = np.hstack( ( C_left, C_right ) )

    nspace = range(-N,N+1)
    C_fix = C_fft[b-A-N:b-A+N+1]
    
    compare = np.array([ C_manual[k] for k in xrange(-N,N+1) ])
    error = C_fix - compare
    
    plt.figure()
    plt.subplot(3,1,1)
    plt.scatter( C_manual.keys(), C_manual.values() )
    plt.subplot(3,1,2)
    C = C_fix
    plt.scatter( nspace, C_fix )
    plt.subplot(3,1,3)
    plt.scatter( nspace, error )
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    


    