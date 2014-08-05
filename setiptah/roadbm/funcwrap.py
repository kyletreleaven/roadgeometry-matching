

import bintrees


class RBTreeLinterp :
    def __init__(self, rbtree ) :
        assert isinstance(rbtree, bintrees.RBTree )
        self.tree = tree
        
    def __call__(self, z ) :
        z1, f1 = self.tree.floor_item(z)
        z2, f2 = self.tree.ceiling_item(z)
        
        if z2 > z1 :
            m = float( f2 - f1 ) / ( z2 - z1 )
            return f1 + m * ( z - z1 )
        else :
            return f1




if __name__ == '__main__' :
    import matplotlib.pyplot as plt
    plt.close('all')
    

    import numpy as np
    
    tree = bintrees.RBTree()
    
    kspace = np.random.rand(5)
    
    tree[0] = np.random.rand()
    tree[1] = np.random.rand()
    for k in kspace :
        tree[k] = np.random.rand()
        
    func = RBTreeLinterp( tree )
    ff = np.vectorize( func )
    
    tspace = np.linspace(0,1,100)
    f = ff( tspace )
    
    plt.plot( tspace, f )
    
    
    
    
    
    
    
    
        
      
    
