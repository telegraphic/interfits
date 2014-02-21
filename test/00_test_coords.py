#! /usr/bin/env python
# encoding: utf-8
from test_main import *
import time                                                

def timeme(method):
    def wrapper(*args, **kw):
        startTime = int(round(time.time() * 1000))
        
        for ii in range(1):
            result = method(*args, **kw)
        endTime = int(round(time.time() * 1000))

        print(endTime - startTime,'ms')
        return result

    return wrapper

@timeme
def test_coords():
    H, d = 1.23, np.deg2rad(4.56)
    xyz = np.array([
        [1, 2, 3],
        [2, 3, 3],
        [2, 1, 4],
        [8, 1, 2],
        [7,1, 8.1]
        ], dtype='float32')
        
    c = coords.computeUVW(xyz, H, d)
    return c

def test_computeBaselineVectors():
    xyz = np.array([
        [1, 2, 3],
        [2, 3, 3],
        [2, 1, 4],
        [8, 1, 2],
        [7,1, 8.1]
        ], dtype='float32')
    
    bls = coords.computeBaselineVectors(xyz, autocorrs=True)
    assert bls.shape == (5 * (5 - 1) / 2 + 5, 3)
    
    bls = coords.computeBaselineVectors(xyz, autocorrs=False)
    assert bls.shape == (5 * (5 - 1) / 2, 3)    
         

if __name__ == "__main__":
    #test_coords()
    
    test_computeBaselineVectors()
    

    
    

