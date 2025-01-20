import numpy as np

def Accuracy(C,X):
    
    NumFixedEffects = X.shape[1]
    C_inv = np.linalg.pinv(C)
    C22_inv = C_inv[NumFixedEffects:, NumFixedEffects: ]
    r2 = np.diag(np.ones((C22_inv.shape[0], C22_inv.shape[1]), dtype = np.int8) - C22_inv*2)
    r2.setflags(write=True)
    r2[r2 < 0] = 0
    r = np.sqrt(r2)
    return r