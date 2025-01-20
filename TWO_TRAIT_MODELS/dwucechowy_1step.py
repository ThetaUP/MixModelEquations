import numpy as np

def MME2(y1,y2,X1,X2,Z1,Z2,A,G,R):

    # creating the y vector
    y = np.concatenate((y1,y2), axis = 0)

    R = np.kron(R, np.eye(len(y1)))

    X1_n, X1_p = X1.shape
    X2_n, X2_p = X2.shape
    Z1_n, Z1_p = Z1.shape
    Z2_n, Z2_p = Z2.shape

    # creating X matrix
    under_X1 = np.zeros((X2_n, X1_p), dtype = np.float16)
    above_X2 = np.zeros((X1_n, X2_p), dtype = np.float16)
    top = np.concatenate((X1, above_X2), axis = 1)
    bottom = np.concatenate((under_X1, X2), axis = 1)
    X = np.concatenate((top, bottom), axis = 0)

    # creating Z matrix
    under_Z1 = np.zeros((Z2_n, Z1_p), dtype = np.float16)
    above_Z2 = np.zeros((Z1_n, Z2_p), dtype = np.float16)
    top = np.concatenate((Z1, above_Z2), axis = 1)
    bottom = np.concatenate((under_Z1, Z2), axis = 1)
    Z = np.concatenate((top, bottom), axis = 0)


    # creating the C matrix
    top_left = np.transpose(X)@np.linalg.pinv(R)@X
    top_right = np.transpose(X)@np.linalg.pinv(R)@Z
    bottom_left = np.transpose(Z)@np.linalg.pinv(R)@X
    bottom_right = np.transpose(Z)@np.linalg.pinv(R)@Z + np.kron(np.linalg.pinv(A), np.linalg.pinv(G)) # this is C22
    top = np.concatenate((top_left, top_right), axis = 1)
    bottom = np.concatenate((bottom_left, bottom_right), axis = 1)
    C = np.concatenate((top, bottom), axis = 0)

    # creating the right-hand-side matrix
    top = np.transpose(X)@np.linalg.pinv(R)@y
    bottom = np.transpose(Z)@np.linalg.pinv(R)@y
    R = np.concatenate((top, bottom), axis = 0)

    
    # dokladnosc oszacowania
    C22 = np.linalg.pinv(bottom_right)
    idx = X1.shape[1]
    C22_small = C22[idx:, idx:]
    trait1 = np.diag(C22_small)[:len(y1)+1]
    trait2 = np.diag(C22_small)[len(y1)+2:]
    var_trait1 = G[0,0]
    var_trait2 = G[1,1]
    r2_trait1 = []
    r2_trait2 = []
    arg_trait1 = (var_trait1 - trait1)/var_trait1
    arg_trait2 = (var_trait2 - trait2)/var_trait2

    r2_trait1 = np.sqrt(abs(var_trait1 - trait1)/var_trait1)
    r2_trait2 = np.sqrt(abs(var_trait2 - trait2)/var_trait2)
    


    # solving the equation
    est = np.linalg.pinv(C)@R

    return (est, C, r2_trait1, r2_trait2)

'''
# testowanie funkcji____________________________________________________________________________________________________________-
y1 = np.array([4.5, 2.9, 3.9, 3.5, 5.0])
y2 = np.array([6.8, 5.0, 6.8, 6.0, 7.5])
G = np.array([[20, 18],
              [18,40]])
R = np.array([[40,11],
             [11,30]])

sex = np.array([1,0,0,1,1])
X1 = np.zeros((5,2), dtype = np.int8)
X1[:,0] = sex
X1[:,1] = 1 - sex
X2 = X1       # x2 WILL JUST BE A COPY OF X1   !!!!!!!!!!
I = np.eye(5)
Z1 = np.zeros((5,8), dtype = np.int8)
Z1[:5, 3:] = I
Z2 = Z1   # Z2 WILL JUST BE A COPY OF Z1 !!!!!!!!!!!!!!!
A = np.array([[1.00,	0.00,	0.00,	0.500,	0.000,	0.50,	0.25,	0.250],
                [0.00,	1.00,	0.00,	0.000,	0.500,	0.50,	0.25,	0.250],
                [0.00,	0.00,	1.00,	0.000,	0.500,	0.00,	0.25,	0.500],
                [0.50,	0.00,	0.00,	1.000,	0.000,	0.25,	0.50,	0.125],
                [0.00,	0.50,	0.50,	0.000,	1.000,	0.25,	0.50,	0.375],
                [0.50,	0.50, 0.00,	0.250,	0.250,	1.00,	0.25,	0.500],
                [0.25,	0.25,	0.25,	0.500,	0.500,	0.25,	1.00,	0.250],
                [0.25,	0.25,	0.50,	0.125,	0.375,	0.50,	0.25,	1.000]])




hope = MME2(y1=y1,
            y2=y2,
            X1=X1,
            X2=X2,
            Z1=Z1,
            Z2=Z2,
            A=A,
            G=G,
            R=R)

print(hope[1])
'''