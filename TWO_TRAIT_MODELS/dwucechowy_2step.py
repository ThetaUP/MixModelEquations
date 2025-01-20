import numpy as np

def MME2G(y1,y2,X1,X2,Z11,Z12,Z21,Z22,A,G,R,M,N):

    # creating the y vector
    y = np.concatenate((y1,y2), axis = 0)

    R = np.kron(R, np.eye(len(y1/2)))

    X1_n, X1_p = X1.shape
    X2_n, X2_p = X2.shape
    Z11_n, Z11_p = Z11.shape
    Z12_n, Z12_p = Z12.shape
    Z21_n, Z21_p = Z21.shape
    Z22_n, Z22_p = Z22.shape

    # creating X matrix
    under_X1 = np.zeros((X2_n, X1_p), dtype = np.float16)
    above_X2 = np.zeros((X1_n, X2_p), dtype = np.float16)
    top = np.concatenate((X1, above_X2), axis = 1)
    bottom = np.concatenate((under_X1, X2), axis = 1)
    X = np.concatenate((top, bottom), axis = 0)

    # creating Z1 matrix
    under_Z11 = np.zeros((Z12_n, Z11_p), dtype = np.float16)
    above_Z12 = np.zeros((Z11_n, Z12_p), dtype = np.float16)
    top = np.concatenate((Z11, above_Z12), axis = 1)
    bottom = np.concatenate((under_Z11, Z12), axis = 1)
    Z1 = np.concatenate((top, bottom), axis = 0)

    # creating Z2 matrix
    under_Z21 = np.zeros((Z22_n, Z21_p), dtype = np.float16)
    above_Z22 = np.zeros((Z21_n, Z22_p), dtype = np.float16)
    top = np.concatenate((Z21, above_Z22), axis = 1)
    bottom = np.concatenate((under_Z21, Z22), axis = 1)
    Z2 = np.concatenate((top, bottom), axis = 0)


    # matrix operations
    X_t = np.transpose(X)
    R_inv = np.linalg.pinv(R)
    Z1_t = np.transpose(Z1)
    Z2_t = np.transpose(Z2)
    A_inv = np.linalg.pinv(A)
    G_inv = np.linalg.pinv(G)

    
    # creating the C matrix
        # row1
    left = X_t@R_inv@X
    middle = X_t@R_inv@Z1
    right = X_t@R_inv@Z2
    row1 = np.concatenate((left, middle, right), axis = 1)

        # row2
    left = Z1_t@R_inv@X
    middle = Z1_t@R_inv@Z1 + np.kron(A_inv, M)
    right = Z1_t@R_inv@Z2
    row2 = np.concatenate((left, middle, right), axis = 1)

        # row3
    left = Z2_t@R_inv@X
    middle = Z2_t@R_inv@Z1
    right = Z2_t@R_inv@Z2 + np.kron(G_inv,N)  # this is for accuracy
    acc_matrix = right
    row3 = np.concatenate((left, middle, right), axis = 1)

    C = np.concatenate((row1, row2, row3), axis = 0)


    # creating the RHSs
    top = X_t@R_inv@y
    middle = Z1_t@R_inv@y
    down = Z2_t@R_inv@y
    RHS = np.concatenate((top, middle, down), axis = 0)

    # estimation
    sol = np.linalg.pinv(C)@RHS

    # accuracy
    C22 = np.linalg.pinv(acc_matrix)
    idx = X1.shape[1]
    C22_small = C22[idx:, idx:]
    trait1 = np.diag(C22_small)[:len(y1)+1]
    trait2 = np.diag(C22_small)[len(y1)+2:]
    var_trait1 = N[0,0]
    var_trait2 = N[1,1]
    r2_trait1 = []
    r2_trait2 = []
    arg_trait1 = (var_trait1 - trait1)/var_trait1
    arg_trait2 = (var_trait2 - trait2)/var_trait2

    r2_trait1 = np.sqrt(abs(var_trait1 - trait1)/var_trait1)
    r2_trait2 = np.sqrt(abs(var_trait2 - trait2)/var_trait2)

    return (sol, C, r2_trait1, r2_trait2)


# test_______________________________________________________________________________________________________-
'''
y1 = np.array([4.5, 2.9, 3.9, 3.5, 5.0])
y2 = np.array([6.8, 5.0, 6.8, 6.0, 7.5])
N = np.array([[20, 18],
              [18,40]])

M = np.array([[10,12],
              [15,7]])
R = np.array([[40,11],
             [11,30]])

sex = np.array([1,0,0,1,1])
X1 = np.zeros((5,2), dtype = np.int8)
X1[:,0] = sex
X1[:,1] = 1 - sex
X2 = X1       # x2 WILL JUST BE A COPY OF X1   !!!!!!!!!!



I = np.eye(5)

# Z1
Z11 = np.zeros((5,8), dtype = np.int8)
Z11[:5, 3:] = I
Z12 = Z11   # Z2 WILL JUST BE A COPY OF Z1 !!!!!!!!!!!!!!!

#Z2
Z21 = np.zeros((5,8), dtype = np.int8)
Z21[:5, 3:] = I
Z22 = Z21   # Z2 WILL JUST BE A COPY OF Z1 !!!!!!!!!!!!!!!


G = np.array([[1.2357915,	-0.3175200,	-0.3494598,	-0.2747769,	-0.2940348],
[-0.3175200,	1.2245186,	-0.2611555,	-0.3297323,	-0.3161109],
[-0.3494598,	-0.2611555,	1.2287459,	-0.3006106,	-0.3175200],
[-0.2747769,	-0.3297323,	-0.3006106,	1.2325035,	-0.3273837],
[-0.2940348,	-0.3161109,	-0.3175200,	-0.3273837,	1.2550493]
])


A = np.array([[1.00,	0.00,	0.00,	0.500,	0.000,	0.50,	0.25,	0.250],
                [0.00,	1.00,	0.00,	0.000,	0.500,	0.50,	0.25,	0.250],
                [0.00,	0.00,	1.00,	0.000,	0.500,	0.00,	0.25,	0.500],
                [0.50,	0.00,	0.00,	1.000,	0.000,	0.25,	0.50,	0.125],
                [0.00,	0.50,	0.50,	0.000,	1.000,	0.25,	0.50,	0.375],
                [0.50,	0.50, 0.00,	0.250,	0.250,	1.00,	0.25,	0.500],
                [0.25,	0.25,	0.25,	0.500,	0.500,	0.25,	1.00,	0.250],
                [0.25,	0.25,	0.50,	0.125,	0.375,	0.50,	0.25,	1.000]])

# we need to resize G so that it has the same dimensions as A
# set the values for all the cows which werent genotyped to 0 (similarly as in Z)
G_big = np.zeros((A.shape[0], A.shape[1]))
rowsisze = A.shape[0] - G.shape[0]
colsize = A.shape[1] - G.shape[1]
G_big[rowsisze:, colsize: ] = G

res = MME2G(y1 = y1,
            y2 = y2,
            X1 = X1,
            X2 = X2,
            Z11 = Z11,
            Z12 = Z12,
            Z21 = Z21,
            Z22 = Z22,
            A = A,
            G = G_big,
            R = R,
            M = M,
            N = N)


print(res[0])
'''