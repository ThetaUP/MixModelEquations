import numpy as np

def MME1(y,X,A,Z,sigma_a,sigma_e):

    alpha = sigma_e / sigma_a
    A_inv = np.linalg.pinv(A)
    X_t = np.transpose(X)
    Z_t = np.transpose(Z)

    TopLeft = X_t@X
    TopRight = X_t@Z
    Top = np.concatenate((TopLeft, TopRight), axis = 1)
    BottomLeft = Z_t@X
    BottomRight = Z_t@Z + A_inv*alpha
    Bottom = np.concatenate((BottomLeft, BottomRight), axis = 1)
    C = np.concatenate((Top, Bottom), axis = 0)
    C_inv = np.linalg.pinv(C)

    rhs = np.concatenate((X_t@y, Z_t@y), axis = 0)
    est = C_inv@rhs

    return (est, C)



# test_____________________________________________________________________________________________________________
"""

A = np.array([[1.00,	0.00,	0.00,	0.500,	0.000,	0.50,	0.25,	0.250],
[0.00,	1.00,	0.00,	0.000,	0.500,	0.50	,0.25,	0.250],
[0.00,	0.00,	1.00,	0.000,	0.500,	0.00,	0.25,	0.500],
[0.50,	0.00,	0.00,	1.000,	0.000,	0.25,	0.50,	0.125],
[0.00,	0.50,	0.50,	0.000,	1.000,	0.25,	0.50,	0.375],
[0.50,	0.50,	0.00,	0.250,	0.250,	1.00,	0.25,	0.500],
[0.25,	0.25,	0.25,	0.500,	0.500,	0.25,	1.00,	0.250],
[0.25,	0.25,	0.50,	0.125,	0.375,	0.50,	0.25,	1.000]])

y = np.array([4.5, 2.9, 3.9, 3.5, 5.0])
sex = np.array([1, 0, 0, 1, 1])
X = np.zeros((5,2), dtype = np.int8)
X[:,0] = sex
X[:,1] = 1-sex

I = np.eye(5)
Z = np.zeros((5,8), dtype = np.int8)
Z[:5, 3:] = I

hope = MME1(y=y, X=X, A=A, Z=Z, sigma_a=20, sigma_e=40)

print(hope[0])
"""