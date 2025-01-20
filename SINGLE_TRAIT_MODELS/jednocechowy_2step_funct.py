import numpy as np

def MME1G(y,X,Z,H_inv,sigma_a,sigma_e):

    alpha = sigma_a/sigma_e
    X_t = np.transpose(X)
    Z_t = np.transpose(Z)
    top = np.concatenate((X_t@X, X_t@Z), axis = 1)
    bottom = np.concatenate((Z_t@X, Z_t@Z + H_inv*alpha), axis = 1)
    C = np.concatenate((top, bottom), axis = 0)
    C_inv = np.linalg.pinv(C)  # pseudo inverse
    rhs = np.concatenate((X_t@y, Z_t@y), axis = 0)
    est = C_inv@rhs
    return (est,C)


'''
# test the function_______________________________________________________________
y = np.loadtxt("TRANSFORMACJE_DANYCH/Y1_FirstLastInsem.txt", dtype = int)
X = np.loadtxt("TRANSFORMACJE_DANYCH/FixedEffects_X.csv", dtype = int, delimiter=',', skiprows=1)
Z = np.loadtxt("TRANSFORMACJE_DANYCH/Z_matrix.csv", dtype = int, delimiter = ',', skiprows=1)
H_inv = np.loadtxt("TRANSFORMACJE_DANYCH/H_inv_matrix.csv", dtype = float, delimiter=',', skiprows=1)

hope = MME1G(y = y, X = X, Z = Z, H_inv = H_inv, sigma_a = 0.5, sigma_e = 0.5)

print(hope[0])
'''