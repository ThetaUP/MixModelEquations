# tylko dla jednocechowego bez genetyki

import jednocechowy as single_step
import numpy as np



def VarEst1(y,X,Z,A,sigma_a,sigma_e):
    n,p = X.shape
    q = A.shape[0]
    n_iter = 1
    tmp = 0.1
    threshold = 0.00001

    while(tmp > threshold):
       
        MME_new = single_step.MME1(y,X,A,Z,sigma_a,sigma_e)
        C_new = np.linalg.pinv(MME_new[1])
        est_new = MME_new[0]
        Ck = C_new[p:(p+q), p:(p+q)]

        # update sigma_a
        a = est_new[p:(p+q)]
        A_inv = np.linalg.pinv(A)
        sigma_a_new = (np.transpose(a)@A_inv@a + np.trace(A_inv@Ck)*sigma_e)/q

        # update sigma_e
        #       get residuals
        eps = (y-X@est_new[:p]) - (Z@est_new[p:(p+q)])
        trace_left = np.concatenate((X,Z), axis = 1)@C_new
        trace_right = np.transpose(np.concatenate((X,Z), axis = 1))
        sigma_e_new = (np.transpose(eps)@eps+np.trace(trace_left@trace_right)*sigma_e)/n

        # update tmp
        sigma_a_diff = np.abs(sigma_a - sigma_a_new)
        sigma_e_diff = np.abs(sigma_e - sigma_e_new)
        tmp = max(sigma_a_diff, sigma_e_diff)

        # update variance components
        sigma_a = sigma_a_new
        sigma_e = sigma_e_new

        n_iter = n_iter + 1


    return (sigma_a, sigma_e, n_iter)


# TEST_____________________________________________________________________________________________________________________________
'''
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


hope = VarEst1(y=y,
               X=X,
               Z=Z,
               A=A,
               sigma_a=500000,
               sigma_e=100000)


print(hope)
'''