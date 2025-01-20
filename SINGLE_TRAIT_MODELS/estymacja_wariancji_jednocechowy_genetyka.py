import numpy as np
import sys
sys.path.append('MODELE_JEDNOCECHOWE/')
import jednocechowy_2step_funct as MixModel

def EM_G(y, X, Z, H_inv, sigma_a, sigma_e):
    n,p = X.shape
    q = H_inv.shape[0]
    n_iter = 1
    tmp = 0.1
    thresh = 0.00001

    while (tmp > thresh):
        mme_new = MixModel.MME1G(y=y, Z=Z, H_inv=H_inv, sigma_a=sigma_a, sigma_e=sigma_e)
        C_new = np.linalg.pinv(mme_new[1])
        Ck = C_new[p:(p+q), p:(p+q)]
        mme2 = mme_new[0]

        a = mme2[p:(p+q)]
        sigma_a_new = (np.transpose(a)@a + np.trace(H_inv@Ck)*sigma_e)/q
        
        if sigma_a_new < 0:
            sigma_a_new = 0.01
        
        res = (y-X@mme2[:p]) - (Z@mme2[p:(p+q)])
        X.tmp1 = np.concatenate((X,Z), axis = 1)@C_new
        X.tmp2 = np.transpose(np.concatenate((X,Z), axis = 1))
        sigma_e_new = (np.transpose(res)@res + np.trace(X.tmp1@X.tmp2)*sigma_e)/n
        
        if sigma_e_new < 0:
            sigma_e_new = 0.01
        

        diff1 = np.abs(sigma_a - sigma_a_new)
        diff2 = np.abs(sigma_e - sigma_e_new)
        tmp = max(diff1, diff2)
        sigma_a = sigma_a_new
        sigma_e = sigma_e_new
        n_iter = n_iter + 1

    return (sigma_a, sigma_e, n_iter)