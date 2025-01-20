import numpy as np
import sys
from matplotlib import pyplot as plt
sys.path.append('MODELE_JEDNOCECHOWE/')
import jednocechowy as MixModel
import estymacja_wariancji_jednocechowy as EM
import dokladnosc_funct as Acc
sys.path.append("MODELE_DWUCECHOWE/")
import dwucechowy_1step as MixModel2


# read data__________________________________________________________________________________________
X = np.loadtxt("TRANSFORMACJE_DANYCH/FixedEffects_X.csv", delimiter=',', skiprows=1)
y1 = np.loadtxt("TRANSFORMACJE_DANYCH/Y1_FirstLastInsem.txt")
y2 = np.loadtxt("TRANSFORMACJE_DANYCH/Y2_LastInsemCalving.txt")
Z = np.loadtxt("TRANSFORMACJE_DANYCH/Z_matrix.csv", delimiter=',', skiprows=1)
A = np.loadtxt("TRANSFORMACJE_DANYCH/A_matrix.csv", delimiter=',', skiprows=1, usecols=range(1,Z.shape[1]+1))



# model jednocechowy bez genetyki____________________________________________________________________________________________________


# estymacja parametrów wariancji
'''
print('estimating variance components...')
var_comp_y1 = EM.VarEst1(y=y1, X=X, Z=Z, A=A, sigma_a = 235.55, sigma_e = 3632.90)
var_comp_y2 = EM.VarEst1(y=y2, X=X, Z=Z, A=A, sigma_a = 412.24, sigma_e = 8134.30)

sigma_a_y1 = var_comp_y1[0]
sigma_e_y1 = var_comp_y1[1]
sigma_a_y2 = var_comp_y2[0]
sigma_e_y2 = var_comp_y2[1]
print('done.')
'''
print('solving equations...')
x = 5
mme1_y1 = MixModel.MME1(y=y1, X=X, A=A, Z=Z, sigma_a = x, sigma_e = x)
mme1_y2 = MixModel.MME1(y=y2, X=X, A=A, Z=Z, sigma_a = x, sigma_e = x)
print('done.')

# wartości hodowlane
sol_y1 = mme1_y1[0]
sol_y2 = mme1_y2[0]
bv_y1 = sol_y1[X.shape[1]:] # remove fixed effect estimates
bv_y2 = sol_y2[X.shape[1]:] # remove fixed effect estimates


# dokładność oceny
C_y1 = mme1_y1[1]
C_y2 = mme1_y2[1]
acc_y1 = Acc.Accuracy(C = C_y1, X = X)
acc_y2 = Acc.Accuracy(C = C_y2, X = X)



# model dwucechowy bez genetyki______________________________________________________________________________________________________

G = np.array([[235.55, -238.38],
              [-283.38, 412.24]])

R = np.array([[3632.90, -1107.40],
              [-1107.40, 8134.30]])

X1 = X
X2 = X
Z1 = Z
Z2 = Z

print('solving equations for 2 trait model...')
mme2 = MixModel2.MME2(y1 = y1, y2 = y2, X1 = X1, X2 = X2, Z1 = Z1, Z2 = Z2, A = A, G = G, R = R)
print('done.')

# wartości hodowlane
sol_dwucechowy = mme2[0]
bv_dwucechowy = sol_dwucechowy[(X1.shape[1]+X2.shape[1]):] # remove fixed effect estimates for both traits
cutoff = int(np.floor(bv_dwucechowy.shape[0]/2))
bv_dwucechowy_y1 = bv_dwucechowy[:cutoff]
bv_dwucechowy_y2 = bv_dwucechowy[cutoff:]

# dokładność oceny
acc_y1_dwucechowy = mme2[2]
acc_y2_dwucechowy = mme2[3]


# porownywanie korelacji wartości hodowlanych________________________________________________________________________
print('plotting...')

#y1
plt.figure()
plt.scatter(bv_y1, bv_dwucechowy_y1)
plt.grid()
plt.xlabel('wartości hodowlane cecha y1 model jednocechowy', fontsize = 8)
plt.ylabel('wartości hodowlane cecha y1 model dwucechowy', fontsize = 8)
plt.title('Korelacja wartości hodowlanych dla cechy y1.')
plt.savefig('no_genotypes_correlation_y1.png', dpi=300, bbox_inches='tight')
plt.close()

#y2
plt.figure()
plt.scatter(bv_y2, bv_dwucechowy_y2)
plt.grid()
plt.xlabel('wartości hodowlane cecha y2 model jednocechowy', fontsize = 8)
plt.ylabel('wartości hodowlane cecha y2 model dwucechowy', fontsize = 8)
plt.title('Korelacja wartości hodowlanych dla cechy y2.')
plt.savefig('no_genotypes_correlation_y2.png', dpi=300, bbox_inches='tight')
plt.close()


# porównywanie dokładności oceny________________________________________________________________________________________________

#y1
plt.figure()
plt.subplot(1,2,1)
plt.boxplot(acc_y1)
plt.grid()
plt.ylim(0,1)
plt.title('Dokładność oceny model jednocechowy, cecha y1', fontsize = 8)
plt.subplot(1,2,2)
plt.boxplot(acc_y1_dwucechowy)
plt.grid()
plt.ylim(0,1)
plt.title('Dokładnośc oceny model dwucechowy, cecha y1', fontsize = 8)
plt.savefig('no_genotypes_accuracy_y1.png', dpi=300, bbox_inches='tight')
plt.close()

#y2
plt.figure()
plt.subplot(1,2,1)
plt.boxplot(acc_y2)
plt.grid()
plt.ylim(0,1)
plt.title('Dokładność oceny model jednocechowy, cecha y2', fontsize = 6)
plt.subplot(1,2,2)
plt.boxplot(acc_y2_dwucechowy)
plt.grid()
plt.ylim(0,1)
plt.title('Dokładnośc oceny model dwucechowy, cecha y2', fontsize = 6)
plt.savefig('no_genotypes_accuracy_y2.png', dpi=300, bbox_inches='tight')
plt.close()

print('all done.')