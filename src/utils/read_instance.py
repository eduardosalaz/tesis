import h5py
import sys
import pickle
import numpy as np

def read_instance(path):
    f = h5py.File(path)
    for key in f.keys():
        print(key)
    dataset = f['instance']
    # dataset[...] this in the repl shows names 
    B = dataset['B']
    S = dataset['S']
    K = dataset['K']
    M = dataset['M']
    P = dataset['P']
    ref_bu_coords = dataset['BU_coords']
    object_bu_coords = f[ref_bu_coords]
    bu_coords = [object_bu_coords[0], object_bu_coords[1]]
    ref_s_coords = dataset['S_coords']
    object_s_coords = f[ref_s_coords]
    s_coords = [object_s_coords[0], object_s_coords[1]]
    ref_D = dataset['D']
    obj_D = f[ref_D]
    cols_D, rows_D = obj_D.shape
    D = [obj_D[col] for col in range(cols_D)]
    D = np.matrix(D).T # fix order of rows and cols
    ref_Sk = dataset['Sk']
    ref_toref_Sk = f[ref_Sk]
    # here's the deal
    # when we point the reference of Sk we actually get another reference
    # since Sk is a vector of vectors, we need to get 2 refs to access each vector
    # in hindsight this was quite logic but it wasn't so easy at first
    Sk = []
    for k in range(K):
        Sk_row = f[ref_toref_Sk[k]]
        Sk.append(Sk_row[...])
    
    ref_Lk = dataset['Lk']
    obj_Lk = f[ref_Lk]
    Lk = obj_Lk[...]

    ref_Uk = dataset['Uk']
    obj_Uk = f[ref_Uk]
    Uk = obj_Uk[...]
    # again, the same goes for V as V is also a vector of vectors
    ref_V = dataset['V']
    ref_toref_V = f[ref_V]
    V = []
    for m in range(M):
        V_row = f[ref_toref_V[m]]
        V.append(V_row[...])
    # same for mu
    ref_mu = dataset['μ']
    ref_toref_mu = f[ref_mu]
    mu = []
    for m in range(M):
        mu_row = f[ref_toref_mu[m]]
        mu.append(mu_row[...])
    ref_T = dataset['T']
    obj_T = f[ref_T]
    T = obj_T[...]
    ref_R = dataset['R']
    obj_R = f[ref_R]
    R = obj_R[...]
    ref_beta = dataset['β']
    obj_beta = f[ref_beta]
    beta = obj_beta[...]
    print(beta)
    print(Sk)
    # then we should create a class and pickle it
    f.close()


if __name__ == '__main__':
    path = sys.argv[1]
    read_instance(path)
    