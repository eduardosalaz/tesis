import h5py
import sys

if __name__ == '__main__':
    path = sys.argv[1]
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
    print(f"B: {B}")
    print(f"S: {S}")
    print(f"K: {K}")
    print(f"M: {M}")
    print(f"P: {P}")
    print("BU coordinates: ")
    print(bu_coords)
    print("S coordinates: ")
    print(s_coords)
    # and so on and so forth
    f.close()