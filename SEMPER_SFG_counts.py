import pickle
import numpy as np

def SEMPER_SFG_counts():
    file='tools/number_counts_SFGs_SEMPER.pkl'
    with open(file, "rb") as file2:
        data = pickle.load(file2)

    new_grid = np.linspace(-3.9, 5.8, 10000) 
    counts=data(new_grid)
    return 10**new_grid, 10**counts