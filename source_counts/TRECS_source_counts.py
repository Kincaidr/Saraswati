import pickle
import numpy as np

def TRECS_counts():
    file= 'TRECS_func.pkl'
    with open(file, 'rb') as file:
        data = pickle.load(file)

    new_grid = np.linspace(-7, 1, 10000) 
    counts=data(new_grid)
    return (10**new_grid)*1e3, counts
