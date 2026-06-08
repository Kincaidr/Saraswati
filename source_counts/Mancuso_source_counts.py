import pickle
import numpy as np
import matplotlib.pyplot as plt

def Mancuso_counts():
    file= 'number_counts_total_mancuso17.pkl'
    with open(file, 'rb') as file:
        data = pickle.load(file)

    new_grid = np.linspace(-3.9, 5.8, 10000) 
    counts=data(new_grid)
    return 10**new_grid, 10**counts

if __name__ == "__main__":
     Mancuso_counts()