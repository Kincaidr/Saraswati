import pickle
import numpy as np
import matplotlib.pyplot as plt


def source_counts(file):
   
    with open(file, 'rb') as file:
        data = pickle.load(file)

    new_grid = np.linspace(-3.9, 5.8, 10000) 
    counts=data(new_grid)
    return new_grid, counts

if __name__ == "__main__":
     file= 'number_counts_total_mancuso17.pkl'
