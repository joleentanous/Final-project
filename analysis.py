import sys
import numpy as np
import pandas as pd



def error():
    print("An Error Has Occurred")
    exit()


def read(data):
    df = pd.read_csv(data, header=None)
    np_array = df.to_numpy()
    return np_array

def main():
    args = sys.argv
    if (len(args) != 3):
        error()
    k = (args[1])
    file_name = (args[2])
    try:
        k = int(k)
        file_name = str(file_name)
    except:
        error()
    matrix = read(file_name)

