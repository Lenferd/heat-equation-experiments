import numpy as np

if __name__ == '__main__':
    nx = 17
    data = np.fromfile("function2-129m17.txt", float, sep=' ')
    data.shape = (nx, nx, nx)

    result_data = data[::2, ::2, ::2]

    np.ndarray.tofile(result_data, "function2-129m9.txt", sep=' ', format="%s\n")
