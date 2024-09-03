import sys
import numpy as np
import symnmf

np.random.seed(1234)

def isInt(arg):
    try:
        return 0, int(arg)
    except:
        try:
            float_arg = float(arg)
            int_arg = int(arg[0: arg.find(".")])
            if (float_arg) == int_arg:
                return 0, int_arg
            else:
                return 1, 0
        except:
            return 1, 0


def similarity_matrix(data):
    N = data.shape[0]
    sim_matrix = np.zeros((N, N))
    for i in range(N):
        for j in range(N):
            if i == j:
                sim_matrix[i][j] = 0
                continue
            distance = np.linalg.norm(data[i] - data[j])
            sim_matrix = np.exp(-(distance**2)/2)
    return sim_matrix


def main():
    if len(sys.argv) != 4:
        print("Error: Incorrect number of arguments")
        return
    err_K, K = isInt(sys.argv[1])
    if err_K:
        print("Error: K must be an integer")
        return
    
    goal = sys.argv[2]
    if goal not in ["symnmf", "sym", "ddg", "norm"]:
        print("Error: goal must be one of symnmf, sym, ddg, norm")
        return
    
    filepath = sys.argv[3]
    if (filepath[-4:] != ".txt"):
        print("Error: file must be a .txt file")
        return
    
    data = np.genfromtxt(filepath, delimiter=',', dtype=float)
    N, D = data.shape
    if (K > N):
        print("Error: K must be less than the number of rows in the data ")
        return

    A = similarity_matrix(data)
    #D = np.diag(np.sum(A, axis=0))
    D_inv_sqrt = np.diag(1/np.sqrt(np.sum(A, axis=0)))
    W = D_inv_sqrt @ A @ D_inv_sqrt
    
    initial_H = np.random.uniform(low=0, high=(2*np.sqrt(np.average(W)/K)), size=(N, K))

    if goal == "symnmf":
        res = symnmf.symnmf(initial_H, W)
    elif goal == "sym":
        res = symnmf.sym(data)
    elif goal == "ddg":
        res = symnmf.ddg(data)
    else:
        res = symnmf.norm(data)

    for row in res:
        target = ",".join([f"{x:.4f}" for x in row])
        print(target)
    return 0


if __name__ == "__main__":
    main()
