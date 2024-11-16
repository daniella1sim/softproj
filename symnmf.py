import sys
import numpy as np
import symnmfmodule as sy
from tools import isInt, error

np.random.seed(1234)


"""
Parses the data from the command line
@type args: List
@param args: A list of command line arguments
@rtype: Tuple(int, string, string)
@returns: A tuple of the parsed data containing K, goal, and filepath to data
"""
def parse_data(args):
    if len(args) != 4:
        error()
    K = isInt(args[1])
    goal = args[2]
    if goal not in ["symnmf", "sym", "ddg", "norm"]:
        error()
    filepath = sys.argv[3]
    if (filepath[-4:] != ".txt"):
        error()
    return K, goal, filepath


"""
Main function that runs the program
@rtype: int
@returns: 0 if program runs successfully else 1
"""
def main():
    K, goal, filepath = parse_data(sys.argv)
    data = np.genfromtxt(filepath, delimiter=',', dtype=float)
    if data.ndim == 1:
        N, D = len(data), 1
        data_list = [[i] for i in data.tolist()]
    else:
        N, D = data.shape
        data_list = data.tolist()
    if (K > N):
        error()
    try:
        if goal == "symnmf":
            W_list = sy.norm(data_list)
            W = np.array(W_list)
            initial_H = np.random.uniform(low=0, high=(2*np.sqrt(np.average(W)/K)), size=(N, K))
            initial_H_list = initial_H.tolist()
            res = sy.symnmf(initial_H_list, W_list)
        elif goal == "sym":
            res = sy.sym(data_list)
        elif goal == "ddg":
            res = sy.ddg(data_list)
        else:
            res = sy.norm(data_list)
    except Exception as e:
        error()
    
    for row in res:
        target = ",".join([f"{x:.4f}" for x in row])
        print(target)
    return 0


if __name__ == "__main__":
    main()
