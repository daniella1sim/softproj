import sys
import numpy as np
import symnmf as sy

np.random.seed(1234)


""""
Verifies argument is an int
@type arg: String
@param arg: A string of input from user to verify
@rtype: Touple(int, int)
@returns: First int return 1 if arg is not an integer and 1 if it is, second returns the arg in integer form.
"""
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


"""
Parses the data from the command line
@type args: List
@param args: A list of command line arguments
@rtype: Tuple(int, string, string)
@returns: A tuple of the parsed data containing K, goal, and filepath to data
"""
def parse_data(args):
    if len(args) != 4:
        print("Error: Incorrect number of arguments")
        return -1, 0, 0
    err_K, K = isInt(args[1])
    if err_K:
        print("Error: K must be an integer")
        return -1, 0, 0
    goal = args[2]
    if goal not in ["symnmf", "sym", "ddg", "norm"]:
        print("Error: goal must be one of symnmf, sym, ddg, norm")
        return -1, 0, 0
    filepath = sys.argv[3]
    if (filepath[-4:] != ".txt"):
        print("Error: file must be a .txt file")
        return -1, 0, 0
    return K, goal, filepath


"""
Main function that runs the program
@rtype: int
@returns: 0 if program runs successfully else 1
"""
def main():
    K, goal, filepath = parse_data(sys.argv)
    if K == -1:
        return 1
    data = np.genfromtxt(filepath, delimiter=',', dtype=float)
    N, D = data.shape
    if (K > N):
        print("An Error Has Occured!")
        return 1

    data_list = data.tolist()

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
    except:
        print("An Error Has Occured!")
        return 1
    
    for row in res:
        target = ",".join([f"{x:.4f}" for x in row])
        print(target)
    return 0


if __name__ == "__main__":
    main()
