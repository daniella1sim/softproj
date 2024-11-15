
"""
Prints an error message and exits the program.

@rtype: None
@returns: Exits the program with a status of 1.
"""
def error():
    print("An Error Has Occured!")
    exit(1)


""""
Verifies argument is an int
@type arg: String
@param arg: A string of input from user to verify
@rtype: Touple(int, int)
@returns: arg in integer form.
"""
def isInt(arg):
    try:
        return int(arg)
    except:
        try:
            float_arg = float(arg)
            int_arg = int(arg[0: arg.find(".")])
            if (float_arg) == int_arg:
                return int_arg
            else:
                error()
        except:
            error()
