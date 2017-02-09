from algorithm.parameters import params

class pymax:
    """Py-max is a max-style problem where the goal is to generate a
    function which outputs a large number. In the standard GP Max
    problem this function can only use the constant (0.5) and
    functions (+, *). The Py-max problem allows more programming:
    numerical expressions, assignment statements and loops. See
    pymax.pybnf.

    """


    maximise = True # True as it ever was

    def __call__(self, ind):
        # ind.phenotype will be a string incl fn defns etc. when we
        # exec it will create a value XXX_output_XXX, but we exec
        # inside an empty dict for safety.
        p = ind.phenotype
        print(p)
        print("")
        d = {}
        exec(p, d)
        s = d['XXX_output_XXX'] # this is the program's output: a number
        return s
