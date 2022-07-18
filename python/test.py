# Main solver functions_classes

def setup():
    import numpy as np
    p  = {'n': 1.0, 'beta': np.arange(0,10,1)}
    return p


def run(p):
    n = p['n']
    beta = p['beta']
    model_results = n*beta
    return model_results
