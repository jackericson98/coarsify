from System.sys_funcs.calcs import calc_dist, calc_com
from System.sys_objs.chain import Ball
from System.cg_designations import proteins, nucleobases, ions, solvents

def coarsify_primo(sys, therm_cush=0.0):
    """
    primo coarse graining model
    :param sys:
    :param therm_cush:
    :return:
    """
    sys.balls = []
    return