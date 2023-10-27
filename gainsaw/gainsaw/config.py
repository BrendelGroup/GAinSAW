'''
config.py - from the gainsaw package

Scope
    a little wrapper to get some parameters and shortcuts into gainsaw

Usage
    from gainsaw import gsconf
    gsconf.read('test.gainsaw.conf')

Authors
    Wenshu Chen, email: wenschen@indiana.edu
    Volker Brendel, email: vbrendel@indiana.edu

Homepage 
    https://github.com/BrendelGroup/GAainSAW"

Licensed under MIT license.
'''

import sys
from os import path
import configparser
from collections import namedtuple


gsconf = configparser.ConfigParser()

# Read the package default config file:
#
dir_path = path.dirname(path.realpath(__file__))
gsconf.read(dir_path + '/' + 'gainsaw.conf')

# Note: A gsconf.read() call in your script will overwrite the default config
#       file and use the configuration parameters from your own config file.


scrs = namedtuple("scoring", "default blastn lastz",
                  defaults = [(-100, -2, "BLASTN"), (-7, -2, "BLASTN"), (-400, -30, "HOXD70")])
gsscrs = scrs()
params = namedtuple("parameters", "qlabel tlabel scoring slop_size mismatch_rate setpfilter usepfilter", defaults = ['unsp', 'unsp', gsscrs.default, 25, 10, 0, 1])
gsparams = params()


def cfcheck(configfile):
    if path.exists(configfile):
        print("... loading existing gainsaw config file ", configfile)
        return 1
    else:
        print("... config file ", configfile, " not found. Exiting.\n")
        sys.exit()
