'''

Licensed under MIT license.
'''
from os import path
import sys
import pybedtools

from .config import gsconf


class BedWrap:
    def __init__(self, interval_file, data_dir='.'):
        '''
        BedWrap arguments:
          interval_file		any valid interval file input accepted by
				 pybedtools.BedTool
          data_dir	(optional) the directory where to look for the input
                         file.
        '''
        if gsconf.has_option(interval_file,'ann') and path.exists(f := gsconf[interval_file]['ann']):
            print("... loading annotation file", f)
            self.bed = pybedtools.BedTool(f)
        elif path.exists(f := data_dir+'/'+interval_file):
            print("... loading annotation file", f)
            self.bed = pybedtools.BedTool(f)
        else:
            sys.stderr.write(f"Problem: specified annotation file {f} not found.")
            sys.exit()

    def check_bed(self):
        '''
        Just checking ...
        '''
        print(self.bed)

    def write_bed(self,filename):
        '''
        Write the BedTool object to the specified file.
        '''
        self.bed.saveas(filename)
