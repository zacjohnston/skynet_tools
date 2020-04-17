"""This module contains the class used to inspect a simulation output
directory. It searchs for the logfile, datfile and checkpoint files.
"""
#TODO: Could also look for plotfiles (_hdf5_plt_%04d), but I never used them.

import sys
import os
from glob import glob

class FlashOutputFiles(object):
    """name container for all simulation files in an output directory

    These files are checkpointfiles, the logfile and the datfile.

    Args:
      pathAndBaseName (str): path to the output directory with the run's basename.

    Attributes:
      path2base (str): path to the output directory with the run's basename
      basename (str): the run's basename
      dirname (str): name of the output directory (not the path!)
      logfile (str or None): path to logfile if present
      datfile (str or None): path to datfile if present
      chknumbers (list or None): list of checkpoint file numbers if any is present

    Methods:
      reload_files: Searches again for all files using privat _search methods.
      chkfiles: Returns a generator giving the paths to chkfiles in ascending order.

    """

    def __init__(self, pathAndBaseName):
        """class constructor taking the Args"""
        dirname = os.path.dirname(pathAndBaseName)
        if not os.path.exists(dirname):
            errorMessage = (
                         "Output directory '{directory}' "
                         "does not exist").format(directory=dirname)
            raise IOError(errorMessage)
        self.path2base = pathAndBaseName
        self.basename = os.path.basename(pathAndBaseName)
        self.dirname = os.path.basename(dirname)
        self.logfile = self._search_file('.log')
        self.datfile = self._search_file('.dat')
        self.chkfile = self.path2base + '_hdf5_chk_{i:04d}'
        self.chknumbers = self._search_chkfiles()


    def _search_file(self,extension):
        """Return path to file if it exists, else return None.
        Argument is the filename extension ('.log','.dat')
        """
        file_ = self.path2base + extension
        if not os.path.isfile(file_):
            file_ = None
        return file_


    def _search_chkfiles(self):
        """Return a list of checkpoint file indices found in the output
        directory.
        """
        chkwildcard = self.path2base + '_hdf5_chk_*'
        chkfiles = glob(chkwildcard)
        if chkfiles == []:
            chknumbers = None
        else:
            chknumbers = [int(chkfile[-4:]) for chkfile in chkfiles]
            chknumbers.sort()

        return chknumbers


    def reload_files(self):
        """Search again for all files."""
        self.logfile = self._search_file('.log')
        self.datfile = self._search_file('.dat')
        self.chknumbers = self._search_chkfiles()


    def chkFilePaths(self):
        """Return a generator for checkpoint file paths if any is present"""
        if self.chknumbers is None:
            errorMessage = (
                    "There are no checkpoint files in "
                    "directory '{directory}' to the given "
                    "basename '{base}'").format(directory=self.dirname,
                                                base=self.basename)
            raise IOError(errorMessage)

        path2chkfile = self.path2base + '_hdf5_chk_{n:04d}'
        chknumbers = self.chknumbers
        for number in chknumbers:
            path = path2chkfile.format(n=number)
            yield path

#TODO: Try unittest package!
def main():
    """function to test this module"""
    pathToOutputDir = 'test/run'

    myFiles = FlashOutputFiles(pathToOutputDir)

    print (myFiles.dirname)
    myChkfiles = myFiles.chkFilePaths()
    files = [chkfiles for chkfiles in myChkfiles]
    print ("first checkpoint:", files[0])
    print ("last checkpoint:", files[-1])
    print ("datfile:", myFiles.datfile)
    print ("\nTest passed!")


if __name__ == '__main__':
    main()

