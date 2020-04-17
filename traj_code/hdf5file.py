"""This module contains the file reader class for FLASH output files.
It is based on the h5py package.
"""

import numpy as np
import h5py
import inspect
import customVariables
import sys

class FLASH1dFile(h5py.File):
    """FLASH 1d file object 
    
    Args:
      pathToFile (str): path to the FLASH hdf5 output file
    
    Attributes:
      definedVars (list): list of defined mesh variables inside the 
                          file or user composed mesh variable in 
                          customVariables.
      simulationTime (float): simulation time where time zero is defined
                              with the start of the FLASH run
      realScalars (dict): dictronary for real scalars
      integerScalars (dict): dictronary for integer scalars
      cellsPerBlock (int): number of cells in each block
      leafBlocks (lsit): indices of blocks reconstructing the domain
                              with highest refinement level
      numBlocks (int): number of leaf blocks
      arrayLength (int): lenght of the output data arrays
    
    Method:
      getData: Returns a numpy array with mesh data for the given variable
    
    Note:    
      For some reason defined variables inside the FLASH file have 
      only four characters.
      
    Example: 
      Imagine we want the density data and we want to know what variables are
      defined for the getData method.
      
      with FLASH1dFile('path/to/file_hdf5_chk_1234') as myFLASH1dFile:
          dens = myFLASH1dFile.getData('dens') # returns the density profile
          print (myFLASH1dFile.definedVars # prints all defined variables)
      
    """
    
    def __init__(self, pathToFile):
        """class constructor taking Args to produce a file instance"""
        super(FLASH1dFile, self).__init__(pathToFile,'r')
        
        self.realScalars = self.loadField('/real scalars')
        self.integerScalars = self.loadField('/integer scalars')
        
        self.simulationTime = self.realScalars['time']
        self.cellsPerBlock = self.integerScalars['nxb']
        self.leafBlocks = self.loadBlockIndices()
        self.numBlocks = len(self.leafBlocks)
        self.arrayLength = self.cellsPerBlock*self.numBlocks
        
        self.definedVars = self.loadDefindedVariables()
    
    
    def __repr__(self):
        """Retrun the string representation of an instance."""
        string = ("FLASH 1d simulation file\n"
                  "file path: {path}\n"
                  "simulation time: {time:f}").format(path=self.filename,
                                                       time=self.simulationTime)
        return string
    
    
    def loadField(self, name):
        """Return a dictionary for a given field name and trims the key strings
        
        Example: 
          realScalars = loadField['/real scalars']
        """
        field = {k.strip():v 
                 for k, v in dict(self[name].value).iteritems()}
        return field
        
    def loadBlockIndices(self):
        """Retrun the list of indices to find the blocks holding the data at
        highest refinement."""
        nodeType = self['/node type'].value
        leafBlocks, = np.where(nodeType==1)
        
        return leafBlocks
        
    def loadDefindedVariables(self):
        """Create the 'get data' method dictionary and return the list 
        of defined variables.
        The field '/unknown names' contains all mesh variables defined inside
        the file. User composed variables are added from module customVariables"""
        self.simVars = [str(name) for name in self['/unknown names'].value.T[0]]
        self.usrMethods = dict(inspect.getmembers(customVariables,
                                        predicate=inspect.isroutine))
        self.usrVars = self.usrMethods.keys()
        definedVars = self.simVars + self.usrVars
        definedVars.sort()
        
        return definedVars
        
    def getMeshData(self,name):
        """Return the array of data collected from the leaf blocks for a
        variable name listed in simVars."""
        dataMatrix = np.empty([self.numBlocks,self.cellsPerBlock])
        fileData = self['/'+name].value[:,0,0,:]
        for i, block in enumerate(self.leafBlocks):
            dataMatrix[i] = fileData[block]
        dataArray = dataMatrix.ravel()
                    
        return dataArray
        
    def getData(self,varName):
        """Return the data array for a given name if name is in definedVars."""
        if varName in self.simVars:
            dataArray = self.getMeshData(varName)
        elif varName in self.usrVars:
            dataArray = self.usrMethods[varName](self)
        else:
            errorMessage = ("Variable '{name}' is neither defined inside"
                    "file '{filename}' nor in module 'customVariables'.\n"
                    "defined variables are: {variables}.").format(name=varName,
                    filename=self.filename, variables=self.definedVars)
            raise KeyError(errorMessage)
            
        return dataArray


def main():
    """function to test this module!"""
    pathToFile = 'test/run_hdf5_chk_0500'
    notPassed = []
    with FLASH1dFile(pathToFile) as file0500:
        defVars = file0500.definedVars
        print ("defined variables:")
        print (defVars)
        for var in defVars:
            try:
                varData = file0500.getData(var)
                isFloat = isinstance(varData[0],np.float64)
                isInt = isinstance(varData[0],np.int32)
                if not (isFloat or isInt):
                    raise
            except:
                notPassed.append(var)
    if notPassed == []:
        print ('\nOk, it works!')
    else:
        print ("something is wrong with {variables}".format(variables=notPassed))


if __name__ == '__main__':
    main()  


