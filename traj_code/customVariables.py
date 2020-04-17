"""Module containing all getData methods for user customized variables

Note:
  The name of each method will be taken as name to call the data in 
  FLASH1dFile's getData method. The method should return a numpy array 
  of lenght FLASH1dFile.arrayLength and the argument should be an 
  instance of FLASH1dFile.
  
Example:
  Let's assume we want the product of entropy and temperature as output 
  and define it as heat.
  
  def heat(simFile):
      (documentation string usually with three double quotes) 
      '''heat defined as entropy times temperature'''
      
      tempData = simFile.getData('temp')
      entrData = simFile.getData('entr')
      heatData = tempData * entrData
      
      return heatData
      
  Now, you can call it with simFile.getData('heat').

"""
import sys
import numpy as np

def abar(simFile):
    """return mean atomic mass number defined as 1./sumY"""
    sumYData = simFile.getData('sumy')
    abarData = 1./sumYData
    
    return abarData 


def radi(simFile):
    """return radial position of the cell center"""
    nb = simFile.numBlocks
    cpb = simFile.cellsPerBlock
    leafBlocks = simFile.leafBlocks
    bbox = simFile['/bounding box'].value.T[:,0]

    # it's easier to itterate through a matrix and flatten it later
    radiMatrix = np.empty([nb,cpb])
    for i, block in enumerate(leafBlocks):
        blockLow = bbox[0,block]
        blockUp = bbox[1,block]
        drHalf = 0.5*(blockUp-blockLow)/cpb
        blockRadii = np.linspace(blockLow + drHalf, 
                                 blockUp - drHalf, cpb)
        radiMatrix[i] = blockRadii
    radiData = radiMatrix.flatten()
    
    return radiData


def mass(simFile):
    """return enclosed mass from cell center"""
    
    solarMass = 1.9891e33   # [g]
    vol = (4./3.)*np.pi
    
    leafBlocks = simFile.leafBlocks
    cpb = simFile.cellsPerBlock
    arrayLength = simFile.arrayLength
    bbox = simFile['/bounding box'].value.T[:,0].T
    
    lRad = np.empty(arrayLength)
    radi = np.empty(arrayLength)
    for i,block in zip(xrange(0,arrayLength,cpb),leafBlocks):
        # maybe the matix way looks more elegant
        blockLow = bbox[block,0]
        blockUp = bbox[block,1]
        dxB = (blockUp-blockLow)/cpb
        lRad[i:i+cpb] = np.arange(blockLow,blockUp,dxB)
        radi[i:i+cpb] = np.linspace(blockLow+dxB/2.,blockUp-dxB/2.,cpb)
    
    # This should be the right way to calculate the enclosed mass.
    # |  x1  |    x2    | => mass(x2) = mass(x1) + mass(x1  |    x2)
    #    |dm1|dm2 |       => mass(x2) = mass(x1) + dm1 + dm2
    densData = simFile.getData('dens')
    massData = np.empty(arrayLength)
    dm1 = vol*(lRad[1:]**3 - radi[:-1]**3)*densData[:-1]
    dm2 = vol*(radi**3 - lRad**3)*densData
    
    massData[0] = dm2[0]
    for i in xrange(1,arrayLength):
        massData[i] = massData[i-1] + dm1[i-1] + dm2[i]
    
    massData /= solarMass
    
    return massData


def rlvl(simFile):
    """return refinement level for each cell"""
    cpb = simFile.cellsPerBlock
    leafBlocks = simFile.leafBlocks
    rlevel = simFile['/refine level'].value
    
    # similar to matrix itteration in 'radi', but more compressed.
    rlvlData = np.array([[rlevel[block]]*cpb
                         for block in leafBlocks]).flatten()
    return rlvlData
    

def dRad(simFile):
    """return the cell size"""
    cpb = simFile.cellsPerBlock
    leafBlocks = simFile.leafBlocks
    blockSize = simFile['/block size'].value.T[0]
    
    # similar to matrix itteration in 'radi', but more compressed.
    dRadData = np.array([[blockSize[block]/cpb]*cpb
                         for block in leafBlocks]).flatten()
    return dRadData
    

