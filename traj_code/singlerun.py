"""This module contains the class collecting and holding all metadata 
for a FLASH CCSN 1d run.
"""

import sys
sys.path.append('../')
import os
import numpy as np
import mathRecipes
from simfiles import FlashOutputFiles
from hdf5file import FLASH1dFile

class CCSN1dRun(FlashOutputFiles):
    """database for all simulation run information
    
    Args:
      pathAndBaseName (str): path to the output directory with the run's basename
    
    Attributes:
      bouncetime (float): time of bounce taken from the last file
                          if not occured it is equal to 0 and None if not found.
    
      inheritance:
      path2base (str): path to the output directory with the run's basename
      basename (str): the run's basename
      dirname (str): name of the output directory (not the path!)
      logfile (str or None): path to logfile if present
      datfile (str or None): path to datfile if present 
      chknumbers (list or None): list of checkpoint file numbers
                                 if any is present
      
                                      
    Methods:
      getBounceTime: Returns the time of core bounce.
      chkfileloop: Loops through all checkpoint files and collects data.
      getChkAtTime: Returns the file path to the checkpoint closest to a given time.
      
      inheritance:
      reload_files: Searches again for all files.
      chkfiles: Returns a generator giving the paths to chkfiles in
                ascending order.
                     
    """
    
    def __init__(self, pathAndBaseName):
        """constructor"""
        super(CCSN1dRun, self).__init__(pathAndBaseName)
        self.bouncetime = self.getBounceTime()
        self._time2Num = None # privat variable used in getChkAtTime()
        
        
        
    def getBounceTime(self):
        """Return the time of bounce listed in the last checkpoint file"""
        try:
            lastFileIndex = self.chknumbers[-1]
            path2lastFile = self.path2base + '_hdf5_chk_{n:04d}'.format(n=lastFileIndex)
            with FLASH1dFile(path2lastFile) as lastFile:
                bouncetime = lastFile.realScalars['bouncetime']
            self.bounce = True
        except:
            # If there is no chk file, it is still posible to extract the
            # bouncetime from the logfile.
            # TODO: write different getBounceTime routines and pass them
            #       to an evaluation function doing all the try, except business.
            bouncetime = 0
            self.bounce = False
        
        return bouncetime
        
        
    def chkfileloop(self,func,processing=None):
        """Loop over all chkfiles and apply function 'func' to each file while 
        collecting the output.
        Optionally, apply function 'processing' to the collected output
        
        Args:
          func (function): function to apply on each checkpoint file.
          processing (function,optional): function for further processing
                                          of the collected output
        
        Returns:
          if processing is None
            list: the collected return values for each func call.
          otherwise
            the return of function processing.
        
        Note:
          The argument of func should be an instance of FLASH1dFile and
          the argument of processing should be a list of collected data
          
          The list of collected results is attached to func, so they are 
          available for recursions (e.g. last shock position to determine 
          the current one).
          
        Example:
          Imagine you want to collect the simulation times and store them
          in a numpy array.
          (numpy imported as np "as always")
          
          def simTime(simFile):
              time = simFile.simulationTime
              return time
              
          def toNumpyArray(myList):
              myArray = np.array(myList)
              return myArray
              
          timeArray = chkfileloop(simTime,toNumpyArray)
          
          # Alternatively, using lambda functions 
          # which is more convenient for such short tasks.
          
          simTime = (lambda simFile: simFile.simulationTime)
          toNumpyArray = (lambda x: np.array(x))
          
          timeArray = chkfileloop(simTime,toNumpyArray)
          
        """
        if self._time2Num is None:
            chkTimes = []
        bouncetime = self.bouncetime
        
        nMax = len(self.chknumbers) # number of chk files
        pathGenerator = self.chkFilePaths()
        result = []
        for n, path in enumerate(pathGenerator,start=1):
            with FLASH1dFile(path) as simFile:
                filetime = simFile.simulationTime - bouncetime
                func.result = result # attach former results to func, so
                                     # you can use it inside func
                func.filetime = filetime # time relative to bounce
                funcResult = func(simFile)
                result.append(funcResult)
                if self._time2Num is None:
                    chkTimes.append(filetime)
                
            # print status
            percent = float(100*n)/nMax
            message = "\rcollecting data: {p:3.1f}%".format(p=percent)
            sys.stdout.write(message)
            sys.stdout.flush()
        print " done"
        
        if self._time2Num is None:
            self.chkTimes = np.array(chkTimes)
            self._time2Num = dict(zip(self.chkTimes,self.chknumbers))
        
        if processing is not None:
            result = processing(result)
            
        return result
        
        
    def getChkAtTime(self,time):
        """Retruns the path to the checkpoint file closest to the given time.
        
        Arg:
          time (float): time to look for
        
        Note:
            Time is considered relative to bouncetime. If bouncetime is not
            available, it uses the simulation time. 
            This method creates the dictionary _time2Num which maps times to 
            checkpoint numbers when called for the first time.
            If the time is out of bounds return an error message
        """
        # TODO: make getChkAtTime list able
        
        if (not self.bounce):
            print "Time is considered as simulation time."
        bouncetime = self.bouncetime
        
        if self._time2Num is None:
            getTime = (lambda simFile: simFile.simulationTime)
            toNumpyArray = (lambda x: np.array(x))
            chkTimes = self.chkfileloop(getTime,toNumpyArray)
            chkTimes -= bouncetime
            self.chkTimes = chkTimes
            self._time2Num = dict(zip(chkTimes,self.chknumbers))
        else:
            chkTimes = self.chkTimes
        
        minTime = chkTimes[0]
        maxTime = chkTimes[-1]
        if time < minTime or time > maxTime:
            errorMessage = ("The given time is out of bounds!\n"
                "time given: {time}\n"
                "[min,max]: [{min_},{max_}]").format(time=time,
                                                        min_=minTime,
                                                        max_=maxTime)
            raise ValueError(errorMessage)
        
        nearestTime = mathRecipes.find_nearest(chkTimes,time)
        chkNum = self._time2Num[nearestTime]
        path2chkfile = self.path2base + '_hdf5_chk_{n:04d}'.format(n=chkNum)
        
        return path2chkfile
        
                
        
def main():
    """function to test this module"""
    import time
    
    run1 = CCSN1dRun('test/run')
    print "time of bounce:",  run1.bouncetime
    
    getTime = (lambda simFile: simFile.simulationTime)
    chkNums = run1.chknumbers
    times = run1.chkfileloop(getTime)
    n2t = dict(zip(chkNums,times))
    print "checkpoint number to simulation time dictionary:"
    print n2t
    
    print "checkpoint closest to bounce:"
    print run1.getChkAtTime(0.)
    
    print "\nOk, it works!"
        
if __name__ == '__main__':
    main()  
        

