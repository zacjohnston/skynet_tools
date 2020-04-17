import numpy as np
import matplotlib.pyplot as plt
import sys
from scipy.interpolate import interp1d

# script, first, second = sys.argv


def read_file( read_list, filename ):

    for each in read_list:
        first = each
        with open(filename+first+'.dat','r') as f:
            r, t = [], []
            i = 1
            for line in f:
                if i < 3:
                    print(line)
                    i += 1
                    continue
                t.append(float(line.split()[0]))
                r.append(float(line.split()[3]))
        
        plt.plot(t,r)
    return

def main():
    #   1.930028464363739094e+00
    match = 1.971338528028059089e+00
    for i in range(198):
        x = np.loadtxt('mass2_0.txt')
        x_1 = np.loadtxt('mass2_'+str(i)+'.txt')
        x_index = np.where(x == match)[0][0]
        x1_index = np.where(x_1 == match)[0][0]
        y = np.loadtxt('temp2_0.txt')
        y_1 = np.loadtxt('temp2_'+str(i)+'.txt')
        x = x[x_index:]
        x_1 = x_1[x1_index:]
        y = y[x_index:]
        y_1 = y_1[x1_index:]
        sub = y-y_1
        nonzero = sub[np.nonzero(y-y_1) ]
        print(i, np.nonzero(y-y_1),nonzero, len(nonzero))
        #xnew = np.loadtxt('mass_grid.txt')
        #f = interp1d(x,y)
        #ynew = f(xnew)
        #ynew_test = np.loadtxt('temphold.txt')
    
        #plt.plot(x,y,'o',xnew,ynew,'-')
        #plt.show()        
        #print(ynew - ynew_test)

if __name__ == "__main__":
    main()
