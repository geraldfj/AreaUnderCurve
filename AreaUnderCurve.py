#The working model for hard coding to estimate area under the curve
#For C30ab and C30bb (IS) hopane


#C30ab          -- 41.996 mins
#C30bb (IS)     -- 44.011 mins

import csv
import sys
import math
import matplotlib
matplotlib.use("TkAgg")
import numpy as np
from scipy.integrate import simps
from numpy import trapz


from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg  #NavigationToolbar to go back, forward, zoom etc
from matplotlib.figure import Figure
from matplotlib import pyplot as plt

import tkinter as tk
from tkinter import ttk

fig1 = plt.figure(1)


#Subprogram to find nearest coordinates on mouse click
def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return array[idx]

#Simple mouse click function to store coordinates
def onclick1(event):
    global ix, iy
    ix, iy = event.xdata, event.ydata

    #assign global variable to access outside of function
    global coords
    coords.append((ix, iy))
    
    #Disconnect after 4 clicks
    if len(coords) == 4:
        fig1.canvas.mpl_disconnect(cid1)
        plt.close(1)
    return


#Selecting the data points between the selected retention times
def data_update(C30ab,C30abRT,C30abRT_left,C30abRT_right):
    global begin, end
    global C30abRTupdate,C30abupdate
    C30abRTupdate = []
    C30abupdate = []
    for i in range(len(C30ab)):
        if C30abRT[i] == C30abRT_left:
            begin = i
        if C30abRT[i] == C30abRT_right:
            end = i

    for i in range(len(C30ab)):
        if (i>= begin) and (i<=end):
            C30abRTupdate.append(float(C30abRT[i]))
            C30abupdate.append(float(C30ab[i]))
    return C30abRTupdate,C30abupdate
            
#Normalizing ydata to remove background
def normal(C30ab,C30abRT_left,C30ab_left,C30abRT_right,C30ab_right):
    global C30abnormal
    C30abnormal = []

    slope = (C30ab_right-C30ab_left)/(C30abRT_right-C30abRT_left)
    constant = C30ab_left - slope*C30abRT_left

    for i in range(len(C30ab)):
        y_ = slope*C30abRT[i] + constant
        C30abnormal.append(float(C30ab[i]-y_))
              
    return C30abnormal

    

#trapezoidal rule area under the curve
def trapez(C30ab):
    h = 0.006
    fxsum = 0
    
    for i in range(1,(len(C30ab)-2)):
        fxsum += C30ab[i]

    integration = h*(C30ab[0] + 2*fxsum + C30ab[(len(C30ab)-1)])*0.5

    return integration
    

def main():
    n = int(input('Enter number of files: '))
    
    for i in range(n):
        #path is data from Auburn GC-QQQ on F1 fraction
        path = 'data'+str(i+1)+'.csv'
        file = open(path,newline='')
        reader = csv.reader(file,delimiter=';')

        header = next(reader)  #First line is the header
        rentime = []
        hopane = []

        #Retention time, what time the compound exits the column
        C30abmin = 41.996
        C30ISmin = 44.011
            
        C30abRT = []
        C30ab = []

        C30ISRT = []
        C30IS = []


        for column in reader:
            rentime.append(i,round(float(column[1]),3))   #Reads the retention time
            hopane.append(i,float(column[137]))           #Reads the abundance of hopane at that retention time
               
        for j in range(len(rentime)):
        #C30ab hopane
            if rentime[i][j] == C30abmin:
                for j in range(i-20,i+30):
                    C30abRT.append(float(rentime[i][j]))
                    C30ab.append(float(hopane[i][j]))
                #C30IS hopane
            if rentime[i][j] == C30ISmin:
                for j in range(i-15,i+30):
                    C30ISRT.append(float(rentime[i][j]))
                    C30IS.append(float(hopane[i][j]))
          

        ax1 = fig1.add_subplot(211)
        ax2 = fig1.add_subplot(212)
        ax1.set_ylabel("Abundance")
        ax2.set_ylabel("Abundance")
        ax2.set_xlabel("Time (mins)")
        ax1.grid(True)
        ax2.grid(True)
        ax1.plot(C30abRT,C30ab, color = 'r')
        ax2.plot(C30ISRT,C30IS, color = 'b')
        ax1.fill_between(C30abRT,C30ab,0,color = 'r',alpha=0.2)
        ax2.fill_between(C30ISRT,C30IS,0,color = 'b',alpha=0.2)
        ax1.plot([C30abmin,C30abmin],[0, 1.2*max(C30ab)],linestyle = '--',linewidth = 0.5,color = 'r')
        ax2.plot([C30ISmin,C30ISmin],[0, 1.2*max(C30IS)],linestyle = '--',linewidth = 0.5,color = 'b')

        #Call click function
        coords = []
        cid1 = fig1.canvas.mpl_connect('button_press_event', onclick1)
        plt.show(1)

        #Selecting x_left,y_left & x_right, y_right datapoints closest to the points selected on graph
        C30abRT_left = np.array(find_nearest(C30abRT, coords[0][0]))
        C30ab_left = np.array(find_nearest(C30ab, coords[0][1]))
        C30abRT_right = np.array(find_nearest(C30abRT, coords[1][0]))
        C30ab_right = np.array(find_nearest(C30ab, coords[1][1]))

        C30ISRT_left = np.array(find_nearest(C30ISRT, coords[2][0]))
        C30IS_left = np.array(find_nearest(C30IS, coords[2][1]))
        C30ISRT_right = np.array(find_nearest(C30ISRT, coords[3][0]))
        C30IS_right = np.array(find_nearest(C30IS, coords[3][1]))



        #Selecting the data points withing the x_left,y_left & x_right, y_right
        C30abRT,C30ab = data_update(C30ab,C30abRT,C30abRT_left,C30abRT_right)
        C30ISRT,C30IS = data_update(C30IS,C30ISRT,C30ISRT_left,C30ISRT_right)

        #Normalizing to C30ab minimum - This is to remove the background abundance level
        C30ab = normal(C30ab,C30abRT_left,C30ab_left,C30abRT_right,C30ab_right)
        C30IS = normal(C30IS,C30ISRT_left,C30IS_left,C30ISRT_right,C30IS_right)

        #Area under the curve
        C30ab_int = trapez(C30ab)
        C30IS_int = trapez(C30IS)
        print(round((C30ab_int),3))
        print(round((C30IS_int),3))


if __name__ == '__main__':
    main()




