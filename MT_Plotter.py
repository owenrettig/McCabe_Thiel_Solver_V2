import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import root_scalar
import os
os.system('cls')

# Import data from excel
data = pd.read_excel("VLE_Data.xlsx", sheet_name=1)
data.iloc[0, 1] = 1e-10

# Assign data from excel to values
xcalc = np.linspace(0,1,5000)
xdata = data["x"].values
ydata = data["y"].values
f = data.iloc[0,3]
z = data.iloc[1,3]
reflux_ratio = data.iloc[2,3]
boilup_ratio = data.iloc[3,3]
xD = data.iloc[4,3]
xB = data.iloc[5,3]
q = data.iloc[6,3]



######################################

# Calculate D, B, L, V, Lbar, Vbar
B = (z * f - xD * f)/(xB - xD)
D = f - B
L = reflux_ratio * D
if q == 1:
    V = L + D
    Vbar = V
    Lbar = L + f
elif q != 0:
    V = L + D
    Vbar = V + (1 - q) * f
    Lbar = L + q * f
else:
    Vbar = L + D
    V = Vbar + f
    Lbar = L
    
# regress vle data
coeffs = np.polyfit(xdata,ydata, deg = 8)
y_pred = np.polyval(coeffs,xdata)
p = np.poly1d(coeffs)
xfit = np.linspace(1e-3,1,500)
yfit = p(xfit)

######################################

# create the feed line equation 
if q == 1:
    xq = np.array([z,z])
    yq = np.array([0,1])
    fintercept = z
elif q == 0:
    xq = np.array([0,1])
    yq = np.array([z,z])
    fintercept = z
else:
    fslope = q / (q - 1)
    fintercept = z / (1 - q)
    yq = fslope * xcalc + fintercept

######################################

# operating line equation
if reflux_ratio != 0:    
    # top slope
    tslope = L / V
    tintercept = (1 - tslope) * xD
    # handle situation were q is equal to 0
    if q == 0:
        h = lambda x: tslope * x + tintercept - z
        result = root_scalar(h , x0 = z)
        fintercept = result.root
    # handle situation where q is not 1 
    if q != 1 and q !=0 :
        h = lambda x: (tslope - fslope) * x + (tintercept - fintercept)
        result = root_scalar(h , x0 = z)
        fintercept = result.root
    # bottom slope
    bslope = ((tslope * fintercept + tintercept) - xB) / (fintercept - xB)
    bintercept = -(bslope - 1) * xB
else:
    tslope = 1

######################################

# combined line
def combined_line(x):
    x = np.asarray(x)
    return np.where(x > fintercept,
                      tslope * x + tintercept,
                      bslope * x + bintercept)

yop = combined_line(xcalc)

######################################

# create stage lines 
def stage_stepper():
    max_stages = 1000
    xs = np.zeros(max_stages)
    ys = np.zeros(max_stages)
    xs[0] = xD
    ys[0] = xD
    i = 1

    while True:
        g = lambda x: p(x) - (ys[i-1])
        xresult = root_scalar(g, x0 = xs[i-1], x1 = xs[i-1] + 0.1)
        xs[i] = xresult.root
        ys[i] = p(xs[i])
        xs[i+1] = xs[i]
        ys[i+1] = combined_line(xs[i+1])
        if ys[i+1] < xB:
            break
        if i > 996:
            print("Separation not possible")
            break
        i += 2
    return xs[:i+2], ys[:i+2]

xs, ys = stage_stepper()


######################################

#plot the chart
def plot_chart():
    plt.figure("Plotted Line")
    plt.plot(xdata,ydata,"r*")
    plt.plot(xfit,yfit,"k")
    plt.plot(xcalc,yop)
    plt.plot(xs,ys)
    plt.plot([0,1],[0,1])
    if q != 1 and q !=0 :
        plt.plot(xcalc,yq)
    else:
        plt.plot(xq,yq)
    plt.xlim([0,1])
    plt.ylim([0,1])
    plt.xlabel("x_1")
    plt.ylabel("y_1")
    plt.grid(True)
    plt.show()
   
plot_chart()