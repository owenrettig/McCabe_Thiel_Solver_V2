import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import root_scalar
import os
os.system('cls')

# Import data from excel
data = pd.read_excel("VLE_Data.xlsx", sheet_name=1)

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
if reflux_ratio != 0:
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
else:
    B = (z * f - xD * f)/(xB - xD)
    D = f - B
    Vbar = boilup_ratio * D
    if q == 1:
        Lbar = Vbar + B
        V = Vbar 
        L = Lbar - f
    elif q == 0:
        Lbar = Vbar + B
        L = Lbar
        V = Vbar + f
    else:
        Lbar = Vbar + B
        L = Lbar - q*f
        V = Vbar + (1-q)*f

######################################
   
# regress vle data
coeffs = np.polyfit(xdata,ydata, deg = 10)
y_pred = np.polyval(coeffs,xdata)
p = np.poly1d(coeffs)
xfit = np.linspace(0,1,5000)
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
    elif q != 1:
        h = lambda x: (tslope - fslope) * x + (tintercept - fintercept)
        result = root_scalar(h , x0 = z)
        fintercept = result.root
    # bottom slope
    bslope = ((tslope * fintercept + tintercept) - xB) / (fintercept - xB)
    bintercept = -(bslope - 1) * xB
else:
    # top slope
    bslope = Lbar / Vbar 
    bintercept = -(bslope - 1) * xB
    #hande situation where q is equal to 0
    if q == 0:
        k = lambda x: bslope * x + bintercept - z
        result = root_scalar(k, x0 = z)
        fintercept = result.root
    # handle situation where q is not 1
    if q != 1:
        k = lambda x: (bslope - fslope) * x + (bintercept - fintercept)
        result = root_scalar(k, x0 = z)
        fintercept = result.root
    # top slope
    tslope = ((bslope * fintercept + bintercept) - xD) / (fintercept - xD)
    tintercept = -(tslope - 1) * xD

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
        xresult = root_scalar(g, x0 = xs[i-1], x1 = xs[i-1] - 0.1)
        xs[i] = xresult.root
        ys[i] = p(xs[i])
        xs[i+1] = xs[i]
        ys[i+1] = combined_line(xs[i+1])
        if ys[i+1] < xB:
            n = (i+1)/2
            print(f"{n} Stages required for separation")
            break
        if i > 996:
            print("Separation not possible")
            break
        i += 2
    return xs[:i+2], ys[:i+2]

xs, ys = stage_stepper()


#plot the chart
def plot_chart():
    fig = plt.figure("Plotted Line")
    ax = plt.gca()
    
    # Set background colors (MATLAB dark mode)
    fig.patch.set_facecolor('#121212')  # outer figure background
    ax.set_facecolor("#1c1c1c")  # plot area background
    
    # MATLAB default color cycle (same colors work in dark mode)
    matlab_colors = ['#0072BD', '#D95319', '#EDB120', '#7E2F8E', '#77AC30', '#4DBEEE', '#A2142F','#ffffff']
    
    plt.plot(xfit, yfit, color=matlab_colors[0], linewidth=1.5)
    plt.plot(xcalc, yop, color=matlab_colors[3], linewidth=1.5)
    plt.plot(xs, ys, color=matlab_colors[4], linewidth=1)
    plt.plot([0,1], [0,1], color=matlab_colors[2], linewidth=1)
    plt.plot(xdata, ydata, "*", color=matlab_colors[7], markersize=3, alpha = 0.3)
    
    if q != 1 and q != 0:
        plt.plot(xcalc, yq, color=matlab_colors[1], linewidth=1)
    else:
        plt.plot(xq, yq, color=matlab_colors[1], linewidth=1)
    
    plt.xlim([0,1])
    plt.ylim([0,1])
    plt.xlabel("$x_1$", fontsize=11, color='white')
    plt.ylabel("$y_1$", fontsize=11, color='white')
    
    # Style axes and ticks (MATLAB dark mode)
    ax.tick_params(colors='white', which='both', labelsize=10)
    ax.spines['top'].set_color('#4d4d4d')
    ax.spines['bottom'].set_color('#4d4d4d')
    ax.spines['left'].set_color('#4d4d4d')
    ax.spines['right'].set_color('#4d4d4d')
    
    plt.grid(True, color='#404040', linestyle='-', linewidth=0.5)
    plt.show()
   
plot_chart()