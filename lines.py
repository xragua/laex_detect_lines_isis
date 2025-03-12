import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
from scipy import signal
from scipy.optimize import curve_fit
from scipy.signal import savgol_filter, find_peaks
from scipy.stats import mode
import sys

percentile_amp = 0;# int(input("Enter value for prominence percentile (a value from 0 to 100): "))


min_rsq = -10000

if len(sys.argv) > 1:
    name = sys.argv[1]
    print(f"Received name: {name}")
else:
    name = ""
    print("No name provided; using default.")

# HELPER FUNCTIONS
#................................................................................................
def moving_average(data, window_size):
    pad_width = window_size // 2

    # Pad the data array with the average of available values at both ends to handle edge cases
    padded_data = np.pad(data, pad_width, mode='edge')

    # Compute the moving average using convolution
    weights = np.ones(window_size) / window_size
    moving_avg = np.convolve(padded_data, weights, mode='valid')

    # Check if the length of the moving average is less than the length of the input data
    if len(moving_avg) < len(data):
        moving_avg = np.append(moving_avg, moving_avg[-1])

    if len(moving_avg) > len(data):
        moving_avg = moving_avg[0:-1]

    return moving_avg
#................................................................................................
def base_calculator(y):

    s=[]
    le=[]

    for i in range(3,int(len(y)/3)):
        s.append(moving_average(y, i))
        le.append(len(moving_average(y, i)))

    st=np.transpose(s)
    base=[]
    for i in range(len(st)):
        base.append(min(st[i]))
        
    return base
    
def find_peaks_new(t, x):
    
    min_distance = min(np.diff(t))  # or provide the correct array for minimum distance calculation

    dips = []
    prominences = []
    widths = []
    width_heights = []
    left_ips = []
    right_ips = []
 
    pro = mode(abs(np.diff(x)), keepdims=False)[0]  # min prominence 3* more frequent count diff

    peaks, properties = find_peaks(x,  prominence=0 , width=0)

    dips.append(peaks)
    prominences.append(properties['prominences'])
    widths.append(properties['widths'])
    width_heights.append(properties['width_heights'])
    left_ips.append(properties['left_ips'])
    right_ips.append(properties['right_ips'])

    dips = np.concatenate(dips)
    prominences = np.concatenate(prominences)
    widths = np.concatenate(widths)
    width_heights = -np.concatenate(width_heights)
    left_ips = np.concatenate(left_ips)
    right_ips = np.concatenate(right_ips)

    pose = [int(divmod(value, 1)[0]) for value in right_ips]
    reste = divmod(right_ips, 1)[1]

    posi = [int(divmod(value, 1)[0]) for value in left_ips]
    resti = divmod(left_ips, 1)[1]

    diff_num = np.insert(np.diff(t), 0, 0, axis=0)

    ti = t[posi] + diff_num[posi] * resti
    te = t[pose] + diff_num[pose] * reste

    dur = te - ti

    t_dip = t[dips]

    sorted_indices = np.argsort(dips)

    sdips = dips[sorted_indices]
    sprominences = prominences[sorted_indices]
    swidths = widths[sorted_indices]
    swidth_heights = width_heights[sorted_indices]
    sleft_ips = left_ips[sorted_indices]
    sright_ips = right_ips[sorted_indices]
    t_dip = t_dip[sorted_indices]
    ti = ti[sorted_indices]
    te = te[sorted_indices]
    sdur = dur[sorted_indices]

    data = {
        'position': sdips,
        'prominences': sprominences,
        'widths': swidths,
        'width_heights': swidth_heights,
        'left_ips': sleft_ips,
        'right_ips': sright_ips,
        'energy': t_dip,
        'ienergy': ti,
        'eenergy': te,
        'twidth': sdur,
        #'diff_ord': diff_ord
    }

    dips_raw = pd.DataFrame(data)
    dips = dips_raw

    return dips.reset_index(drop=True)
#................................................................................................
def gaussian(x, amplitude, center, sigma):
    """
    """
    return amplitude * np.exp(-(x - center)**2 / (2 * sigma**2))
#................................................................................................
def n_gaussian(x, *params):
    y = np.zeros_like(x)
    
    for i in range(0, len(params), 3):
        amplitude,center,sigma = params[i:i+3]
        y += gaussian(x, amplitude, center, sigma)
    return y
#................................................................................................
def p0_generator(x,y,good_peaks_dataframe):
    
    p0=[]
    bound_low=[]
    bound_high=[]

    for i in range(len(good_peaks_dataframe)):

        p0.append(y[good_peaks_dataframe.position.loc[i]])
        p0.append(good_peaks_dataframe.energy.loc[i])
        
        if ((good_peaks_dataframe.twidth.loc[i]<0.01) & (good_peaks_dataframe.twidth.loc[i]>0)):
            p0.append(good_peaks_dataframe.twidth.loc[i])
        else:
            p0.append(0.005)
            
        bound_low.append(y[good_peaks_dataframe.position.loc[i]]*0.2)
        bound_low.append(good_peaks_dataframe.energy.loc[i]*0.99)
        bound_low.append(0)

        bound_high.append(y[good_peaks_dataframe.position.loc[i]]*2)
        bound_high.append(good_peaks_dataframe.energy.loc[i]*1.01)
        bound_high.append(0.1)

    return p0,(bound_low, bound_high)
#...............................................................................................
    
#################################################################################################


spectra = pd.read_csv(r"spec_0.dat", sep="\\s+",  comment="#", header=None, engine='python')


#FIND LINE CANDIDATES
x = np.array((spectra[0]+spectra[1])/2)
y = np.array(spectra[2])
sy = np.array(spectra[3])

sorted_indices = np.argsort(x)

x = x[sorted_indices]
y = y[sorted_indices]
sy = sy[sorted_indices]

base = base_calculator(y)

new_y = np.where(y < np.array(base)*1, base, y)
ylines = new_y-base
#................................................................................................
#FIND BLOCKS: LINE CANDIDATES #######################################################################
sorted_indices = np.argsort(x)

x = x[sorted_indices]
y = y[sorted_indices]
sy = sy[sorted_indices]

base = base_calculator(y)

new_y = np.where(y < np.array(base)*1, base, y)
ylines = new_y-base


# Find indices where y is non-zero
nonzero_indices = np.nonzero(ylines)[0]

# Initialize dictionary for storing blocks
blocks = {}
xblocks = {}
yblocks = {}
syblocks = {}
contblocks = {}

# Initialize variables to track block start and end indices

start_index = nonzero_indices[0]

n=0

for i in range(len(nonzero_indices)):
        
    if nonzero_indices[i] - nonzero_indices[i-1] == 1:  # Continue current block
        
        end_index = nonzero_indices[i]

        continue
        
        
    else:  # End current block and start a new one
        
        end_index = nonzero_indices[i]
        
        block = np.array(ylines[start_index-1:end_index])
        xblock = np.array(x[start_index-1:end_index ])
        yblock = np.array(y[start_index-1:end_index ])
        syblock = np.array(sy[start_index-1:end_index ])
        
        contblock = np.array(base[start_index-1:end_index])
        
        #plt.plot(xblock,block,"k")
        #plt.plot(xblock,yblock,"r")
        
        blocks[n] = block
        xblocks[n] = xblock
        yblocks[n] = yblock
        syblocks[n] = syblock
        contblocks[n] = contblock
        
        start_index = nonzero_indices[i]
        n=n+1


# Iterate over the blocks and remove consecutive zeros
for i in range(len(blocks)):
    
    for j in range(1, len(blocks[i])):

        if ((blocks[i][j] == 0) and (blocks[i][j-1] == 0) and (len(blocks[i]) > 2)):

            blocks[i] = blocks[i][:j]
            
            xblocks[i] = xblocks[i][:j]
            yblocks[i] = yblocks[i][:j]
            syblocks[i] = syblocks[i][:j]

            
            contblocks[i] = contblocks[i][:j]
            
            break  # Exit the inner loop once consecutive zeros are removed

###################################################################################

nrows = int(len(blocks) / 6) # Calculate number of rows based on 6 columns
ncols = 6  # Number of columns per row


fitted_lines_and_errors_g = pd.DataFrame(columns=['amplitude','center','sigma',
                                                  'eamplitude','ecenter','esigma',
                                                  'rsq'])
z=0

for j in range(len(xblocks)):

    xblock = xblocks[j]  # Access corresponding xblock and yblock elements outside the loop
    yblock = yblocks[j]
    syblock = syblocks[j]
    contblock = contblocks[j]
    
    if (len(yblock) > 3) and (max(yblock) > 0):
        
        res_dif = 0.001

        p = find_peaks_new(xblock, yblock)
        
        
        if len(p)>0:

            the_good_list=[]

            for i in range(len(p)):
            
                prominence_ratio = p.prominences[i]/max(p.prominences)

                idx1 =(abs((p.energy[i]-p.energy))<res_dif)

                if (all((p.prominences[i]-p.prominences[idx1])>=0) & (prominence_ratio>0.1)):

                    the_good_list.append(i)
                    
            good_peaks = p.loc[the_good_list].sort_values(by="prominences").reset_index(drop=True)
            
            max_peaks = int(np.floor(len(xblock)/4))
            good_peaks = good_peaks[0:max(1,max_peaks)]
            
            
            if len(good_peaks)>0:
            
                p0,bounds = p0_generator(xblock,yblock,good_peaks)

                try:
                    popt, pcov = curve_fit(n_gaussian, xblock, yblock, p0=p0, bounds=bounds, maxfev=100000)
                   
                    errors = np.sqrt(np.diag(pcov))
                    yfit = n_gaussian(xblock,*popt)
                    rsq = 1 - np.sum((yblock - yfit)**2) / np.sum((yblock - np.mean(yblock))**2)
              
                    popt_ = np.reshape(popt, (-1, 3))
                    errors_ = np.reshape(errors, (-1, 3))

                    for k in range(len(good_peaks)):
                        
                        new_line = np.concatenate([popt_[k],errors_[k],[rsq]])

                        fitted_lines_and_errors_g.loc[z] = np.concatenate([new_line])
                        z=z+1
                    
                except Exception as e:

                    print(f"Error fitting block {j}: {e}")
                    
#########################################################################################################
min_diff_positions = []

for i in range(len(fitted_lines_and_errors_g)):
    min_diff_index = np.argmin(np.abs(x - fitted_lines_and_errors_g.center.loc[i]))
    min_diff_positions.append(min_diff_index)

base_on_line = [base[pos] for pos in min_diff_positions]
value_on_line = [y[pos-3:pos+3] for pos in min_diff_positions]
fitted_lines_and_errors_g["base_on_line"] = base_on_line
fitted_lines_and_errors_g["value_on_line"] = value_on_line
fitted_lines_and_errors_g["relative_amplitude"] = fitted_lines_and_errors_g.value_on_line-fitted_lines_and_errors_g.base_on_line
    
#########################################################################################################


amplitude = fitted_lines_and_errors_g.relative_amplitude
samplitude = fitted_lines_and_errors_g.eamplitude
center = fitted_lines_and_errors_g.center
scenter = fitted_lines_and_errors_g.ecenter
sigma = fitted_lines_and_errors_g.sigma
ssigma = fitted_lines_and_errors_g.esigma


# Calculating area under the curve
sqrt_2pi = np.sqrt(2 * np.pi)
fitted_lines_and_errors_g["area"] = sqrt_2pi * amplitude * sigma

# Partial derivatives for error propagation
error_amplitude_partial = sqrt_2pi * sigma * samplitude
error_center_partial = sqrt_2pi * amplitude * scenter
error_sigma_partial = sqrt_2pi * amplitude / (2 * sigma)

# Calculating error in area using error propagation formula
fitted_lines_and_errors_g["delta_area"] = ((error_amplitude_partial * samplitude)**2
                                                  + (error_center_partial * scenter)**2
                                                  + (error_sigma_partial * ssigma)**2)**0.5

#########################################################################################################
fitted_lines_and_errors_g=fitted_lines_and_errors_g.reset_index(drop=True)
#########################################################################################################
min_amp = 0.001



clean_lines = fitted_lines_and_errors_g#[ (fitted_lines_and_errors_g.relative_amplitude > min_amp)
                                        #& (fitted_lines_and_errors_g.rsq > min_rsq)].reset_index(drop=True)

clean_lines = clean_lines[['center','ecenter','area','delta_area', 'sigma','esigma','amplitude' ,'eamplitude','relative_amplitude', 'rsq']]
#clean_lines.to_csv('clean_lines.csv', index=False)


###################### OUTPUT FOR ISIS ######################################################################


all_egauss_lines = all_egauss_lines = "+".join([f"egauss({i}) \n" for i in range(1, len(clean_lines) + 1)])

with open(f"set_line_model{name}.sl", "w") as file:
    # Write the definition of the 'lines' function with all egauss terms
    file.write(f"public define linemodel(){{\n\t{all_egauss_lines};\n}}\n\n")

with open(f"set_line_parameters{name}.sl", "w") as file:

    # Write the set_par commands for center
    for i in range(len(clean_lines)):
        file.write(f"set_par(\"egauss({i+1}).center\", {clean_lines['center'].iloc[i]}, 1, "
                   f"{max(clean_lines['center'].iloc[i] - 2*clean_lines['ecenter'].iloc[i], 0)}, "
                   f"{clean_lines['center'].iloc[i] + 2*clean_lines['ecenter'].iloc[i]});\n")

    file.write("\n")

    # Write the set_par commands for area
    for i in range(len(clean_lines)):
        #file.write(f"set_par(\"egauss({i+1}).area\", {clean_lines['area'].iloc[i]}, 0, "
        file.write(f"set_par(\"egauss({i+1}).area\", 0, 1, "
                   f"{0}, "
                   f"{100000000});\n")

    file.write("\n")

    # Write the set_par commands for sigma
    for i in range(len(clean_lines)):
        #file.write(f"set_par(\"egauss({i+1}).sigma\", {clean_lines['sigma'].iloc[i]}, 0, "
        file.write(f"set_par(\"egauss({i+1}).sigma\", 0.0001, 1, "
                   #f"set_par(\"egauss({i+1}).sigma\", {clean_lines['sigma'].iloc[i]}, 0, "
                   f"{0}, "
                   f"{0.005});\n")
                   #f"{min(clean_lines['sigma'].iloc[i] + clean_lines['esigma'].iloc[i],0.15)});\n")

