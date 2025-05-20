import numpy as np 
from matplotlib import pyplot as plt
from matplotlib import rc
import glob as glob

plt.rcParams["axes.grid"] = False

rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

def generate_specific_rows(filePath, userows=[]):
    with open(filePath) as f:
        for i, line in enumerate(f):
            if i in userows:
                yield line

def meanAndError(val, count):
    if sum(count) > 0:
        mean = sum(val*count)/sum(count)
        error = np.sqrt( (sum(count* ((val-mean)**2) ))/sum(count) )
    else:
        mean = 0.
        error = 0.
    return [mean, error]

def AverageFrom3DHisto(path_total):
    #Function that takes data from a 2D histogram with x and y axis, and count in z
    #and calculate the average and error of the y values for each x bin separately
    #Error bands (eb) are calulates for:
    #   - The error comparing between groups, eb_1 (Similar to jackknife error, though just np.std between groups)
    #   - The total errror of all data, ignoring groups, eb_2
    group_counter = 0 #Counts how many different groups with the same random_seed_end exists
    temp_mean = []
    temp_error = []
    x = np.array([])
    y = np.array([])
    y_err = np.array([])
    y_err_2 = np.array([])
    shift_x = 0.
    shift_y = 0.
    #Read all in path with corresponding random_seed_end
    for np_name in glob.glob(path_total):
        a = [0]
        gen = generate_specific_rows(np_name, userows=a)
        data = np.loadtxt(gen, unpack='true')
 
        N1 = int(data[2])
        N2 = int(data[5])

        shift_x = (data[4]-data[3])/(2.*N1)
        shift_y = (data[7]-data[6])/(2.*N2)


        data = np.loadtxt(np_name, skiprows=1)
        x = data[0::N2,0]
        x2 = data[0:N2,1]

        z = np.transpose(np.array(data[:,2]).reshape(N1, N2))
        
        #Save mean y for each x, for this group only
        temp_mean_vec = []
        temp_error_vec = []
        for column in z.T:
            me_temp = meanAndError(x2, column)
            temp_mean_vec.append(me_temp[0])
            temp_error_vec.append(me_temp[1])
        temp_mean.append(np.asarray(temp_mean_vec))
        temp_error.append(np.asarray(temp_error_vec))
        group_counter = group_counter + 1
    
    #Find mean and error across all groups
    if group_counter > 0 and temp_mean:
        temp_mean = np.vstack(temp_mean)
        temp_error = np.vstack(temp_error)

        y = np.mean(temp_mean, axis=0)
        y_err = np.std(temp_mean, axis=0)
        y_err_2 = np.mean(temp_error, axis=0)  

    #Create error bands
    eb_1_bottom = y - y_err
    eb_1_bottom[eb_1_bottom < 0] = 0.
    eb_1_top = y + y_err
    eb_1 = [eb_1_bottom, eb_1_top]

    eb_2_bottom = y - y_err_2
    eb_2_bottom[eb_2_bottom < 0] = 0.
    eb_2_top = y + y_err_2
    eb_2 = [eb_2_bottom, eb_2_top]

    return np.asarray(x)+shift_x, np.asarray(y)+shift_y, np.asarray(y_err), np.asarray(y_err_2), eb_1, eb_2, group_counter


def Plot2D(path_total):
    #Function to collect the data from 1D histograms and return as x (bins) and y (count per bin)
    #and also calculate the error and error bands, with error_1 
    group_counter = 0
    temp_mean_vec = []
    x_shift = 0
    x = []
    y = []
    y_tot = []
    for np_name in glob.glob(path_total):
        data = np.loadtxt(np_name, skiprows=1)
        if group_counter == 0:
            x = data[:,0]
            x_shift = np.abs(x[1] - x[0])/2.
            y_tot = data[:,1]
            y_temp = data[:,1]
            group_counter = 1
        else:
            y_tot = y_tot + data[:,1]
            y_temp = data[:,1]
        temp_me = meanAndError(x,y_temp)
        temp_mean_vec.append(temp_me[0])
    
    
    mean_err_tot = meanAndError(x,y_tot)
    y_mean = mean_err_tot[0]
    y_err_1 = np.std(temp_mean_vec)
    y_err_2 = mean_err_tot[1]
 
    return np.asarray(x), y_tot, y_mean, y_err_1, y_err_2, group_counter


