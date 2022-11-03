import McsPy
import h5py
import matplotlib.pyplot as plt
from McsPy import ureg, Q_
import numpy as np
import pandas as pd
import os
import math

plt.style.use('default')
channels = [12, 13, 21, 22, 23, 24, 31, 32, 33, 34, 42, 43]


# Calculates Synchronization Index Between 2 Channels
def synch(ch1, ch2,t):
# Select Channels to compare
# Change to channels of interest (e.g. ch12, ch21, etc.)
    s1 = ch1
    s2 = ch2

# Set up FFT and convolution
# Create complex Morlet Wavelet
    srate = len(ch1)
    center_freq = 5                                         # filter frequency
    time = np.linspace(-1, 1, srate)

  # Create complex sine wave
    sine_wave = np.exp((1j*2) * np.pi * center_freq * time)

  # Create Gaussian window
    s = 7 / (2 * np.pi * center_freq)
    gaus_win = np.exp((-np.power(time,2) / (2 * np.power(s,2))))

  # Create Morlet Wavelet
    cmw = np.multiply(sine_wave, gaus_win)
    half_wavN = (len(time)-1)/2

  #FFT of wavelet
    n_wavelet = len(cmw)
    n_data    = len(t)
    n_conv    = n_wavelet + (n_data-1)

  #FFT of wavelet
    waveletX = np.fft.fft(cmw)
    max_waveletX = np.amax(waveletX)
    waveletXf = np.divide(waveletX, max_waveletX)

  # compute Hz for plotting
    hz_min = 0.0
    hz_max = np.floor((n_conv/2)+1)
    hz_half = np.multiply(srate, 0.5)
    hz = np.linspace(hz_min, hz_max, len(t))


  ## Calculate phase angle for each channel
  # analytic signal of channel 1
    fft_data = np.fft.fft(s1)                                   #fft of channel 1
    ch1_conv = np.multiply(waveletX, fft_data)                  #convolution of channel 1
    an_s1 = np.fft.ifft(ch1_conv)                               #inverse fft to return analog signal

  # collect real and phase data
    phase_data1 = np.angle(an_s1)
    real_data1 = np.real(an_s1)

      # analytic signal of channel 2
    fft_data = np.fft.fft(s2)                                   #fft of channel 2
    ch2_conv = np.multiply(waveletX, fft_data)                  #convolution of channel 2
    an_s2 = np.fft.ifft(ch2_conv)                               #inverse fft to return analog signal

      # # collect real and phase data
    phase_data2 = np.angle(an_s2)
    real_data2 = np.real(an_s2)

      ### CALCULATE PHASE ANGLE DIFFERENCE
    phase_angle_differences = phase_data2-phase_data1               #phase angle difference
    euler_phase_differences = np.exp(1j*phase_angle_differences)    #euler representation of angles
    mean_complex_vector = np.mean(euler_phase_differences)          #mean vector (in complex space)
    phase_synchronization = np.absolute(mean_complex_vector)        #length of mean vector (M from Me^ik)
    return phase_synchronization

#Given Data array, Channels in Data, and the number of measurements for each timepoint (t), this function creates a dictionary that contains all combinations of synchrony indices
def synch_data(data, channels,t):
  synchron_data = {}
  for i in channels:
    temp = {} #new dictionary for each channel
    sum = 0
    for j in channels:
      if i != j: #iterates through all non-equal channel pairs
        key = "ch"+str(j)
        ch1 = data["ch"+str(i)]
        ch2 = data["ch"+str(j)]
        temp[key] = synch(ch1, ch2,t) #calculates and stores synchrony index
        sum += synch(ch1, ch2,t)
    temp["avg"] = sum / (len(channels) - 1) # uses Sum to calculate Average
    synchron_data[i] = temp # Adds channel's synchrony data to the master dictionary
  return synchron_data

# Combination function based on factorials
def nCr(n, r):
    f = math.factorial
    return f(n) // f(r) // f(n - r)

#given 2 points, p1 and p2, and a weight of the line, w, with optional tags, creates a line of thickness w between p1 and p2
def connectpoints(p1, p2, w, tag='k-'):
    x1, x2 = p1[0], p2[0]
    y1, y2 = p1[1], p2[1]
    plt.plot([x1, x2], [y1, y2], tag, linewidth=1-(w**2))

#finds the maximum value of a list of values
def max_val(values):
    k = 0.000001
    for name in values:
        if values[name] > k:
            k = values[name]
    return k

for well in os.listdir(r'E:\Archived Data\Venkat data\Blast48 20220627\CSV_Files\CSV Files'): #iterates through each well folder
    print(well)
    completed_wells = ['Well_0_Blast', "Well_0_Control", 'Well_10_Blast', 'Well_10_Control', 'Well_11_Blast', 'Well_11_Control', 'Well_1_Blast', 'Well_1_Control', 'Well_2_Blast', "Well_3_Control", 'Well_3_Blast', "Well_2_Control", 'Well_4_Blast', "Well_4_Control", 'Well_5_Blast', "Well_5_Control", 'Well_6_Blast', "Well_6_Control", 'Well_7_Blast', "Well_7_Control", "Well_9_Blast", "Well_9_Control", "Well_8_Blast"] #note all completed wells in case of having to re-run
    if str(well) not in completed_wells:
        well_path = os.path.join(r'E:\Archived Data\Venkat data\Blast48 20220627\CSV_Files\CSV Files', well)
        clus_coeff_well = []
        x = []
        n_nodes = []
        sig_node_well = []
        for num_w,file in enumerate(os.listdir(well_path)): #iterates through each data file (phase)
            print(file)

            #creates pandas dataframe from the file
            data = pd.read_csv(os.path.join(well_path, file))
            print(file[:-4])
            if file[:-4].split(sep = '_')[-1] == 'Control':
                x.append(0)
            else:
                x.append(int(file[:-4].split(sep = '_')[-1]))
            r_dict = {'Unnamed: 0': 'time'}
            for i in range(12):
                r_dict[str(i)] = 'ch' + str(channels[i])
            data = data.rename(r_dict, axis='columns')
            print(data.keys())
            t = data['time']

            for i in range(12):
                r_dict[str(i)] = 'ch' + str(channels[i])
            data = data.rename(r_dict, axis='columns')


            sum_val = 0
            synchron_data = synch_data(data, channels,t) #extracts synchrony index data
            print(synchron_data)
            for channel in synchron_data:
              sum_val += (synchron_data[channel]["avg"])
            threshold = sum_val/len(synchron_data.keys()) #sets connection threshold as average (two neuronal masses/channels are considered to be connected if their synchrony index is above this threshold)

            for channel in synchron_data:
              for chan in synchron_data[channel]:
                synchron_data[channel][chan] = [synchron_data[channel][chan], synchron_data[channel][chan] >= threshold]

            # creates a dictionary containing all connections
            connections = {}
            for chan in synchron_data:
              temp = []
              for channel in synchron_data[chan]:
                if channel != 'avg':
                  if synchron_data[chan][channel][1]:
                    temp.append(int(channel[2:]))
              connections[chan] = temp

            #Generates Clustering Coefficient data
            cluster_coeff = {}
            print(connections)
            for channel in connections:
                sum_connections = 0
                for chan in connections[channel]:
                    for chan2 in connections[channel]:
                        if chan2 != chan:
                            if chan in connections[chan2]:
                                sum_connections += 1
                try:
                    cluster_coeff[channel] = sum_connections / 2 / nCr(len(connections[channel]), 2) # calculates clustering coefficients for each channel/node
                except:
                    cluster_coeff[channel] = 0 #if 0 connections exist, the above equation will result in a division by 0, so this sets the CC to 0

            avg_cc = np.mean(list(cluster_coeff.values()))
            clus_coeff_well.append(avg_cc) # calculates and records average clustering coefficient
            cc_threshold = avg_cc

            # calculates and records numbers of significant nodes (nodes with CC above average CC)
            n_node = 0
            sig_nodes = cluster_coeff
            sig_node_well.append(sig_nodes)
            for channel in channels:
                if cluster_coeff[channel] > avg_cc:
                    n_node += 1
            n_nodes.append(n_node)

            #Creates coordinates
            coords = []
            for i in range(4):
                for j in range(4):
                    coords.append([.2 * (i + 1), .2 * (5 - (j + 1))])
            coords.remove([.2, .2])
            coords.remove([.2, .8])
            coords.remove([.8, .2])
            coords.remove([.8, .8])

            k = max_val(cluster_coeff)

            # Assigns channels to their respective coordinates
            channel_coord = {}
            for i in range(len(channels)):
                channel_coord[channels[i]] = coords[i]

            # Generates template for map
            figure, axes = plt.subplots()
            axes.set_aspect(1)
            for channel in channel_coord:
                draw_circle = plt.Circle(tuple(channel_coord[channel]), 0.03,
                                         color=str(1 - (cluster_coeff[channel] ** 5 + 1 - k)))
                axes.add_artist(draw_circle)
            axes.add_artist(draw_circle)
            for i in channels:
                for j in channels:
                    if i != j:
                        ch1 = data["ch" + str(i)]
                        ch2 = data["ch" + str(j)]
                        connectpoints(channel_coord[i], channel_coord[j], synch(ch1, ch2,t))
            phat = os.path.join(well_path, file[:-4]+'_Synchrony_Map')
            plt.savefig(phat)

        CC_channels = {'Phase': x, 'Average Clustering Coefficient': clus_coeff_well, 'Number of Nodes': n_nodes}
        for channel in channels:
            CC_channels[channel] = []
            for ph in sig_node_well:
                CC_channels[channel].append(ph[channel])

        #Plots and saves Average Clustering Coefficient, Number of Significant Nodes by phase, Saves CC values of each channel for each phase
        CC_data = pd.DataFrame(data = CC_channels)
        CC_data = CC_data.sort_values(by='Phase')
        CC_data.to_csv(os.path.join(well_path, str(well)+'_Clustering Coefficient Data.csv'))
        plt.subplot(1, 2, 1)
        plt.plot(CC_data['Phase'].tolist(), CC_data['Average Clustering Coefficient'].tolist(), label='Average Clustering Coefficient')
        plt.title('Average Clustering Coefficient')
        plt.subplot(1, 2, 2)
        plt.plot(CC_data['Phase'].tolist(), CC_data['Number of Nodes'].tolist(), label='Number of Significant Nodes')
        plt.title('Nodes')
        plt.savefig(os.path.join(well_path, str(well)+'_Clustering Coefficient Plots'))

