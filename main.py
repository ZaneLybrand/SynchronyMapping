import McsPy
import h5py
import matplotlib.pyplot as plt
from McsPy import ureg, Q_
import numpy as np
import pandas as pd
import os
import math


dest_dir = 'E:\Archived Data\Venkat data\CSV_Files'  # Directory to export csv files to
origin_dir = 'D:\data\Venkat data\H5_Files'  # Directory (Folder) to get h5 files from
image_dir = 'D:\data\Venkat data\Channel_plots'

def extract_streams(recording):
    streams = {}
    if recording.analog_streams is not None:
        an=[]
        for s in recording.analog_streams:
            an.append(recording.analog_streams[s])
        streams['Analog'] = an
    return streams


'''for filename in os.listdir(origin_dir):
    f = os.path.join(origin_dir, filename)
    file_name_list = filename.split('_')
    data = McsPy.McsData.RawData(f)
    rec = data.recordings[0]
    stream = extract_streams(rec)
    ana_str = stream['Analog'][0]
    channels = list(ana_str.channel_infos.keys())
    ch_data = {}
    for i in range(12):
        temp_data = {}
        for j in range(12):
            data = McsPy.McsData.RawData(f)
            rec = data.recordings[0]
            stream = extract_streams(rec)
            ana_str = stream['Analog'][0]
            temp_data[j] = list(ana_str.get_channel(channels[i*12+j])[0])
        fname = 'Well_'+str(i)
        if 'Control_Blast' in filename:
            fname+= '_Control'
        else:
            fname+= '_Blast'
        if 'Control_mwd' in filename:
            fname += '_Control'
        else:
            fname += '_Phase_'+file_name_list[file_name_list.index('Phase') + 1]
        fname += '.csv'
        destination = os.path.join(dest_dir, fname)
        pd.DataFrame.from_dict(temp_data).to_csv(destination)'''
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


def synch_data(data, channels,t):
  synchron_data = {}
  for i in channels:
    temp = {}
    sum = 0
    for j in channels:
      if i != j:
        key = "ch"+str(j)
        ch1 = data["ch"+str(i)]
        ch2 = data["ch"+str(j)]
        temp[key] = synch(ch1, ch2,t)
        sum += synch(ch1, ch2,t)
    temp["avg"] = sum / (len(channels) - 1) #Average
    synchron_data[i] = temp
  return synchron_data


def nCr(n, r):
    f = math.factorial
    return f(n) // f(r) // f(n - r)


def connectpoints(p1, p2, w, tag='k-'):
    x1, x2 = p1[0], p2[0]
    y1, y2 = p1[1], p2[1]
    plt.plot([x1, x2], [y1, y2], tag, linewidth=1-(w**2))

def max_val(values):
    k = 0.000001
    for name in values:
        if values[name] > k:
            k = values[name]
    return k


for num_w,file in enumerate(os.listdir(dest_dir)):
    if os.path.isfile(os.path.join(dest_dir, file)) and 'Well_0_' not in file and 'Well_10_' not in file and os.path.splitext(file)[1] == '.csv':
        data = pd.read_csv(os.path.join(dest_dir, file))
        print(file[:-4])
        image_dir = os.path.join(dest_dir, file[:-4]+'_images')
        os.mkdir(image_dir)
        r_dict = {'Unnamed: 0': 'time'}
        for i in range(12):
            r_dict[str(i)] = 'ch' + str(channels[i])
        data = data.rename(r_dict, axis='columns')
        print(data.keys())
        t = data['time']

        for i in range(12):
            r_dict[str(i)] = 'ch' + str(channels[i])
        data = data.rename(r_dict, axis='columns')
        ch_len = 100000
        fig, ax = plt.subplots()
        for chan1 in channels:
            for chan2 in channels:
                x = []
                y = []
                if chan1 < chan2:
                    for i in range(int(len(data['ch12']) / ch_len)):
                        x.append(i)
                        y.append((synch([data['ch' + str(chan1)][i + j] for j in range(ch_len)],
                                        [data['ch' + str(chan2)][i + j] for j in range(ch_len)],t)))
                    ax.set_ylabel('Synchronization Index')
                    ax.set_xlabel('Iterations')
                    ax.plot(x,y)
                    name = 'Channel ' + str(chan1) + ' and Channel ' + str(chan2) + '.png'
                    title = 'Channel ' + str(chan1) + ' and Channel ' + str(chan2) + ":"
                    ax.set_title(title)
                    plt.savefig(os.path.join(image_dir,name))
                    ax.clear()
        sum_val = 0
        synchron_data = synch_data(data, channels,t)
        print(synchron_data)
        for channel in synchron_data:
          sum_val += (synchron_data[channel]["avg"])
        threshold = sum_val/len(synchron_data.keys()) - 0.2
        for channel in synchron_data:
          for chan in synchron_data[channel]:
            synchron_data[channel][chan] = [synchron_data[channel][chan], synchron_data[channel][chan] >= threshold]

        connections = {}
        for chan in synchron_data:
          temp = []
          for channel in synchron_data[chan]:
            if channel != 'avg':
              if synchron_data[chan][channel][1]:
                temp.append(int(channel[2:]))
          connections[chan] = temp
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
                cluster_coeff[channel] = sum_connections / 2 / nCr(len(connections[channel]), 2)
            except:
                cluster_coeff[channel] = 0
        avg_cc = np.mean(list(cluster_coeff.values()))
        cc_threshold = avg_cc
        sig_nodes = [for channel in cluster_coeff if cluster_coeff[channel]]
        coords = []
        for i in range(4):
            for j in range(4):
                coords.append([.2 * (i + 1), .2 * (5 - (j + 1))])
        coords.remove([.2, .2])
        coords.remove([.2, .8])
        coords.remove([.8, .2])
        coords.remove([.8, .8])

        k = max_val(cluster_coeff)

        channel_coord = {}
        for i in range(len(channels)):
            channel_coord[channels[i]] = coords[i]

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
        phat = os.path.join(dest_dir, file[:-4]+'_Synchrony_Map')
        plt.savefig(phat)
