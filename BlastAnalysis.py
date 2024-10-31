import os
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import json
import math
import seaborn as sns
import McsPy
import shutil
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
  return synchron_data #returns a dictionary of dictionaries

# Combination function based on factorials
def nCr(n, r):
    f = math.factorial
    return f(n) // f(r) // f(n - r)

def str_to_dict(s):
    s=s.replace('\'', '\"')
    return json.loads(s)


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

# Extracts raw data streams
def extract_streams(recording):
    streams = {}
    if recording.analog_streams is not None:
        an = []
        for s in recording.analog_streams:
            an.append(recording.analog_streams[s])
        streams['Analog'] = an
    return streams

data_dir = r"G:\Shared drives\Lybrand MEA Drive\Opto blast\hd5 files\opto blast 8w\OptoPharm stim"
or_dir = r"G:\Shared drives\Lybrand MEA Drive\Blast48 CSV Files"
dest_dir = r"G:\Shared drives\Lybrand MEA Drive\NEW Blast Analysis"

for filename in os.listdir(data_dir):
    print(filename)
    f = os.path.join(data_dir, filename)
    file_name_list = filename.split('_')
    if '.h5' in file_name_list[-1]:
        data = McsPy.McsData.RawData(f)
        rec = data.recordings[0]
        stream = extract_streams(rec)
        ana_str = stream['Analog'][0]
        channels = list(ana_str.channel_infos.keys())
        print(len(channels))
        ch_data = {}
        for i in range(4):
            temp_data = {}
            for j in range(12):
                data = McsPy.McsData.RawData(f)
                rec = data.recordings[0]
                stream = extract_streams(rec)
                ana_str = stream['Analog'][0]
                temp_data[j] = list(ana_str.get_channel(channels[i*12+j])[0])
            fname = 'Well_'+str(i)
            if 'Gabazine' in filename:
                fname+= r"_GABA"
            elif 'lido' in filename:
                fname+= r"_Lidocaine"
            elif "CNQX_AP5" in filename:
                fname+= r"_Glutamate"
            else:
                fname+= r"_Baseline"

            if 'Control_8w' in filename:
                fname+= r"_Control"
            else:
                fname+= r"_Blast"
            fname += r".csv"
            destination = os.path.join(or_dir, fname)
            pd.DataFrame.from_dict(temp_data).to_csv(destination)

for file in os.listdir(or_dir):
    well = str(file)
    if well in os.listdir(dest_dir) and 'Well_' in well:
        print(file)
        or_well = os.path.join(or_dir,file)
        dest_well = os.path.join(dest_dir, well)
        csv_path = os.path.join(dest_well, 'CSV_Files')
        '''
        for csv_file in os.listdir(or_well):
            if "Clustering" not in csv_file and "Synchrony" not in csv_file:
                os.rename(os.path.join(or_well, csv_file), os.path.join(csv_path, csv_file))
        '''
        synch_map_path = os.path.join(dest_well, 'Synchrony Maps')
        try:
            os.makedirs(synch_map_path)
        except:
            continue
        comp_synch = os.path.join(synch_map_path, 'Complete Synchrony Maps')
        norm_synch = os.path.join(synch_map_path, 'Normal Synchrony Maps')
        simp_synch = os.path.join(synch_map_path, 'Simplified Synchrony Maps')
        os.makedirs(comp_synch)
        os.makedirs(norm_synch)
        os.makedirs(simp_synch)
        synch_cc_data = os.path.join(dest_well, 'Other Data')
        try:
            os.makedirs(synch_cc_data)
        except:
            continue
        heat_dir = os.path.join(dest_well, 'Heatmaps')
        try:
            os.makedirs(heat_dir)
        except:
            continue
        phases = []
        clus_coeff = []
        n_nodes = []
        total_conn = []
        sig_node_well = []
        syn_data = {'threshold': []}
        for channel in channels:
            syn_data[channel] = []
        for p, file in enumerate(os.listdir(csv_path)):

            #reads the csv file
            data = pd.read_csv(os.path.join(csv_path, file))
            print(file[:-4])

            #adds the phase
            if file[:-4].split(sep='_')[-1] == 'Control':
                phases.append(0)
            else:
                phases.append(int(file[:-4].split(sep='_')[-1]))

            #rename the index to the time and the channels to their appropriate channel number
            r_dict = {'Unnamed: 0': 'time'}
            for i in range(12):
                r_dict[str(i)] = 'ch' + str(channels[i])
            data = data.rename(r_dict, axis='columns')
            print(data.keys())
            t = data['time']
            for i in range(12):
                r_dict[str(i)] = 'ch' + str(channels[i])
            data = data.rename(r_dict, axis='columns')

            #calculate the synchrony index values
            synchron_data = synch_data(data, channels, t)  # extracts synchrony index data
            print(synchron_data)

            #generate heatmaps
            l = list(synchron_data.values())
            for ch in range(len(l)):
                l[ch]=list(l[ch].values())
            hm = np.random.rand(12,12)
            for i in range((12)):
                hm[i][i] = 1
                for j in range(12):
                    if j < i:
                        hm[i][j] = l[i][j]
                    elif j>i:
                        hm[i][j] = l[i][j-1]
            plt.figure(figsize=(12, 12))
            heat_map = sns.heatmap(hm, linewidth=1, annot=True)
            plt.title(well + ' Phase ' + str(p))
            plt.savefig(os.path.join(heat_dir, str(well) + ' Phase ' + str(p) + ' Heatmap'))
            plt.close('all')

            #calculates threshold
            sum_val = 0
            for channel in synchron_data:
                sum_val += (synchron_data[channel]["avg"])
            threshold = sum_val / len(synchron_data.keys())  # sets connection threshold as average (two neuronal masses/channels are considered to be connected if their synchrony index is above this threshold)

            #generates synchrony data
            for channel in channels:
                syn_data[channel].append(synchron_data[channel].copy())
            syn_data['threshold'].append(threshold)

            #calculates connections
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

            # Generates Clustering Coefficient data
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
                    cluster_coeff[channel] = sum_connections / 2 / nCr(len(connections[channel])  # calculates clustering coefficients for each channel/node
                except:
                    cluster_coeff[channel] = 0  # if 0 connections exist, the above equation will result in a division by 0, so this sets the CC to 0

            avg_cc = np.mean(list(cluster_coeff.values()))
            clus_coeff.append(avg_cc)  # calculates and records average clustering coefficient
            cc_threshold = avg_cc

            # calculates and records numbers of significant nodes (nodes with CC above average CC)
            n_node = 0
            sig_nodes = cluster_coeff
            sig_node_well.append(sig_nodes)
            for channel in channels:
                if cluster_coeff[channel] > avg_cc:
                    n_node += 1
            n_nodes.append(n_node)
        #Complete Synchrony Maps
            # Creates coordinates
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
            plt.xlim(0, 1)
            plt.ylim(0, 1)
            for channel in channel_coord:
                draw_circle = plt.Circle(tuple(channel_coord[channel]), 0.01*(len(connections[channel])),
                                         color=str(1 - (cluster_coeff[channel] ** 5 + 1 - k)))
                axes.add_artist(draw_circle)
            axes.add_artist(draw_circle)
            for i in channels:
                for j in channels:
                    if i != j:
                        ch1 = data["ch" + str(i)]
                        ch2 = data["ch" + str(j)]
                        connectpoints(channel_coord[i], channel_coord[j], synch(ch1, ch2, t))
            phat = os.path.join(comp_synch, file[:-4] + '_Complete_Synchrony_Map')
            plt.savefig(phat)
            plt.close('all')

        #Normal Synchrony Maps:
            # Creates coordinates
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
            plt.xlim(0, 1)
            plt.ylim(0, 1)
            for channel in channel_coord:
                draw_circle = plt.Circle(tuple(channel_coord[channel]), 0.01 * (len(connections[channel])),
                                         color=str(1 - (cluster_coeff[channel] ** 5 + 1 - k)))
                axes.add_artist(draw_circle)
            axes.add_artist(draw_circle)
            for i in channels:
                for j in channels:
                    if i != j:
                        ch1 = data["ch" + str(i)]
                        ch2 = data["ch" + str(j)]
                        if synch(ch1, ch2, t) >= threshold:
                            connectpoints(channel_coord[i], channel_coord[j], 0)

            phat = os.path.join(norm_synch, file[:-4] + '_Normal_Synchrony_Map')
            plt.savefig(phat)
            plt.close('all')
        #Simplified Synchrony Maps
            # Creates coordinates
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
            plt.xlim(0,1)
            plt.ylim(0,1)
            for channel in channel_coord:
                if cluster_coeff[channel] == 1:
                    draw_circle = plt.Circle(tuple(channel_coord[channel]), 0.01 * (len(connections[channel])),
                                         color=str(0))
                else:
                    draw_circle = plt.Circle(tuple(channel_coord[channel]), 0.02,
                                         color=str(0.3))

                axes.add_artist(draw_circle)
            axes.add_artist(draw_circle)
            for i in channels:
                for j in channels:
                    if i != j:
                        ch1 = data["ch" + str(i)]
                        ch2 = data["ch" + str(j)]
                        if (cluster_coeff[i] == 1 or cluster_coeff[j] == 1) and synch(ch1, ch2, t) > threshold:
                            connectpoints(channel_coord[i], channel_coord[j], 0)
            phat = os.path.join(simp_synch, file[:-4] + '_Simplified_Synchrony_Map')
            plt.savefig(phat)
            plt.close('all')


        CC_channels = {'Phase': phases, 'Average Clustering Coefficient': clus_coeff, 'Number of Nodes': n_nodes}
        for channel in channels:
            CC_channels[channel] = []
            for ph in sig_node_well:
                CC_channels[channel].append(ph[channel])
        print(CC_channels)
        # Plots and saves Average Clustering Coefficient, Number of Significant Nodes by phase, Saves CC values of each channel for each phase
        CC_data = pd.DataFrame.from_dict(CC_channels)
        CC_data = CC_data.sort_values(by='Phase')
        CC_data.to_csv(os.path.join(synch_cc_data, str(well) + '_Clustering Coefficient Data.csv'), index = False)
        plt.subplot(1, 2, 1)
        plt.plot(CC_data['Phase'].tolist(), CC_data['Average Clustering Coefficient'].tolist(),
                 label='Average Clustering Coefficient')
        plt.title('Average Clustering Coefficient')
        plt.subplot(1, 2, 2)
        plt.plot(CC_data['Phase'].tolist(), CC_data['Number of Nodes'].tolist(), label='Number of Significant Nodes')
        plt.title('Nodes')
        plt.savefig(os.path.join(synch_cc_data, str(well) + '_Clustering Coefficient Plots'))

        syn_data['Phases']=phases
        pd.DataFrame.from_dict(syn_data).to_csv(os.path.join(synch_cc_data, str(well) + '_Synchrony Data.csv'), index = False)
        plt.close('all')
for well in os.listdir(dest_dir):
    if "Stim CSV" not in str(well):
        print(well)
        dest_well = os.path.join(dest_dir, well)
        synch_cc_data = os.path.join(dest_well, 'Other Data')
        CDataWell = os.path.join(synch_cc_data, str(well) + '_Clustering Coefficient Data.csv')
        SDataWell = os.path.join(synch_cc_data, str(well) + '_Synchrony Data.csv')
        coefficient_data = pd.read_csv(CDataWell)
        synch_data = pd.read_csv(SDataWell)
        coeffic_data = coefficient_data.values.tolist()
        coefficient_data = []
        sig_nodes = []
        for i in coeffic_data:
            coefficient_data.append(i[4:])
            sig_nodes.append(i[3])
        s_data = synch_data.values.tolist()
        synch_data = []
        conn_data = []
        average = []
        total = []
        phases = []
        per_nodes = []
        ind_connections = []
        channels = ['ch12', 'ch13', 'ch21', 'ch22', 'ch23', 'ch24', 'ch31', 'ch32', 'ch33', 'ch34', 'ch42', 'ch43']
        #iterates by phase (i is the phase number)
        for i in range(len(s_data)):
            print("Phase:" + str(i))
            hub_nodes = {}
            d=s_data[i]
            #l is a list of dictionaries that shows the sychrony indexes between all channels. i.e. [{ch13: 0.1, ch21: .2, ...}, {ch12:0.1, ch21:.3, ...}, ...]
            l = []
            for j in range(11):
                l.append(str_to_dict(d[j+2]))
            c = coefficient_data[i]
            cluster_coeff = {}
            for j in range(len(c)):
                if c[j] == 1:
                    hub_nodes[channels[j]] = l[j]
            for x in range(len(c)):
                cluster_coeff[channels[x]] = c[x]
            threshold = 0
            for dic in l:
                threshold+= dic['avg']
            threshold = threshold/12
            for hub in hub_nodes:
                conn = []
                for chan in hub_nodes[hub]:
                    if hub_nodes[hub][chan] >= threshold:
                        conn.append(chan)
                if 'avg' in conn:
                    conn.remove('avg')
                hub_nodes[hub] = conn
            hubs = list(hub_nodes.keys())
            print("Per Nodes: " + str(len(hubs)))
            print(hubs)
            if len(hubs) > 1:
                networks = [hubs[0]]
                for hub in hubs[1:]:
                    new = True
                    for h in networks:
                        if hub not in hub_nodes[h]:
                            new = False
                    if new:
                        networks.append(hub)
                ind_con = len(networks)
            elif len(hubs) == 1:
                networks = [hubs[0]]
                ind_con = 1
            else:
                networks = []
                ind_con = 0
            print(networks)
            print("Ind Con: " + str(ind_con))

            ind_connections.append(ind_con)
            per_nodes.append(len(hubs))
            phases.append(i)

        pd.DataFrame({'Phases':phases, 'Independent Networks':ind_connections, 'Perfect Nodes':per_nodes}).to_csv(os.path.join(synch_cc_data, str(well) + '_Independent_Networks.csv'))