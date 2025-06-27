class Analyzer:
    #Global Field Variables used in every method of this class:
    #For SBND specifically
    numchannels = 1280 # 11265 #Number of wires in the TPC
    N = 2128 # 3400 #Minimum number of time ticks present in SBND data arrays we're analyzing.
    PSDlength = 1065 # 1701 #the length of the PSD of the minwvfm
    Seconds_Per_Timetick = 0.5e-6 #0.5 microseconds per tick
    #Path_to_Validation = '/exp/dune/data/users/mking/ICEBERG_Noise_Ar_39/iceberg_noise/Sim_Validation_04192024/'
    Path_to_Validation = '/exp/dune/data/users/mking/ICEBERG_Noise_Ar_39/iceberg_noise/Update_Noise_Model_09292024/'


    def LoadWaveformFile(Waveform_File_Path='Waveforms.root',returnIndex=0):
        import uproot
        import awkward as awkward
        import numpy as np

        print ("Opening file: ",Waveform_File_Path)
        infile = uproot.open(Waveform_File_Path)
        key = infile.keys()[returnIndex] #only take daq channel
        H, xedges, yedges = infile[key].to_numpy()
        xcenters = Analyzer.EdgesToCenters(xedges)
        timeticks = yedges[:-1]

        return H, xedges, yedges, xcenters, timeticks, key

    def EdgesToCenters(edges_array):
        return (edges_array[:-1] + edges_array[1:]) / 2

    def CentersToEdges(centers_array,mode=0):
        import numpy as np

        if (mode==0):
            edges_array = np.zeros(len(centers_array)+1)
            bin_radius = (centers_array[1]-centers_array[0])/2
            edges_array[:-1] = centers_array - bin_radius*np.ones(len(centers_array))
            edges_array[-1] = edges_array[-2]+bin_radius*2
            return edges_array

        if (mode==1):
            return np.append(centers_array,centers_array[-1]+centers_array[1])

    ##Method takes 20-30 minutes to run for 100 events
    ##Method takes in a waveform histogram and returns the FFT for that event
    def FFTEvent(H,key):
        import numpy as np
        from pathlib import Path

        FFTEvent = np.zeros((Analyzer.numchannels,Analyzer.PSDlength),dtype=float)

        StartingChannelIndex = 0

        #If in-progress file exists, read from that.
        my_file = Path(Analyzer.Path_to_Validation+'FFTs/'+key+'_FFT_In_Progress.npz')
        if my_file.is_file():
            ExistingFFTEvent = np.load(Analyzer.Path_to_Validation+'FFTs/'+key+'_FFT_In_Progress.npz')
            FFTEvent = ExistingFFTEvent['FFT']
            StartingChannelIndex=ExistingFFTEvent['num']

        if StartingChannelIndex == Analyzer.numchannels:
            return FFTEvent

        for iChannel in np.arange(StartingChannelIndex,np.shape(H)[0],1):
            FFTEvent[iChannel] = Analyzer.getFFT(H[iChannel])
            if iChannel%200 == 0:
                np.savez(Analyzer.Path_to_Validation+'FFTs/'+key+'_FFT_In_Progress.npz',FFT = FFTEvent,num=iChannel)
        
        np.savez(Analyzer.Path_to_Validation+'FFTs/'+key+'_FFT_In_Progress.npz',FFT = FFTEvent,num=Analyzer.numchannels)

        return FFTEvent
    
    def FFTEventFast(H,key):
        import numpy as np
        from pathlib import Path

        my_file = Path(Analyzer.Path_to_Validation+'FFTs/'+key+'_FFT_In_Progress.npz')
        if my_file.is_file():
            ExistingFFTEvent = np.load(Analyzer.Path_to_Validation+'FFTs/'+key+'_FFT_In_Progress.npz')
            FFTEvent = ExistingFFTEvent['FFT']
            StartingChannelIndex=ExistingFFTEvent['num']
            if StartingChannelIndex == Analyzer.numchannels:
                return FFTEvent

        FFT = np.abs(np.fft.rfft(H)) #Take the absolute value of the FFT to get the power spectrum
        np.savez(Analyzer.Path_to_Validation+'FFTs/'+key+'_FFT_In_Progress.npz',FFT = FFT,num=Analyzer.numchannels)
        
        return FFT

    def AverageFFT(OutputFFT_File_Path = 'AverageFFT.npz', Waveform_File_Path='waveform.root',nEvents = 100):
        import uproot
        import awkward as awkward
        import numpy as np

        AverageFFT = np.zeros((Analyzer.numchannels,Analyzer.PSDlength),dtype=float)
        Events_Used_Count=0

        for iEvent in range(nEvents):
            print(iEvent)

            #CosmicEvents = [2,4,6,14,21,26,29,34,35,39,43,44,47,59]
            #If iEvent+1 is equal to an element in CosmicEvents, skip that event
            #if iEvent+1 in CosmicEvents:
            #    continue

            Output = Analyzer.LoadWaveformFile(Waveform_File_Path=Waveform_File_Path,returnIndex = iEvent)
            H = Output[0]
            key = Output[5]
            FFT = Analyzer.FFTEventFast(H,key)
            ## Zero out the first bin of the FFT for each channel to remove the DC component
            FFT[:,0] = np.zeros(Analyzer.numchannels)
            AverageFFT += FFT
            Events_Used_Count += 1
            del Output, FFT
        print(Events_Used_Count)
        AverageFFT /= Events_Used_Count
        np.savez(OutputFFT_File_Path,FFT = AverageFFT,nEvents=Events_Used_Count)
        return AverageFFT
    
    def UpdateAverageFFT(AverageFFT_File_Path = 'AverageFFT.npz', OutputFFT_File_Path = 'AverageFFT.npz', Waveform_File_Path='waveform.root',nStart = 0,nEvents = 100):
        import numpy as np

        SumFFT = np.load(AverageFFT_File_Path)['nEvents']*np.load(AverageFFT_File_Path)['FFT']

        for iEvent in np.arange(nStart,nEvents):
            print(iEvent)
            Output = Analyzer.LoadWaveformFile(Waveform_File_Path=Waveform_File_Path,returnIndex = iEvent)
            H = Output[0]
            key = Output[5]
            FFT = Analyzer.FFTEventFast(H,key)
            SumFFT += FFT
            del Output, H, key, FFT
        AverageFFT = SumFFT/nEvents

        np.savez(OutputFFT_File_Path,FFT = AverageFFT,nEvents=nEvents)
        return AverageFFT
    
    def AverageFromFFTFiles(directory_in_str = 'FFTs/'):
        import os
        import numpy as np

        AverageFFT = np.zeros((Analyzer.numchannels,Analyzer.PSDlength),dtype=float)
        Count = 0

        directory = os.fsencode(directory_in_str)
    
        for file in os.listdir(directory):
            filename = os.fsdecode(file)
            if filename.endswith(".npz"):
                FFTFile = np.load(filename)
                if FFTFile['num'] == Analyzer.numchannels:
                        AverageFFT += FFTFile['FFT']
                        Count += 1
            else:
                continue

    def Plot2DWaveform(H,xedges,yedges):
        import matplotlib.pyplot as plt

        fig2,ax2 = plt.subplots(1,figsize=(20,10))
        pc2 = ax2.pcolorfast(xedges,yedges,H.T)

        plt.show()

    def getFFT(Waveform):
        import numpy as np

        ADCFFT = np.abs(np.fft.rfft(Waveform))
        ADCFFT[0]=0

        return ADCFFT
    
    def getMedian(Function,Window_Size = 30):
        import numpy as np

        if len(Function)< 2*Window_Size:
            return np.median(Function)*np.ones(len(Function))
        median = Function.copy()

        for i in range(Window_Size):
            median[i] = np.median(Function[0:i+Window_Size])
        for i in range(Window_Size,len(Function)-Window_Size):
            median[i] = np.median(Function[i-Window_Size:i+Window_Size])
        for i in range(len(Function)-Window_Size,len(Function)):
            median[i] = np.median(Function[i-Window_Size:len(Function)-1])

        return median

    def Plot2DFFT(FFT,vmin = -1, vmax = 7,log=True):
        import numpy as np
        import matplotlib.pyplot as plt
        import matplotlib as mpl

        debug = True

        freq = np.fft.rfftfreq(Analyzer.N,Analyzer.Seconds_Per_Timetick)
        freq_edges = Analyzer.CentersToEdges(freq,mode=1)
        channel_edges = Analyzer.CentersToEdges(np.arange(Analyzer.numchannels))

        fig,ax = plt.subplots(1,figsize=(20,10))
        #pc = ax.pcolorfast(channel_edges,freq_edges,FFT.T)

        cmap = mpl.cm.get_cmap('viridis')
        cmap.set_under('white')
        cmap.set_over('black')

        #debug
        if debug:
            print(np.shape(channel_edges))
            print(np.shape(freq_edges))
            print(np.shape(FFT))

        if log:
            FFT = np.log10(FFT)

        ax.pcolormesh(channel_edges,freq_edges,FFT.T,cmap = cmap,shading='flat',vmin = vmin, vmax = vmax) #default -1,7
        fig.colorbar(ax.pcolormesh(channel_edges,freq_edges,FFT.T,cmap = cmap,shading='flat',vmin = vmin, vmax = vmax) ) #default -1,7

        plt.xlabel('Channel')
        plt.ylabel('Frequency [Hz]')

        plt.show()
    


