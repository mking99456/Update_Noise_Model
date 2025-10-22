To look at ICEBERG Noise Data:
1. Set up a build of dunesw that can read raw digits from a .root file and can execute ROOT macros.
2. Run ICEBERG_Waveform_Creator.C on a decoded ICEBERG root data file to create a root file consisting of TH2D histograms. These histograms' bin content consists of an event ADC value as a function of wire number and time tick.
3. Warm_Noise_Analysis.ipynb looks at the baseline, rms, and correlations for a noise event by reading in a root file consisting of event data stored in a TH2D histogram.
