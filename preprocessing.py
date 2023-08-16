import matplotlib.pyplot as plt
import numpy as np
import uproot
from array import array
from scipy.signal import find_peaks
import sys

def waveform_baseline(waveform, N):
    # Find baseline of waveform using the first N samples.
    # @parameters: waveform and N
    # @return    : average of the first N samples
    if N>0: return np.average(waveform[0:N+1])
    if N<0: return np.average(waveform[N:])

def pulse_charge(waveform, centers, width):
    # Find area (charge) in region center - width to center + width.
    # @parameters: waveform, baseline, center and width
    # @return    : area
    charges = np.array([], dtype="float")
    for center in centers:
        charge = 0.
        kmin = center - width 
        kmax = center + width 
        if center - width < 0 : kmin = 0
        if center + width > len(waveform) : kmax = len(waveform)
        for k in range(kmin, kmax): charge += waveform[k]
        charges = np.append(charges, [charge], axis=None)
    return charges

def pulse_amplitude(waveform, pulses):
    # Find the amplitude of each entry in pulses
    # @params: _waveform_ values and positions of minima in _pulses_
    # @return: numpy array with _waveform_[_pulses_[:]]
    amps = np.array([], dtype="float")
    for pulse in pulses: amps = np.append(amps, [waveform[pulse]], axis=None)
    return amps
    # Future: return np.array([waveform[pulse] for pulse in pulses])

def find_pulses(waveform, threshold):
    # PMT Test Facilities (PTF) algorithm to find pulses above _threshold_ in a waveform
    # @params: _waveform_ values and _threhsold_ value
    # @return: position of pulses in the waveform
    ymin = 99999999.
    imin = 0
    inside_pulse = False 
    pulses = np.array([], dtype="int")
    for k in range(0, len(waveform)):
        if (waveform[k] < ymin):
            ymin = waveform[k]
            imin = k 
        if ((waveform[k] < threshold) and (not inside_pulse) ):
            inside_pulse = True
        if ((waveform[k] >= threshold) and inside_pulse):
            inside_pulse = False
            pulses = np.append(pulses, [imin], axis=None)
            imin = 0
            ymin = 99999999.
    return pulses

def get_run_number(root_file_path):
    # Split path of files named "root_run_XXXXXX.root" to find XXXXXX
    # @params: _root_file_path_ containing the Linux/Mac path to the root file
    # @return: string with XXXXXX
    split1 = root_file_path.split("_")
    run_number = split1[-1].replace(".root","")
    return run_number

def walk_forward(waveform, reference):
    # Walk forward on the _waveform_ samples to find sample _j_ below _reference_ value
    # @params: _waveform_ values and _reference_ value
    # @return: _j_ or 10000 if exception occur
    try:
        sample_value = waveform[0]
        j = 0
        while ((sample_value > reference) and (j < len(waveform)-1)):
            j += 1
            sample_value = waveform[j]
        return j-1
    except IndexError:
        return 10000

def walk_backward(waveform, reference, start_point):
    # Walk forward on the _waveform_ samples to find sample _j_ below _reference_ value
    # @params: _waveform_ values, _reference_ threshold value and _start_point_
    # @return: _j+1_ from where the waveform crosses _reference_
    j = start_point
    while ((waveform[j] < reference) and (j > 0)):
        j -= 1
    return j+1

def linear_interpolation(x0, y0, x1, y1):
    # Find the linear _a_ and angular _b_ parameters
    # @params: Initial (_x0_,_y0) and final (_x1_, _y1_) coordinates
    # @return: _a_ and _b_
    b = (y1-y0)/(x1-x0)
    a = y0-x0*b
    return a, b

def CFD_timing_extrapolation(waveform, pulse_positions, percentage, start_rise_time=0.1, end_rise_time=0.9):
    
    # Find pulse timing as a _percentage_ of rise-time amplitude. 
    # @params: _waveform_, _pulse_positions_, _percentage_ or rise-time amplitude
    # @params: optional _start_rise_time_ and _end_rise_time_ as percentages of pulse amplitude
    # @return: numpy array with

    CFD_samples = np.array([], dtype="int")
    
    for pulse in pulse_positions:

        # Find max amplitude, and start and end of rise-time amplitude
        ymax    = waveform[pulse]
        y_start = start_rise_time*ymax
        y_end   = end_rise_time*ymax

        # Find positions of start and end of rise time
        i_start = walk_backward(waveform, y_start, pulse)
        i_end   = walk_backward(waveform, y_end, pulse)

        # Make linear interpolation between 
        a, b    = linear_interpolation(i_start, y_start, i_end, y_end)

        # If b be nan or inf, take i_start as timing
        if np.isinf(b) or np.isnan(b):
            CFD_samples = np.append(CFD_samples, [i_start], axis=None)

        # Otherwise, apply extrapolation through the baseline
        else: 
            rise_amplitude = (end_rise_time-start_rise_time)*ymax
            CFD_time = (percentage*rise_amplitude-a)/b
            CFD_samples = np.append(CFD_samples, [CFD_time], axis=None)

    return CFD_samples

def main():
    
    # Branches, channels and file to analyze
    branches = ["midas_data_D300", "midas_data_D301", "midas_data_D302"]
    channels = [f"Channel{i}" for i in range(0,8)]
    root_file_path = sys.argv[1]
    
    # Output npz file
    run_number = get_run_number(root_file_path)
    output_folder_path = "."
    final_file = f"{output_folder_path}/pulse_information-{run_number}"

    # Set threshold
    V_per_unit = 1./16384.*2.*1.e3
    threshold = 20. #mV
    threshold = threshold/V_per_unit

    # Get number of entries 
    entries = []
    for branch_index in np.arange(0, len(branches)):
        with uproot.open(root_file_path) as root_file:
            entries.append(root_file[branches[branch_index]].num_entries)

    # Output arrays
    max_pulses = 100
    baselines    = np.empty((entries[0], 32), dtype="float")
    event_number = np.empty((entries[0]), dtype="int")    
    count        = np.empty((entries[0], 32), dtype="int")
    pulses       = np.empty((entries[0], 32, max_pulses), dtype="int")
    charges      = np.empty((entries[0], 32, max_pulses), dtype="float")
    amps         = np.empty((entries[0], 32, max_pulses), dtype="float")
    CFD_timing   = np.empty((entries[0], 32, max_pulses), dtype="float")
    pulses_scipy  = np.empty((entries[0], 32, max_pulses), dtype="int")
    charges_scipy = np.empty((entries[0], 32, max_pulses), dtype="float")
    amps_scipy    = np.empty((entries[0], 32, max_pulses), dtype="float")
    count_scipy   = np.empty((entries[0], 32), dtype="int")

    # Start analysis
    print("Analyzing run: ", run_number)
    for branch_index in np.arange(0, len(branches)):
        for channel_index in np.arange(len(channels)):
            with uproot.open(root_file_path) as root_file:
                waveforms = root_file[branches[branch_index]].arrays([channels[channel_index], "eventNumber"], library="np")
            
            channel_number = channel_index+8*branch_index
            print(f"Analyzing channel {channel_number}:", end="\r")
            for k in range(0, int(0.001*entries[branch_index])):
                waveform = waveforms[channels[channel_index]][k]
                if len(waveform) > 0:
                    
                    baseline = waveform_baseline(waveform, -20)
                    waveform -= baseline
                    pulse  = find_pulses(waveform, -1.*threshold)
                    charge = pulse_charge(waveform, pulse, 5)
                    amp    = pulse_amplitude(waveform, pulse)
                    
                    pulse_sci, _ = find_peaks(-1.*waveform, height=threshold)
                    charge_sci   = pulse_charge(waveform, pulse_sci, 5)
                    amp_sci      = pulse_amplitude(waveform, pulse_sci)

                    CFD_times = CFD_timing_extrapolation(waveform, pulse_sci, 0., 0.4, 0.8)

                    # print("\n")
                    # print(k) 
                    # print(pulse_sci) 
                    # print(CFD_times)

                    baselines[k][channel_number] = baseline
                    count[k][channel_number] = len(pulse)
                    event_number[k] = waveforms["eventNumber"][k]
                    
                    if len(pulse) > max_pulses: to_store = max_pulses
                    if len(pulse) <= max_pulses: to_store = len(pulse)
                    for l in range(0, to_store): 
                        # print("PTF: ", len(pulse), k, l, to_store, float(charge[l]))
                        pulses[k][channel_number][l] = int(pulse[l])
                        charges[k][channel_number][l] = float(charge[l])
                        amps[k][channel_number][l] = float(amp[l])

                    count_scipy[k][channel_number] = len(pulse_sci)
                    if len(pulse_sci) > max_pulses: to_store = max_pulses
                    if len(pulse_sci) <= max_pulses: to_store = len(pulse_sci)
                    for l in range(0, to_store): 
                        # print("Scipy: ", len(pulse_sci), k, l, to_store, float(charge_sci[l]))
                        pulses_scipy[k][channel_number][l] = int(pulse_sci[l])
                        charges_scipy[k][channel_number][l] = float(charge_sci[l])
                        amps_scipy[k][channel_number][l] = float(amp_sci[l])
                        CFD_timing[k][channel_number][l] = float(CFD_times[l])

                    if k%100==0:
                        print(f"Analyzing {channel_number}: ", k, "/", len(waveforms[channels[channel_index]]), end="\r")
                    else:
                        if (k==len(waveforms[channels[channel_index]])-1):
                            print(f"Analyzing {channel_number}: ", k+1, "/", len(waveforms[channels[channel_index]]), end="\n")
    
    # Output data to compressed files
    print("Saving data to numpy file...")
    np.savez_compressed(f"{final_file}", baseline=baselines, event_number=event_number, count=count, pulse=pulses, charge=charges, amplitude=amps)
    np.savez_compressed(f"{final_file}_scipy", baseline=baselines, event_number=event_number, count_scipy=count_scipy, pulse_scipy=pulses_scipy, CFD_timing=CFD_timing, charge_scipy=charges_scipy, amp_scipy=amps_scipy)
    print(f"Pulse information saved to: {final_file}.npz")

if __name__ == "__main__":
    main()
