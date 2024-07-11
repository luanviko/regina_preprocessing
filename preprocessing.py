import uproot, sys, progressbar, warnings
from scipy.signal import find_peaks
import matplotlib.pyplot as plt
import numpy as np
from array import array

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

def get_run_number(root_file_path):
    # Split path of files named "root_run_XXXXXX.root" to find XXXXXX
    # @params: _root_file_path_ containing the Linux/Mac path to the root file
    # @return: string with XXXXXX
    split1 = root_file_path.split("_")
    run_number = split1[-1].replace(".root","")
    return run_number

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

def CFD_timing_extrapolation(waveform, pulse_positions, percentage, start_rise_time, end_rise_time):
    # Extrapolate the linear fit to extract the timing of the pulse at the baseline.
    # @parameters: waveform, position of pulses, percentage, begin and end of rise time.
    # @return    : the timing of each pulse
    CFD_samples = np.array([], dtype="int")
    for pulse in pulse_positions:
        # Find max amplitude, and start and end of rise-time amplitude
        ymax    = waveform[pulse]
        y_start = start_rise_time*ymax
        y_end   = end_rise_time*ymax
        # Find positions of start and end of rise time
        i_start = walk_backward(waveform, y_start, pulse)
        i_end   = walk_backward(waveform, y_end, pulse)
        # Make linear interpolation between points
        a, b    = linear_interpolation(i_start, y_start, i_end, y_end)
        # Extrapolate to the baseline
        rise_amplitude = (end_rise_time-start_rise_time)*ymax
        CFD_time = (percentage*ymax-a)/b 
        # If b be nan or inf, take i_start as timing
        if np.isnan(CFD_time) or np.isinf(CFD_time):
            CFD_samples = np.append(CFD_samples, [i_start], axis=None)
        # Otherwise, apply extrapolation through the baseline
        else:
            CFD_samples = np.append(CFD_samples, [CFD_time], axis=None)
    return CFD_samples


def by_threshold(waveform, pulses, threshold):
    by_threshold_timing = np.array([], dtype="int")
    by_threshold_value  = np.array([], dtype="float")
    for pulse in pulses:
        i_threshold = walk_backward(waveform, threshold, pulse)
        by_threshold_timing = np.append(by_threshold_timing, [i_threshold], axis=None)
        by_threshold_value = np.append(by_threshold_value, [waveform[i_threshold]], axis=None)
    return by_threshold_timing

def to_root(file_name, data):
    # Save data dictionary to a tree in a root file
    # @params: _file_name_ with path to root file and _data_ dict.
    print(f"Saving to: {file_name}.root")
    with uproot.recreate(file_name+".root") as output_file:
        output_file["pulse_information"] = data

def to_npz(file_name, data):
    # Save data dictionary to an npz file
    # @params: _file_name_ with path to npz file and _data_ dict.   
    print(f"Saving to: {file_name}")
    np.savez_compressed(file_name, **data)

def main():

    # Comment this filter if debugging
    warnings.simplefilter("ignore", category=RuntimeWarning)
    
    # Branches, channels and file to analyze
    branches = ["midas_data_D300", "midas_data_D301", "midas_data_D302"]
    channels = [f"Channel{i}" for i in range(0,8)]
    root_file_path = sys.argv[1]
    total_channels = int(len(branches)*len(channels))
    
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

    # Output fixed-size arrays to be filled with a single value per waveform
    baselines    = np.empty((entries[0], total_channels), dtype="float")
    event_number = np.empty((entries[0], total_channels), dtype="int")
    count        = np.empty((entries[0], total_channels), dtype="int")

    # Output variable-size arrays to be filled multiple values per waveform
    max_pulses = 30
    pulses     = np.zeros((entries[0], total_channels, max_pulses), dtype="int")
    charges    = np.zeros((entries[0], total_channels, max_pulses), dtype="float")
    amps       = np.zeros((entries[0], total_channels, max_pulses), dtype="float")
    CFD_timing = np.zeros((entries[0], total_channels, max_pulses), dtype="float")
    threshold_timings = np.zeros((entries[0], total_channels, max_pulses), dtype="float")

    # Start analysis
    for branch_index in np.arange(0, len(branches)):
        for channel_index in np.arange(len(channels)):
            with uproot.open(root_file_path) as root_file:
                waveforms = root_file[branches[branch_index]].arrays([channels[channel_index], "eventNumber"], library="np")
            
            # Loop over every channel in every branch
            channel_number = channel_index+8*branch_index
            widgets=[f'Run {run_number}. Channel {channel_number}. Preprocessing: ', progressbar.Percentage(), progressbar.Bar('\u2587', '|', '|'), ' ', progressbar.Timer()]
            
            bar = progressbar.ProgressBar(widgets=widgets, maxval=int(1.*entries[branch_index])).start()
            for k in range(0, int(1.*entries[branch_index])):
                waveform = waveforms[channels[channel_index]][k]
                bar.update(k+1)
                if len(waveform) > 0:
                    
                    # Apply baseline correction
                    baseline = waveform_baseline(waveform, -15)
                    waveform -= baseline
                    
                    # Find pulse information
                    pulse, _          = find_peaks(-1.*waveform, height=threshold)
                    charge            = pulse_charge(waveform, pulse, 2)
                    amp               = pulse_amplitude(waveform, pulse)
                    CFD_times         = CFD_timing_extrapolation(waveform, pulse, 0., 0.4, 0.8)
                    threshold_timing  = by_threshold(waveform, pulse, -15)

                    # # Save fixed-size information
                    baselines[k][channel_number] = baseline
                    count[k][channel_number]     = len(pulse)
                    event_number[k]              = waveforms["eventNumber"][k]
                    
                    # Save variable-size information
                    for l in range(0, len(amp)):
                        if len(amp) < max_pulses:
                            pulses[k][channel_number][l]            = pulse[l]
                            charges[k][channel_number][l]           = charge[l]
                            amps[k][channel_number][l]              = amp[l]
                            CFD_timing[k][channel_number][l]        = CFD_times[l]
                            threshold_timings[k][channel_number][l] = threshold_timing[l]


    data = {"baseline":baselines, "event_number":event_number, "count":count, "pulse":pulses, "charge":charges, "amplitude":amps, "CFD_timing":CFD_timing, "threshold_timing":threshold_timings}
    to_root(final_file, data)
    to_npz(final_file, data)

if __name__ == "__main__":
    main()
