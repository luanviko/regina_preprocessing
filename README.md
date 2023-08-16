# regina_preprocessing

(Any typo is on purpose)

Basic nomenclature:
* Baseline: Average of N samples at the start or end of a waveform.
* Pulse position: Sample with minimum of pulse.
* Pulse amplitude: Distance from pulse minimum to baseline.
* Pulse charge: Sum of sample amplitude times sample sum over a width around pulse position.
* Rise time: Linear portion of pulse, the distance between 10% to 90% of pulse amplitude.
* Pulse timing: Sample where the pulse starts (CFD algorithm) or reaches max. amplitude (STT algorithm). 

# CFD Timing Algorithm for Multi-Pulse Analysis

Explaining the plot:
* Blue crosses: Start and end of rise time being considered for timing determination.
* Red crosses: Timing of the pulse as the extrapolated line crosses the baseline.

Basic algorithm:
* Find the position and amplitude of a pulse.
* Determine rise time.
* Walk backward from the pulse position, and find the end of rise-time.
* Continue to walk backward and find the start of rise-time.
* Find linear and angular coefficients of the line between the start and end of rise-time.
* Extrapolate this line to the baseline.
* Save the point where the extrapolated line crosses the baseline as the timing of the pulse.

