# regina_preprocessing

(Any typo is on purpose)

# CFD Timing Algorithm for Multi-Pulse Analysis

Basic nomenclature:
* Pulse position:
* Pulse amplitude:
* Rise time:
* Pulse timing:

Basic algorithm logic:
* Find the position and amplitude of a pulse.
* Determine rise time.
* Walk backward from the pulse position, and find the end of rise-time.
* Continue to walk backward and find the start of rise-time.
* Find linear and angular coefficients of the line between the start and end of rise-time.
* Extrapolate this line to the baseline.
* Save the point where the extrapolated line crosses the baseline as the timing of the pulse.

Explaining the plot:
* Blue crosses: Start and end of rise time being considered for timing determination.
* Red crosses: Timing of the pulse as the extrapolated line crosses the baseline.
