%INITPROPER Initializes the proper library
% Author: Timothy Setterfield (Timothy.P.Setterfield@jpl.nasa.gov)

addpath('./thirdparty/proper/');

% Initialize proper
propcommon

% Set default values of input parameters
do_table              = 0;    % do not print table
print_it              = 0;    % do not print intermediate msgs, surface labels
print_total_intensity = 0;    % do not print total intensity
prop_phase_offset     = 0;    % do not apply phase offset
prop_verbose          = 0;    % do not print informational messages


