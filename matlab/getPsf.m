function [ psf, dx ] = getPsf( D, f, zO, zD, N, l, Wl, mask, printFigs )
%GETPSF Get the PSF at the focal plane array for a given mask
%   Inputs:
%       D           Aperture diameter [m]
%       f           Focal length [m]
%       zO          Object position (-ve) [m]
%       zD          Detector position (+ve) [m]
%       N           Beam sampling (power of 2 preferred) [px]
%       l           nLx1 Vector of wavelengths [m]
%       Wl          nLx1 Vector of weights for each wavelength [-]
%       mask        NxN Mask to use
%       printFigs   Whether or not to print progress figures (default 0)
%   Outputs:   
%       psf         NxN Point spread function [m]
%       dx          The distance between pixels of the psf [m]
%   Author: Timothy Setterfield (Timothy.P.Setterfield@jpl.nasa.gov)

% Check inputs
if sign(zO) == 1
    error('Object distance must be negative');
end
if sign(zD) == -1
    error('Detector distance must be positive');
end
if mod(log2(N),1) ~= 0
    disp('Warning: Powers of 2 are preferred for beam sampling');
end
if length(l) ~= length(Wl)
    error('Number of weights must equal the number of wavelengths');
end

% Define beam ratio
beamRatio = 0.1;

if printFigs
    figure(1); clf;
end

% Create the wavefront object and point spread function for each wavelength
psf = zeros(N,N);
for i = 1:length(l)
    
    tic;
    wfo = prop_begin(D,l(i),N,beamRatio); % define beam parameters
    wfo = prop_multiply(wfo,mask);        % apply mask at EP = Stop = XP
    wfo = prop_define_entrance(wfo);      % normalize
    wfo = prop_lens(wfo,zO);              % curve incoming wf from pt src
    wfo = prop_lens(wfo,f);               % curve outgoing wf using lens
    wfo = prop_propagate(wfo,zD);         % propagate wf to detector
    disp(['wfo formed with proper in ' num2str(toc) ...
                's, dx=' num2str(wfo.dx)]);
       
    % Add wavelength-weighted electrical field energy to point spread func
    psf = psf + Wl(i)*prop_get_amplitude(wfo).^2;
    
    if printFigs
        figure(1); subplot(ceil(sqrt(length(l))),ceil(sqrt(length(l))),i);
        imagesc(Wl(i)*prop_get_amplitude(wfo).^2); axis equal;
        title(['PSF \lambda = ' num2str(l(i)*1e9) ' nm']);
        drawnow;
    end
    
end

% Normalize to a total transmission of 1
psf = psf/sum(psf(:)); 

if printFigs
    figure;
    subplot(1,2,1); imagesc(mask); axis equal; title('Mask');   
    subplot(1,2,2); imagesc(psf); axis equal; 
    title('PSF (weighted sum of \lambda{s})');
    drawnow;
end   

% Get the distance between pixels of the psf
dx = wfo.dx;

end

