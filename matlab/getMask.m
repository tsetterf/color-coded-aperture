function [ mask, maskObj ] = getMask( maskName, Nim, radFrac )
%GETMASK Gets the Nim x Nim mask with the specified name
%   Inputs:
%       maskName    Name of the mask to generate (options on each line
%                     do not overlap)
%                       'tridispred','tridispgreen','tridispblue'
%                       'circle'
%       Nim         Size of mask
%       radFrac     Optional argument for 'autodispgreen' mask with radius 
%                   of circle as a fraction of beam 
%   Outputs:
%       mask        Nim x Nim mask
%       maskObj     A mask object with properties:
%           color       Color of the mask
%           area        Area of the mask in units of radius squared
%           shapes      Cell array of shape objects describing the mask
%               type        Type of shape 'circle' or 'polygon'
%               transmit    Light is tranmitted 'inside' or 'outside' shape
%               data        Polygon points [[x1;y1] [x2;y2] ...] [R]
%           transmissivity  Fraction of light transmitted (between 0 and 1)
%   Author: Timothy Setterfield (Timothy.P.Setterfield@jpl.nasa.gov)

addpath('./thirdparty/proper/');

beamRatio = 0.1;                                % beam ratio [ - ]
D         = 0.05/1.4;                           % aperture diameter [m]
R         = D/2;                                % aperture radius [m]
wfo       = prop_begin(D,500e-9,Nim,beamRatio); % define wfo for sampling

maskObj.shapes{1}.type = 'circle';              % type 'circle'/'polygon'
maskObj.shapes{1}.transmit = 'inside';          % 'inside'/'outside'
maskObj.shapes{1}.data = getCircle(0, 0, 1);    % [x,y,r/R]

switch maskName
    case 'tridispred'
        tR = 0.7;
        mask = prop_ellipse(wfo,0.12,0.12, ...
                            'xc',+0.6495,'yc',-0.375,'norm',1) * tR;
    case 'tridispgreen'
        tG = 0.7;
        mask = prop_ellipse(wfo,0.12,0.12, ...
                            'xc',0,'yc',0.75,'norm',1) * tG;
    case 'tridispblue'
        tB = 0.7;
        mask = prop_ellipse(wfo,0.12,0.12, ...
                            'xc',-0.6495,'yc',-0.375,'norm',1) * tB;
    otherwise
        error(['Mask name ' maskName ' not supported']);
end

% Arrange all polygon data in a counterclockwise order for polybool func
for i = 1:length(maskObj.shapes)
    [maskObj.shapes{i}.data(1,:),maskObj.shapes{i}.data(2,:)] = ...
        poly2cw(maskObj.shapes{i}.data(1,:),maskObj.shapes{i}.data(2,:));
end


end

% Get the polygon points for a circle centered at (x,y), with radius r
function data = getCircle(x, y, r)

    theta     = linspace(0,2*pi,100);
    data(1,:) = x + r*cos(theta);       % x coordinates
    data(2,:) = y + r*sin(theta);       % y coordinates

end

