function [ qc ] = quatconj( q)
%QCONJ Returns the conjugate quaternion
    
    qc = [-q(1) -q(2) -q(3) q(4)]';

end

