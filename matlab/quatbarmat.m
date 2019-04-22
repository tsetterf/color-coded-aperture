function [ Qbar ] = quatbarmat( q )
%QUATBARMAT Creates the matrix required for quaternion multiplication in the reverse 
% of conventional order (but the same order as rotation matrices)
% i.e. qAtoC = qAtoB * qBtoC = quatbarmat(qBtoC) * qAtoB

Qbar = [  q(4)   q(3)  -q(2)  q(1);
         -q(3)   q(4)   q(1)  q(2);
          q(2)  -q(1)   q(4)  q(3);
         -q(1)  -q(2)  -q(3)  q(4)  ];

end

