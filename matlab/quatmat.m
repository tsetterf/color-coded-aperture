function [ Q ] = quatmat( q )
%QUATMAT Creates the matrix needed for conventional matrix multiplication
% i.e. qA2C = qA2B * qB2C = quatmat(qA2B) * qB2C

Q = [  q(4)  -q(3)   q(2)  q(1);
       q(3)   q(4)  -q(1)  q(2);
      -q(2)   q(1)   q(4)  q(3);
      -q(1)  -q(2)  -q(3)  q(4)  ];

end

