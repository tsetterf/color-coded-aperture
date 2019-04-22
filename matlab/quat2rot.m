function R = quat2rot(q)
%QUAT2ROT Get a rotation matrix equivalent to the quaternion

% Note that Eigen in C++ will give you the transpose of this

    q = q/norm(q); % ensure that quaternion is normalized

    q1 = q(1);
    q2 = q(2);
    q3 = q(3);
    q4 = q(4);

    R(1,1) = q1^2 - q2^2 - q3^2 + q4^2;
    R(1,2) = 2*(q1*q2 + q3*q4);
    R(1,3) = 2*(q1*q3 - q2*q4);

    R(2,1) = 2*(q1*q2 - q3*q4);
    R(2,2) = -q1^2 + q2^2 - q3^2 + q4^2;
    R(2,3) = 2*(q2*q3 + q1*q4);

    R(3,1) = 2*(q1*q3 + q2*q4);
    R(3,2) = 2*(q2*q3 - q1*q4);
    R(3,3) = -q1^2 - q2^2 + q3^2 + q4^2;

end
