clc
clear all

% Cayley Transform
% C = (I-Q)/(I+Q)
% Q = (I-C)/(I+C)

% Q -> skew-symmetric
% Q = [0, -q(3), q(2);
%     q(3), 0, -q(1);
%     -q(2), q(1), 0]

% C -> proper orthogonal matrix
% I -> identity NxN

%% N = 3
I = eye(3);
C = [0.813797, 0.296198, -0.5;
    0.235888, 0.617945, 0.75;
    0.531121, -0.728292, 0.433012];

Q3 = (I-C)/(I+C)

% CRPs
q3 = [Q3(3,2); Q3(1,3); Q3(2,1)]

%% N = 4
I = eye(4);
D = [0.505111, -0.503201, -0.215658, 0.667191;
    0.563106, -0.034033, -0.538395, -0.626006;
    0.560111, 0.748062, 0.272979, 0.228387
    -0.337714, 0.431315, -0.767532, 0.332884];

Q4 = (I-D)/(I+D)