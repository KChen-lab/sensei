%% global parameters
rng(2020);

M1_list = 6:2:12;
M0_list = 6:2:12;

alpha = 0.05;

%% parameters for 1

m1 = 0.5;
v1 = 0.1;

v1 = v1 * v1;
a1 = ((1 - m1) / v1 - 1 / m1) * m1 * m1;
b1 = a1 * (1 / m1 - 1);

N1 = 1000;


%% parameters for 0

m0 = 0.3;
v0 = 0.15;

v0 = v0 * v0;
a0 = ((1 - m0) / v0 - 1 / m0) * m0 * m0;
b0 = a0 * (1 / m0 - 1);

N0 = 1000;


%% parameters for 1

m1 = 0.05;
v1 = 0.012;

v1 = v1 * v1;
a1 = ((1 - m1) / v1 - 1 / m1) * m1 * m1;
b1 = a1 * (1 / m1 - 1);

N1 = 1000;


%% parameters for 0

m0 = 0.03;
v0 = 0.015;

v0 = v0 * v0;
a0 = ((1 - m0) / v0 - 1 / m0) * m0 * m0;
b0 = a0 * (1 / m0 - 1);

N0 = 1000;
