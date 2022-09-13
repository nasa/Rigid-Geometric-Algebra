%% Test Vectors
syms a0 a1 a2 a3 a23 a31 a12 a321 b0 b1 b2 b3 b23 b31 b12 b321 real
%a = ga3([a0 a1 a2 a3 a23 a31 a12 a321])
%b = ga3([b0 b1 b2 b3 b23 b31 b12 b321])
syms a4 a43 a42 a41 a412 a431 a423 a1234 b4 b43 b42 b41 b412 b431 b423 b1234 real
a = rga([a0 a1 a2 a3 a4 a23 a31 a12 a43 a42 a41 a321 a412 a431 a423 a1234]);
b = rga([b0 b1 b2 b3 b4 b23 b31 b12 b43 b42 b41 b321 b412 b431 b423 b1234]);
%x = ga3(randn(8,1));
%y = ga3(randn(8,1));
x = rga(randn(16,1));
y = rga(randn(16,1));
%[e0, e1, e2, e3, e23, e31, e12, e321] = ga3bases
[e0,e1,e2,e3,e4,e23,e31,e12,e43,e42,e41,e321,e412,e431,e423,e1234] = rgabases;
m1 = e0 + e1 + e2 + e3 + 2*e23 + 2*e31 + 2*e12 + e321;