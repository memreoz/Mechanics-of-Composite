close all;clear all;clc
%% Composite material properties
% E-Glass/Epoxy

mat.E1 = 53.78E9; % [Pa]
mat.E2 = 17.93E9;
mat.v12 = 0.25;
mat.G12 = 8.62E9;
mat.v21 = (mat.E2/mat.E1)*mat.v12;

% Stress limits
lim.Xt = 1035e6;
lim.Xc = 1035e6;
lim.Yt = 27.6e6;
lim.Yc = 138e6;
lim.S = 41.4e6;

% Thermal expansion coefficients
a1 = 6.3e-6; %C
a2 = 20.52e-6;
% Loading conditions
Load = [1 0 0 0 0 0]';

% Laminate Geometrical Features
N = 3; % number of ply
theta = [0 90 0];
Ply_Thickness = [0.127e-3 1.27e-3 0.127e-3];
h = Ply_Thickness(1)+Ply_Thickness(2)+Ply_Thickness(3); %Total thickness
ml = h/2; %midline

z0 = -0.762e-3;
m1 = -0.6985e-3;
z1 = -0.635e-3;
m2 = 0;
z2 = 0.635e-3;
m3 = 0.6989e-3;
z3 = 0.762e-3;
z = [z0 z1 z2 z3];

%% Örnek 5.3 için TEST 
% Graphite/Epoxy

% mat.E1 = 181E9; % [Pa]
% mat.E2 = 10.30E9;
% mat.v12 = 0.28;
% mat.G12 = 7.17E9;
% mat.v21 = (mat.E2/mat.E1)*mat.v12;
% 
% % Stress limits
% lim.Xt = 1500e6;
% lim.Xc = 1500e6;
% lim.Yt = 40e6;
% lim.Yc = 246e6;
% lim.S = 68e6;
% 
% Load = [1 0 0 0 0 0]';
% 
% theta = [0 90 0];
% Ply_Thickness = [5e-3 5e-3 5e-3];
% h = Ply_Thickness(1)+Ply_Thickness(2)+Ply_Thickness(3); %Total thickness
% ml = h/2; %midline
% 
% z0 = -0.0075;
% m1 = -0.005;
% z1 = -0.0025;
% m2 = 0;
% z2 = 0.0025;
% m3 = 0.005;
% z3 = 0.0075;
% z = [z0 z1 z2 z3];

%% Reduced Stiffness Matrices
Q11 = mat.E1 / (1-mat.v12*mat.v21);
Q12 = (mat.v12*mat.E2) / (1-(mat.v12*mat.v21));
Q22 = mat.E2 / (1-mat.v12*mat.v21);
Q66 = mat.G12;

Qij1 = Qij(Q11, Q12, Q22, Q66, theta(1));
Qij2 = Qij(Q11, Q12, Q22, Q66, theta(2));
Qij3 = Qij(Q11, Q12, Q22, Q66, theta(3));

% Qij1 = Qij(0, 0, 0, 0, theta(1));
% Qij2 = Qij(0, 0, 0, 0, theta(2));
% Qij3 = Qij(0, 0, 0, 0, theta(3));

%% Transformation matrices
T1 = T(theta(1));
T2 = T(theta(2));
T3 = T(theta(3));

%% Stiffness matrix
Stiff = Stiff(Qij1, Qij2, Qij3, z);
Stiffinverse = inv(Stiff);
%%
sk = inv(Stiff)*Load;
% Global strains
ex = sk(1);
ey = sk(2);
gxy = sk(3);
eg = [ex ey gxy]';
% Curvatures kappa
kx = sk(4);
ky = sk(5);
kxy = sk(6);
k = [kx ky kxy]';
%% Global strains includes [ex ey gxy] for each location
% 1st layer global strain matrix
eg_1top = eg + z0*k;
eg_1mid = eg + m1*k;
eg_1bot = eg + z1*k;
% 2nd layer strain matrix
eg_2top = eg + z1*k;
eg_2mid = eg + m2*k;
eg_2bot = eg + z2*k;
% 3th layer strain matrix
eg_3top = eg + z2*k;
eg_3mid = eg + m3*k;
eg_3bot = eg + z3*k;
%% Reuter matrix
R = [1 0 0;0 1 0;0 0 2];

%% Local strains for ply layers outer plies(1,3) inner ply(2)
% 1st layer local strain matrix
el_1top = R * T1 * inv(R) * eg_1top;
el_1mid = R * T1 * inv(R) * eg_1mid;
el_1bot = R * T1 * inv(R) * eg_1bot;
% 2nd layer local strain matrix
el_2top = R * T2 * inv(R) * eg_2top;
el_2mid = R * T2 * inv(R) * eg_2mid;
el_2bot = R * T2 * inv(R) * eg_2bot;
% 3th layer local strain matrix
el_3top = R * T3 * inv(R) * eg_3top;
el_3mid = R * T3 * inv(R) * eg_3mid;
el_3bot = R * T3 * inv(R) * eg_3bot;

%% Global stresses includes [sx sy zxy]
% 1st layer global stress matrix
sg_1top = Qij1 * eg_1top;
sg_1mid = Qij1 * eg_1mid;
sg_1bot = Qij1 * eg_1bot;
% 2nd layer global stress matrix
sg_2top = Qij2 * eg_2top;
sg_2mid = Qij2 * eg_2mid;
sg_2bot = Qij2 * eg_2bot;
% 3th layer global stress matrix
sg_3top = Qij3 * eg_3top;
sg_3mid = Qij3 * eg_3mid;
sg_3bot = Qij3 * eg_3bot;

%% Local stresses includes [s1 s2 zxy]
% 1st layer global stress matrix
sl_1top = T1 * sg_1top;
sl_1mid = T1 * sg_1mid;
sl_1bot = T1 * sg_1bot;
% 2nd layer global stress matrix
sl_2top = T2 * sg_2top;
sl_2mid = T2 * sg_2mid;
sl_2bot = T2 * sg_2bot;
% 3th layer global stress matrix
sl_3top = T3 * sg_3top;
sl_3mid = T3 * sg_3mid;
sl_3bot = T3 * sg_3bot;

%% Load distribution
% Potion of load Nx taken by ply 1 : ex * t_k
Nx1 = sg_1mid(1) * Ply_Thickness(1);
% Potion of load Nx taken by ply 2
Nx2 = sg_2mid(1) * Ply_Thickness(2);
% Potion of load Nx taken by ply 3
Nx3 = sg_3mid(1) * Ply_Thickness(2);

% Percentage of load Nx taken by ply 1
Nxratio1 = (Nx1 / Load(1)) *100;
% Percentage of load Nx taken by ply 2
Nxratio2 = (Nx2 / Load(1)) *100;
% Percentage of load Nx taken by ply 3
Nxratio3 = (Nx3 / Load(1)) *100;

%% Table
pos = {'Top';'Middle';'Bottom';'Top';'Middle';'Bottom';'Top';'Middle';'Bottom'};
ply = {'1(0)';'2(90)';'3(0)';'1(0)';'2(90)';'3(0)';'1(0)';'2(90)';'3(0)'};
label_eg = {'Ply No','Position','ex','ey','gxy'};
label_el = {'Ply No','Position','e1','e2','g12'};
label_sg = {'Ply No','Position','sx','sy','zxy'};
label_sl = {'Ply No','Position','s1','s2','z12'};
eglobal = [eg_1top';eg_1mid';eg_1bot';eg_2top';eg_2mid';eg_2bot';eg_3top';eg_3mid';eg_3bot'];
elocal = [el_1top';el_1mid';el_1bot';el_2top';el_2mid';el_2bot';el_3top';el_3mid';el_3bot'];
sglobal = [sg_1top';sg_1mid';sg_1bot';sg_2top';sg_2mid';sg_2bot';sg_3top';sg_3mid';sg_3bot'];
slocal = [sl_1top';sl_1mid';sl_1bot';sl_2top';sl_2mid';sl_2bot';sl_3top';sl_3mid';sl_3bot'];

%% Failure Consideration
% Tsai-Wu
H1 = 1/lim.Xt - 1/lim.Xc;
H11 = 1/(lim.Xt*lim.Xc);
H2 = 1/lim.Yt - 1/lim.Yc;
H22 = 1/(lim.Yt*lim.Yc);
H6 = 0;
H66 = 1/lim.S^2;
% H12 = -1/(2*lim.Xt^2); %Tsai-Hill
% H12 = -1/(2*lim.Xt*lim.Xc);  %Hoffman
H12 = -(1/2)*sqrt(1/(lim.Xt*lim.Xc*lim.Yt*lim.Yc)); % Mises-Hencky


for i=1:9
    SR = lim.Xt/1;
    %Wu(i) = H1*slocal(i,1)*SR+H2*slocal(i,2)*SR+H6*slocal(i,3)+H11*slocal(i,1)^2*SR^2+H22*slocal(i,2)^2*SR^2+H66*slocal(i,3)^2+2*H12*slocal(i,1)*slocal(i,2)*SR^2;
    Wu(i) = H1*slocal(i,1)+H2*slocal(i,2)+H6*slocal(i,3)+H11*slocal(i,1)^2+H22*slocal(i,2)^2+H66*slocal(i,3)^2+2*H12*slocal(i,1)*slocal(i,2);
    
%     s1 = slocal(i,1);
%     s2 = slocal(i,2);
%     s3 = slocal(i,3);
%     
%     Wu(i) = H1*s1 + H2*s2 + H6*s3 + H11*s1^2 + H22*s2^2 + H66*s3^2 + 2*H12*s1*s2;
    %Tsai-Hill
    Hill(i) = (slocal(i,1)/lim.Xt)^2 - ((slocal(i,1)*slocal(i,2))/lim.Xt^2) + (slocal(i,2)/lim.Yt)^2 + (slocal(i,3)/lim.S)^2  ;
    
end
Wu = Wu';
Hill = Hill';

%mid-ply failure
% 
% Qij2 = [  0,0,0;
%           0,0,0;
%           0,0,0];
% Qij2 = Qij(0, 0, 0, 0, theta(2));





