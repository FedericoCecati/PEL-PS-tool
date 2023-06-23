clear all; clc;

syms Kp_PLL Ki_PLL Kp Ki Kp_DC Ki_DC R_c C_c Kp_AC L_g R_g v_dc_ref v_ac_ref wcc C_dc L_f real;
syms Rl_12 Ll_12 Rl_23 Ll_23 R_g L_g  Lf      real
parametersSYM = [ Kp_PLL Ki_PLL Kp Ki Kp_DC Ki_DC R_c C_c Kp_AC L_g R_g v_dc_ref v_ac_ref wcc C_dc L_f ];

v_dc_ref = 1100;
v_ac_ref = 563;

omega  = 2*pi*50;           % reference omega, f = 50 Hz
jw =    [0      -omega;     % cross coupiling dq matrix
         omega   0 ];       % cross coupiling dq matrix


%% Definition of the variables as symbolic

% State Variables of VSC1: x1
i_c1  = sym('i_c1',[2, 1],'real');
v_dc1 = sym('v_dc1','real');
phi_dc1 = sym('phi_dc1','real');
phi_ac1 = sym('phi_ac1','real');

delta1 = sym('delta1','real');
phi_q1 = sym('phi_q1','real');

v_cc1 = sym('v_cc1',[2, 1],'real');

x1 = [ i_c1 ; v_dc1 ; phi_dc1 ; delta1 ; phi_q1 ; v_cc1 ];

% State Variables of VSC2: x2
i_c2  = sym('i_c2',[2, 1],'real');
v_dc2 = sym('v_dc2','real');
phi_dc2 = sym('phi_dc2','real');
phi_ac2 = sym('phi_ac2','real');

delta2 = sym('delta2','real');
phi_q2 = sym('phi_q2','real');

v_cc2 = sym('v_cc2',[2, 1],'real');

x2 = [ i_c2 ; v_dc2 ; phi_dc2 ; delta2 ; phi_q2 ; v_cc2 ];

% State Variables of VSC3: x3
i_c3  = sym('i_c3',[2, 1],'real');
v_dc3 = sym('v_dc3','real');
phi_dc3 = sym('phi_dc3','real');
phi_ac3 = sym('phi_ac3','real');

delta3 = sym('delta3','real');
phi_q3 = sym('phi_q3','real');

v_cc3 = sym('v_cc3',[2, 1],'real');

x3 = [ i_c3 ; v_dc3 ; phi_dc3 ; delta3 ; phi_q3 ; v_cc3 ];

% State Variables: Power Network
i_g1DQ   = sym('i_g1DQ',[2, 1],'real');
i_g2DQ   = sym('i_g2DQ',[2, 1],'real');
i_g3DQ   = sym('i_g3DQ',[2, 1],'real');
il_12DQ  = sym('il_12DQ',[2, 1],'real');
il_23DQ  = sym('il_23DQ',[2, 1],'real');
i_pccDQ  = sym('i_pccDQ',[2, 1],'real');

xPN = [ i_pccDQ ; il_12DQ ; il_23DQ ];

% State vector: x
x = [ x1 ; x2 ; x3; xPN ];
% Disturbance vector: d
Power1 = sym('Power1','real');
Power2 = sym('Power2','real');
Power3 = sym('Power3','real');
%
e = sym('e',[2, 1],'real');
d = [Power1 ; Power2 ; Power3 ; e ];

%% Power network nodal equations, according to the Kirkhoff laws

i_g1DQ = i_pccDQ - il_12DQ;
i_g2DQ = il_12DQ - il_23DQ;
i_g3DQ = il_23DQ;

%% Equations Converter 1

iq_ref1 = 0;
Tdelta1 = [cos(delta1) sin(delta1); -sin(delta1) cos(delta1)]; % from DQ to dq

v_g1 = Kp*(i_c1-Tdelta1*i_g1DQ) + v_cc1;

i_ref1 =  [  Kp_DC*(v_dc1 - v_dc_ref)  + Ki_DC*phi_dc1 ;
             ( iq_ref1 + Kp_AC*(v_ac_ref - sqrt(v_g1'*v_g1)) ) ];   %%%%%%%%  Voltage controller
        
wcc = Kp / Lf;
Tdelta1 = [cos(delta1) sin(delta1); -sin(delta1) cos(delta1)]; % from DQ to dq
        
di_c1    = -wcc *i_c1 + wcc *i_ref1;
dv_dc1   = -3/2 * 1/C_dc * v_g1'*Tdelta1*i_g1DQ/v_dc1   +   1/C_dc*Power1 / v_dc1;
dphi_dc1 = v_dc1 - v_dc_ref;
ddelta1  = [0 Kp_PLL]*v_g1   +   Ki_PLL*phi_q1;
dphi_q1  = [0 1] * v_g1;
dv_cc1   = -Ki*Tdelta1*i_g1DQ + Ki*i_c1 ;

v_g1DQ   = Tdelta1'*v_g1;

dx1 = [ di_c1 ; dv_dc1 ; dphi_dc1 ; ddelta1 ; dphi_q1 ; dv_cc1 ];

   
%% Equations Converter 2

iq_ref2 = 0;
Tdelta2 = [cos(delta2) sin(delta2); -sin(delta2) cos(delta2)]; % from DQ to dq

v_g2 = Kp*(i_c2-Tdelta2*i_g2DQ) + v_cc2;

i_ref2 =  [    Kp_DC*(v_dc2 - v_dc_ref)  + Ki_DC*phi_dc2 ;
               ( iq_ref2 + Kp_AC*(v_ac_ref - sqrt(v_g2'*v_g2)) ) ];   %%%%%%%% % Voltage controller


        
di_c2    = -wcc *i_c2 + wcc *i_ref2;
dv_dc2   = -3/2 * 1/C_dc * v_g2'*Tdelta2*i_g2DQ/v_dc2   +   1/C_dc*Power2 / v_dc2;
dphi_dc2 = v_dc2 - v_dc_ref;
ddelta2  = [0 Kp_PLL]*v_g2   +   Ki_PLL*phi_q2;
dphi_q2  = [0 1] * v_g2;
dv_cc2   = -Ki*Tdelta2*i_g2DQ + Ki*i_c2;

v_g2DQ = Tdelta2'*v_g2;

dx2 = [ di_c2 ; dv_dc2 ; dphi_dc2 ; ddelta2 ; dphi_q2 ; dv_cc2 ];

%% Equations Converter 3

iq_ref3 = 0;
Tdelta3 = [cos(delta3) sin(delta3); -sin(delta3) cos(delta3)]; % from DQ to dq

v_g3 = Kp*(i_c3-Tdelta3*i_g3DQ) + v_cc3;

i_ref3 =  [    Kp_DC*(v_dc3 - v_dc_ref)  + Ki_DC*phi_dc3 ;
               ( iq_ref3 + Kp_AC*(v_ac_ref - sqrt(v_g3'*v_g3)) )  ];   %%%%%%%% % Voltage controller
        
        
di_c3    = -wcc *i_c3 + wcc *i_ref3;
dv_dc3   = -3/2 * 1/C_dc * v_g3'*Tdelta3*i_g3DQ/v_dc3   +   1/C_dc*Power3 / v_dc3;
dphi_dc3 = v_dc3 - v_dc_ref;
ddelta3  = [0 Kp_PLL]*v_g3   +   Ki_PLL*phi_q3;
dphi_q3  = [0 1] * v_g3;
dv_cc3 = -Ki*Tdelta3*i_g3DQ + Ki*i_c3 ;

v_g3DQ = Tdelta3'*v_g3;


dx3 = [ di_c3 ; dv_dc3 ; dphi_dc3 ; ddelta3 ; dphi_q3 ; dv_cc3 ];
   
%% Power network and grid equations, according to the Kirkhoff laws

di_pccDQ  = -(R_g)/(L_g)*i_pccDQ  - jw*i_pccDQ  + 1/(L_g) * (v_g1DQ - e); % Grid 
dil_12DQ  = -(Rl_12)/(Ll_12)*il_12DQ - jw*il_12DQ + 1/(Ll_12) * (v_g2DQ - v_g1DQ); % Transmission line 1-2
dil_23DQ  = -(Rl_23)/(Ll_23)*il_23DQ - jw*il_23DQ + 1/(Ll_23) * (v_g3DQ - v_g2DQ); % Transmission line 2-3

dxPN = [ di_pccDQ  ; dil_12DQ ; dil_23DQ ];  % states concatenation


%% Parameters declaration

omega  = 2*pi*50;           % reference omega, f = 50 Hz
jw =    [0      -omega;     % cross coupiling dq matrix
         omega   0 ];       % cross coupiling dq matrix

fs = 2000;
Ts = 1/fs;

Power1 = 1e6;  % 1 MW
Power2 = 1e6;  % 1 MW
Power3 = 1e6;  % 1 MW
Power  = Power1 + Power2 + Power3;


v_dc_ref = 1100;
C_dc = 22e-3;

% Distribution lines %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ll_12 = 0.01 / (2*pi*50);    
Rl_12 = 2*pi*50 * Ll_12 * 2.5;

Ll_23 = 0.02 / (2*pi*50);   
Rl_23 = 2*pi*50 * Ll_23 * 1.5;

% PLL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PLL_D = 1/sqrt(2);
PLL_Tsettle = 0.2;  % original 0.2
Vamp = 563;

Kp_PLL = 2*4.6 /PLL_Tsettle;
Ti_PLL = PLL_Tsettle*(PLL_D^2) / 2.3;
Ki_PLL = Kp_PLL/Ti_PLL;


Kp_PLL = Kp_PLL / Vamp;
Ki_PLL = Ki_PLL / Vamp;

% Current Loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Lf = 0.1e-3;    % L Filter
Rf = 0.33*2*pi*50*Lf;

Kp  = Lf*1.8/(3*Ts);
Ti  = Lf/Rf;
Ki  = Kp/Ti;
T_c = Kp/Lf;

% Dc and AC voltage Loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Kp_DC = 0.12 * C_dc/Ts * 1.8;   % DC Voltage loop
Ti_DC = 17*Ts;                  % DC Voltage loop
Ki_DC = Kp_DC/Ti_DC / 4;        % DC Voltage loop

Kp_AC = -25; % AC Voltage droop


% Grid Impedance

SCR_des = 2.4; 
RXratio = 0.3;

Z_g = 690^2 / (Power * SCR_des);
X_g = sqrt(   Z_g^2 / ( RXratio^2 + 1 )   );

L_g = X_g / omega;
R_g = RXratio * X_g;


%% Linearized system

fs = [ dx1 ; dx2 ; dx3 ; dxPN ];
As = jacobian(fs,x); % linearization of the model

%% Static Power Flow analysis - Steady state computation

ys = [ v_g1(1) / v_ac_ref ; delta1 /(2*pi) * 180   ;   3/2* v_g1DQ'*i_g1DQ / 1e6  ;  3/2* v_g1DQ'*[0 -1; 1 0] * i_g1DQ  / 1e6 ;
       v_g2(1) / v_ac_ref ; delta2 /(2*pi) * 180   ;   3/2* v_g2DQ'*i_g2DQ / 1e6 ;  3/2* v_g2DQ'*[0 -1; 1 0] * i_g2DQ  / 1e6 ;
       v_g3(1) / v_ac_ref ; delta3 /(2*pi) * 180   ;   3/2* v_g3DQ'*i_g3DQ / 1e6 ;  3/2* v_g3DQ'*[0 -1; 1 0] * i_g3DQ   / 1e6  ];   %% THESE ARE THE VARIABLES DISPLAYED IN THE POWER FLOW ANALYSIS. EVERY VARIABLE OF INTEREST FOR THE POWER FLOW CAN BE ADDED: CURRENTS, POWER, DC VOLTAGES, ...

parametersNUM = eval(parametersSYM);

x01 =  [  Power1/v_ac_ref  0    1100    0    0    0    v_ac_ref   0  ]'; % flat start
x02 =  [  Power2/v_ac_ref  0    1100    0    0    0    v_ac_ref   0  ]'; % flat start
x03 =  [  Power3/v_ac_ref  0    1100    0    0    0    v_ac_ref   0  ]'; % flat start
x0l12 = [ 0  0 ]'; % flat start
x0l23 = [ 0  0 ]'; % flat start
x0g   = [ 0  0 ]'; % flat start
x0 = [ x01 ; x02 ; x03 ; x0l12 ; x0l23 ; x0g ]; % flat start

xe = x0;
de = [ Power1   ; Power2   ; Power3   ;  [ 563 0 ]'   ];  %%%%%%%%%%%%%%%
ue = [ 0 ; 0]; %%%%%%%%%%%%%%%

%% Newton Raphson iterative algorithm
disp('Executing Newton-Raphson algorithm for static analysis... ')
f_norm = 1e6;
tolerance = 1e-4;
iteration = 0;

fNR = subs(fs,d,de);

ANR = subs(As,d,de);

while(f_norm > tolerance)
    
    f = subs(fNR,x,xe);
    f = eval(f);
    A = subs(ANR,x,xe);
    A = eval(A);
tic
    xe = xe - A\f;
    f_norm = norm(f,1);
    iteration = iteration + 1;
toc    
end
disp('Finished! ')
disp(' ')
disp( [ 'Number of Newton-Raphson iterations: ' , num2str(iteration) ] );

y = subs(ys,x,xe);
y = eval(y);

disp(' ')
disp('################################')
disp('Results of the static power flow analysis')
disp('################################')
disp('Steady-state equilibrium point of each VSC, in the format: ')
disp(' --- ');
disp('|v| (p.u.) ')
disp('delta (deg) ')
disp('P (p.u.) ')
disp('Q (p.u.) ')
disp(' --- ');
disp(' ')
disp(' ')

format shortG
disp('Equilibrium point of VSC1')
disp(y(1:4)); % show the static equilibrium point
disp([ 'Power factor: ', num2str(y(3)/sqrt(y(3)^2 + y(4)^2)) ])
disp(' --- ');
disp('Equilibrium point of VSC2')
disp(y(5:8)); % show the static equilibrium point
disp([ 'Power factor: ', num2str(y(7)/sqrt(y(7)^2 + y(8)^2)) ])
disp(' --- ');
disp('Equilibrium point of VSC3')
disp(y(9:12)); % show the static equilibrium point
disp([ 'Power factor: ', num2str(y(11)/sqrt(y(11)^2 + y(12)^2)) ])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(364); plot(eig(A)/(2*pi),'*', 'MarkerSize',12,'Linewidth',2);  sgrid; hold on;   % eigenvalue analysis

disp(' ')
disp(' ')
disp('################################')
disp('Results of the small-signal dynamic eigenvalue analysis: Damping and natural frequency of the grid oscillatory modes')
disp('################################')

damp(eig(A))

disp(' ')
disp(' ')


