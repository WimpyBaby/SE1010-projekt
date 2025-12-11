% Initialize global variables

global c A v p a1 a2 g m L R df db h h1 b1 dh rd rb bb bd D d a gamma rh

c = 0.3;
A = 0.45;
v = 120/3.6;
p = 1.2;
a1 = 4;
a2 = -10;
g = 9.81;
m = 120;
L = 1.1;
R = 14;
df = 0.8;
db = 0.3;
h = 0.44;
h1 = 0.2;
b1 = 0.16;
dh = 0.33;
rd = 0.11;
rb = 0.09;
bb = 0.28;
bd = 0.0;
D = 0.028;
d = 0.6*D;
a = 0;
gamma = 0.25;
rh = 0.5*dh;

%%
% Calculating Nb for Case 1 

Fl = 0.5*p*c*A*v^2;
Fd = Fl;
Nb = (Fd*h+Fl*h1+m*g*df)/(df+db);

disp("Nb: " + Nb)

%%
% Case 1 calculating the Fk value
Fl = 0.5*p*c*A*v^2;
Fd = Fl;
Fk = (Fd*rh)/(rd);

disp("Fk is: " + Fk + ", Fd is: " + Fd);

%%
% Case 1 xz-plane Correct
syms R1z R2z

b1 = 0.16;
L = 1.1;

kraft1xz = Nb + R1z + R2z == 0;
Moment1xz = -R1z.*b1 - R2z.*(L-b1) - (Nb.*L)./2 == 0;

sol1xz = solve([kraft1xz, Moment1xz], [R1z, R2z]);

disp("Case 1 xz-plane:")
disp("R1z: " + string(vpa(sol1xz.R1z)))
disp("R2z: " + string(vpa(sol1xz.R2z)))


%%
% Case 1 xy-plane 
syms R1y R2y

Fl = 0.5*p*c*A*v^2;
Fd = Fl;

kraft = Fd + R1y + R2y + Fk == 0;
moment = -R1y*b1 - Fk *(L/2) - R2y*(L-b1) - (Fd/2)*L == 0;

solve_Ry = solve([kraft, moment], [R1y, R2y]);

disp("Reaktion forces in Y direction");
disp("R1y: " + string(vpa(solve_Ry.R1y)));
disp("R2y: " + string(vpa(solve_Ry.R2y)));


%%
% Case 2 calculating the Fk value 

Fl = 0;
Fd = Fl + m*a1;
Fk = (Fd*rh)/(rd);

disp("Fk is: " + Fk + ", Fd is: " + Fd);

%%
% Calculating Nb for case 2

Fl = 0.5*p*c*A*a^2;
Fd = m*a1 + Fl;
Nb = (Fd*h+Fl*h1+m*g*df)/(df+db);

disp("Nb: " + Nb)

%%
% Case 2 xz-plane 

syms R1z R2z

kraft = Nb + R1z + R2z == 0;
moment = -R1z*b1 - R2z*(L-b1) - (Nb/2)*L == 0;

solve_Rz = solve([kraft, moment], [R1z, R2z]);

disp("Reaktion forces in z direction case 2");
disp("R1z: " + string(vpa(solve_Rz.R1z)));
disp("R2z: " + string(vpa(solve_Rz.R2z)));


%%
% Case 2 xy-plane

syms R1y R2y

kraft = Fd + R1y + Fk + R2y == 0;
moment = -R1y*b1 - Fk*(L/2) - R2y*(L-b1) - (Fd/2)*L == 0;

solve2_Ry = solve([kraft, moment], [R1y, R2y]);

disp("Reaktion forces in y direction case 2");
disp("R1y: " + string(vpa(solve2_Ry.R1y)));
disp("R2y: " + string(vpa(solve2_Ry.R2y)));

%%
% Case 3 calculating Fb

Fl = 0.5*p*c*A*v^2;
Fd = Fl + m*a2;
Fb = (Fd*rh)/rb;

disp("Fb is: " + Fb + ", Fd is: " + Fd);

%%
% Calculating Nb for case 3

Fl = 0.5*p*c*A*v^2;
Fd = m*a2 + Fl;
Nb = (Fd*h+Fl*h1+m*g*df)/(df+db);

disp("Nb: " + Nb)
%%
% Case 3 xz-plane

syms R1z R2z 

kraft = Nb + R1z + R2z - (1/sqrt(2))*Fb == 0;
moment = -R1z*b1 + (0.5/sqrt(2))*Fb*(0.5*L-bb) + (0.5/sqrt(2))*Fb*(0.5*L+bb) - R2z*(L-b1) - Nb*(L/2) == 0;

solve3_Rz = solve([kraft, moment], [R1z, R2z]);

disp("Reaktion forces in z direction case 3");
disp("R1z: " + string(vpa(solve3_Rz.R1z)));
disp("R2z: " + string(vpa(solve3_Rz.R2z)));

%%
% Case 3 xy-plane

syms R1y R2y 

kraft = Fd + R1y + R2y + (1/sqrt(2))*Fb == 0;
moment = -R1y*b1 - (1/(2*sqrt(2)))*Fb*(0.5*L-bb) - (1/(2*sqrt(2)))*Fb*(0.5*L+bb) - R2y*(L-b1) - Fd*(L/2) == 0;

solve3_Ry = solve([kraft, moment], [R1y, R2y]);

disp("Reaktion forces in z direction case 3");
disp("R1y: " + string(vpa(solve3_Ry.R1y)));
disp("R2y: " + string(vpa(solve3_Ry.R2y)));

%%
% Case 4 calculating Fd and Fk

% calculate Nb and Nf for constant speed
Fl = 0.5*p*c*A*(gamma*v)^2;
Fd = Fl;
Fk = (Fd*rh)/(rd);

disp("Fk is: " + Fk + ", Fd is: " + Fd);

%%
% Calculating V and H
Nb = (Fd*h+Fl*h1+m*g*df)/(df+db);
Nf = m*g-Nb;

syms Vfi Vfy Vby Vbi Hbi Hfi Hby Hfy

momentyz = Vfi*(df+db) + Vfy*(df+db) + Fl*(h1+h) - m*g*db == 0;
z = Vbi + Vby + Vfy + Vfi - m*g == 0;
x = ((m*(gamma*v)^2)/R) - Hbi - Hfi - Hby - Hfy == 0;
momentxz1 = 0.5*L*Vbi + 0.5*L*Vfi + (h*m*(gamma*v)^2)/R - 0.5*L*Vby - 0.5*L*Vfy == 0;
momentxz2 = (db*m*(gamma*v)^2)/R - Hfi*(df+db) - Hfy*(df+db) == 0;
samband1 = Vfi == (Nf/Nb)*Vbi;
samband2 = Vfy == (Nf/Nb)*Vby;
samband3 = Hfi == (Vfi/Vfy)*Hfy;
samband4 = Hbi == (Vbi/Vby)*Hby;

solve_VH = solve([momentyz, z, x, momentxz1, momentxz2, samband1, samband2, samband3, samband4], [Vfi, Vfy, Vby, Vbi, Hbi, Hfi, Hby, Hfy]);

% disp("Vbi:" + string(vpa(solve_VH.Vbi)) + " Vby:" + string(vpa(solve_VH.Vby)) + " Vfi:" + string(vpa(solve_VH.Vfi)) + " Vfy:" + string(vpa(solve_VH.Vfy)));
% 
% disp("Hfy:" + string(vpa(solve_VH.Hfy)) + " Hfi:" + string(vpa(solve_VH.Hfi)) + " Hby:" + string(vpa(solve_VH.Hby)) + " Hbi:" + string(vpa(solve_VH.Hbi)));

variableNames = fieldnames(solve_VH);

% Loop through each field and assign it to a variable in the workspace
for i = 1:length(variableNames)
    varName = variableNames{i};
    varValue = solve_VH.(varName);
    assignin('base', varName, varValue);
end

% Now all variables Vfi, Vfy, Vby, Vbi, Hbi, Hfi, Hby, Hfy are in the workspace
% You can display or use them directly
disp("Vbi: " + string(vpa(Vbi)) + " Vby: " + string(vpa(Vby)) + " Vfi: " + string(vpa(Vfi)) + " Vfy: " + string(vpa(Vfy)));
disp("Hfy: " + string(vpa(Hfy)) + " Hfi: " + string(vpa(Hfi)) + " Hby: " + string(vpa(Hby)) + " Hbi: " + string(vpa(Hbi)));


%%
% Case 4 xz-plane

syms R1z R2z Rx

kraft_v = Vbi + R1z + R2z + Vby == 0;
kraft_h = Rx - Hbi - Hby == 0;
moment = 0.5*L*Vbi + R1z*(0.5*L-b1) - R2z*(0.5*L-b1) - 0.5*L*Vby + Hbi*rh + Hby*rh == 0; 

solve_xz = solve([kraft_v, kraft_h, moment], [R1z, R2z, Rx]);

disp("R1z:" + string(vpa(solve_xz.R1z)) + " R2z:" + string(vpa(solve_xz.R2z)) + " Rx:" + string(vpa(solve_xz.Rx)));

%%
% Case 4 xy-plane

syms R1y R2y Rx

kraft_h = Rx - Hbi - Hby == 0;
kraft_v = Fd + R1y + R2y + Fk == 0;
moment = 0.5^2*Fd*L-R1y*(0.5*L-b1) + R2y*(0.5*L-b1) - 0.5^2*Fd*L == 0;

solve_xy = solve([kraft_v, kraft_h, moment], [R1y, R2y, Rx]);

disp("R1y:" + string(vpa(solve_xy.R1y)) + " R2y:" + string(vpa(solve_xy.R2y)) + " Rx:" + string(vpa(solve_xy.Rx)));

%% Max gamma

numer = g*L*R;
denom = 2*h*v^2;

gamma_max = sqrt(numer/denom);

disp(gamma_max)
