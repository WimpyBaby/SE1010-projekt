% Kod för beräkning av NB som används senare

Nf = 1;
Nb = 2;
mg = 140*9.81;
Fl =0.5*0.3*0.5*1.21*27.7778.^2;
df = 0.8;
db = 0.3;
h1 = 0.2;
h = 0.44;


Nf = (mg.*db - Fl.*(h1+h))./(df+db);
Nb = mg - Nf;

%%
% Case 1 xz-plane
syms R1z R2z

b1 = 0.16;
L = 1.1;

kraft1xz = Nb + R1z + R2z == 0;
Moment1xz = -R1z.*b1 - R2z.*(L-b1) - (Nb.*L)./2 == 0;

sol1xz = solve([kraft1xz, Moment1xz], [R1z, R2z]);

vpa(sol1xz.R1z)
vpa(sol1xz.R2z)

%%
% Case 1 xy-plane
syms R1y R2y

Fd = 70;
Fk = 105;

kraft1xy = Fd - R1y -R2y + Fk == 0;
Moment1xy = R1y*b1 - 0.5*Fk*L + R2y*(L-b1) - 0.5*Fd*L == 0;

sol1xy = solve([kraft1xy, Moment1xy], [R1y, R2y]);

vpa(sol1xy.R1y)
vpa(sol1xy.R2y)


%%
% Case 1 yz-plane
syms R1y R2y

Fd = 70;
Fk = 105;
rd = 0.11;
rh = 0.5*0.33;

kraft1yz = Fd - R1y -R2y + Fk == 0;
Moment1yz = -Fd*rh - Fk*rd == 0;

sol1yz = solve([kraft1yz, Moment1yz], [R1y, R2y]);

vpa(sol1yz.R1y)
vpa(sol1yz.R2y)

%%
% Case 2 xz-plane
syms R1z R2z

b1 = 0.16;
L = 1.1;

kraft2xz = R1z + R2z + Nb == 0;
Moment2xz = -(R1z*b1 + R2z*(L-b1) + 0.5*Nb*L);

sol2xz = solve([kraft2xz, Moment2xz], [R1z, R2z]);

vpa(sol2xz.R1z)
vpa(sol2xz.R2z)

%%
% Case 2 xy-plane
syms R1y R2y

m = 140;
Fd = 630;
Fk = 945;
a = 4.5;
Fa =m*a;
b1 = 0.16;
L = 1.1;

kraft2xy = Fd - R1y - R2y + Fk == -Fa;
Moment2xy = R1y*b1 - 0.5*Fk*L + R2y*(L-b1) - 0.5*Fd*L == 0;

sol2xy = solve([kraft2xy, Moment2xy], [R1y, R2y]);

vpa(sol2xy.R1y)
vpa(sol2xy.R2y)


%%
% Case 2 yz-plane
syms R1y R2y

Fd = 70;
Fk = 105;
rd = 0.11;
rh = 0.5*0.33;
Fa =m*a;
m = 140;
a = 4.5;

kraft2yz = R1y + R2y - Fd + Fk == Fa;
Moment2yz = -Fd*rh - Fk*rd == 0;

sol2yz = solve([kraft2yz, Moment2yz], [R1y, R2y]);

vpa(sol2yz.R1y)
vpa(sol2yz.R2y)

%%
% Case 3 xy-plane
syms R1z R2z

m = 140;
Fd = 630;
Fk = 945;
a = 4.5;
Fa =m*a;
b1 = 0.16;
L = 1.1;

kraft2xy = Nb + R1z + R2z + sqrt(2)*Fb == 0;
Moment2xy = R1y*-b1 + Fb*(-(L/2)-bb) + fb == 0;

sol2xy = solve([kraft2xy, Moment2xy], [R1y, R2y]);

vpa(sol2xy.R1y)
vpa(sol2xy.R2y)