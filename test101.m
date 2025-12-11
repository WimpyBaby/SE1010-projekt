clear
clc
% Reservera minne för plottning
x1 = zeros(1,145);
y1 = zeros(1,145);
y2 = zeros(1,145);
y3 = zeros(1,145);
%indata
h0 = 1.5; % Plåttjocklek innan valsning, mm
%h1 = 0.5; % Plåttjocklek efter valsning, samma startvärde, mm
R = 94/2; % Valsens radie, valsdiameter = 94 mm i labbvalsverk
sigmaf = 50; % Materialets flytspänning i N/mm2
b0 = 27; % Plåtens bredd före valsning, 27 mm
m = 0.5; % Friktionsfaktor, varierar mellan 0 och 1
n = 141; % Valsarnas varvtal, varv/min
verkningsgrad = 0.8; % Drivlinans verkningsgrad, 80 %
omega = 2*pi*n/60; % Vinkelhastighet, rad/s
% Formler, se kursmaterialet, finns utrymme för beräkningseffektivisering

target_Pmotor = 2000;
Pmotor = target_Pmotor + 1; % Initialize Pmotor to a value greater than target_Pmotor

while Pmotor >= target_Pmotor
    h1 = 0.5;

    for ii = 1:1:10
        x1(ii) = h1;
        deltah = h0 - h1;
        hmedel = (h0 + h1) / 2;
        L = sqrt(R * deltah);
        y2(ii) = L;
        F = 1/sqrt(3) * sigmaf * L * b0 * (2 + m * L / hmedel);
        y3(ii) = F;
        Mv = F * L;
        Pvals = Mv * omega * 0.001; % Omvandla Nmm till Nm för effektberäkning
        Pmotor = Pvals / verkningsgrad;
        y1(ii) = Pmotor;
        R = R - 0.0001;
    end
end

subplot(1,3,1);
plot(x1,y1);
grid;
title('Motoreffekt');
xlabel('h1 [mm]'); ylabel('Motoreffekt [W]');

subplot(1,3,2);
plot(x1,y2);
grid;
title('Kontaktlängd');
xlabel('h1 [mm]'); ylabel('Kontaktlängd [mm]');

subplot(1,3,3);
plot(x1,y3);
grid;
title('Valskraft');
xlabel('h1 [mm]'); ylabel('Valskraft [N]');