%% Initialization
clear;
close all;

%Using Alyssa's excel spreadsheet, these should be the only quantities you
%have to change

hfig = openfig("pour8.fig");
pour = 8;

disp("Taking data for trial #" + pour);

%Plot the original figure for reference
title("Original Figure")

%Define gravitational constant
g = 9.81; %m/s^2

%Get the properties of the pour from the excel spreadsheet
data = readmatrix('PourProperties.xlsx');
rho=data(pour,3); %fluid density, kg/m^3
mu=data(pour,4); %viscosity, Pa*s
Q=data(pour,6); %Volume, m^3

disp("Nominal/Recorded Values");
disp("\rho = " + rho);
disp("\mu = " + mu);
disp("Q = " + Q * 10^6 + "mL");

%% Plot theoretical and measured radii

%Get the x and y of our measured data
axsobj = hfig.Children;
dataobjs = axsobj.Children;

%Apply a time shift and convert from diameter to radius
t = dataobjs(2).XData-1;
R = (dataobjs(2).YData/2)/1000;

%Plot it in a more aesthetic way
figure;
plot(t,R*1000,'x')
hold on;

%Generate a theoretical plot of r and dr/dt (in mm and mm/s, respectively)
t_theo = 0:.1:max(t);
R_theo = 0.779 .* (g .* rho .* Q .^3 .* t_theo ./ mu) .^ (1/8);
dR_theo = (0.779/8) .* (g .* rho .* Q .^3 ./ mu) .^ (1/8) .* (t_theo.^(-7/8));

%Plot the theoretical next to the measured data
plot(t_theo,R_theo*1000,'--');

title("Radius vs. Time for Trial " + pour);
xlabel("t (s)");
ylabel("R (mm)");

legend(["Measured", "Theoretical"]);

%% Plot U = dR/dt theoretical and measured

figure;

%Measure the t interval
t_int = t(2)-t(1);

%Shift the x to be between intervals, and truncate
deriv_t = t + t_int/2;
deriv_t = deriv_t(1:end-1);

%Initialize derivative array
dR = zeros(1, numel(R)-1);

for i = 1:(numel(R)-1)
    dR(i) = (R(i+1) - R(i)) / (t_int);
end

plot(deriv_t, dR*1000, 'x')
hold on;
plot(t_theo, dR_theo*1000, '--')

title("Radius Expansion Speed for Trial " + pour);
xlabel("t (s)");
ylabel("U=dR/dt (mm/s)");
legend(["Measured", "Theoretical"]);

%% plot height vs. time, assuming cylindrical spread
figure;

h = Q ./ ( R .^2 .* pi);

plot(t, h*1000);
title("Cylindrical Height for Trial " + pour);
ylabel("H (mm)")
xlabel("t (s)");

%% find value of mu from data
figure;
mu_data = (0.779./R).^8 .* g .* rho .* Q.^3 .* t;

plot(t,mu_data);
title("Calculated Mu vs. Time");

mu_avg = sum(mu_data)/numel(mu_data);
disp("Average mu: " + mu_avg);


%% Re and Bond

sig = 0.073; %N/m
U = dR;
L = h;
L = L(2:end);

re = U .* rho .* L ./ mu;
figure
plot(t(2:end), re);
title("Re in time");
avg_re = sum(re)/numel(re);
max_re = max(re);
disp("Reynold Number Average: " + avg_re);
disp("Reynold Number Maximum: " + max_re);
xlabel("t (s)");
ylabel("Re");

bo = rho .* g .* L.^2 ./ sig;
%disp("Bond number: " + bo);
figure
plot(t(2:end),bo);
title("Bo in time");
xlabel("t (s)");
ylabel("Bo");

avg_bo = sum(bo)/numel(bo);
min_bo = min(bo);
disp("Bond Number Average: " + avg_bo);
disp("Bond Number Minimum: " + min_bo);


%% Plot adjusted theoretical prediction
figure;


plot(t, R*1000, "x");

tr = 0.779 .* (g .* rho .* Q .^3 .* t_theo ./ mu_avg) .^ (1/8);

hold on;

plot(t_theo, tr*1000, "--");
title("Adjusted Theoretical and Measured for Trial " + pour);
legend(["Measured", "Theoretical"]);
xlabel("t (s)");
ylabel("R (mm)");

%% calculate dimensionless quantities for a set of pours
%This part of the code is separate from the rest. It takes in a set of
%Pours and plots them dimensionlessly against each other along with the 
%Theoretical

clear;

pours = [4,5,6,7,8];
data = readmatrix('PourProperties.xlsx');
g = 9.81;

figure;

for i = pours
    
    fig_i = openfig("pour" + i + ".fig");
    
    %Get the x and y of our measured data
    axsobj_i = fig_i.Children;
    dataobjs_i = axsobj_i.Children;

    %Apply a time shift and convert from diameter to radius
    t_i = dataobjs_i(2).XData-1;
    R_i = (dataobjs_i(2).YData/2)/1000;
    
    rho_i=data(i,3); %fluid density, kg/m^3
    mu_i=data(i,4); %viscosity, Pa*s
    Q_i=data(i,6); %Volume, m^3
    
    mu_data_i = (0.779./R_i).^8 .* g .* rho_i .* Q_i.^3 .* t_i;
    mu_avg_i = sum(mu_data_i)/numel(mu_data_i);
    
    close(fig_i);
    
    Rstar = R_i ./ (Q_i.^(1/3));
    tstar = t_i .* rho_i .* g .* Q_i.^(1/3) ./ mu_avg_i;

    plot(tstar, Rstar, "x");
    hold on;
end

%These are for the theoretical calculation. Their value doesn't matter, I'm
%just setting them this way to get the correct scaling. Changing them won't
%change the shape of the plot. 
rho_t = 1;
Q_t = 1;
mu_t = .1;

t_dimless = 0:0.1:20;
R_dimless = 0.779 .* (g .* rho_t .* Q_t .^3 .* t_dimless ./ mu_t) .^ (1/8);

t_dimless = t_dimless .* rho_t .* g .* Q_t .^ (1/3) ./ mu_t;
R_dimless = R_dimless ./ (Q_t .^ (1/3));

plot(t_dimless, R_dimless, "k-");

title("Dimensionless Radial Growth");
xlabel("t*");
ylabel("R*");

legend(["Trial " + pours, "Theoretical"]);

