%Using Alyssa's excel spreadsheet, these should be the only quantities you
%have to change
hfig = openfig("pour6.fig");
pour = 6;

title("Original Figure")

axsobj = hfig.Children;
dataobjs = axsobj.Children;

x = dataobjs(2).XData-1
y = dataobjs(2).YData/2

figure;
plot(x,y,'x')
hold on;

g = 9.81; %m/s^2

data = readmatrix('PourProperties.xlsx');
rho=data(pour,3); %fluid density, kg/m3
mu=data(pour,4); %viscosity, Pa s
Q=data(pour,6); %Volume (m^3)

tx = 0:.1:max(x);

r = 0.779 .* (g .* rho .* Q .^3 .* tx ./ mu) .^ (1/8);
dr = (0.779/8) .* (g .* rho .* Q .^3 ./ mu) .^ (1/8) .* (tx.^(-7/8));

plot(tx,r*1000,'--');

title("Radius vs. Time for Pour " + pour);
xlabel("t (s)");
ylabel("R (mm)");

figure;

x_int = x(2)-x(1)

deriv_x = x + x_int/2;

deriv_y = zeros(size(y));

for i = 1:numel(x)-1
    deriv_y(i) = (y(i+1) - y(i)) / (x_int);
end

plot(deriv_x, deriv_y, 'x')
hold on;
plot(tx, dr*1000, '--')

title("Radius Expansion Speed for Pour " + pour);
xlabel("t (s)");
ylabel("U (mm/s)");

sig = 0.073; %N/m
U = 0;
L = 0;

re = U * rho * L / mu;
disp("Reynolds number: " + re);

bo = rho * g * L^2 / sig;
disp("Bond number: " + bo);

%calculate dimensionless quantities
figure;

Rstar = y ./ (Q.^(1/3));
tstar = x .* rho .* g .* Q.^(1/3) ./ mu;

plot(tstar, Rstar, "-");
title("Dimensionless Radial Growth");
xlabel("R*");
ylabel("t*");

%plot height vs. time, assuming cylindrical spread
figure;

h = Q ./ ((y./1000).^2 .* pi);

plot(x, h*1000);
title("Cylindrical Height");
ylabel("H (mm)")
xlabel("t (s)");

