function T = TransferFunction(args, params, freq)


c1 = params(2);
ro1 = params(3);

alfa0 = args(1);
freq0 = args(2);
c2 = args(3);
n = args(4);
ro2 = args(5);
h = args(6);

fsampl = params(1);
alf = alfa0 * (abs(freq) ./ freq0) .^ n;  % slopinimo koeficientas
C = c2;
k = 2 * pi * (freq) ./ C + 1i * alf; % randamas kompleksinis banginis skaicius vertinant slopinima
%veloc = 2 * pi * (freq) ./ (k);  % randamas sluosknio ultragarso greitis
veloc = 2 * pi * (freq) ./ real(k);

Z1 = c1 * ro1; % pirmojo sluoksnio impedansas
Z2 = veloc * ro2; % antro sluoksnio impedansas
T = (4 * Z1 * Z2) ./ ...
    ( (Z1+Z2).^2 .* exp(-1i*k*h) - (Z1-Z2).^2.*exp(1i*k*h) ); %Brekhovskikh & Godin 1990 (2.4.18)
T(1) = 0;
end
