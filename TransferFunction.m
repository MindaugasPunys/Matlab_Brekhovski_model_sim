function T=TransferFunction(args,params,freq)
        alfa0=args(1);
        c2=args(2);
        n=args(3);
        ro2=args(4);
        h=args(5);
        
        fsampl=params(1);
        ro1=params(3);
        freq0=params(4);
        
             Temperature=args(7);
             Pressure=args(8);
             c1=331.38*sqrt(1+Temperature/273.15); % velocity in air
             %gama=1.4; % air adiabatic constant 
             %Ks=gama*Pressure; % coefficient of stiffness
             ro1=1.4*Pressure/(c1^2);  % Ks/(c1^2); % air density
        
        alf=alfa0*(abs(freq)./freq0).^n;  % slopinimo koeficientas 
        C=c2;
        k=2*pi*(freq)./C+1i*alf; % randamas kompleksinis banginis skaicius vertinant slopinima
        %veloc=2*pi*(freq)./(k);  % randamas sluosknio ultragarso greitis 
        veloc=2*pi*(freq)./real(k);

Z1=c1*ro1; % pirmojo sluoksnio impedansas
Z2=veloc*ro2; % antro sluoksnio impedansas
T=(4*Z1*Z2)./((Z1+Z2).^2.*exp(-1i*k*h)-(Z1-Z2).^2.*exp(1i*k*h)); %Brekhovskikh & Godin 1990 (2.4.18)
T(1)=0;
end
