function T = TransferFunction(args, params, freq)
DEBUG = 0;
% fsampl = params(1);
c1 = params(2);
ro1 = params(3);

alfa0 = args(1);
freq0 = args(2);
c2 = args(3);
n = args(4);
ro2 = args(5);
h = args(6);

alf = alfa0 * (abs(freq) ./ freq0) .^ n;  % slopinimo koeficientas
C = c2;
k = 2 * pi * (freq) ./ C + 1i * alf; % randamas kompleksinis banginis skaicius vertinant slopinima
veloc = 2 * pi * (freq) ./ real(k);

Z1 = c1 * ro1; % pirmojo sluoksnio impedansas
Z2 = veloc * ro2; % antro sluoksnio impedansas
T = (2 * Z1 * Z2) ./ ...
    ( (Z1+Z2).^2 .* exp(-1i*k*h) - (Z1-Z2).^2.*exp(1i*k*h) ); %Brekhovskikh & Godin 1990 (2.4.18)
T(1) = 0;

if DEBUG
    %% Optimization for HLS
    % Euler's Cosine Identity replacement
    % https://proofwiki.org/wiki/Cosine_of_Complex_Number
    cos_kh = cos(real(k.*h)) .* cosh(imag(k.*h)) - 1i .* sin(real(k.*h)) .* sinh(imag(k.*h));
    % Euler's Sine Identity replacement
    % https://proofwiki.org/wiki/Sine_of_Complex_Number
    sin_kh = sin(real(k.*h)) .* cosh(imag(k.*h)) + 1i .* cos(real(k.*h)) .* sinh(imag(k.*h));

    % T_a = (2 * Z1 * Z2) ./ ...
    %         ((2 * Z1 * Z2) .* cos(k.*h) - 1i .* sin(k.*h) .* (Z1.^2 + Z2.^2));
    numerator = Z1 * Z2;
    denominator = (2 * Z1 * Z2) .* cos_kh - 1i .* sin_kh .* (Z1.^2 + Z2.^2);
    T_a = (Z1 * Z2) ./ ...
        ((2 * Z1 * Z2) .* cos_kh - 1i .* sin_kh .* (Z1.^2 + Z2.^2));
    T_b = numerator ./ denominator;

    % figure(10)
    % plot(real(T_a)); hold on;
    % plot(real(T_b), '--'); hold off;
    % figure(11)
    % plot(imag(T_a)); hold on;
    % plot(imag(T_b), '--'); hold off;

    T_a(1) = 0;
    T_b(1) = 0;

    %% TEST
    % dir = 'C:/GIT/VITIS/Model_sim/solution/csim/report/TransferFunction_csim.log';
    % dir = 'C:/GIT/VITIS/Model_sim/solution/csim/report/WaveSynthesis_csim.log';
    dir = 'C:/GIT/VITIS/Model_sim/solution/csim/report/pso_process_csim.log';

    freq = freq(2:end);
    denominator = denominator(2:end) / 1e9;
    cos_kh = cos_kh(2:end);
    sin_kh = sin_kh(2:end);

    cos_kh_HLS = plot_xilinx_data(dir, 'cos_kh');
    cos_kh_HLS = cos_kh_HLS(:,1)' + 1i * cos_kh_HLS(:,2)';
    cos_kh_MSE = mean(abs(cos_kh_HLS - cos_kh).^2)

    [freq_plot, cos_kh_HLS] = plot_with_freq(freq, cos_kh_HLS);
    [~, cos_kh] = plot_with_freq(freq, cos_kh);

    plot(freq_plot / 1e6, real(cos_kh_HLS), freq_plot / 1e6, imag(cos_kh_HLS), ...
         freq_plot / 1e6, real(cos_kh), '--', freq_plot / 1e6, imag(cos_kh), '--');
    title('cos(kh)'); legend('HLS re', 'HLS im', 'Matlab re', 'Matlab im');
    xlabel('Argumentas, MHz'); ylabel('Amplitudė');
    grid on;

    sin_kh_HLS = plot_xilinx_data(dir, 'sin_kh');
    sin_kh_HLS = sin_kh_HLS(:,1)' + 1i * sin_kh_HLS(:,2)';
    sin_kh_MSE = mean(abs(sin_kh_HLS - sin_kh).^2)

    [freq_plot, sin_kh_HLS] = plot_with_freq(freq, sin_kh_HLS);
    [~, sin_kh] = plot_with_freq(freq, sin_kh);
    plot(freq_plot / 1e6, real(sin_kh_HLS), freq_plot / 1e6, imag(sin_kh_HLS), ...
         freq_plot / 1e6, real(sin_kh), '--', freq_plot / 1e6, imag(sin_kh), '--');
    title('sin(kh)'); legend('HLS re', 'HLS im', 'Matlab re', 'Matlab im');
    xlabel('Argumentas, MHz'); ylabel('Amplitudė');
    grid on;

    denominator_HLS = plot_xilinx_data(dir, 'denominator');
    denominator_HLS = denominator_HLS(:,1)' + 1i * denominator_HLS(:,2)';
    denominator_MSE = mean(abs(denominator_HLS - denominator).^2)

    [freq_plot, denominator_HLS] = plot_with_freq(freq, denominator_HLS);
    [~, denominator] = plot_with_freq(freq, denominator);
    plot(freq_plot / 1e6, real(denominator_HLS), freq_plot / 1e6, imag(denominator_HLS), ...
         freq_plot / 1e6, real(denominator), '--', freq_plot / 1e6, imag(denominator), '--');
    title('Perdavimo funkcijos daliklis'); legend('HLS re', 'HLS im', 'Matlab re', 'Matlab im');
    xlabel('Argumentas, MHz'); ylabel('Amplitudė');
    grid on;
    
    % denominator_diff = (denominator_HLS - denominator);
    % figure();
    % 
    % plot(freq / 1e6, abs(real(denominator_diff) ./ max(real(denominator)) * 100), freq / 1e6, abs(imag(denominator_diff) ./ max(imag(denominator)) * 100));
    % title('Matlab ir HLS daliklio normalizuota paklaida'); legend('Re', 'Im');
    % xlabel('Dažnis, MHz'); ylabel('Paklaida, %');
    % grid on;

    % T_HLS = plot_xilinx_data(dir, 'T');
    % T_HLS = T_HLS(:,1)' + 1i * T_HLS(:,2)';
    % T_MSE = mean(abs(T_HLS - T(2:end)).^2)
    % hold on;
    % plot(real(T), '--');
    % hold on;
    % plot(imag(T), '--');
    % title('T koeficientai'); legend('HLS re', 'HLS im', 'Matlab re', 'Matlab im');
    % xlabel('Indeksas'); ylabel('Amplitudė');

    if ( max(real(T) - real(T_a)) > 0.001 || max(imag(T) - imag(T_a)) > 0.001 )
        disp("difference!---------------------------------------------------");
    end
end

end
