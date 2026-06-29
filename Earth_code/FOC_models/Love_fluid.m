%-------------------------------------------------------------------------%
% --- This is a basic version of the code implemented specifically ------ %
% --- for the Earth. For all the details, please refer to the ----------- %
% --- general code for super-Earths. ------------------------------------ %
%-------------------------------------------------------------------------%
%INTEGRATE FROM NEAR THE CENTRE AND NOT FROM NEAR CMB DUE TO
%INSTABILITIES AND MISSING DEGREES LNs

clear clc

%Select 'Tidal' / 'Load' / 'Shear' type of LNs
Love_type = 'Shear';

N_max = 25;  %max harmonic degree
n = 2:N_max;  %harmonic degree

%result arrays
h_array_inc = zeros(1, N_max-1);
k_array_inc = zeros(1, N_max-1);
l_array_inc = zeros(1, N_max-1);

h_array_c = zeros(1, N_max-1);
k_array_c = zeros(1, N_max-1);
l_array_c = zeros(1, N_max-1);


%LNs calculation
for i = 2:N_max
    [h, k, l] = main_incompressible(i, Love_type);
    h_array_inc(i-1) = h;
    k_array_inc(i-1) = k;
    l_array_inc(i-1) = l;
end

for i = 2:N_max
    [h, k, l] = main_compressible(i, Love_type);
    h_array_c(i-1) = h;
    k_array_c(i-1) = k;
    l_array_c(i-1) = l;
end


%PLOTTING------------------------------------------------------------------
%Incompressible LNs
figure(1)
subplot(1,3,1)
plot(n,h_array_inc, '-b','LineWidth',1,'Marker','o','MarkerSize',2.5,'MarkerFaceColor','k')
title('h (Fluid outer core, Incompressible)')
xlabel('Harmonic degree n')
ylabel('h_n')
grid on

subplot(1,3,2)
plot(n,k_array_inc, '-b','LineWidth',1,'Marker','o','MarkerSize',2.5,'MarkerFaceColor','k')
title('k (Fluid outer core, Incompressible)')
xlabel('Harmonic degree n')
ylabel('k_n')
grid on

subplot(1,3,3)
plot(n,l_array_inc, '-b','LineWidth',1,'Marker','o','MarkerSize',2.5,'MarkerFaceColor','k')
title('l (Fluid outer core, Incompressible)')
xlabel('Harmonic degree n')
ylabel('l_n')
grid on


%Compressible LNs
figure(2)
subplot(1,3,1)
plot(n,h_array_c, '-r','LineWidth',1,'Marker','o','MarkerSize',2.5,'MarkerFaceColor','k')
title('h (Fluid outer core, Compressible)')
xlabel('Harmonic degree n')
ylabel('h_n')
grid on

subplot(1,3,2)
plot(n,k_array_c, '-r','LineWidth',1,'Marker','o','MarkerSize',2.5,'MarkerFaceColor','k')
title('k (Fluid outer core, Compressible)')
xlabel('Harmonic degree n')
ylabel('k_n')
grid on

subplot(1,3,3)
plot(n,l_array_c, '-r','LineWidth',1,'Marker','o','MarkerSize',2.5,'MarkerFaceColor','k')
title('l (Fluid outer core, Compressible)')
xlabel('Harmonic degree n')
ylabel('l_n')
grid on


%Comparison
figure(3)
subplot(1,3,1)
plot(n,h_array_inc, '-b','LineWidth',1,'Marker','o','MarkerSize',2,'MarkerFaceColor','b')
hold on
plot(n,h_array_c, '-r','LineWidth',1,'Marker','o','MarkerSize',2,'MarkerFaceColor','r')
hold off
title('Love number h'' ')
legend('Incompressible', 'Compressible', 'Location','northeast')
xlabel('Harmonic degree n')
ylabel('h_n'' ')
grid on

subplot(1,3,2)
plot(n,k_array_inc, '-b','LineWidth',1,'Marker','o','MarkerSize',2,'MarkerFaceColor','b')
hold on
plot(n,k_array_c, '-r','LineWidth',1,'Marker','o','MarkerSize',2,'MarkerFaceColor','r')
hold off
title('Love number k'' ')
legend('Incompressible', 'Compressible', 'Location','northeast')
xlabel('Harmonic degree n')
ylabel('k_n'' ')
grid on

subplot(1,3,3)
plot(n,l_array_inc, '-b','LineWidth',1,'Marker','o','MarkerSize',2,'MarkerFaceColor','b')
hold on
plot(n,l_array_c, '-r','LineWidth',1,'Marker','o','MarkerSize',2,'MarkerFaceColor','r')
hold off
title('Love number l'' ')
legend('Incompressible', 'Compressible', 'Location','northeast')
xlabel('Harmonic degree n')
ylabel('l_n'' ')
grid on