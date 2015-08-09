function [Tnfn_SAVD] = CalResponseSpectra_SP(fn, dt, acc, damping, G)
%---------------------------------------------------------
% SDOF Response Spectra
% Using Piecewise Exact Method and a given vector of frequency points
%
% Jongwon Lee (c)
% University of Michigan
%
% March 20 2007
% Jan. 2 2009   Modified for sinusoidal excitation
% Sept. 1, 2010 Transformed to a function
% March 17, 2011 PSA and PSV were added.
%
% [Output]
% Tnfn_SAVD = a combined matrix consisting of the followings in column:
%   Tn = natural period of SDOF oscillator (sec)
%   fn = corresponding natural freqeuncy of SDOF oscillator (Hz)
%   PSA = (total) pseudo-spectral acceleration in G
%   SA = (total) spectral acceleration in G
%   PSV = (relative) pseudo-spectral velocity in unit of given G times sec.
%   SV = (relative) spectral velocity in unit of given G times sec.
%   SD = (relative) spectral displacement in unit of given G times sec^2.
%
% [Input]
% fn = a given vector of natural frequency of oscillater (Hz)
% Np = desired number of uniformly-spaced points within a period decade
%      in log scale; e.g., 100 pts. per order of magnitude (e.g., 1 - 10 s)
% dt = time interval. (sec)
% acc = acceleration time history (g)
% damping = damping of SDOF system in decimal
% G = gravitational acceleration; e.g., 9.81 m/s^2; 32.174 ft/s^2
%---------------------------------------------------------

% Generate time history.
t = 0:dt:(length(acc)-1)*dt;

acc_g1 = acc.*G;

wn = 2.*pi.*fn;                             % [rad/sec] Natural freq
beta = damping;                             % damping ratio in decimal
wd = wn*sqrt(1 - beta^2);                   % damped freq.

% Piecewise Exact Method===================================================
u_max = zeros(1,length(wd));
ud_max = zeros(1,length(wd));
udd_max = zeros(1,length(wd));
utdd_max = zeros(1,length(wd));
PSV = zeros(1,length(wd));
PSA = zeros(1,length(wd));

u_exact = zeros(1,length(acc_g1));
ud_exact = zeros(1,length(acc_g1));
udd_exact = zeros(1,length(acc_g1));
utdd_exact = zeros(1,length(acc_g1));


for i = 1:length(wd);
    A = exp(-beta*wn(i)*dt)*(beta/sqrt(1-beta^2)*sin(wd(i)*dt)+cos(wd(i)*dt));
    B = exp(-beta*wn(i)*dt)*(1/wd(i)*sin(wd(i)*dt));
    C = 1/wn(i)^2*(2*beta/wn(i)/dt + exp(-beta*wn(i)*dt)*(((1-2*beta^2)/wd(i)/dt-beta/sqrt(1-beta^2))...
        *sin(wd(i)*dt) - (1+2*beta/wn(i)/dt)*cos(wd(i)*dt)));
    D = 1/wn(i)^2*(1 - 2*beta/wn(i)/dt + exp(-beta*wn(i)*dt)*((2*beta^2-1)/wd(i)/dt*sin(wd(i)*dt)...
        +2*beta/wn(i)/dt*cos(wd(i)*dt)));
    
    A_d = -exp(-beta*wn(i)*dt)*(wn(i)/sqrt(1-beta^2)*sin(wd(i)*dt));
    B_d = exp(-beta*wn(i)*dt)*(cos(wd(i)*dt) - beta/sqrt(1-beta^2)*sin(wd(i)*dt));
    C_d = 1/wn(i)^2*(-1/dt + exp(-beta*wn(i)*dt)*((wn(i)/sqrt(1-beta^2)+beta/dt/sqrt(1-beta^2))*sin(wd(i)*dt)...
        +1/dt*cos(wd(i)*dt)));
    D_d = 1/wn(i)^2/dt*(1 - exp(-beta*wn(i)*dt)*(beta/sqrt(1-beta^2)*sin(wd(i)*dt) + cos(wd(i)*dt)));
    
    
    for j = 1:length(acc_g1)-1
        u_exact(1,j+1) = u_exact(j)*A + ud_exact(j)*B + (-1)*acc_g1(j)*C + (-1)*acc_g1(j+1)*D;
        ud_exact(1,j+1) =u_exact(j)*A_d + ud_exact(j)*B_d + (-1)*acc_g1(j)*C_d + (-1)*acc_g1(j+1)*D_d;
        udd_exact(1,j+1) = -(2*wn(i)*beta*ud_exact(1,j+1)+wn(i)^2*u_exact(1,j+1)+acc_g1(j+1));
        utdd_exact(1, j+1) = udd_exact(1, j+1) + acc_g1(j+1);
    end
    
    u_max(1,i) = max(abs(u_exact));
    ud_max(1,i) = max(abs(ud_exact));
    udd_max(1,i) = max(abs(udd_exact));
    utdd_max(1,i) = max(abs(utdd_exact));
    PSV(1,i) = u_max(1,i)*wn(i);
    PSA(1,i) = PSV(1,i)*wn(i);
end

SA = (utdd_max).';   % (Total or Absolute) Spectral Acceleration in g
SV = ud_max.';          % (Relative) Spectral Velocity
SD = u_max.';           % (Relative) Spectral Displacement
PSA = (PSA).';       % (Total) Pseudo-spectral Acceleration in g
PSV = PSV.';            % (Relative) Pseudo-spectral Velocity

Tn = 1./fn;
Tnfn_SAVD = [Tn fn PSA SA PSV SV SD]; %combine in a matrix.
% Tnfn_SAVD = flipud(Tnfn_SAVD);    %flip up and down.

% figure; % Response Spectra
% num = 410;
% subplot(num+1); plot(t, acc_g1./G, 'r-', 'LineWidth', 1.5);
% xlabel('time(sec)');
% ylabel('Base Accel.(g)'); grid on;
% title('Piece-wise Exact');
%
% % subplot(num+2); semilogx(Tn, udd_max./G, 'b-', 'LineWidth', 1.5);
% % xlabel('Period(sec)');
% % ylabel('Relative PSA(g)'); grid on;xlim([Tmin Tmax]);
%
% subplot(num+2); semilogx(Tn, utdd_max./G, 'b-', 'LineWidth', 1.5);
% xlabel('Period(sec)');
% ylabel('Total PSA(g)'); grid on;
%
% subplot(num+3); semilogx(Tn, ud_max, 'b-', 'LineWidth', 1.5);
% xlabel('Period(sec)');
% ylabel('Relative PSV(Length/sec)'); grid on;xlim([Tmin Tmax]);
%
% subplot(num+4); semilogx(Tn, u_max, 'b-', 'LineWidth', 1.5);
% xlabel('Period(sec)');
% ylabel('SD(m)'); grid on; xlim([Tmin Tmax]);
