%% Beam steering efficiency 
% clc ; close all ; clear all ; 
addpath('G:\내 드라이브\Reserach_Source\A. Simulation\fmm')

%% Parameter
nm = 1e-9 ; um = 1e-6 ; lam0 = 266*nm; 
h = 500*nm ; P = 260*nm; n_m = 1.92; nSiO2 = 1.46 ; nx = 6; ny = 6;

r = [0 linspace(40*nm,73*nm,34)]; % length
%P/2- 20

for i=1:numel(r)
    clear c
    c = fmm ;
    c.setopt('verbose',false,'nvm',true,'basis','lr')
    c.set('lam0',lam0,'ax',P,'ay',P,'nx',nx,'ny',ny,'eta',45,'n1',nSiO2,'n2',1)
    if r(i)~= 0
                c.add('multiptc','d',100e-9,'nh',1,... 
                          'rect',{'n',n_m,'x',P/2,'y',P/2,'xspan',P,'yspan',P,'nvm',true})
    c.add('multiptc','d',h,'nh',1,'circ',{'n',n_m,'x',P/2,'y',P/2,'radius',r(i),'nvm',true});
    end
    c.compute
    t0(i) = c.out.tl(nx+1,ny+1 );
    T0(i) = c.out.Tl(nx+1,ny+1);
    % c.visualize('fft')
end

phase = angle(t0)
phase2 = mod(phase,2*pi)
Amplitude = T0
figure
plot(r*1e9,phase2,'s'); xlim([min(r)*1e9 max(r)*1e9]), ylim([0,2*pi])
xlabel('Radius (nm)'); ylabel('Phase')
figure
plot(r*1e9,Amplitude,'s'); xlim([min(r)*1e9 max(r)*1e9]), ylim([0,1])
xlabel('Radius (nm)'); ylabel('Amplitude')
%r(1)=0;

[r'*1e9, phase']