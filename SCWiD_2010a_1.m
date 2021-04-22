%% Clear Command
clear all; close all; clc;

%% INPUT GUI
% Alerts added to guide the user which file to select
% alert1=msgbox('STEP1:  Select INTERIOR data file','Icon','help');
% waitfor(alert1);
[DataFileName,PathName1] = uigetfile('./*.dat','Select data file');
datafile = [PathName1 DataFileName];
% alert2=msgbox('STEP2:  Select PORTLOC data file','Icon','help');
% waitfor(alert2);
[PortsFileName,PathName2] = uigetfile('./*.dat','Select ports location file');
portsfile = [PathName2 PortsFileName];
% alert3=msgbox('STEP3:  Verify/Edit Input Data ','Icon','help');
% waitfor(alert3);
box=inputdlg({'Reference Pressure (Psi): ','Reference Temperature (R): '...
    ,'Reference Density (lbm/ft^3): ','Gamma: ','Speed of Sound (ft/s):'...
    ,'Theta Opening Time (NON-DIM): '},'Input Data',1,{'14.7','518',...
    '0.073473','1.4','1138.898','0'});

%% Type of Post-processing
S1=' Full Post-Processing';
S2='Basic Post-Processing';
str=[S1;S2];
[SELECT,ok] = listdlg('PromptString','Select Post-Processing Type ',...
                'SelectionMode','single',...
                'ListString',str,'ListSize',[200 100]);
switch(SELECT)
    case(1)
PORTOWALL = questdlg('Do you want to average data over ports or the entire walls?', ...
                         'Select Port/Wall', ...
                         'Port', 'Wall', 'Port');
end

%% Initializing Values
n=1;     %ncells=100;
dt=0;    %tcycle=6200;

%% PROCESSING Message
wb1 = waitbar(0,'Reading Data...');

%% PROCESSING INPUT DATA 
Pref=str2double(box(1)); Tref=str2double(box(2)); Dref=str2double(box(3));
gam=str2double(box(4)); sos=str2double(box(5)); THETOT=str2double(box(6));

%% Loading Ports Data
fid = fopen(datafile,'r');   %load PortsFile
Port=load(portsfile);    %Port = load PortsFile;
 
ncells=50;
tcycle=3140;
% *********************************************
THETAL1=Port(2)+(9*3.14/180);    %Inlet port
THETAL2=Port(3)+(9*3.14/180);
THETAR1=Port(6)+(9*3.14/180);    %Exhaust port
THETAR2=Port(7)+(9*3.14/180);
%% Data check
% make sure data is available in this time range:
kcyci=0;
kcycf=1;
if (kcyci==kcycf);
   k=kcyci-1;
else
   k=0;
end

DX=1./ncells;   
for j=1:ncells+2
   X(j)=(j-1)*DX;
end

for l=1:(kcycf-kcyci)
   cycle(l)=kcyci*tcycle+(l-1)*tcycle;
end
i=1; m=1;
time=fscanf(fid,'%f',[1]);

while(feof(fid)==0)
    Mt=fscanf(fid,'%f%f%f%f%f%f%f',[6,ncells+2]);
   	t(i)=(time+cycle(m))/1000;
   	PR(:,i)=Mt(1,:)'
   	DEN(:,i)=Mt(2,:)'
      VEL(:,i)=Mt(3,:)'
      Z(:,i)=Mt(4,:)'
      ZO(:,i)=Mt(5,:)'
      ZI(:,i)=Mt(6,:)'
%      AREA(:,i)=Mt(7,:)'
%      Twall(:,i)=Mt(5,:)';
      i=i+1 
      waitbar(i/i)
   %end
   time=fscanf(fid,'%f',[1]);
   if time==0
      m=m+1;
	end    
end
close(wb1)
%% Setting Variables for Plotting
j=64;
%TEM=(PR./DEN)*1000;
%Mach.No=(VEL./1000)./((TEM./1000).^0.5);
%FFLX=(PR./1000)./(TEM./1000).*(VEL./1000).*(Z./1000);
%Pstg=PR.*((1+(gam-1)/2.*(Mach.No.^2)).^(gam/(gam-1)));
%Tstg=TEM.*(1+(gam-1)/2.*(Mach.No.^2));
%expS=(TEM./1000).*((PR./1000).^((1-gam)/gam));
TEM=(PR./DEN)*1;
Mach.No=(VEL./1)./((TEM./1).^0.5);
FFLX=(PR./1)./(TEM./1).*(VEL./1).*(Z./1);
Pstg=PR.*((1+(gam-1)/2.*(Mach.No.^2)).^(gam/(gam-1)));
Tstg=TEM.*(1+(gam-1)/2.*(Mach.No.^2));
expS=(TEM./1).*((PR./1).^((1-gam)/gam));
cmmap=clrmap(j);

%% Right Wall calculations (at x=1)
%switch(SELECT)
    %case(1)
%wb2 = waitbar(0,'Processing Data...');
%cc1=1;
%switch (PORTOWALL)
    %case('Port')
%NSTEPS1=100; NDIVS1=10;
%SPOS1=THETAR1; EPOS1=THETAR2;
    %case('Wall')
%NSTEPS1=length(t)-1; NDIVS1=10;
%SPOS1=0; EPOS1=6.283;
%end
%waitbar(cc1/cc1);
 %[Pavg1,Vavg1,Tavg1,PSavg1,TSavg1,Pfluc1,Pprms1,PSfluc1,...
     %PSprms1,Tmax1,Tmin1,TSmax1,TSmin1,Vfluc1,Vprms1]...
     %=AVERAGE(t,VEL,PR,TEM,Pstg,Tstg,ncells+2,NSTEPS1,NDIVS1,THETOT,SPOS1,EPOS1);
%close(wb2)
%end

%% Left Wall Calculations (at x=0)
%switch (SELECT)
    %case(1)
%wb3 = waitbar(0,'Analyzing Data ...');
%switch (PORTOWALL)
    %case('Port')
%NSTEPS2=100; NDIVS2=10;
%SPOS2=THETAL1; EPOS2=THETAL2;
    %case('Wall')
%NSTEPS2=length(t)-1; NDIVS2=10;
%SPOS2=0; EPOS2=6.283;
%end
%waitbar(cc1/cc1);
%[Pavg0,Vavg0,Tavg0,PSavg0,TSavg0,Pfluc0,Pprms0,PSfluc0,...
    %PSprms0,Tmax0,Tmin0,TSmax0,TSmin0,Vfluc0,Vprms0]...
    %=AVERAGE(t,VEL,PR,TEM,Pstg,Tstg,1,NSTEPS2,NDIVS2,THETOT,SPOS2,EPOS2);
%close(wb3)
%end

%% Adjusting Reference Angle 
%%COMMENTED OUT BY PAWAN TO ADJUST PLOTTING REFERENCE ANGLE TO ZERO
%%BELOW MODULE IS USED FOR PURDUE WR REF. ANGLES
% Port=Port.*360./6.283;
% pp=Port(2)+9; %(THETOT*180/2/pi);
% tcycle=tcycle*360/6.283;
% t=t*360/6.283 - pp;
% Port = Port - pp;

%%ADDED BY PAWAN TO ADJUST PLOTTING REFERENCE ANGLE TO ZERO
Port=Port.*360./6.283;
t=t.*(180.0/pi);
% THETOT=0.4762;
THETADD=THETOT/2.0;
THETADD_Deg=THETADD*(180.0/pi);
pp=THETADD_Deg;
tcycle=tcycle*360/6.283;
%pp=Port(2)+9; %(THETOT*180/2/pi);
t=t-pp;
%t=t*360/6.283 - pp;
Port = Port - pp;
%% COMPUTING SPILLAGE PERCENTAGE
%qq=0;
%ww=1;                                 
%while qq>=0
    %qq=qq+1;            
    %if t(qq)>=50 %Port(2)+9
        %tspil(ww)=t(qq);                   
        %Fspil(ww)=FFLX(ncells+2,qq);
        %ww=ww+1;
    %end                  
    %if t(qq)>= Port(7)+9
        %qq=-1;
    %end
%end
%SPILFLUX=trapz(tspil,Fspil); 
%INLTFLUX=trapz(t,FFLX(1,:));
%SPIL_PERCNT=SPILFLUX/INLTFLUX*100

%% ANGULAR POSITION SETUP
jj=0;
for ii=1:length(t)
    jj=jj+1;
    if t(ii)>= -60
        break
    end
end
ang_pos=[t(1,jj:length(t)) t(1,1:jj-1)];
%% OUTPUT Visualizations
clf;
DWG(1)=figure(1);
subplot(131);
colormap jet;
plot((VEL(1,:)./1000),t,'b',(VEL(ncells+2,:)./1000),t,'r--','LineWidth',1.5);
%plot((VEL(1,:)./1),t,'b',(VEL(ncells+2,:)./1),t,'r--','LineWidth',1.5);
% plot((PR(1,:)./1000),t,'b',(PR(ncells+2,:)./1000),t,'r--','LineWidth',1.5);
%pcolor(X,t,log(PR'./TEM'));
%shading flat;
pos0 = get(gca,'Position');
%axis([-0.6 0.6 0 180])
set(gca,'Position',pos0 - [0.027 0 0.05 0])
set(gca,'YTick',[-60:30:300])
%klb0=colorbar('Position',[0.24 0.2 0.012 0.6]);
%xlabel('Position/Length');
ylabel('Angular position (degrees)');
legend('Left', 'Right')
grid;
title('Velocity');
% hold
% XX=[0;0];
% tt=[Port(1,1)+9; Port(2,1)+9];
% plot(XX,tt,'black','LineWidth',1);
% tt=[Port(2,1)+9; Port(3,1)+9];
% plot(XX,tt,'white','LineWidth',2);
% tt=[Port(3,1)+9; Port(4,1)+9];
% plot(XX,tt,'black','LineWidth',1);
% 
% XX=[1;1];
% tt=[Port(5,1)+9; Port(6,1)+9];
% plot(XX,tt,'black','LineWidth',1);
% tt=[Port(6,1)+9; Port(7,1)+9];
% plot(XX,tt,'white','LineWidth',2.5);
% tt=[Port(7,1)+9; Port(8,1)+9];
% plot(XX,tt,'black','LineWidth',1);
% tt=[Port(8,1)+9; Port(9,1)+9];
% plot(XX,tt,'<black','LineWidth',1.5);
% % tt=[Port(9,1)+9; Port(10,1)+9];
% % plot(XX,tt,'black','LineWidth',1);
% % tt=[Port(10,1)+9; Port(11,1)+9];
% % plot(XX,tt,'white','LineWidth',3);
% % tt=[Port(11,1)+9; Port(12,1)+9];
% % plot(XX,tt,'black','LineWidth',1);


subplot(132);
colormap jet;
%colormap(gray);
%pcolor(X,t,TEM'./1000);
pcolor(X,t,TEM'./1);
%pcolor(X,t,expS');
shading flat;
axis([-0.01 1.01 0 180])
pos1 = get(gca,'Position');
set(gca,'Position',pos1 - [0.017 0 0.05 0])
set(gca,'YTick',[-60:30:300])
%caxis([0.3 0.7])
colorbar;
klb1=colorbar('Position',[0.58 0.2 0.012 0.6]);
xlabel('Position/Length');
title('Temperature');

hold
XX=[-0.01;-0.01];
tt=[Port(1,1); Port(2,1)];
plot(XX,tt,'black','LineWidth',4);
tt=[Port(2,1); Port(3,1)];
plot(XX,tt,'white','LineWidth',3.5);
tt=[Port(3,1); Port(4,1)];
plot(XX,tt,'black','LineWidth',4);
tt=[Port(4,1); Port(5,1)];
plot(XX,tt,'white','LineWidth',3.5);
tt=[Port(5,1); Port(6,1)];
plot(XX,tt,'black','LineWidth',4);
tt=[Port(6,1); Port(7,1)];
plot(XX,tt,'white','LineWidth',3.5);
tt=[Port(7,1); Port(8,1)];
plot(XX,tt,'black','LineWidth',4);
tt=[Port(8,1); Port(9,1)];
plot(XX,tt,'white','LineWidth',3.5);
tt=[Port(9,1); Port(10,1)];
plot(XX,tt,'black','LineWidth',4);

XX=[1.01;1.01];
tt=[Port(11,1); Port(12,1)];
plot(XX,tt,'black','LineWidth',4);
tt=[Port(12,1); Port(13,1)];
plot(XX,tt,'white','LineWidth',3.5);
tt=[Port(13,1); Port(14,1)];
plot(XX,tt,'black','LineWidth',4);
tt=[Port(14,1); Port(15,1)];
plot(XX,tt,'white','LineWidth',3.5);
tt=[Port(15,1); Port(16,1)];
plot(XX,tt,'black','LineWidth',4);
tt=[Port(16,1); Port(17,1)];
plot(XX,tt,'white','LineWidth',3.5);
tt=[Port(17,1); Port(18,1)];
plot(XX,tt,'black','LineWidth',4);
tt=[Port(18,1); Port(19,1)];
plot(XX,tt,'white','LineWidth',3.5);
tt=[Port(19,1); Port(10,1)];
plot(XX,tt,'black','LineWidth',4);

subplot(133);
colormap jet;
%colormap(gray);
pcolor(X,t,log10(PR'./1000));
%pcolor(X,t,log10(PR'./1));
%pcolor(X,t,expS');
shading flat;
axis([-0.01 1.01 0 180])
pos2 = get(gca,'Position');
set(gca,'Position',pos2 + [0 0 -0.05 0])
set(gca,'YTick',[-60:30:300])
% caxis([0 0.5])
colorbar;
klb2=colorbar('position',[0.9 0.2 0.012 0.6]);
xlabel('Position/Length');
title('Log Pressure');

hold
XX=[-0.01;-0.01];
tt=[Port(1,1); Port(2,1)];
plot(XX,tt,'black','LineWidth',4);
tt=[Port(2,1); Port(3,1)];
plot(XX,tt,'white','LineWidth',3.5);
tt=[Port(3,1); Port(4,1)];
plot(XX,tt,'black','LineWidth',4);
tt=[Port(4,1); Port(5,1)];
plot(XX,tt,'white','LineWidth',3.5);
tt=[Port(5,1); Port(6,1)];
plot(XX,tt,'black','LineWidth',4);
tt=[Port(6,1); Port(7,1)];
plot(XX,tt,'white','LineWidth',3.5);
tt=[Port(7,1); Port(8,1)];
plot(XX,tt,'black','LineWidth',4);
tt=[Port(8,1); Port(9,1)];
plot(XX,tt,'white','LineWidth',3.5);
tt=[Port(9,1); Port(10,1)];
plot(XX,tt,'black','LineWidth',4);

XX=[1.01;1.01];
tt=[Port(11,1); Port(12,1)];
plot(XX,tt,'black','LineWidth',4);
tt=[Port(12,1); Port(13,1)];
plot(XX,tt,'white','LineWidth',3.5);
tt=[Port(13,1); Port(14,1)];
plot(XX,tt,'black','LineWidth',4);
tt=[Port(14,1); Port(15,1)];
plot(XX,tt,'white','LineWidth',3.5);
tt=[Port(15,1); Port(16,1)];
plot(XX,tt,'black','LineWidth',4);
tt=[Port(16,1); Port(17,1)];
plot(XX,tt,'white','LineWidth',3.5);
tt=[Port(17,1); Port(18,1)];
plot(XX,tt,'black','LineWidth',4);
tt=[Port(18,1); Port(19,1)];
plot(XX,tt,'white','LineWidth',3.5);
tt=[Port(19,1); Port(10,1)];
plot(XX,tt,'black','LineWidth',4);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DWG(2)=figure(2);
% 
% subplot(131);
% colormap(cmmap);
% %colormap(gray);
% pcolor(X,t,Z'./1000);
% %pcolor(X,t,Z'./1);
% shading flat;
% % caxis([0 1]);
% pos3 = get(gca,'Position');
% set(gca,'Position',pos3 + [0.029 0 -0.05 0])
% set(gca,'YTick',[-60:30:300])
% colorbar;
% %klb3=colorbar('position',[0.916 0.2 0.012 0.6]);
% xlabel('Position/Length');
% title('Fuel Concentration');
% hold
% XX=[-0.01;-0.01];
% tt=[Port(1,1); Port(2,1)+9];
% plot(XX,tt,'black','LineWidth',4);
% tt=[Port(2,1)+9; Port(3,1)+9];
% plot(XX,tt,'white','LineWidth',3.5);
% tt=[Port(3,1)+9; Port(4,1)+9];
% plot(XX,tt,'black','LineWidth',4);
% 
% XX=[1.01;1.01];
% tt=[Port(5,1); Port(6,1)+9];
% plot(XX,tt,'black','LineWidth',4);
% tt=[Port(6,1)+9; Port(7,1)+9];
% plot(XX,tt,'white','LineWidth',3.5);
% tt=[Port(7,1)+9; Port(12,1)+9];
% plot(XX,tt,'black','LineWidth',4);
% tt=[Port(8,1)+9; Port(9,1)+9];
% plot(XX,tt,'<black','LineWidth',2);
% % tt=[Port(9,1)+9; Port(10,1)+9];
% plot(XX,tt,'black','LineWidth',1);
% tt=[Port(10,1)+9; Port(11,1)+9];
% plot(XX,tt,'white','LineWidth',3);
% tt=[Port(11,1)+9; Port(12,1)+9];
% plot(XX,tt,'black','LineWidth',1);
% 
% subplot (132);
% colormap(cmmap);
% %colormap(gray);
% pcolor(X,t,ZI'./1000);
% %pcolor(X,t,ZI'./1);
% shading flat;
% % caxis([0 1]);
% pos3 = get(gca,'Position');
% set(gca,'Position',pos3 + [0.029 0 -0.05 0])
% set(gca,'YTick',[-60:30:300])
% colorbar;
% %klb3=colorbar('position',[0.916 0.2 0.012 0.6]);
% xlabel('Position/Length');
% title('Intermediate Concentration');
% hold
% XX=[-0.01;-0.01];
% tt=[Port(1,1); Port(2,1)+9];
% plot(XX,tt,'black','LineWidth',4);
% tt=[Port(2,1)+9; Port(3,1)+9];
% plot(XX,tt,'white','LineWidth',3.5);
% tt=[Port(3,1)+9; Port(4,1)+9];
% plot(XX,tt,'black','LineWidth',4);
% 
% XX=[1.01;1.01];
% tt=[Port(5,1); Port(6,1)+9];
% plot(XX,tt,'black','LineWidth',4);
% tt=[Port(6,1)+9; Port(7,1)+9];
% plot(XX,tt,'white','LineWidth',3.5);
% tt=[Port(7,1)+9; Port(12,1)+9];
% plot(XX,tt,'black','LineWidth',4);
% tt=[Port(8,1)+9; Port(9,1)+9];
% plot(XX,tt,'<black','LineWidth',2);
% % tt=[Port(9,1)+9; Port(10,1)+9];
% plot(XX,tt,'black','LineWidth',1);
% tt=[Port(10,1)+9; Port(11,1)+9];
% plot(XX,tt,'white','LineWidth',3);
% tt=[Port(11,1)+9; Port(12,1)+9];
% plot(XX,tt,'black','LineWidth',1);
% 
% subplot(133)
% colormap(cmmap);
% %colormap(gray);
% pcolor(X,t,ZO'./1000);
% %pcolor(X,t,ZO'./1);
% shading flat;
% caxis([0 1]);
% pos3 = get(gca,'Position');
% set(gca,'Position',pos3 + [0.029 0 -0.05 0])
% set(gca,'YTick',[-60:30:300])
% %colorbar;
% klb3=colorbar('position',[0.916 0.2 0.012 0.6]);
% xlabel('Position/Length');
% title('Oxidant Concentration');
% hold
% XX=[-0.01;-0.01];
% tt=[Port(1,1); Port(2,1)+9];
% plot(XX,tt,'black','LineWidth',4);
% tt=[Port(2,1)+9; Port(3,1)+9];
% plot(XX,tt,'white','LineWidth',3.5);
% tt=[Port(3,1)+9; Port(4,1)+9];
% plot(XX,tt,'black','LineWidth',4);
% 
% XX=[1.01;1.01];
% tt=[Port(5,1); Port(6,1)+9];
% plot(XX,tt,'black','LineWidth',4);
% tt=[Port(6,1)+9; Port(7,1)+9];
% plot(XX,tt,'white','LineWidth',3.5);
% tt=[Port(7,1)+9; Port(12,1)+9];
% plot(XX,tt,'black','LineWidth',4);
% tt=[Port(8,1)+9; Port(9,1)+9];
% plot(XX,tt,'<black','LineWidth',2);
% % tt=[Port(9,1)+9; Port(10,1)+9];
% plot(XX,tt,'black','LineWidth',1);
% tt=[Port(10,1)+9; Port(11,1)+9];
% plot(XX,tt,'white','LineWidth',3);
% tt=[Port(11,1)+9; Port(12,1)+9];
% plot(XX,tt,'black','LineWidth',1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %% Pressure Traces
% cntr=0;
% for j=1:length(t)
%     if t(j) < 0;
%         cntr=cntr+1;
%         t(j)=t(j)+360;
%     end
% end
% time = [t(:,cntr+1:length(t)),t(:,1:cntr)];
% Press = [PR(:,cntr+1:length(t)),PR(:,1:cntr)].*Pref/1000;
% NN=3;
% figure(NN)
% plot(time,Press(17,:),'-k','LineWidth',2)
% grid on
% title('Pressure trace at PT2 (Passage 16)')
% set(gca,'XTick',[0:30:360])
% xlabel('Angular position (degrees)')
% ylabel('Pressure (psia)')
% NN=NN+1;
% figure(NN)
% plot(time,Press(25,:),'-k','LineWidth',2)
% grid on
% set(gca,'XTick',[0:30:360])
% title('Pressure trace at PT3 (Passage 16)')
% xlabel('Angular position (degrees)')
% ylabel('Pressure (psia)')
% NN=NN+1;
% figure(NN)
% plot(time,Press(48,:),'-k','LineWidth',2)
% grid on
% set(gca,'XTick',[0:30:360])
% title('Pressure trace at PT5 (Passage 16)')
% xlabel('Angular position (degrees)')
% ylabel('Pressure (psia)')
% NN=NN+1;
% figure(NN)
% plot(time,Press(78,:),'-k','LineWidth',2)
% grid on
% set(gca,'XTick',[0:30:360])
% title('Pressure trace at PT8 (Passage 16)')
% xlabel('Angular position (degrees)')
% ylabel('Pressure (psia)')

% figure(6)
% for b=550:20:900
%     plot(X,Z(:,b)./1000,'-r','LineWidth',2)
%     grid on
%     xlabel('Position x/L')
%     ylabel('Fuel concentration')
%     hold on
%     pause
% end

% PRESSURE AND TEMPERATURE LINE PLOTS - PAWAN 07/14/2020
%DWG(2)=figure(2);
%subplot(131);
%plot(VEL(1,:)./1000./((TEM(1,:)./1000).^0.5),t,'b',VEL(end,:)./1000./((TEM(end,:)./1000).^0.5),t,'r--','Linewidth',1.5);
%xlabel('Mach No.');
%legend('Left','Right');
%title('Mach No.');
%grid;

% Following subplot (Vel and Tem) has been added just for pocket examination purpose.
% Should be commented out for pocketless config. Above Mach No. plot should
% be used instead. -Pawan
DWG(2)=figure(2);
subplot(131);
plot((VEL(1,:)./1000),t,'b',(TEM(1,:)),t,'r:','Linewidth',1.5);
legend('LWall Vel','LWall Tem');
title('Vel and Tem');
grid;

subplot(132);
plot((TEM(1,:)./1),t,'b', (TEM(ncells+2,:)./1),t,'r--','Linewidth',1.5);
xlabel('Temperature');
legend('Left','Right');
title('Temperature');
grid;

subplot(133);
plot((PR(1,:)./1000),t,'b', (PR(ncells+2,:)./1000),t,'r--','Linewidth',1.5);
xlabel('Pressure');
legend('Left','Right');
title('Pressure');
grid;

%% End of Script
hold off
status=fclose(fid);
