%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% last update 13April2018, lne %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% This program solves the time independante-Schrodinger equations
%%% in 2D on an inhomogeneous grid

clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

saveXY=0;
saveV=0;
savePSI=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n=15;                   %% number of solution asked 
Mass = 0.067;           %% effective mass, constant over all the structure...
Fx=0;%5e7;              %% Electric field [V/m] in the x-direction
Fy=0;%5e7;              %% Electric field [V/m] in the y-direction

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% Potential definition %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Two vectors (x and y) and one matrix V0 must be defined
% They can be non homogeneous (meaning: dx~=cst,dy~=cst)
% x and y [meter]
% V0 [eV]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Nx1=60;                  %% Meshing point in x-direction in the 1st domain
Nx2=10;                  %% Meshing point in x-direction in the 2nd domain
Ny1=15;                  %% Meshing point in y-direction in the 1st domain
Ny2=40;                  %% Meshing point in y-direction in the 2nd domain

Mx=20e-9;               %% map X [m]
My=20e-9;               %% map Y [m]

x1=linspace(-Mx/2,0,Nx1);
x2=linspace(x1(end)+1e-12,Mx/2,Nx2);
x=[x1 x2];

y1=linspace(-My/2,0,Ny1);
y2=linspace(y1(end)+1e-9,My/2,Ny2);
y=[y1 y2];

Nx=length(x);
Ny=length(y);

[X,Y]=meshgrid(x,y);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Choose your between the next 3 potentials or build your own!

% Pot_Rectangular
% Pot_Elliptical
 Pot_Hexagonal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Vb=1.5;                 %% Potential barrier height[eV]
V0=(idx)*0 + (1-idx)*Vb ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% NOTHING TO CHANGE ANYMORE!!! %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

V0=(Fx*X)+V0;        % adding the electric field Fx to the potential in the x-direction
V0=(Fy*Y)+V0;        % adding the electric field Fx to the potential in the y-direction
V0=V0-min(min(V0));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% Selection of the model %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
E1=[];
display('=======================================')

tic
if length(x)*length(y)>1e4
  N=length(x)*length(y);
  display(strcat('Warning: Take care, H=',num2str(N),'x',num2str(N),'elements'))
end
[E,psi] = Schroed2D_FEM_sq_f(x,y,V0,Mass,n);  % m=cste
display(strcat('-> Finite Elements method =',num2str(toc),'sec'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Display Results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

display('=======================================')
display('Results:')
display('=======================================')
display(strcat('E(eV)='))
display(strcat(num2str(E)))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% figures %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('Name','Potential','position',[100 100 1200 400])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(1,2,1)
hold on;grid on;
surf(x*1e9,y*1e9,V0)

colormap(jet)
colorbar
view(30,30)
%shading flat

xlabel('x (nm)')
ylabel('y (nm)')
zlabel('Energy (eV)')
title(strcat('Potential, Fx=',num2str(Fx,'%.1e'),'[V/m], Fy=',num2str(Fy,'%.1e'),'[V/m]'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(1,2,2)
hold on;grid on;
pcolor(x*1e9,y*1e9,V0)

colormap(jet)
colorbar
%shading flat

xlabel('x (nm)')
ylabel('y (nm)')
zlabel('Energy (eV)')
title(strcat('Potential, Fx=',num2str(Fx,'%.1e'),'[V/m], Fy=',num2str(Fy,'%.1e'),'[V/m]'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c=0;
ii=0;
for i=1:n
    if i>45
      break
    end
    if i==1 || i==16 || i==31
      figure('Name','FEM method','position',[100 100 1600 900])
%      figure('Name','FEM method','position',[-1800 100 1600 900])
      c=c+1;
      ii=0;
    end
    ii=ii+1;
    
    subplot(3,5,ii,'fontsize',10)
    hold on
    
    pcolor(x*1e9,y*1e9,(psi(:,:,i)) )  
    contour(x*1e9,y*1e9,V0,1,'linewidth',3,'linecolor','w')
    
    xlabel('x (nm)')
    ylabel('y (nm)')
    title(strcat('E',num2str(i),'=',num2str(E(i,1)*1000,'%.1f'),'meV'))
    %axis equal
    %shading flat
    colormap(jet)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% Data save %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if saveXY==1;
    x=x';y=y';
    save('data_x.txt','x','-ascii')
    save('data_y.txt','y','-ascii')
end

if saveV==1;
    save('data_V.txt','V0','-ascii')
end

if savePSI==1;
    
  for i=1:n
    M1 = psi1(:,:,i);
    save(strcat('data_psi',num2str(i),'_FEM.txt'),'M1','-ascii')
  end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%