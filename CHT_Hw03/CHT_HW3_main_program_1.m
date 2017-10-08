%ME8833KU - Diego Vaca - Main program
%Homework 3
%Discretization: Finite Volume
% Problem 1.
%%input of parameters
%Since the eqn is non-dimensionalized, k and alpha are not present.
clear % Clear variables and functions from memory
clc
L=1;%Normalized size of the domain in x 
M=1;% Normalized size of the domain in y.
Nx=41; %grid size in x. 
%For this homework, the number of cells must be odd since the information of a cell in the center is requested. 
%If grid size is even, the temperature in the center point should be calculated as the average of the 4 neighbors cells.
Ny=41;
dx=L/Nx; %Distance between grid points in x
dy=M/Ny; %Distance between grid points in y
dxe=dx; %for this case, grid is constant
dxw=dx;
dyn=dy;
dys=dy;
dxb=dx/2;
dyb=dy/2;
%boundary conditions. T represent theta in the non-dimensional eqn.
Tleft=0; 
Tright=0;
Ttop=1;
Tbottom=0;
dt=0.01; %dt is time step. Represent dtau in non-dimensional eqn.
ap0=dx*dy/dt;
%%
%Generation of initial field.
%T0 represents field at t=0. T0 also stores the converged field for each
%time step.
T0{1}=zeros(Nx,Ny);
%%
%Generation of vectors with coefficients.
tolerance=1; %tolerance to check it steady state is reached.
k=1;
while tolerance>1e-5
convergence=1;
Tstar=zeros(Nx,Ny); %Tstar is the initial guess for each time step.
while convergence>1e-5
    for i=1:Ny %i iterates rows
        for j=1:Nx %j iterates columns
            if i==1
                if j==1
                    abe(j)=0; %Use of left and top boundaries. Constant T
                    abw(j)=dy/dxb;
                    abn(j)=dx/dyb;
                    abst(j)=0;
                    ae(j)=dy/dxe;
                    aw(j)=0;
                    an(j)=0;
                    as(j)=dx/dys;
                    ap(j)=ae(j)+aw(j)+an(j)+as(j)+abn(j)+abst(j)+abe(j)+abw(j)+ap0; %abst, abn, etc are boundary coefficients.
                    b(j)=as(j)*Tstar(i+1,j)+abw(j)*Tleft+abn(j)*Ttop+ap0*T0{k}(i,j); 
                 elseif j==Nx %Use of right and top boundaries. Constant T
                    abe(j)=dy/dxb;
                    abw(j)=0;
                    abn(j)=dx/dyb;
                    abst(j)=0;
                    ae(j)=0;
                    aw(j)=dy/dxw;
                    an(j)=0;
                    as(j)=dx/dys;
                    ap(j)=ae(j)+aw(j)+an(j)+as(j)+abn(j)+abst(j)+abe(j)+abw(j)+ap0;
                    b(j)=as(j)*Tstar(i+1,j)+abe(j)*Tright+abn(j)*Ttop+ap0*T0{k}(i,j);
                else %Use of top boundary. Constant T
                    abe(j)=0;
                    abw(j)=0;
                    abn(j)=dx/dyb;
                    abst(j)=0;
                    ae(j)=dy/dxe;
                    aw(j)=dy/dxw;
                    an(j)=0;
                    as(j)=dx/dys;
                    ap(j)=ae(j)+aw(j)+an(j)+as(j)+abn(j)+abst(j)+abe(j)+abw(j)+ap0;
                    b(j)=as(j)*Tstar(i+1,j)+abn(j)*Ttop+ap0*T0{k}(i,j);
                end
            elseif i~=1 && i~=Ny %Coefficient for cells between top and bottom boundaries.
                if j==1 %Use of left boundary. Constant T
                    abe(j)=0;
                    abw(j)=dy/dxb;
                    abn(j)=0;
                    abst(j)=0;
                    ae(j)=dy/dxe;
                    aw(j)=0;
                    an(j)=dx/dyn;
                    as(j)=dx/dys;
                    ap(j)=ae(j)+aw(j)+an(j)+as(j)+abn(j)+abst(j)+abe(j)+abw(j)+ap0;
                    b(j)=an(j)*Tstar(i-1,j)+as(j)*Tstar(i+1,j)+abw(j)*Tleft+ap0*T0{k}(i,j);
                elseif j==Nx %Use of right boundary. Constant T
                    abe(j)=dy/dxb;
                    abw(j)=0;
                    abn(j)=0;
                    abst(j)=0;
                    ae(j)=0;
                    aw(j)=dy/dxw;
                    an(j)=dx/dyn;
                    as(j)=dx/dys;
                    ap(j)=ae(j)+aw(j)+an(j)+as(j)+abn(j)+abst(j)+abe(j)+abw(j)+ap0;
                    b(j)=an(j)*Tstar(i-1,j)+as(j)*Tstar(i+1,j)+abe(j)*Tright+ap0*T0{k}(i,j);
                else 
                    abe(j)=0;
                    abw(j)=0;
                    abn(j)=0;
                    abst(j)=0;
                    ae(j)=dy/dxe;
                    aw(j)=dy/dxw;
                    an(j)=dx/dyn;
                    as(j)=dx/dys;
                    ap(j)=ae(j)+aw(j)+an(j)+as(j)+abn(j)+abst(j)+abe(j)+abw(j)+ap0;
                    b(j)=an(j)*Tstar(i-1,j)+as(j)*Tstar(i+1,j)+ap0*T0{k}(i,j);
                end
            elseif i==Ny %Use of left and bottom boundaries. Constant T
                if j==1
                    abe(j)=0;
                    abw(j)=dy/dxb;
                    abn(j)=0;
                    abst(j)=dx/dyb;
                    ae(j)=dy/dxe;
                    aw(j)=0;
                    an(j)=dx/dyn;
                    as(j)=0;
                    ap(j)=ae(j)+aw(j)+an(j)+as(j)+abn(j)+abst(j)+abe(j)+abw(j)+ap0;
                    b(j)=an(j)*Tstar(Ny-1,j)+abw(j)*Tleft+abst(j)*Tbottom+ap0*T0{k}(i,j);
                elseif j==Nx %Use of right and bottom boundaries. Constant T
                    abe(j)=dy/dxb;
                    abw(j)=0;
                    abn(j)=0;
                    abst(j)=dx/dyb;
                    ae(j)=0;
                    aw(j)=dy/dxw;
                    an(j)=dx/dyn;
                    as(j)=0;
                    ap(j)=ae(j)+aw(j)+an(j)+as(j)+abn(j)+abst(j)+abe(j)+abw(j)+ap0;
                    b(j)=an(j)*Tstar(Ny-1,j)+abe(j)*Tright+abst(j)*Tbottom+ap0*T0{k}(i,j);
                else %Use of bottom boundary. Constant T
                    abe(j)=0;
                    abw(j)=0;
                    abn(j)=0;
                    abst(j)=dx/dyb;
                    ae(j)=dy/dxe;
                    aw(j)=dy/dxw;
                    an(j)=dx/dyn;
                    as(j)=0;
                    ap(j)=ae(j)+aw(j)+an(j)+as(j)+abn(j)+abst(j)+abe(j)+abw(j)+ap0;
                    b(j)=an(j)*Tstar(Ny-1,j)+abst(j)*Tbottom+ap0*T0{k}(i,j);
                end    
            end        
        end
[T]=TDMA1(ae,aw,ap,b,Nx); %calls TDMA function to solve eqns.
    for j=1:Nx %Creation of matrix using vectors coming from TDMA1 function.
        Tmatrix(i,j)=T(j);
    end
    end
    convergence=max(abs(Tmatrix-Tstar));
   Tstar=Tmatrix; %Update of guess matrix
end
T0{k+1}=Tmatrix; %Store converged matrix for each time step.
tolerance=max(abs(T0{k+1}-T0{k})); %tolerance to reach steady state
k=k+1;
end

%%
%Plots
cell_c_row=fix((Ny+1)/2);
cell_c_column=fix((Nx+1)/2);
for k=1:size(T0,2)
    Tcenter(k)=T0{k}(cell_c_column,cell_c_row);
    time(k)=(k-1)*dt;
end
figure()
plot(time,Tcenter)
legend('Temperature')
xlabel('\tau')
ylabel('\theta')
str1 = sprintf('Variation of \\theta with \\tau in central cell (row %d, file %d)',cell_c_row,cell_c_column);
title(str1)

%Vector for y coordinates.
for i=1:Ny+2 %yc indicates coordinates of centroids
    if i==1
        yc(i)=0; %top wall
    elseif i==2;
        yc(i)=dy/2;
    elseif i==Ny+2;
        yc(i)=1;
    else
        yc(i)=yc(i-1)+dy;
    end
end

for k=1:size(T0,2)
    for i=1:Ny+2
        if i==1
            Tplot{k}(i)=Ttop;
        elseif i==Ny+2
            Tplot{k}(i)=Tbottom;
        else
            Tplot{k}(i)=T0{k}(i-1,cell_c_column);
        end
    end
end

figure()
plot(yc,Tplot{2},'-d')
hold on
plot(yc,Tplot{4},'-*')
hold on
plot(yc,Tplot{8},'-o')
hold on
plot(yc,Tplot{12},'-s')
hold on
plot(yc,Tplot{16},'-^')
hold on
plot(yc,Tplot{size(T0,2)})
t1=2*dt;t2=4*dt;t3=8*dt;t4=12*dt;t5=16*dt;t6=size(T0,2)*dt;
legend(['\tau=',num2str(t1)],['\tau=',num2str(t2)],['\tau=',num2str(t3)],['\tau=',num2str(t4)],['\tau=',num2str(t5)],['Stady state. \tau=',num2str(t6)])
xlabel('y. y = 0 is top wall')
ylabel('\theta')
title('Variation of \theta with y at central line')