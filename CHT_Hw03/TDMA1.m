function [T] = TDMA1(ae,aw,ap,b,Nx)
%This function solves a system of eqns using TDMA algorithm found in
%Numerical Heat Transfer and Fluid FLow by Patankar.
%ai = ap, bi = ae, ci = aw, di = b
        for j=1:Nx %j iterates columns
            if j==1 %Calculation of P1 and Q1
                P(j)=ae(j)/ap(j);
                Q(j)=b(j)/ap(j);
            else %Forward calulation Pi and Qi from i=2 to i=N
                P(j)=ae(j)/(ap(j)-aw(j)*P(j-1));
                Q(j)=(b(j)+aw(j)*Q(j-1))/(ap(j)-aw(j)*P(j-1));
            end
        end
    T(Nx)=Q(Nx);%Assign to TN = QN
    % Back calculation Ti from TN-1 to T1
        for j=Nx-1:-1:1 
            T(j)=P(j)*T(j+1)+Q(j);
        end
end

