clear all
clc

d = 3; %Physical dimension of the states in the ensemble
r = 3; %Dimensionality of the classical simulation
m = 6; %Number of states in the ensemble

id = eye(d); %Idenity in physical dimension
omega = exp(2*pi*i/d);

nb = 3; %Number of outcomes for each measurement
ny = 2; %Number of measurements

%MUB measurements

%Computational basis
for l = 0 : d-1
    mub{1,l+1} = id(:,l+1);
end

%Other MUBs
for j = 2 : d+1
    for l = 0 : d-1
        mub{j,l+1} = 0;
        for k = 0 : d-1
            mub{j,l+1} = mub{j,l+1} + 1/sqrt(d)*omega^(k*(l + (j-2)*k))*id(:,k+1);
        end
    end
end

%Measurement operators
for b = 1 : nb
    for y = 1 : ny
        M{b,y} = mub{y,b}*mub{y,b}';
    end
end

%Coefficients for the witness operator
for x = 1 : m
    L = toSeveralBases(x-1,[nb ny])+1;
    for b = 1 : nb
        for y = 1 : ny
            coef(b,x,y) = [L(1) == b]*[L(2) == y];
        end
    end
end

%Computation of an upper bound on the witness
witmax = classical_witness(d,r,m,ny,nb,coef,M)
