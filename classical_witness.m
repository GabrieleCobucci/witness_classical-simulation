function witmax = classical_witness(d,r,m,ny,nb,coef,M)

%-------------------------------------------------------------------------%
%This function computes an upper bound on the witness for classical
%simulability of a quantum ensemble,
%i.e. W = \sum_{b,x,y} c_{b,x,y} Tr(\rho_x M_{b|y}).

%Inputs:
% - d: physical dimension of the states in the ensemble
% - r: dimensionality of the classical simulation
% - m: number of states in the ensemble
% - ny: number of the possible measurements to perform
% - nb: number of the possible outputs for each measurement
% - coef(b,x,y): coefficients c_{b,x,y} of the witness operator
% - M{b,y}: measurement operators M_{b|y} of the witness operator

%Output:
% - witmax: upper bound of the witness operator over all deterministic
% strategies for an r-dimensional classically simulable ensemble
%-------------------------------------------------------------------------%

%Identity in physical dimension
id = eye(d);

%Number of deterministic strategies
ndet = r^m;

%Deterministic strategies
for l = 0 : ndet-1
    ds{l+1} = [toSeveralBases(l,r*ones(1,m))+1];
end

%Ox operators
for x = 1 : m
    O{x} = 0;
    for y = 1 : ny
        for b = 1 : nb
            O{x} = O{x} + coef(b,x,y)*M{b,y};
        end
    end
end

wit = 0;
for mu = 1 : ndet

    %Deterministic strategies
    ds = toSeveralBases(mu-1,r*ones(1,m))+1;

    %Constraints
    C = [];

    %Variables
    for j = 1 : r
        E{j} = sdpvar(d,d,'hermitian','complex');
        C = [C, E{j} >= 0, trace(E{j}) == 1];
    end

    %Orthonormality relaxation
    s = 0;
    for j = 1 : r
        s = s + E{j};
    end

    C = [C, s <= id];

    %Objective function (witness)
    w = 0;
    for i = 1 : r
        for x = 1 : m
            w = w + [ds(x) == i]*trace(O{x}*E{i});
        end
    end

    %SolveSDP
    disp('Options')
    ops=sdpsettings('solver','sedumi', 'cachesolvers', 1);
    diagnostic=solvesdp(C,-real(w),ops)

    wit(mu)=double(w);

    %Maximum over deterministic strategies
    witmax = max(wit)

end