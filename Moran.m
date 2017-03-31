function [Moran_I, Z_I] = Moran(y, Adj)
%% Implementing Moran's Test
% Input:    y is the vector of observations at j zones
%           Adj is the weight matrix of spatial effect
% Output:   Moran_I is the value of Moran's I
%           Z_I is the Z-statistic of I
% Use:      [Moran_I, Z_I] = Moran(y, A)
% Author:   Sun Chenshuo
% Link:     http://www.sunchenshuo.com/as-a-phd-student/expertise/methodology/spatial-autocorrelation/

%% Calculate Moran' I
    ym = mean(y);
    n = length(Adj);
    MIa = 0;
    MIb1 = 0;
    MIb2 = 0;
    for i = 1:n
        for j = 1:n
            MIa = MIa + Adj(i, j) * (y(i) - ym) * (y(j) - ym);
            MIb1 = MIb1 + Adj(i, j);
        end
    end
    for i = 1:n
        MIb2 = MIb2 + (y(i) - ym)^2;
    end
    MIb2 = MIb2 / n;
    Moran_I = MIa / (MIb1 * MIb2);
    clear MIa MIb1 MIb2 i j

%% Calculate Z(I)
    S0 = 0;
    S1 = 0;
    S2 = 0;
    S = 0;
    k = 0;
    for i = 1:n
        S = S + (y(i) - ym)^2;
        k = k + (y(i) - ym)^4;
    end
    k = k / n;
    S = S / n;
    k = k / S^2;
    for i = 1:n
        for j = 1:n
            S0 = S0 + Adj(i, j);
            S1 = S1 + 1/2 * (Adj(i, j) + Adj(j, i))^2;
        end
        S2 = S2 + (sum(Adj(i,:)) + sum(Adj(:,i)))^2;
    end
    
    EI = -1/(n - 1);
    VarI = (n*((n^2 - 3*n  + 3) * S1 - n * S2 + 3*S0)-k*((n^2 - n) * S1 - 2*n*S2 + 6*S0^2))/((n - 1)*(n - 2)*(n - 3)*S0^2) - EI^2;

    Z_I = abs((Moran_I - EI)/sqrt(VarI)); 
