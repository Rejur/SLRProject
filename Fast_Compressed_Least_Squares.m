function [M, time] = Fast_Compressed_Least_Squares(X, f, cu, cv, alpha, beta, gamma)
    %%Fast Compressed Least-Squares.
    %%Get M.

    %% init parameters.
    tic;
    k = 20;
    Pi = [1, 0, 0 ; 0, 1, 0 ; 0, 0, 1 ];
    K = [f, 0, cu ; 0, f, cv; 0, 0, 1, 0];
    % invC = inv(Pi * K);
    % [f * x + cu * z; y * f + cv * z; z; 1] = [f * x / z + cu; y * f / z +
    % cv / z; 1]
    lenN = size(X, 2);
    zero3 = ones(3, 1);
    Al = [];
    Ar = [];
    for i = 1:lenN
        x_i_l = X(1:2,i);
        x_i_r = X(3:4,i);
        x_ip1_l = X(5:6,i);
        x_ip1_r = X(7:8,i);
        X_i_l = (Pi * K) \ [x_i_l; 1];
        Ail = [zero3', f * X_i_l', (cv - x_ip1_l(2)) * X_i_l', 0, f, (cv - x_ip1_l(2)), 0;
            -f * X_i_l', zero3', (x_ip1_l(1) - cu) * X_i_l', -f, 0, (x_ip1_l(1) - cu), 0;
            f * x_ip1_l(2) * X_i_l', -f * x_ip1_l(1) * X_i_l', (cu * x_ip1_l(2) - cv * x_ip1_l(1)) * X_i_l', f * x_ip1_l(2), -f * x_ip1_l(1), (cu * x_ip1_l(2) - cv * x_ip1_l(1)), 0];
        Air = [zero3', f * X_i_l', (cv - x_ip1_r(2)) * X_i_l', 0, f, (cv - x_ip1_r(2)), alpha;
            -f * X_i_l', zero3', (x_ip1_r(1) - cu) * X_i_l', -f, 0, (x_ip1_r(1) - cu), beta;
            f * x_ip1_r(2) * X_i_l', -f * x_ip1_r(1) * X_i_l', (cu * x_ip1_r(2) - cv * x_ip1_r(1)) * X_i_l', f * x_ip1_r(2), -f * x_ip1_r(1), (cu * x_ip1_r(2) - cv * x_ip1_r(1)), gamma];
        
        Al = [Al; Ail];
        Ar = [Ar; Air];
    end
    Gamma = Al' * Al + Ar' * Ar;
    fun = @(w) EM(w, Gamma);
    options = optimoptions(@lsqnonlin,'Algorithm','levenberg-marquardt');
    x0 = [0; 0; 0; 0; 0; 0];
    x = lsqnonlin(fun,x0,[],[],options);
    
    
    
    %% Do some final work.
    trash, M = EM(x, Gamma);
    time = toc;
end

function [value, EXP] = EM(w, Gamma)
    Z = [0 -w(3) -w(2) w(4); w(3) 0 -w(1) w(5); -w(2) w(1) 0 w(6); 0 0 0 1];
    theta = sqrt(w(1)^2 + w(2)^2 + w(3)^2);
    I = eye(4);
    EXP = I + Z + ((1 - cos(theta)) .* (Z^2)) ./ (theta ^ 2) + ((theta - sin(theta)) .* (Z^3)) ./ (theta ^ 3);
    tmp = [EXP(1:3, 1) ; EXP(1:3, 2); EXP(1:3, 3); EXP(1:4, 4)];
    value = tmp' * Gamma * tmp;
end
