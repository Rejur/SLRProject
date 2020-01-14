function [W] = Generate_data(Nm, Nc, P)
    INF = 10000;
    Pi = [1 0 0 0; 0 1 0 0; 0 0 1 0];
    K = [2 0 2 0; 0 2 2 0; 0 0 1 0; 0 0 0 1];
    Il = [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1];
    Ir = [1 0 0 -2; 0 1 0 0; 0 0 1 0; 0 0 0 1];
    
    F = [0 0 0; 0 0 1; 0 -1 0];
    W = zeros(Nm, 8, Nc);
    %display(Nc * P);
    for j = 1:Nm
        M = [1 0 0 j; 0 1 0 j; 0 0 1 j; 0 0 0 1];
        Wt = [];
        for i = 1:Nc
            if i < Nc * P
                delta = 1.5;
            else
                delta = 100; 
            end
            % display(delta);
            X = [randi(INF, [3, 1]); 1];
            tmp = Pi * K * Il * X;
            tmp = tmp ./ tmp(3);
            xl = tmp(1:2,1) + normrnd(0, delta, 2, 1);
            tmp = Pi * K * Ir * X;
            tmp = tmp ./ tmp(3);
            xr = tmp(1:2,1) + normrnd(0, delta, 2, 1);
            %display([xr;1]'*F*[xl;1]);
            %display(xr);
            %display(xl);
            tmp = Pi * K * Il * M * X;
            tmp = tmp ./ tmp(3);
            xnl = tmp(1:2,1) + normrnd(0, delta, 2, 1);
            tmp = Pi * K * Ir * M * X;
            tmp = tmp ./ tmp(3);
            xnr = tmp(1:2,1) + normrnd(0, delta, 2, 1);
            %display(xnr);
            %display(xnl);
            %display([xnr;1]'*F*[xnl;1]);
            Wt = [Wt, [xl; xr; xnl; xnr]];
        end
        % display(size(Wt));
        % display(size(W(j, :, :)));
        W(j, :, :) = Wt;
    end
    % display(size(W));
end