function [W] = Generate_data(Nm, Nc)
    INF = 1000000000;
    Pi = [1 0 0 0; 0 1 0 0; 0 0 1 0];
    K = [2 0 2 0; 0 2 2 0; 0 0 1 0; 0 0 0 1];
    Il = [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1];
    Ir = [1 0 0 -2; 0 1 0 0; 0 0 1 0; 0 0 0 1];
    
    
    W = [];
    for j = 1:Nm
        M = [1 0 0 j; 0 1 0 j; 0 0 1 j; 0 0 0 1];
        Wt = [];
        for i = 1:Nc
            X = randi(INF, [3, 1]);
            tmp = Pi * K * Il * X;
            tmp = tmp ./ tmp(3);
            xl = tmp(1:2,1);
            tmp = Pi * K * Ir * X;
            tmp = tmp ./ tmp(3);
            xr = tmp(1:2,1);
            tmp = Pi * K * Il * M * X;
            tmp = tmp ./ tmp(3);
            xnl = tmp(1:2,1);
            tmp = Pi * K * Ir * M * X;
            xnr = tmp(1:2,1);
            Wt = [Wt, [xl; xr; xnl; xnr]];
        end
        W = [W; Wt];
    end
end