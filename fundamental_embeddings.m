function [X_tilde, X_tildeSF, X_tilde_1, X_tilde_2] = fundamental_embeddings(xt1, xt2, xnt1, xnt2)
    len = size(xt1, 2);
    X_tilde_1 = zeros([8, len]);
    X_tilde_2 = zeros([8, len]);
    X_tilde = zeros([16, len]);
    X_tildeSF = zeros([8, 2 * len]);

    
    for index = 1:len
        X_tilde_1(:, index) = [xt1(1, index) * xt2(1, index); xt1(2, index) * xt2(1, index); xt2(1, index); xt1(1, index) * xt2(2, index); xt1(2, index) * xt2(2, index); xt2(2, index); xt1(1, index); xt1(2, index)];
        X_tilde_2(:, index) = [xnt1(1, index) * xnt2(1, index); xnt1(2, index) * xnt2(1, index); xnt2(1, index); xnt1(1, index) * xnt2(2, index); xnt1(2, index) * xnt2(2, index); xnt2(2, index); xnt1(1, index); xnt1(2, index)];
        X_tilde(:, index) = [X_tilde_1(:, index) ; X_tilde_2(:, index)];
        X_tildeSF(:, 2 * (index - 1) + 1) = X_tilde_1(:, index);
        X_tildeSF(:, 2 * (index - 1) + 2) = X_tilde_2(:, index);
    end
end
