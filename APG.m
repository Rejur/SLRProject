function [Xin] = APG(W)

lambda = 0.0001;
len = size(W, 1);
Xin = cell(1, len);
delta = 10^-5;

iter = 10000;
for k = 1:len
    tmpWi = W(k, :, :);
    Nc = size(W, 3);
    Wi = zeros(8, Nc);
    Wi(:,:) = tmpWi(1, :, :);
    Akp1 = zeros(8, Nc); Ak = zeros(8, Nc); Apre = zeros(8, Nc);
    Ekp1 = zeros(8, Nc); Ek = zeros(8, Nc); Epre = zeros(8, Nc);
    tkp1 = 1; tk = 1; tpre = 1;
    %[tx, ty] = eig(Wi);
    %[tu, ts, tv] = svds(Wi, 1);
    %eigenvalue=ts(1, 1) * ts(1, 1);
    miuk = 0.9 * norm(Wi, 2);
    miukp1 = 0;
    miubar = 10^-5 * miuk;
    tmpXin = [];
    %%
    while(true)
    %for i = 1:200
        %disp(Ak);
        %disp(Ek);
        YA = Ak + ((tpre - 1) / tk) .* (Ak - Apre);
        YE = Ek + ((tpre - 1) / tk) .* (Ek - Epre);
        GA = YA - (YA + YE - Wi) ./ 2;
        %disp(GA);
        [U, S, V] = svd(GA);
        %disp(S);
        %disp(Sh(S, miuk / 2));
        %disp(miuk / 2);
        Akp1 = U * Sh(S, miuk / 2) * V';
        %disp(Ak);
        GE = YE - (YA + YE - Wi) ./ 2;
        Ekp1 = Sh(GE, lambda * miuk / 2);
        tkp1 = (1 + sqrt(4 * tk * tk + 1)) / 2;
        miukp1 = max(0.9 * miuk, miubar);
        
        Apre = Ak; Ak = Akp1;
        Epre = Ek; Ek = Ekp1;
        tpre = tk; tk = tkp1;
        miuk = miukp1;
        SA = 2 * (YA - Ak) + (Ak + Ek - YA - YE);
        SE = 2 * (YE - Ek) + (Ak + Ek - YA - YE);
        %disp(norm(SA, 'fro'));
        if sqrt(norm(SA, 'fro') ^ 2 + norm(SE, 'fro')^2) < 0.001
           break; 
        end
    end
    %%
    A = Ak;
    E = Ek;
    disp(A);
    disp(Wi);
    disp(E);
    %disp(Wi);
    for i = 1:Nc
        if norm(Wi(:,i), 1) <= min(0.5, norm(E(:, i), 1) * 1.1)
           tmpXin = [tmpXin, Wi(:, i)]; 
        end
    end

    display(size(tmpXin));
    Xin{i} = tmpXin;
end

%display(Xin);

end


function [r] = Sh(S, uk)
    len1 = size(S, 1);
    len2 = size(S, 2);
    r = zeros(size(S));
    for i = 1:len1
        for j = 1:len2
            if S(i, j) > uk
                r(i, j) = S(i, j) - uk;
            elseif S(i, j) < -uk
                r(i, j) = S(i, j) + uk;
            else
                r(i, j) = 0;
            end
        end
    end
end