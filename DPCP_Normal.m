function [Xin, precision, recal] = DPCP_Normal(W, p)
%% init parameters.
%display(size(W));
len = size(W, 1);
delta = 10^(-9);
T = 1000;
epsilon_J = 10^(-6);
Xin = cell(1, len);
yf = [0; 0; 0; 0; 0; 1; 0; -1; 0];
precision = 0;
recal = 0;
for i = 1:len
    tmpWi = W(i, :, :);
    Nc = size(W, 3);
    fenjie = Nc * p;
    Wi = zeros(8, Nc);
    Wi(:,:) = tmpWi(1, :, :);
    X = Wi;
    [f, distance, time] = DPCP_IRLS_modified(X, delta, -1, epsilon_J, 2);
    tSum = norm(distance, 1);
    %% plot
    %display(normc(f));
    pX = 1:Nc;
    %display(pX);
    %display(size(normc(f)'*X));
    %figure; subplot(1,1,1); stem(pX, distance');
    %title('normal-subspace distance for each embedding to Span(h)^\perp');
    %%
    tmpXin = [];
    %display(size(distance));
    %display((distance));
   % display(f);
    % display(distance(1));
    
    TP = 0;
    lq = max(((tSum / (2 * Nc)) * 0.3), 2);
    for j = 1:Nc
        if distance(j) <= lq
            if j < fenjie
                TP = TP + 1;
            end
            tmpXin = [tmpXin , Wi(:,j)];
        end
    end
    precision = precision + (TP) / size(tmpXin, 2);
    recal = recal + (TP) / (fenjie - 1);
    %display(size(tmpXin));
    Xin{1} = tmpXin;
end
    %display(Xin);

precision = precision / len;
recal = recal / len;

end