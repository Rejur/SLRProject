function [Xin] = DPCP(W)
%% init parameters.
%display(size(W));
len = size(W, 1);
delta = 10^(-9);
T = 1000;
epsilon_J = 10^(-6);
Xin = cell(1, len);
yf = [0; 0; 0; 0; 0; 1; 0; -1; 0];
for i = 1:len
    tmpWi = W(i, :, :);
    Nc = size(W, 3);
    Wi = zeros(8, Nc);
    Wi(:,:) = tmpWi(1, :, :);
    
    [trash, X, trash, trash] = fundamental_embeddings(Wi(1:2,:), Wi(3:4,:), Wi(5:6,:), Wi(7:8,:));
    [f, distance, time] = DPCP_IRLS_modified(X, delta, -1, epsilon_J, 1);
    tSum = 0.0;
    for x = 1:size(X,2)
        tSum = tSum + norm(X(:,x)' * f);
    end
    %display(tSum);
    %% plot
    figure; subplot(1,1,1); stem(abs(normc(f)'*X));
title('fundenmental-subspace distance for each embedding to Span(h)^\perp');
    %%
    tmpXin = [];
   %display((distance));
   % display(f);
    %iSum = 0;
    %for x = 1:1000
        %iSum = iSum + norm(X(:,x)' * f);
    %end
    %display(iSum);
    %display(iSum / 1000);
    %iSum = 0;
    %for x = 1000:2000
%        iSum = iSum + norm(X(:,x)' * normc(f));
    %end
    qwe = zeros(1,2 * Nc);
    for m = 1:2 * Nc
        qwe(m) = norm(X(:,m)' * normc(f));
    end
    %display(iSum);
    %display(iSum / 1000);
    % display(distance(1));
    %display(tSum / (2 * Nc));
    lq = max(((tSum / (2 * Nc)) * 0.3), 5);
    display(lq);
    for j = 1:Nc
        if qwe(2 * (j - 1) + 1) <= lq && qwe(2 * (j - 1) + 2) <= lq
            tmpXin = [tmpXin , Wi(:,j)];
        end
    end
    %display(size(tmpXin));
    Xin{i} = tmpXin;
end
    display(Xin);

end