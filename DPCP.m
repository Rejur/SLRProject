function [Xin, precision, recal] = DPCP(W, p)
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
    
    [trash, X, trash, trash] = fundamental_embeddings(Wi(1:2,:), Wi(3:4,:), Wi(5:6,:), Wi(7:8,:));
    [f, distance, time] = DPCP_IRLS_modified(X, delta, -1, epsilon_J, 1);
    tSum = 0.0;
    for x = 1:size(X,2)
        tSum = tSum + norm(X(:,x)' * f);
    end
    %display(tSum);
    %% plot
    %figure; subplot(1,1,1); stem(abs(normc(f)'*X));
%title('fundenmental-subspace distance for each embedding to Span(h)^\perp');
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
    %display(lq);
    TP = 0;
    for j = 1:Nc
        if qwe(2 * (j - 1) + 1) <= lq && qwe(2 * (j - 1) + 2) <= lq
            if j < fenjie
                TP = TP + 1;
            end
            tmpXin = [tmpXin , Wi(:,j)];
        end
    end
    precision = precision + (TP) / size(tmpXin, 2);
    recal = recal + (TP) / (fenjie - 1);
    %display(size(tmpXin));
    Xin{i} = tmpXin;
end
precision = precision / len;
recal = recal / len;
    %display(Xin);

end