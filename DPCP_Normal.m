function [Xin] = DPCP_Normal(W)
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
    X = Wi;
    [f, distance, time] = DPCP_IRLS_modified(X, delta, -1, epsilon_J, 2);
    %% plot
    %display(normc(f));
    pX = 1:Nc;
    %display(pX);
    %display(size(normc(f)'*X));
    figure; subplot(1,1,1); stem(pX, distance');
    title('normal-subspace distance for each embedding to Span(h)^\perp');
    %%
    tmpXin = [];
    %display(size(distance));
    %display((distance));
   % display(f);
    % display(distance(1));
    iSum = 0;
    for x = 1:500
        iSum = iSum + norm(X(:,x)' * f);
    end
    display(iSum);
    display(iSum / 1000);
    iSum = 0;
    for x = 500:1000
        iSum = iSum + norm(X(:,x)' * f);
    end
    display(iSum);
    display(iSum / 1000);
    for j = 1:Nc
        if distance(j) <= 1.5
            tmpXin = [tmpXin , Wi(:,j)];
        end
    end
    %display(size(tmpXin));
    Xin{1} = tmpXin;
end
    display(Xin);



end