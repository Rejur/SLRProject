function [Xin] = RANSAC(W)

len = size(W, 1);
Xin = cell(1, len);
sigma = 15;
iter = 1000;
for k = 1:len
    tmpWi = W(k, :, :);
    Nc = size(W, 3);
    Wi = zeros(8, Nc);
    Wi(:,:) = tmpWi(1, :, :);
    Nc = size(W, 3);
    tmpXin = [];
    maxtotal = 0;
    for i = 1:iter
        idx = randperm(Nc, 6);
        sample = Wi(:, idx);
        [U, S, V] = svd(sample);
        b = U(:, 8 - 1:8);


        distance = zeros(1, Nc);
        total = 0;
        for j = 1:Nc
            distance(j) = norm(b' * Wi(:, j));
            if distance < sigma
                total = total + 1;
            end
        end
        
        if maxtotal < total
            maxtotal = total;
            maxb = b;
        end
    end

    for j = 1:Nc
        distance(j) = norm(b'* Wi(:, j));
        if distance(j) < sigma
            tmpXin = [tmpXin, Wi(:, j)];
        end
    end

    %disp(distance);
    %display(maxtotal);
    display(233);
    display(size(tmpXin));
    Xin{i} = tmpXin;
end

%display(Xin);

end