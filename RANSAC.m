function [Xin, precision, recal] = RANSAC(W, p)

len = size(W, 1);
Xin = cell(1, len);
sigma = 15;
iter = 1000;
precision = 0;
recal = 0;
for k = 1:len
    tmpWi = W(k, :, :);
    Nc = size(W, 3);
    fenjie = Nc * p;
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
    TP = 0;
    for j = 1:Nc
        distance(j) = norm(b'* Wi(:, j));
        if distance(j) < sigma
            if j < fenjie
                TP = TP + 1;
            end
            tmpXin = [tmpXin, Wi(:, j)];
        end
    end
    precision = precision + (TP) / size(tmpXin, 2);
    recal = recal + (TP) / (fenjie - 1);
    %disp(distance);
    %display(maxtotal);
    %display(233);
    %display(size(tmpXin));
    Xin{i} = tmpXin;
end
precision = precision / len;
recal = recal / len;
%display(Xin);

end