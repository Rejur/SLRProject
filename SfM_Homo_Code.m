 load('SfM_Homo_Data.mat');
% load('My_data_400_inliers.mat');
% A set of L=250 correspondences are given, i.e., x1(:, i) <-> x2(:, i) for
% i = 1, ..., L. Points are represented in homogeneous coordinates. As a
% ground-truth information, the first 100 coorespondences are inliers while
% the rest 150 of them are outliers. While you may use the ground-truth for
% debugging purpose, your output should not depend on the ground-truth
% information that which ones are inliers/outliers. We will run your code
% on different test data to ensure the correctness of your result.


% YOUR CODE: Compute a matrix X_tilde that contains embeddings for all the
% correspondences. Remember that in homography matrix estimation, each
% correspondence gives 3 embeddings. Therefore, X_tilde is expected to have
% size 9*3L.
% iL = 250;
% x1 = x1(:,1:iL);
% x2 = x2(:,1:iL);

X_tilde = homographic_embeddings(x1, x2);

p1 = null(X_tilde','r');
% display(p1);
% Normalize the embeddings to have the unit norm.
X_tilde = normc(X_tilde);

% YOUR CODE: Learn a subspace from X_tilde, using your modification of
% DPCP-IRLS as described in problem 3. 
delta = 10^(-9);
T = 1000;
epsilon_J = 10^(-6);
[distance, h, num_iter, time] = DPCP_IRLS_modified(X_tilde, delta, T, epsilon_J);

cnt = 0;
len = size(X_tilde, 2);
% display(distance);
for index = 1:len
    if distance(:,index) < 10^(-6)
        cnt = cnt + 1;
    end
end
display(cnt);
   
% Visualize the subspace distance for each embedding to Span(h)^\perp
figure; subplot(1,2,1); stem(abs(normc(h)'*X_tilde));
title('subspace distance for each embedding to Span(h)^\perp');
subplot(1,2,2);
data = [cnt, len - cnt];
b = bar(data);
ch = get(b,'children');
set(gca,'XTickLabel',{'inlier','outlier'});
title('Number of the inliers and outliers');
% For debugging purpose, you may also want to learn a subspace using the
% original DPCP-IRLS algorithm, and perform the same visualization. What do
% you find?



