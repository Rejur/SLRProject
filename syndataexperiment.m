

for i = 1:1
%% Generate data
W = Generate_data(1, i * 100, 1.1);


XinDPCP = DPCP(W);

%XinDPCPn = DPCP_Normal(W);

XinRANSAC = RANSAC(W);

XinAPG = APG(W);

%%


end