
figure;
p = 0.4;
x = [];
y = [];
for j = 1:5
    Nc = [100, 200, 400, 700, 1000]
    psDPCP = 0; psDPCPN = 0;psRANSAC = 0; psAPG = 0;
    rsDPCP = 0; rsDPCPN = 0; rsRANSAC = 0; rsAPG = 0;
    parfor i = 1:30
    %% Generate data
    W = Generate_data(2, Nc(j), p);
    
    disp("Execute DPCP");
    tic;
    [XinDPCP, precisionDPCP, recalDPCP] = DPCP(W, p);
    DPCPt = toc;
    disp("Execute DPCPN");   
    tic;
    [XinDPCPn, precisionDPCPN, recalDPCPN] = DPCP_Normal(W, p);
    DPCPNt = toc;
    disp("Execute RANSAC"); 
    tic;
    [XinRANSAC, precisionRANSAC, recalRANSAC] = RANSAC(W, p);
    RANSACt = toc;
    disp("Execute APG"); 
    tic;
    [XinAPG, precisionAPG, recalAPG] = APG(W, p);
    APGt = toc;
    psDPCP = psDPCP + precisionDPCP / 30;
    psDPCPN = psDPCPN + precisionDPCPN / 30;
    psRANSAC = psRANSAC + precisionRANSAC / 30;
    psAPG = psAPG + precisionAPG / 30;
    rsDPCP = rsDPCP + recalDPCP / 30;
    rsDPCPN = rsDPCPN + recalDPCPN / 30;
    rsRANSAC = rsRANSAC + recalRANSAC / 30;
    rsAPG = rsAPG + recalAPG / 30;
    %%
    %% Estimation M
    %B = [2, 0, 0];
    %disp(XinDPCP);
    %M = Fast_Compressed_Least_Squares(XinDPCP{1}, 2, 2, 2, B);
    %disp(M);

    end
    %% plot figure
    x = [x;psDPCP, psDPCPN, psRANSAC, psAPG];
    colormap(cool);
    y = [y;rsDPCP, rsDPCPN, rsRANSAC, rsAPG];
    colormap(cool);
    subplot(1,2,1);
    bar(x);
    subplot(1,2,2);
    bar(y);
    colormap(cool);
end