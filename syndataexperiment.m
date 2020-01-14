
figure;
p = 0.4;
x = [];
y = [];
for j = 1:5
    Nc = [100, 200, 400, 700, 1000]
    psDPCP = 0; psDPCPN = 0;psRANSAC = 0; psAPG = 0;
    rsDPCP = 0; rsDPCPN = 0; rsRANSAC = 0; rsAPG = 0;
    for i = 1:30
    %% Generate data
    W = Generate_data(2, Nc(j), p);
    
    disp("Execute DPCP");
    [XinDPCP, precisionDPCP, recalDPCP] = DPCP(W, p);
    disp("Execute DPCPN");   
    [XinDPCPn, precisionDPCPN, recalDPCPN] = DPCP_Normal(W, p);
    disp("Execute RANSAC"); 
    [XinRANSAC, precisionRANSAC, recalRANSAC] = RANSAC(W, p);
    disp("Execute APG"); 
    [XinAPG, precisionAPG, recalAPG] = APG(W, p);
    psDPCP = psDPCP + precisionDPCP / 30;
    psDPCPN = psDPCPN + precisionDPCPN / 30;
    psRANSAC = psRANSAC + precisionRANSAC / 30;
    psAPG = psAPG + precisionAPG / 30;
    rsDPCP = rsDPCP + precisionDPCP / 30;
    rsDPCPN = rsDPCPN + precisionDPCPN / 30;
    rsRANSAC = rsRANSAC + precisionRANSAC / 30;
    rsAPG = rsAPG + precisionAPG / 30;
    %%
    %% Estimation M
    %B = [2, 0, 0];
    %disp(XinDPCP);
    %M = Fast_Compressed_Least_Squares(XinDPCP{1}, 2, 2, 2, B);
    %disp(M);

    end
    %% plot figure
    x = [x;precisionDPCP, precisionDPCPN, precisionRANSAC, precisionAPG];
    colormap(cool);
    y = [y;recalDPCP, recalDPCPN, recalRANSAC, recalAPG];
    colormap(cool);
    subplot(1,2,1);
    bar(x);
    subplot(1,2,2);
    bar(y);
    colormap(cool);
end