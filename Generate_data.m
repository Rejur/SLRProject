L = 400;
tmp = ones(1, L);
H = [ 1 2 3
    4 5 6
    7 8 9];
x1 = randn(2, L);
x1 = [x1; tmp];
x1t = H * x1;
x2 = x1t .* 2;
