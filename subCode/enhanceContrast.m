function pMap = enhanceContrast(I, b)
t = 0.5*max(I(:));
v1 = mean(I(I>=t));
v2 = mean(I(I<t));
pMap = 1./(1+exp(-b*(I-0.5*(v1+v2))));