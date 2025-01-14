results = ones(3,10,3);

T1 = [0 0 0 0 0.041 0.454 0.607 0.694 0.816 0.978;
 0.472 0.098 0.033 0.019 0.023 0.039 0.055 0.068 0.073 0.016;
 0.298 0.344 0.148 0.092 0.118 0.178 0.193 0.170 0.094 0.006];

T2 = [0 0 0 0 0.267 0.594 0.712 0.791 0.917 nan;
 0.395 0.073 0.025 0.016 0.018 0.031 0.043 0.048 0.035 nan;
 0.465 0.459 0.196 0.126 0.141 0.195 0.189 0.145 0.049 nan];

T3 = [0 0 0 0.053 0.388 0.674 0.779 0.854 0.982 nan;
 0.345 0.06 0.021 0.014 0.017 0.03 0.039 0.039 0.005 nan;
 0.55 0.477 0.199 0.126 0.133 0.172 0.153 0.103 0.013 nan];

results(:,:,1) = T1;
results(:,:,2) = T2;
results(:,:,3) = T3;

results_mat = matfile('ResultsFromLit.mat','Writable',true);
save('ResultsFromLit.mat','results')