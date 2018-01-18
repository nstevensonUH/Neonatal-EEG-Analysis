function prd1 = predict_svr(mdl, dat, mu1, sig1);

svl = mdl.totalSV; % number of support vectors
alpha = mdl.sv_coef; % (alpha_i^* - alpha_i;)
sv = mdl.SVs; % support vectors
gamma = mdl.Parameters(4); % RBF kernel parameter
b = mdl.rho; % bias value
val = (dat-mu1)./sig1; % normalise data (zscore)

%%%%%%%%%%%%%%%%%%%%%%%%
 prd = zeros(1,svl);
        for ii = 1:svl
            prd(ii) = alpha(ii).*exp(-gamma.*norm(val-sv(ii,:)).^2);
        end
    %%%%%%%%%%%%%%%%%%%%%%%%
prd1 = sum(prd)-b; % prediction
