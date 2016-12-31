function clipVal = EstimateDynamicParas(adjcMatrix, colDistM)
[meanMin1, meanTop, meanMin2] = GetMeanMinAndMeanTop(adjcMatrix, colDistM, 0.01);
clipVal = meanMin2;
