%NetworkCalculations_TP.mat This code conducts stepwise multiple linear
%regression modeling, as described in Larsen, Newman, Saunders, and Harvey,
%"Complex networks of functional connectivity in an isolated wetland
%reconnected to its floodplain." Inputs are specified in the User Input
%section below. Outputs are printed to the screen. Recommended use is to
%execute the code part by part rather than all at once, as this helps in
%keeping track of which multiple linear regression model is being run. Once
%the critical variables in the multiple linear regression model are
%determined, their associated correlation coefficients can be located by
%examining (double-clicking) on the rpTPdiff and rpTP arrays. We selected
%the Pearson correlation value in row 1. Column headings for these arrays
%are given by TotalDiffLabels and TotalLabels, respectively. The code
%requires the statistics toolbox for Matlab.

%This version of the code conducts analyses for total phosphorus and is the
%fully commented verison.

% Copyright 2017 Laurel Larsen
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
% 
% 1. Redistributions of source code must retain the above copyright notice,
% this list of conditions and the following disclaimer.
% 
% 2. Redistributions in binary form must reproduce the above copyright
% notice, this list of conditions and the following disclaimer in the
% documentation and/or other materials provided with the distribution.
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
% IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
% THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
% PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
% CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
% EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
% PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
% PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
% LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
% NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

%Clear the memory, start with a clean slate. :)
clear
clc
close all

%USER INPUT SECTION
A = importdata('/File path/InputSS.csv'); %Enter the path and file name of the input spreadsheet.
b = A.textdata;
c = A.data;
%Below, enter the columns for the response variables at the site of
%interest (residuals and the residuals of the "diff model"), the value of the
%analyte for the sites other than the site of interest, the values of the
%residuals at all sites other than the site of interest, and the value of
%the "diff" at all sites other than the site of interest
ResidlogTPdiff = c(:,383);
ResidlogTP = c(:,279);
SitelogTP = c(:,[461:472]);
logTPdiffResids = c(:,[384,386,387,388,390,391,392,393,394,395]);
logTPResids = c(:, [280,282,283,284,286,287,288,289,290,291]);
logTPdiff = c(:,345:356);

%LABELS FOR SITES FOR RESPONSE AND TARGET VARIABLES
sites = {'C1', 'DB1', 'DB2', 'DB3', 'RS1', 'RS2', 'S1', 'UB1', 'UB2', 'UB3'};
allsites = {'C1','C2','DB1', 'DB2', 'DB3','L67A', 'RS1', 'RS2', 'S1', 'UB1', 'UB2', 'UB3'};
Residlabels = {'Resid_C1', 'Resid_DB1', 'Resid_DB2', 'Resid_DB3', 'Resid_RS1', 'Resid_RS2', 'Resid_S1', 'Resid_UB1', 'Resid_UB2', 'Resid_UB3'};
Residlabels2 = {'Resid_C1', 'Resid_C2', 'Resid_DB1', 'Resid_DB2', 'Resid_DB3', 'Resid_L67A', 'Resid_RS1', 'Resid_RS2', 'Resid_S1', 'Resid_UB1', 'Resid_UB2', 'Resid_UB3'};
Difflabels = {'C1_diff','C2_diff', 'DB1_diff', 'DB2_diff', 'DB3_diff','L67A_diff', 'RS1_diff', 'RS2_diff', 'S1_diff', 'UB1_diff', 'UB2_diff', 'UB3_diff'};
Previouslabels = {'C1_previous','C2_previous', 'DB1_previous', 'DB2_previous', 'DB3_previous','L67A_previous', 'RS1_previous', 'RS2_previous', 'S1_previous', 'UB1_previous', 'UB2_previous', 'UB3_previous'};
sampledsites = b(2:end,2);

%MODEL EXECUTION
for ii = 1:10; %Loop over all sites
    %Find the selected site and corresponding inputs.
    myrows = find(ismember(sampledsites, sites{ii}));
    mycols = find(ismember(allsites, sites{ii}));
    TPdiffInput = [ResidlogTPdiff(myrows), SitelogTP(myrows, setdiff(1:12, mycols)), logTPdiffResids(myrows, setdiff(1:10, ii)), logTPdiff(myrows, setdiff(1:12, mycols))];
    TPInput = [ResidlogTP(myrows), SitelogTP(myrows, setdiff(1:12, mycols)), logTPResids(myrows, setdiff(1:10, ii))]; 
    TotalDiffLabels = [allsites(setdiff(1:12, mycols)), Residlabels(setdiff(1:10,ii)), Difflabels(setdiff(1:12, mycols))];
    TotalLabels = [allsites(setdiff(1:12, mycols)), Residlabels(setdiff(1:10, ii))]; 
    
    rpTPdiff = NaN(3, size(TPdiffInput,2)-1);
    for jj = 1:size(TPdiffInput, 2)-1 %Loop over all sites that may be contributing information to the site of interest.
        TPdiffIn = TPdiffInput(:, [1, jj+1]);
        TPdiffIn = TPdiffIn(~isnan(sum(TPdiffIn, 2)), :); %Remove rows with NaNs
        TPdiffIn = TPdiffIn-repmat(mean(TPdiffIn,1), size(TPdiffIn,1), 1);
        [r, p] = corrcoef(TPdiffIn); %Use a Pearson's approach to generate correlation coefficients. Store the correlation (r) and the significance (p) as the first two indices of the rpTPDiff variable. 
        rpTPdiff(1,jj) = r(1,2); %Store value
        rMonte = NaN(1,1000);
        for kk = 1:1000 %Use a Monte Carlo approach to generate correlations, stored as the third index of the rpTPDiff variable.
            r = corrcoef([randsample(TPdiffIn(:,1), size(TPdiffIn,1)), randsample(TPdiffIn(:,2), size(TPdiffIn,1))]);
            rMonte(kk) = r(1,2);
        end
        pMonte = 1-numel(find(abs(rMonte)<=abs(rpTPdiff(1,jj))))./1000;
        rpTPdiff(2,jj) = p(1,2);
        rpTPdiff(3,jj) = pMonte;
    end
    sigs = find(rpTPdiff(3,:)<=0.05); %Find significant correlations using an alpha = 0.05 significance criterion. 

    d1 = [TotalDiffLabels; num2cell(rpTPdiff)];
    sigs2 = find(rpTPdiff(3,:)<=0.1); %Find significant correlations using an alpha = 0.1 significance criterion. 
    sigs3 = find(ismember(sigs2, sigs));
    T = zeros(numel(sigs)+1, numel(sigs2)+1); %This is the list of the variables that the backwards regression model will be initiated with.
    for ll = 1:numel(sigs)
        T(1+ll, sigs3(ll)) = 1;
    end
    TPdiffInput2 = TPdiffInput(:,[1,sigs2+1]); %The variables with significant correlations are inputs to the stepwise MLR model. 
    TPdiffInput2 = TPdiffInput2(~isnan(sum(TPdiffInput2,2)),:); %Take out any rows with NaNs.
    %Read in the data and conduct the forward and backwards regression
    %model for the "diff" data.
    readin = array2table([TPdiffInput2(:,2:end), TPdiffInput2(:,1)], 'VariableNames', [TotalDiffLabels(sigs2), 'ResidTPdiff']);
    mdlDdb = stepwiselm(readin,  T, 'ResponseVar', 'ResidTPdiff','Upper', 'linear') %Backwards, starting with all significant correlations that conform to the alpha = 0.05 criterion.
    mdlDdb.ModelCriterion.AICc %Print the AICc value.
    mdlDdf = stepwiselm(readin, 'constant', 'ResponseVar', 'ResidTPdiff', 'Upper', 'linear') %Forwards.
    mdlDdf.ModelCriterion.AICc %Print the AICc value.
    
    %NOW REPEAT THE SAME PROCESS FOR THE "NON-DIFF" DATA. THIS PART FOLLOWS
    %THE SAME STRUCTURE AS ABOVE. SEE ABOVE FOR DETAILED COMMENTS.
    %First, find the basic correlations.
    rpTP = NaN(3, size(TPInput,2)-1);
    for jj = 1:size(TPInput, 2)-1
        TPIn = TPInput(:, [1, jj+1]);
        TPIn = TPIn(~isnan(sum(TPIn, 2)), :); %Remove rows with NaNs
        TPIn = TPIn-repmat(mean(TPIn,1), size(TPIn,1), 1);
        [r, p] = corrcoef(TPIn);
        rpTP(1,jj) = r(1,2);
        rMonte = NaN(1,1000);
        for kk = 1:1000
            r = corrcoef([randsample(TPIn(:,1), size(TPIn,1)), randsample(TPIn(:,2), size(TPIn,1))]);
            rMonte(kk) = r(1,2);
        end
        pMonte = 1-numel(find(abs(rMonte)<=abs(rpTP(1,jj))))./1000;
        rpTP(2,jj) = p(1,2);
        rpTP(3,jj) = pMonte;
    end
    sigs = find(rpTP(3,:)<=0.05);

    d3 = [TotalLabels; num2cell(rpTP)];
    sigs2 = find(rpTP(3,:)<=0.1);
    sigs3 = find(ismember(sigs2, sigs));
    T = zeros(numel(sigs)+1, numel(sigs2)+1); %How to initialize the backwards regression
    for ll = 1:numel(sigs)
        T(1+ll, sigs3(ll)) = 1;
    end
    %Read in the inputs to the stepwise regression models.
    TPInput2 = TPInput(:,[1,sigs2+1]);
    TPInput2 = TPInput2(~isnan(sum(TPInput2,2)),:);
    readin = array2table([TPInput2(:,2:end), TPInput2(:,1)], 'VariableNames', [TotalLabels(sigs2), 'TP']);
    %Run the regressions and print the outputs to the screen. 
    mdlDb = stepwiselm(readin,  T, 'ResponseVar', 'TP','Upper', 'linear')
    mdlDb.ModelCriterion.AICc
    mdlDf = stepwiselm(readin, 'constant', 'ResponseVar', 'TP', 'Upper', 'linear')
    mdlDf.ModelCriterion.AICc

end
    