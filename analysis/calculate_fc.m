%% functional connectivity
%% Jessica Jesser
%% January 2015, modification October 2015

%% September 2022, edited by Tianlu Wang
% Description: calculate the whole brain FC and homotopic FC from full
% connectivity matrices
%
% Input: 3D matrix, 3rd dimension containing subjects
%
% Output: lists of whole brain FC and homotopic fc, size n_subjects x 1 
%
%% AAL ROIs
% 1:45 left cerebrum
% 46:54 left cerebellum
% 55:62 vermis
% 63:71 right cerebellum
% 72:116 right cerebrum (in reverse order)
% List of AAL ROIs from REST toolbox:
% AALTag = {'PreCG.L ','SFGdor.L','ORBsup.L ', 'MFG.L ','ORBmid.L ','IFGoperc.L ','IFGtriang.L ','ORBinf.L ','ROL.L ','SMA.L ', 'OLF.L ','SFGmed.L ','ORBsupmed.L ','REC.L ','INS.L ','ACG.L ','DCG.L ','PCG.L ','HIP.L ','PHG.L ','AMYG.L ','CAL.L ','CUN.L ','LING.L ','SOG.L ','MOG.L ','IOG.L ','FFG.L ','PoCG.L ','SPG.L ','IPL.L ','SMG.L ','ANG.L ','PCUN.L ','PCL.L ','CAU.L ','PUT.L ','PAL.L ','THA.L ','HES.L ','STG.L ','TPOsup.L ','MTG.L ','TPOmid.L ','ITG.L ','Crbl crus1.L','Crbl crus2.L','Crbl 3.L','Crbl 4 5.L','Crbl 6.L','Crbl  7b.L','Crbl 8.L','Crbl 9.L','Crbl 10.L','Vermis 1 2','Vermis 3','Vermis 4 5','Vermis 6','Vermis 7','Vermis 8','Vermis 9','Vermis 10','Crbl 10.R','Crbl 9.R','Crbl 8.R','Crbl  7b.R','Crbl 6.R','Crbl 4 5.R','Crbl 3.R','Crbl crus2.R','Crbl crus1.R','ITG.R ','TPOmid.R ','MTG.R ','TPOsup.R ','STG.R ','HES.R ','THA.R ','PAL.R ','PUT.R ','CAU.R ','PCL.R ','PCUN.R ','ANG.R ','SMG.R ','IPL.R ','SPG.R ','PoCG.R ','FFG.R ','IOG.R ','MOG.R ','SOG.R ','LING.R ','CUN.R ','CAL.R ','AMYG.R ','PHG.R ','HIP.R ','PCG.R ','DCG.R ','ACG.R ','INS.R ','REC.R ','ORBsupmed.R ','SFGmed.R ','OLF.R ','SMA.R ','ROL.R ','ORBinf.R ','IFGtriang.R ','IFGoperc.R ','ORBmid.R ','MFG.R ','ORBsup.R ','SFGdor.R ','PreCG.R '};

%%

function [fcwb, fch] = calculate_fc(X)

% FC of entire matrix
fcwb = zeros(size(X,3),1);
for idx = 1:size(X,3)
    fcwb(idx) = mean(nonzeros(tril(X(:,:,idx),-1)));
end

% FC between homotopic regions
fch = zeros(size(X,3),1);
for idx = 1:size(X,3)
    subj_fch = zeros(54,1);
    for idxr = 1:54
        subj_fch(idxr) = X(idxr,117-idxr,idx);
    end
    fch(idx) = mean(subj_fch);
end

