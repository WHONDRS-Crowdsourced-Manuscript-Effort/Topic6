% function runCrowdsourcedAnalysis

clear all; close all; clc;

%% Input files

fticrdataFilename = "ExperimentalData/Processed_S19S_Sediments_Water_2-2_newcode.csv";
metadataFilename = "ExperimentalData/WHONDRS_S19S_Metadata_v3.csv";
sedrespdataFilename = "ExperimentalData/WHONDRS_S19S_Sediment_Incubations_"+...
    "Respiration_Rates.csv";

% -------------------------------------------------------------------------
% Naming convention:
% -------------------------------------------------------------------------
% S19S_####; #### spans 0001 to 0100
% -------------------------------------------------------------------------
% FTICR
% Sed sample ID: S19S_####_Sed_Field_ICR.$_XXX; $ = D,M,U are replicates
% Surf sample ID: S19S_####_ICR.$_XXX; $ = 1,2,3 are replicates
% -------------------------------------------------------------------------
% Sediment respiration
% Sed sample ID: S19S_####_SED_INC-$; $ = D,M,U are replicates
% Column identifiers: 
% Sample_ID
% rate_mg_per_L_per_h
% -------------------------------------------------------------------------
% Metadata
% *** Remove second row ***
% Sample ID: S19S_####
% Column identifiers: 
% Sample_ID
% Stream_Order (number, Not_Provided) 
% General_Vegetation (*grass**, **shrub**, **tree**, not vegetated) + combo
% Intermittent_or_Perennial (intermittent, perennial, Not_Provided)
% SW_pH (**number**)
% -------------------------------------------------------------------------

%% Read files

tbl_fticr = readtable(fticrdataFilename);
tbl_meta = readtable(metadataFilename);
tbl_resp = readtable(sedrespdataFilename);

%% Preprocessing

% Remove non carbon sources
tbl_fticr(tbl_fticr.C==0,:) = [];

%% Run lambda

dataDescrp = "fullData";
phspan = unique(tbl_meta.pH);
wrt = 'n';
tbl_OutputMaster = runLambda_v0p31_crowdsourced(tbl_fticr,phspan,wrt,dataDescrp);
assignin('base',"tbl_OutputMaster",tbl_OutputMaster)

%% Metadata classification

sampCol = 39;       % column # where sample name starts
samp = tbl_fticr.Properties.VariableNames(39:end)';

% Classify sediment and surface-water data
samp_sed = samp(contains(samp,"sed",'IgnoreCase',true)); 
idx_samp_sed = find(contains(samp,samp_sed));
samp_sw = samp(~contains(samp,"sed",'IgnoreCase',true)); 
idx_samp_sw = find(contains(samp,samp_sw));

idx_comp_sed = [];
for iSamp = 1:length(idx_samp_sed)
    idx_comp_sed = [idx_comp_sed; 
                    find(tbl_fticr{:,idx_samp_sed(iSamp)+sampCol-1})];
end
idx_comp_sed = unique(idx_comp_sed);

idx_comp_sw = [];
for iSamp = 1:length(idx_samp_sw)
    idx_comp_sw = [idx_comp_sw; 
                    find(tbl_fticr{:,idx_samp_sw(iSamp)+sampCol-1})];
end
idx_comp_sw = unique(idx_comp_sw);

% Classify intermittent and perennial data
samp_int = tbl_meta.Sample_ID(contains...
    (tbl_meta.Intermittent_or_Perennial,"intermittent",'IgnoreCase',true));
idx_samp_int = find(contains(extractBefore(samp,10),samp_int));
samp_int = tbl_fticr.Properties.VariableNames(idx_samp_int+sampCol-1)';
samp_pern = tbl_meta.Sample_ID(contains...
    (tbl_meta.Intermittent_or_Perennial,"perennial",'IgnoreCase',true));
idx_samp_pern = find(contains(extractBefore(samp,10),samp_pern));
samp_pern = tbl_fticr.Properties.VariableNames(idx_samp_pern+sampCol-1)';

idx_comp_int = [];
for iSamp = 1:length(idx_samp_int)
    idx_comp_int = [idx_comp_int; 
                    find(tbl_fticr{:,idx_samp_int(iSamp)+sampCol-1})];
end
idx_comp_int = unique(idx_comp_int);

idx_comp_pern = [];
for iSamp = 1:length(idx_samp_pern)
    idx_comp_pern = [idx_comp_pern; 
                    find(tbl_fticr{:,idx_samp_pern(iSamp)+sampCol-1})];
end
idx_comp_pern = unique(idx_comp_pern);

n_sed = length(samp_sed);
n_sw = length(samp_sw);
n_int = length(samp_int);
n_pern = length(samp_pern);
n_sed_int = length(intersect(samp_sed,samp_int));
n_sed_pern = length(intersect(samp_sed,samp_pern));
n_sw_int = length(intersect(samp_sw,samp_int));
n_sw_pern = length(intersect(samp_sw,samp_pern));

%% Sample pH-specific results

dat = zeros(size(tbl_fticr{:,sampCol:end}))';
for iSamp = 1:size(dat,1)
    idx_meta = find(contains(tbl_meta.Sample_ID,extractBefore(samp(iSamp),10)));
    idx_pH = find(phspan==tbl_meta.pH(idx_meta));
    tblTemp = tbl_OutputMaster.tblOut{idx_pH};
    idx = find(tbl_fticr{:,sampCol+iSamp-1});
    if ~isempty(idx)
%         dat(iSamp,idx) = tblTemp.lambda(idx);
        dat(iSamp,idx) = tblTemp.delGcox(idx);
    end
end

tbl_OutpH = [table(tbl_fticr.MolForm,'VariableNames',"Compound") array2table(dat','VariableNames',samp)];

%% PCA grouping

% grouping = cell(size(dat,1),3);
% grouping(:,:) = cellstr("na");
% grouping(idx_samp_sed,1) = cellstr("Sediment");
% grouping(idx_samp_sw,1) = cellstr("Surface water");
% grouping(idx_samp_int,2) = cellstr("Intermittent");
% grouping(idx_samp_pern,2) = cellstr("Perennial");
% grouping(intersect(idx_samp_sed,idx_samp_int),3) = cellstr("Sediment + Intermittent");
% grouping(intersect(idx_samp_sed,idx_samp_pern),3) = cellstr("Sediment + Perennial");
% grouping(intersect(idx_samp_sw,idx_samp_int),3) = cellstr("Surface water + Intermittent");
% grouping(intersect(idx_samp_sw,idx_samp_pern),3) = cellstr("Surface water + Perennial");
% 
% tbl_grouping = cell2table(grouping);

% writetable(tbl_OutpH,"PCAdata_lambda.csv")
% writetable(tbl_OutpH,"PCAdata_delGcox.csv")
% writetable(tbl_grouping,"PCAgrouping.csv")

%% Respiration rates

nOC = 1;
nO2 = 9;
tblTemp = tbl_OutputMaster.tblOut{phspan==7};
idxExcl = find(isnan(tblTemp.lambda)|isinf(tblTemp.lambda));
tblTemp(idxExcl,:) = [];

[VhOC,VhO2] = meshgrid(linspace(0.01,2,50),linspace(0.1,10,50));
muRel_avg = exp(-abs(mean(tblTemp.stoichMet(:,nOC)))./VhOC).*exp(-abs(mean(tblTemp.stoichMet(:,nO2)))./VhO2);

figure(1)
surf(VhOC,VhO2,muRel_avg)
xlabel("VhOC")
ylabel("VhO2")
zlabel("Averaged \mu_{rel}")
colorbar

VhOC_lim = 0.05;
VhOC_ab = 2;
VhO2_lim = 0.5;
VhO2_ab = 10;

muRel_OC_lim = exp(-abs(tblTemp.stoichMet(:,nOC))./VhOC_lim).*exp(-abs(tblTemp.stoichMet(:,nO2))./VhO2_ab);
muRel_O2_lim = exp(-abs(tblTemp.stoichMet(:,nOC))./VhOC_ab).*exp(-abs(tblTemp.stoichMet(:,nO2))./VhO2_lim);
muRel_OCO2_lim = exp(-abs(tblTemp.stoichMet(:,nOC))./VhOC_lim).*exp(-abs(tblTemp.stoichMet(:,nO2))./VhO2_lim);

resp_OC_lim = zeros(size(tblTemp,1),length(samp));
resp_O2_lim = zeros(size(tblTemp,1),length(samp));
resp_OCO2_lim = zeros(size(tblTemp,1),length(samp));

tblfticrTemp = tbl_fticr;
tblfticrTemp(idxExcl,:) = [];
for iSamp = 1:size(resp_OC_lim,2)
    idx = find(tblfticrTemp{:,sampCol+iSamp-1});
    if ~isempty(idx)
        resp_OC_lim(idx,iSamp) = muRel_OC_lim(idx);
        resp_O2_lim(idx,iSamp) = muRel_O2_lim(idx);
        resp_OCO2_lim(idx,iSamp) = muRel_OCO2_lim(idx);
    end
end
%%
idx_samp_resp = [];
for iSamp = 1:length(samp)
    idx_samp_resp = [idx_samp_resp; find(contains(extractBefore(tbl_resp.Sample_ID,10),extractBefore(samp(iSamp),10)))];
end


return

%% Distirbution of thermodynamic properties

% figure(f)
% subplot(1,2,1)
% hLambda = histogram(tblOut.lambda);
% xlabel("\lambda")
% ylabel("Frequency")
% subplot(1,2,2)
% hDelG = histogram(tblOut.delGcox);
% xlabel("\DeltaG_{cox}")
% ylabel("Frequency")
% 
% figure(f)
% subplot(1,2,1)
% ecdf(tblOut.lambda)
% h = get(gca,'children');
% set(h,'LineWidth',2)
% xlabel("\lambda")
% ylabel("Frequency")
% subplot(1,2,2)
% ecdf(tblOut.delGcox)
% h = get(gca,'children');
% set(h,'LineWidth',2)
% xlabel("\DeltaG_{cox}")
% ylabel("Frequency")

%% SINDy

return

outp = goSindy(depVar,Theta,lambdaGuess,nlambda,lambda,plt);

%% Postprocessing


% end