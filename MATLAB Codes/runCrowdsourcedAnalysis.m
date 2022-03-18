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
phspan = 7;
wrt = 'n';
tblOut = runLambda_v0p31_crowdsourced(tbl_fticr,phspan,wrt,dataDescrp);
assignin('base',"tblOut",tblOut)

idxExcl = find(isnan(tblOut.lambda)|isinf(tblOut.lambda));
tblOut(idxExcl,:) = [];
tbl_fticr(idxExcl,:) = [];

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
samp_ = char(samp); samp_ = string(samp_(:,1:9)); 
idx_samp_int = find(contains(samp_,samp_int));
samp_int = tbl_fticr.Properties.VariableNames(idx_samp_int+sampCol-1)';
samp_pern = tbl_meta.Sample_ID(contains...
    (tbl_meta.Intermittent_or_Perennial,"perennial",'IgnoreCase',true));
idx_samp_pern = find(contains(samp_,samp_pern));
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

%% PCA

dat = zeros(size(tbl_fticr{:,sampCol:end}))';

for iSamp = 1:size(dat,1)
    idx = find(tbl_fticr{:,sampCol+iSamp-1});
    if ~isempty(idx)
        dat(iSamp,idx) = tblOut.lambda(idx);
    end
end

grouping = cell(size(dat,1),3);
grouping(:,:) = cellstr("na");
grouping(idx_samp_sed,1) = cellstr("Sediment");
grouping(idx_samp_sw,1) = cellstr("Surface water");
grouping(idx_samp_int,2) = cellstr("Intermittent");
grouping(idx_samp_pern,2) = cellstr("Perennial");
grouping(intersect(idx_samp_sed,idx_samp_int),3) = cellstr("Sediment + Intermittent");
grouping(intersect(idx_samp_sed,idx_samp_pern),3) = cellstr("Sediment + Perennial");
grouping(intersect(idx_samp_sw,idx_samp_int),3) = cellstr("Surface water + Intermittent");
grouping(intersect(idx_samp_sw,idx_samp_pern),3) = cellstr("Surface water + Perennial");

T1 = array2table(dat','VariableNames',samp,'RowNames',cellstr("cpd_"+string(1:size(dat,2))));
T2 = cell2table(grouping);

% writetable(T1,"data.csv")
% writetable(T2,"grouping.csv")

%% Distirbution of thermodynamic properties



%% SINDy

outp = goSindy(depVar,Theta,lambdaGuess,nlambda,lambda,plt);

%% Postprocessing


% end