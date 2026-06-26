close all
clear all
clc

addpath(genpath(fullfile('T1D_VPP', 'Single Hormone Population')));

nPats = 99;
rng(42)
seeds = randi([0, 1000], 99, 1);
samples = 30*288;

bolusModulations = 0.5:0.1:1.5;
basalModulations = 0.5:0.1:1.5;
mealModulations = 0:0.2:2.0;

glucose = zeros(samples,nPats);
bolus = zeros(samples,nPats);
iob = zeros(samples,nPats);
basal = zeros(samples,nPats);
meal = zeros(samples,nPats);
weights = zeros(nPats,1);
crs = zeros(nPats,1);
cfs = zeros(nPats,1);
gts = zeros(nPats,1);

%%
names = {'baseline'};

for n = 1:length(names)

    disp(['Creating DatasetPJ_' names{n} '...']);

    for nn = 1:nPats

        [G, CHO, IB, Ib, BW, CR, CF, GT] = runSimulation(nn,0,nan,1,1, seeds(nn));
        glucose(:,nn) = G(1:samples);
        bolus(:,nn) = IB(1:samples);    % u/min
        basal(:,nn) = Ib(1:samples);    % u/min
        meal(:,nn) = CHO(1:samples);
        weights(nn) = BW;
        crs(nn) = CR;
        cfs(nn) = CF;
        gts(nn) = GT;

    end

    saveResults(glucose, bolus, basal, meal, iob, weights, crs, cfs, gts, ['DatasetPJ_' names{n}])

end

bolus_baseline = bolus;
writematrix(bolus_baseline, 'bolus_baseline.csv')

%%
names = {'bolusMinus50','bolusMinus40','bolusMinus30','bolusMinus20','bolusMinus10',...
        'baseline_bolus',...
        'bolusPlus10','bolusPlus20','bolusPlus30','bolusPlus40','bolusPlus50'};

for n = 1:length(names)

    disp(['Creating DatasetPJ_' names{n} '...']);

    for nn = 1:nPats

        bolus_input = bolus_baseline(:,nn)'*bolusModulations(n);

        [G, CHO, IB, Ib, BW, CR, CF, GT] = runSimulation(nn,1,bolus_input,1,1, seeds(nn));
        glucose(:,nn) = G(1:samples);
        bolus(:,nn) = IB(1:samples);
        basal(:,nn) = Ib(1:samples);
        meal(:,nn) = CHO(1:samples);
        weights(nn) = BW;
        crs(nn) = CR;
        cfs(nn) = CF;
        gts(nn) = GT;

    end

    saveResults(glucose, bolus, basal, meal, iob, weights, crs, cfs, gts, ['DatasetPJ_' names{n}])

end

names = {'basalMinus50','basalMinus40','basalMinus30','basalMinus20','basalMinus10',...
        'baseline_basal',...
        'basalPlus10','basalPlus20','basalPlus30','basalPlus40','basalPlus50'};

for n = 1:length(names)

    disp(['Creating DatasetPJ_' names{n} '...']);

    for nn = 1:nPats

        [G, CHO, IB, Ib, BW, CR, CF, GT] = runSimulation(nn,1,bolus_baseline(:,nn)',1,basalModulations(n), seeds(nn));
        glucose(:,nn) = G(1:samples);
        bolus(:,nn) = IB(1:samples);
        basal(:,nn) = Ib(1:samples);
        meal(:,nn) = CHO(1:samples);
        weights(nn) = BW;
        crs(nn) = CR;
        cfs(nn) = CF;
        gts(nn) = GT;

    end
    saveResults(glucose, bolus, basal, meal, iob, weights, crs, cfs, gts, ['DatasetPJ_' names{n}])
end

names = {'mealMinus100','mealMinus80','mealMinus60','mealMinus40','mealMinus20',...
        'baseline_meal',...
        'mealPlus20','mealPlus40','mealPlus60','mealPlus80','mealPlus100'};

for n = 1:length(names)

    disp(['Creating DatasetPJ_' names{n} '...']);

    for nn = 1:nPats

        [G, CHO, IB, Ib, BW, CR, CF, GT] = runSimulation(nn,1,bolus_baseline(:,nn)',mealModulations(n),1, seeds(nn));
        glucose(:,nn) = G(1:samples);
        bolus(:,nn) = IB(1:samples);
        basal(:,nn) = Ib(1:samples);
        meal(:,nn) = CHO(1:samples);
        weights(nn) = BW;
        crs(nn) = CR;
        cfs(nn) = CF;
        gts(nn) = GT;

    end
    saveResults(glucose, bolus, basal, meal, iob, weights, crs, cfs, gts, ['DatasetPJ_' names{n}])
end
% 

names = {'hypotreatment25'};
bolus_baseline = readmatrix('bolus_baseline.csv');

for n = 1:length(names)

    disp(['Creating DatasetPJ_' names{n} '...']);

    for nn = 1:nPats

        [G, CHO, IB, Ib, BW, CR, CF, GT] = runSimulation(nn,1,bolus_baseline(:,nn)',1,1,seeds(nn));
        glucose(:,nn) = G(1:samples);
        bolus(:,nn) = IB(1:samples);
        basal(:,nn) = Ib(1:samples);
        meal(:,nn) = CHO(1:samples);
        weights(nn) = BW;
        crs(nn) = CR;
        cfs(nn) = CF;
        gts(nn) = GT;

    end
    saveResults(glucose, bolus, basal, meal, iob, weights, crs, cfs, gts, ['DatasetPJ_' names{n}])
end


names = {'correctionBolus50CF'};


for n = 1:length(names)

    disp(['Creating DatasetPJ_' names{n} '...']);

    for nn = 1:nPats

        bolus_input = bolus_baseline(:, nn);
        bolus_input(9.5*12+1) = 50/cfs(nn)/5;

        [G, CHO, IB, Ib, BW, CR, CF, GT] = runSimulation(nn,1,bolus_input',1,1,seeds(nn));
        glucose(:,nn) = G(1:samples);
        bolus(:,nn) = IB(1:samples);
        basal(:,nn) = Ib(1:samples);
        meal(:,nn) = CHO(1:samples);
        weights(nn) = BW;
        crs(nn) = CR;
        cfs(nn) = CF;
        gts(nn) = GT;

    end
    saveResults(glucose, bolus, basal, meal, iob, weights, crs, cfs, gts, ['DatasetPJ_' names{n}])
end



