function saveResults(glucose, bolus, basal, meal, iob, weights, crs, cfs, gts, label)

    saveFolderName = fullfile('results', label);
    if ~exist(saveFolderName, 'dir')
        mkdir(saveFolderName)
    end
    
    % Cut
    startSample = 1;
    glucoseAll = glucose(startSample:end,:);
    basalAll = basal(startSample:end,:);
    bolusAll = bolus(startSample:end,:);
    mealAll = meal(startSample:end,:);
    iobAll = iob(startSample:end,:);
    
    for p = 1:size(glucose,2)

        t = datetime(2000,1,1,6,0,0):minutes(5):(datetime(2000,1,1,6,0,0)+minutes(5*size(glucoseAll,1)));
        t = t(1:end-1)';

        glucose = glucoseAll(:,p);
        glucose(glucose<0) = 10^(-5); 
        cho = mealAll(:,p);
        bolus = bolusAll(:,p);
        basal = basalAll(:,p);
        iob = iobAll(:,p);

        data = table(t,glucose,bolus,basal,cho,iob);

        writetable(data, fullfile(saveFolderName,['patient_' num2str(p) '.csv']));

    end

    patient = (1:99)';
    bw = weights;
    
    data = table(patient, weights, crs, cfs, gts);
    writetable(data, fullfile(saveFolderName,['patInfo.csv']));
    
end