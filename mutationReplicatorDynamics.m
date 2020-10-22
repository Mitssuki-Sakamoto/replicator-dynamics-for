

datetimestr =  datestr(datetime('now'), "yyyy-mm-dd-HH-MM-SS");
logdir = "outcomes/" + datetimestr + '/';
mkdir(logdir);
saveCsvsDir = logdir+"csvs/";
saveImagesDir = logdir+"images/";
mkdir(saveCsvsDir);
mkdir(saveImagesDir);
payoffCsvFile =  "payoff_0.90_0.10_0.90.csv";
% 少々雑
payoffName = split(payoffCsvFile,".csv");
payoffName = payoffName(1);
payoffMatrix = csvread(payoffCsvFile,1,1);
mutaionValue = 0.01;
dt = 0.2;

mutationRates = ones(length(payoffMatrix));
nStragtegies = length(payoffMatrix);
for i = 1:nStragtegies
    for j = 1:nStragtegies
        if i == j
            mutationRates(i,j) = 1 - mutaionValue;
        else
            mutationRates(i,j) = mutaionValue/(nStragtegies-1);
        end
    end
end

rd = {
    @(populations) mutationRD1(payoffMatrix, populations, mutationRates);
    @(populations) mutationRD2(payoffMatrix, populations, mutaionValue);
    @(populations) mutationRD3(payoffMatrix, populations, mutationRates);
    };

for ir = 1:length(rd)
    populations = ones(length(payoffMatrix),1) ./ length(payoffMatrix);
    populationsHistories = [populations];
    count = 0;
    tic;
    while 1
        count = count + dt;
        dx = rd{ir}(populations);
        populations = populations + (dx * dt);
        populationsHistories = [populationsHistories, populations];
        if max(reshape(dx,1,[])) < 0.00001 || count > 1000
            disp(count);
            break;
        end
    end
    toc;
    fileName = "rd"+ ir + "_" + payoffName +".csv";
    csvwrite(saveCsvsDir + fileName, populationsHistories.')
    figure;
    plot(populationsHistories.');
    ylim([0 1]);
    f = gcf;
    fileName =  "rd"+ ir + "_"  + payoffName + ".png";
    exportgraphics(f, saveImagesDir + fileName);
    disp(populations);
end



% 強化学習のサーベイ論文
function dv = mutationRD1(payoffMatrix, populations, mutationRates)
    strategyAveragePayoffs = payoffMatrix * populations;
    populationAveragePayoff =  populations.' * strategyAveragePayoffs;
    dv = populations.*((strategyAveragePayoffs - populationAveragePayoff)) ...
        + sum(mutationRates .* (repmat(populations,1,length(populations))-repmat(populations.',length(populations), 1))).';  
end

% 新レプリーター
function dv = mutationRD2(payoffMatrix, populations, mutationValue)
    strategyAveragePayoffs = payoffMatrix * populations;
    populationAveragePayoff =  populations.' * strategyAveragePayoffs;
    dv = populations.*((strategyAveragePayoffs - populationAveragePayoff)) ...
        + mutationValue .* (repmat(1/length(populations), length(populations), 1) - populations);  
end

% 旧レプリーター
function dv = mutationRD3(payoffMatrix, populations, mutationRates)
    strategyAveragePayoffs = payoffMatrix  * populations;
    populationAveragePayoff =  populations.' * strategyAveragePayoffs;
    dv = sum(repmat(populations,1,length(populations)) .* repmat(strategyAveragePayoffs,1,length(strategyAveragePayoffs)) .* mutationRates).' - populations .* populationAveragePayoff;
end