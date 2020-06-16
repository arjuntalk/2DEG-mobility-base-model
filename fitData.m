% Load the csv
data = readmatrix('input_data.xlsx');
darkData = data(:,1:2);
%lightData = data(:,10:11);

datan2D = darkData(:,1);
datamob = darkData(:,2);

% datan2D = lightData(:,1);
% datamob = lightData(:,2);

% getMobility takes in density in units of 1E10/cm^2
% and return mobility in units of 1E6cm^2/Vs
x = datan2D;
y = datamob; 

% Specify bounds for each parameter
% First index is startpoint, 2nd is lower bound, 3rd is upper bound
deltabnds = [log10(0.11e-9),log10(0.11e-9),log10(0.11e-9)]; %roughness amplitude (m)
droughbnds = [log10(14e-9),log10(14e-9),log10(14e-9)]; %corelation length (m)
Nddbnds = [log10(1.5e15),log10(1.5e15),log10(1.5e15)]; %surface charge (/m2)
Nbkbnds = [log10(1.67e20),log10(1.67e20),log10(1.67e20)]; %background concentration (/m3)


ft = fittype( 'getMobility( x, logNdd, logNbk, logDelta, logDrough)' );
% Apparently fit doesn't keep track of the order you gave the coefficients
% for it to fit... So we need to check and rearrange everything in the that
% matlab is storing the coefficients
coeffOrder = coeffnames(ft);
stPts = zeros(1,length(coeffOrder));
upBnds = zeros(1,length(coeffOrder));
lwBnds = zeros(1,length(coeffOrder));
for ii = 1:length(coeffnames(ft))
    if strcmp(coeffOrder{ii},'logNdd')
        stPts(ii) = Nddbnds(1); lwBnds(ii) = Nddbnds(2); upBnds(ii) = Nddbnds(3);
    elseif strcmp(coeffOrder{ii},'logNbk')
        stPts(ii) = Nbkbnds(1); lwBnds(ii) = Nbkbnds(2); upBnds(ii) = Nbkbnds(3);
    elseif strcmp(coeffOrder{ii},'logDelta')
        stPts(ii) = deltabnds(1); lwBnds(ii) = deltabnds(2); upBnds(ii) = deltabnds(3);
    elseif strcmp(coeffOrder{ii},'logDrough')
        stPts(ii) = droughbnds(1); lwBnds(ii) = droughbnds(2); upBnds(ii) = droughbnds(3);
    end
end
[f,gof] = fit( x, y, ft, 'StartPoint', stPts, 'Lower', lwBnds,...
    'Upper', upBnds);
Rsquared = gof.rmse^2;
fprintf(1,'Done fitting!\n');
fprintf(1,'Generating figure comparing theory and data...\n');

% Make a figure of the model and data to compare
modeln2D = linspace(0.1,30,200);
modelmob = getMobility(modeln2D, f.logNdd, f.logNbk, f.logDelta, f.logDrough);

figure('Color','white');
hold on;
plot(modeln2D, modelmob,'Linestyle','-','Color','r','Linewidth',2);
plot(datan2D,datamob,'o','MarkerSize',10,'MarkerEdgeColor','k','Linewidth',2);
set(gca,'TickLabelInterpreter','latex','Fontsize',14,'YScale','log');
xlabel('$n_{2D}$ $[1{\rm E}10/{\rm cm}^2]$','Fontsize',20,'Interpreter','latex');
ylabel('$\mu$ $[10{\rm E}6{\rm \,cm}^2{\rm V}^{-1}{\rm s}^{-1}]$','Fontsize',20,'Interpreter','latex');
grid on;
lgd = legend({'Theory','Data'});
lgd.Location = 'best';
lgd.Interpreter = 'latex';
lgd.FontSize = 14;

fprintf(1,'Done!\n');


%Added by Arjun to write fit parameters to output file 'output.csv'
fileID = fopen('output.csv','w');
fprintf(fileID, 'Roughness(delta)(nm),Roughness(correlation length)(nm),Delta doping(surface charge)(/m2),Background doping(/m2)\n');
fprintf(fileID, '%e,%e,%e,%e\n',(10^(f.logDelta))*(1e9),(10^(f.logDrough))*(1e9),10^(f.logNdd),10^(f.logNbk));
%fprintf(fileID, '%e',10^(f.logNdd));