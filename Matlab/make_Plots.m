function make_Plots(p,oc,t,n)
close all;
%% figure1
% plotting amuont of species
figure(1)

% interesting timePoint, which should get plotted
timePoint=[1 5 p.timeStep];

for i=1:3
    % species to plot i=[1=H2, 2=N2, 3=NH3]
    species=[i:p.n:p.N*p.n];
    
    % amount of  each species to different times stored in one matrix
    plotMat=flipud(n(timePoint,species)'); % flipud so compartment 150 is up and 1 is down
    % row: time
    % column: compartment
    
    % compute minimum and maximum value
    nMin=min(plotMat(:));
    nMax=max(plotMat(:));
    
    % each species gets a new subplot
    subplot(1,3,i)
    plot(plotMat',[1:p.N],'+')
    title(['species: ',p.speciesNames{i}])
    % y-axis
    yticks([1 25 50 75 100 125 150])
    yticklabels({'150','125','100','75','50','25','1'})
    ylabel('compartment [1]')
    % x-axis
    xticks([nMin mean([nMin,nMax]) nMax]);
    xticklabels({[num2str(nMin)],[num2str(mean([nMin,nMax]))],[num2str(nMax)]})
    xlabel('amount of species [kmol]')
    % legend
    legend([num2str(t(timePoint(1))),'h'],...
        [num2str(t(timePoint(2))),'h'],...
        [num2str(t(timePoint(3))),'h']);
    
end % for


%% figure2
% plotting temperature
figure(2)
% calculating temperature for each compartment at each time
A=zeros(p.n*p.N,p.N);
for i=1:p.N
    A(p.n*(i-1)+1:p.n*i,i)=[1;1;1;1];
end % for
% sum over all species in one compartment
nSum=n*A;
% Temperature
T=(oc.p * p.V * p.epsilon)./( p.N * nSum * p.R); % [12]
T=(T'); % transpose and switch upside down of T matrix
imagesc(T); % plot
colorbar; % show colorbar
title('Temperature [K]')

% y-axis
yticks([1 25 50 75 100 125 150])
yticklabels({'1','25','50','75','100','125','150'})
ylabel('compartment [1]')
% x-axis
for i=1:p.timeStep
labelxTicks{i}=[num2str(t(i)),'h'];
end % for
xticks([1:p.timeStep])
xticklabels(labelxTicks)
xlabel('time t [h]')



%% figure3
% plotting amount as heat map
figure(3)
desiredSpec=1;
nSpec=zeros(p.timeStep,p.N);
for i=1:p.N
nSpec(:,i)=n(:,desiredSpec+p.n*(i-1)); % 
end % for
% row: point in time
% column: compartment

imagesc(nSpec'); % plot
colorbar; % show colorbar
title(['amount of ',p.speciesNames{desiredSpec},' [kmol]'])

% y-axis
yticks([1 25 50 75 100 125 150])
yticklabels({'1','25','50','75','100','125','150'})
ylabel('compartment [1]')
% x-axis
for i=1:p.timeStep
labelxTicks{i}=[num2str(t(i)),'h'];
end % for
xticks([1:p.timeStep])
xticklabels(labelxTicks)
xlabel('time t [h]')

end % function
