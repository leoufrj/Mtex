%define os graos
grains = calcGrains(ebsd,'angle',10*degree)

ebsd('notIndexed').color = str2rgb('red');
ebsd('aust').color = str2rgb('yellow');

%plota o band contrast
figure(1);
plot(ebsd,ebsd.bc);
colormap gray;
hold on
plot(grains.boundary,'linewidth',0.5)
plot(ebsd('notIndexed'))
legend off
hold off

%plota as regioes indexadas e nao indexadas
figure(2);
plot(ebsd);

figure(3);
plot(grains.boundary)
hold on
plot(ebsd,'FaceAlpha',0.5);
hold off

%%
%preencher os locais nao indexados
[grains,ebsd.grainId,ebsd.mis2mean] = calcGrains(ebsd);
notIndexed = grains('notIndexed')

figure(4);
plot(notIndexed,log(notIndexed.grainSize ./ notIndexed.boundarySize))
mtexColorbar

% the "not indexed grains" we want to remove
toRemove = notIndexed(log(notIndexed.grainSize ./ notIndexed.boundarySize)<-0.8);

% now we remove the corresponding EBSD measurements
ebsd(toRemove) = [];

% and perform grain reconstruction with the reduces EBSD data set
figure(5)
[grains,ebsd.grainId,ebsd.mis2mean] = calcGrains(ebsd);
plot(grains,'lineWidth',1)

figure(6)
plot(grains('aust'),grains('aust').area)
mtexColorbar
hold on
plot(grains.boundary,'lineWidth',1)
plot(ebsd('notIndexed'))
legend off
hold off

%% plota o KAM

ebsd_g = ebsd.gridify;

kam = ebsd_g.KAM / degree;

figure(7)
plot(ebsd_g,kam)
caxis([0,5])
mtexColorbar
mtexColorMap Parula
hold on
plot(grains.boundary,'lineWidth',1)
plot(ebsd('notIndexed'))
legend off
hold off

%% plota a ipf

ebsd('notIndexed').color = str2rgb('black');
figure(8)
plot(ebsd('aust'),ebsd('aust').orientations)
hold on
plot(grains.boundary,'lineWidth',1)
plot(ebsd('notIndexed'))
legend off
hold off

ipfKey = ipfColorKey(ebsd('aust'));

%plota o triangulo estenográfico com os graos considerando aqueles com
%diametro maior que 6.

crit_area_1 = 10
crit_area_2 = 400
figure(9)
plot(ipfKey)
hold on
plot(ebsd('aust').orientations(grains('aust').area > crit_area_1),'Marker','+','MarkerFaceColor','none','MarkerEdgeColor','k','MarkerSize',10)
plot(ebsd('aust').orientations(grains('aust').area > crit_area_2),'Marker','^','MarkerFaceColor','w','MarkerEdgeColor','k','MarkerSize',8)
hold off

figure(10)
plot(grains('aust'),grains('aust').area > crit_area_1)
largeGrains = grains(grains.area > crit_area_2);
text(largeGrains,largeGrains.area)



%%
sS = slipSystem.fcc(ebsd('aust').CS)
sS = sS.symmetrise;

% rotate slip systems into specimen coordinates
sSLocal = grains('aust').meanOrientation * sS

% compute Schmid factor
sigma = stressTensor.uniaxial(vector3d.Z)
SF = sSLocal.SchmidFactor(sigma);

% take the maxium allong the rows
[SFMax,active] = max(SF,[],2);

% plot the maximum Schmid factor
figure(11)
plot(grains('aust'),SFMax,'micronbar','off','linewidth',1)
hold on
plot(ebsd('notIndexed'))
legend off
mtexColorbar
mtexColorMap Parula
hold off

%% Calculo do GOS (grain Orientation Spread)
%https://mtex-toolbox.github.io/GrainOrientationParameters.html

ori = ebsd(ebsd.grainId == 4).orientations
mis2mean = inv(grains(4).meanOrientation) .* ori
mis2mean = calcGROD(ebsd('aust'), grains)

GOS = ebsd('aust').grainMean(mis2mean.angle);

% plot it
figure(12)
plot(grains, GOS ./ degree)
caxis([0,2])
mtexColorbar('title','GOS in degree')
hold on
legend off
mtexColorMap Parula
hold off

%calcula o GAM
gam = ebsd('aust').grainMean(ebsd('aust').KAM);

figure(13)
plot(grains,gam./degree)
mtexColorbar('title','GAM in degree')
setColorRange([0,1])
hold on
plot(ebsd('notIndexed'))
legend off
mtexColorMap Parula
hold off

%% outra forma dos GOS - baseado no GROD
%https://mtex-toolbox.github.io/EBSDGROD.html

grod = ebsd('aust').calcGROD(grains);

figure(14)
plot(ebsd('aust'),grod.angle./degree,'micronbar','off')
mtexColorbar('title','misorientation angle to meanorientation in degree')
mtexColorMap Parula
hold on
plot(ebsd('notIndexed'))
plot(grains.boundary,'lineWidth',1)
legend off
mtexColorMap Parula
hold off

GOS2 = grainMean(ebsd('aust'), grod.angle);

figure(15)
plot(grains, GOS ./ degree)
mtexColorbar('title','GOS in degree')
mtexColorMap Parula
hold on
plot(ebsd('notIndexed'))
setColorRange([0,3])
plot(grains.boundary,'lineWidth',1)
legend off
mtexColorMap Parula
hold off

%% plotagens em geral

gB = grains.boundary('aust','aust');

figure(17)
plotAngleDistribution(gB.misorientation)

figure(16)
plot(grains('aust'),grains('aust').diameter > 7)

%%
figure(18);
plot(grains('aust'),grains('aust').meanOrientation,'FaceAlpha',0.15)
hold on
plot(grains.boundary('aust','aust'),grains.boundary('aust','aust').misorientation.angle./degree,'linewidth',2.5)
colormap plasma;
plot(ebsd('notIndexed'))
legend off
setColorRange([10,60])
mtexColorbar('title','misorientation angle')
hold off

%%