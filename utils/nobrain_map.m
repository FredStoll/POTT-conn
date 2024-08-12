function nobrain_map(gr,area2test,scales,col,fig,dir)

if nargin<3
    scales=[];
end
if nargin<4
    col = 'Blues';
end
if nargin<5
    fig = figure("Color",[1 1 1]);
end
if nargin<6
    dir = [];
end

if isobject(gr)

    %- normalization for line sizes
    norm = @(data,minW,maxW) (1-0.4)*(  (data-minW)/(maxW-minW)) +0.4;

    %- extract the size for the various lines
    if ~isempty(scales)
        min_w = scales(1); max_w = scales(2);
    else
        min_w = min(gr.Edges.Weight);
        max_w = max(gr.Edges.Weight);
    end
    gr.Edges.Weight(gr.Edges.Weight<0)=0;
    LWidths = 0.01+(6*((gr.Edges.Weight-min_w)/(max_w-min_w)));

    %- mean weight for every node
    meanWeight = [];
    for i = 1 : length(area2test)
        meanWeight(i) = mean(gr.Edges.Weight(sum(ismember(gr.Edges.EndNodes,area2test(i))')==1));
    end

else
    meanWeight = gr(:,1);
end

%-load the brain images (only boundaries + area colored)
brain = imread('Orbital_view_brain_lines_8areas.bmp');
brain_col = imread('Orbital_view_brain_8areas.bmp');

%- area location and color on the brain map
areamapping =  {'a12r'      'a12m'    'a12o'   'a12l'     'a11ml'    'a13l'      'a13m'      'LAI'      };
colormapping = [200 200 0 ; 0 255 0 ; 255 0 9 ; 0 0 255 ; 0 200 200 ; 0 100 0  ; 100 0 100 ; 0 0 150 ];
locXY =        [867 630 ;  1073 946 ;1313 1226; 1281 944;   739 742 ; 1070 1293; 832 1225  ; 995 1528];

%- transform value to get a white background and a range that match the measure to plot
if ~isempty(scales)
    minScale = scales(1); maxScale = scales(2);
else
    minScale = min(meanWeight)-((max(meanWeight)- min(meanWeight))*0.1); %0.1
    maxScale =  max(meanWeight)+((max(meanWeight)- min(meanWeight))*0.1); %0.7
end
brainsize = size(brain_col);
backgrd=((maxScale- minScale)*0.05);
Perf_map = (minScale+backgrd)*ones(brainsize(1:2));

for i = 1 : length(area2test)
    b = ismember(areamapping,area2test{i}); % match area with brain
    ar1 = (brain_col(:,:,1)==colormapping(b,1) & brain_col(:,:,2)==colormapping(b,2) & brain_col(:,:,3)==colormapping(b,3));

    Perf_map(ar1)=meanWeight(i);
end
boun = (brain(:,:,1)<100 & brain(:,:,2)<100  & brain(:,:,3)<100 );

Perf_map(boun)=minScale;
Perf_map(:)=0.01; %- if you comment that line, you can see brain pic

colorsJ = cbrewer('seq', col, 98);
colorsJ = [0 0 0 ; 1 1 1 ; 1 1 1 ; 1 1 1 ; 1 1 1 ; 1 1 1 ; 1 1 1 ; colorsJ];

colorsValues = [minScale:(maxScale-minScale)/(length(colorsJ)-1):maxScale];

%axes1 = axes('Parent',fig);
imagesc(Perf_map(5:end-5,5:end-5),[minScale maxScale]);colormap(colorsJ); %- 5 to remove some black border
axis off;
colorbar
set(gca,'FontSize',16);

if isobject(gr) & isempty(dir) %- only for graph plot with no arrows
    %- add lines, with size = connection strenght
    hold on
    LWidths_color = 1-repmat(norm(gr.Edges.Weight,minScale,maxScale),1,3);
    for i = 1 : length(LWidths)
        el_1 = find(ismember(areamapping,gr.Edges.EndNodes{i,1}));
        el_2 = find(ismember(areamapping,gr.Edges.EndNodes{i,2}));
        vect = [locXY(el_1,:) ; locXY(el_2,:)  ];
        line(vect(:,1),vect(:,2),'LineWidth',LWidths(i),'Color',LWidths_color(i,:));hold on
    end
    for ar = 1 : length(locXY)
        [minDistance, indexOfMin] = min(abs(colorsValues-meanWeight(ar)));
        plot(locXY(ar,1),locXY(ar,2),'.','MarkerSize',50,'Color',colorsJ(indexOfMin,:)) %- [66 146 200]/255
    end
elseif isobject(gr) & ~isempty(dir) %- only for graph plot with arrows

    %- add lines, with size = connection strenght
    hold on
    for i = 1 : length(LWidths)
        el_1 = find(ismember(areamapping,gr.Edges.EndNodes{i,1}));
        el_2 = find(ismember(areamapping,gr.Edges.EndNodes{i,2}));
        vect = [locXY(el_1,:) ; locXY(el_2,:)  ];
        line(vect(:,1),vect(:,2),'LineWidth',LWidths(i),'Color',[.6 .6 .6]);hold on
    end
    for ar = 1 : length(locXY)
        plot(locXY(ar,1),locXY(ar,2),'.','MarkerSize',50,'Color',[.8 .8 .8]) %- [66 146 200]/255
    end

    for i = 1 : length(dir(:,1))
        if dir(i,3)<0.05
            arrow_col = colorsJ(end-5,:);
        else
            arrow_col = colorsJ(15,:);
        end
        if dir(i,4)<0
            vect = [locXY(dir(i,1),:) ; locXY(dir(i,2),:)  ];
            quiver(vect(1,1),vect(1,2),vect(2,1)-vect(1,1),vect(2,2)-vect(1,2),1,'LineWidth',2,'Color',arrow_col);hold on
        else
            vect = [locXY(dir(i,2),:) ; locXY(dir(i,1),:)  ];
            quiver(vect(1,1),vect(1,2),vect(2,1)-vect(1,1),vect(2,2)-vect(1,2),1,'LineWidth',2,'Color',arrow_col);hold on
        end
    end
end

