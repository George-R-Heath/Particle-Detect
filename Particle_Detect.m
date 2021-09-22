%A simple single particle height analysis tool to measure particle heights
%and positions

%see varibale 'stats' for compiled: frame, x, y, heights
%heights are background subtracted based on gaussian fitting to baseline

new = 1;         %set to 1 to load new tiff image file, 0 to keep current image
gus = 2;         %set level of gaussian filtering (for detection only!)
thresh = 0.1;    %set threshold height 
% or use auto threshold:
auto_thresh = 3; %set to 0 for no auto threshold otherwise set value to number of sd away from baseline
min_peak_sep = 2;  %reduce close neighbours with lower heights 
%only particles shown in the final plot with triangles 
%%

if new == 1
    clearvars -except gus thresh new auto_thresh min_peak_sep
[f,path] = uigetfile('*.tif');
f = fullfile(path,f); %tif (in nm) filename
info = imfinfo(f);  n = numel(info);
else
    new = 0;
end 

edg =2;
for i = 1:n
A = imread(f,i);
A = im2double(A);

b=find(A(:)~=0);% get the locations of nonzero. 

    if (gus~=0)
    Ag = imgaussfilt(A,gus);
    else 
    Ag = A;    
    end
 fo = fitoptions('method','NonlinearLeastSquares');
 [hy, x] =  hist(Ag(:),500);
 fsq  = fit(x',hy','gauss2',fo); 
    if auto_thresh >0
 thresh = fsq.b1 + fsq.c1*auto_thresh;
    else
    end
 
    At=Ag.*(Ag>thresh);
    if (gus~=0)
      At = imgaussfilt(At,gus);
    else 
    end
if any(At(:)) 
             
                sd=size(At);
                [x y]=find(At(edg:sd(1)-edg,edg:sd(2)-edg));
                
                % initialize outputs
                cent{i}=[];
                cent_h{i}=[];
                cent_map{i}=zeros(sd);
                cent_map_p_top{i} =zeros(sd);

                x=x+edg-1;
                y=y+edg-1;
                for j=1:length(y)
                    if (At(x(j),y(j))>=At(x(j)-1,y(j)-1 )) &&...   %find peaks
                            (At(x(j),y(j))>=At(x(j)-1,y(j))) &&...
                            (At(x(j),y(j))>=At(x(j)-1,y(j)+1)) &&...
                            (At(x(j),y(j))>=At(x(j),y(j)-1)) && ...
                            (At(x(j),y(j))>=At(x(j),y(j)+1)) && ...
                            (At(x(j),y(j))>=At(x(j)+1,y(j)-1)) && ...
                            (At(x(j),y(j))>=At(x(j)+1,y(j))) && ...
                            (At(x(j),y(j))>=At(x(j)+1,y(j)+1))
                       cent{i} = [cent{i} ;  y(j) ; x(j)];
                       cent_h{i} = [cent_h{i} ;A(x(j),y(j))]; %read height of unfiletered
                    end
                end
x = cent{i}(1:2:end);
y = cent{i}(2:2:end);
xy{i} = [x y];
 keep = [];
 
    %remove nearest neighbours
    if any(xy{i})
    for j = 1:numel(xy{i}(:,1))
    [ne,d] = knnsearch(xy{i},xy{i}(j,:),'k',20);
     pos = d<min_peak_sep;
        if sum(pos)>1
        [del_pos del_h] = max(cent_h{i}(ne(pos)));
        keep(j) = (ne(del_h));
        else
        keep(j) = j;
        end
    end
  
keep = unique(keep'.');
xy_clean{i} = xy{i}(keep,:);
height_clean{i} = cent_h{i}(keep,:);

     end
end
end
stats=[];heights=[];xy_all=[];
for i = 1:numel(height_clean)
 frame = [];frame(1:numel(height_clean{i}),1) = i;
xy_all = [xy_all; xy_clean{i}];
heights = [heights; height_clean{i}-fsq.b1];
stats = [stats; frame, xy_clean{i}, height_clean{i}-fsq.b1];
end

tiledlayout(2,3, 'Padding', 'none', 'TileSpacing', 'none');
nexttile
imagesc(A)
nexttile
imagesc(Ag)
nexttile

hist(Ag(:),500)
    hold on
    plot(fsq)
plot([thresh thresh], [0 max(hy)],'k','LineWidth',2)

nexttile
imagesc(At)
hold on
plot(cent{i}(1:2:end),cent{i}(2:2:end),'rx','MarkerSize',6)
plot(xy_clean{i}(:,1),xy_clean{i}(:,2),'ko')
nexttile
imagesc(A)
hold on
plot(cent{i}(1:2:end),cent{i}(2:2:end),'rx','MarkerSize',6)
plot(xy_clean{i}(:,1),xy_clean{i}(:,2),'ko')
colormap(jet)
nexttile
histogram(heights,30)


