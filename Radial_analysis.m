%code used to perform profile height analysis in x, y, xy and x-y

pleng = 34;     %set total profile length (even value in pixels)
loctor = 0.15;  %use peak located at central loctor fraction of total (max value = 1)


%%
clear p_x p_y p_fy
for j =1:numel(xy_all(:,1))
    %pleng = round(heights(j)).*4;
for jj=1:pleng
p_x{j}(jj,1)= xy_clean{i}(j,1)-(jj-pleng/2);
p_y{j}(jj,1)= xy_clean{i}(j,2)-(jj-pleng/2);
end
p_yf{j}= flip(p_y{j});
   
for jj=1:pleng
    if p_x{j}(jj,1)>0 && p_x{j}(jj,1)< numel(A(1,:))         
    prof_x(j,jj) = A(xy_clean{i}(j,2),p_x{j}(jj,1));
    else
    prof_x(j,jj) = 0;
    end
            
    if p_y{j}(jj,1)>0 && p_y{j}(jj,1)<numel(A(:,1))
    prof_y(j,jj) = A(p_y{j}(jj,1),xy_clean{i}(j,1));
    else
    prof_y(j,jj) = 0;
    end
            
    if p_x{j}(jj,1)>0 && p_x{j}(jj,1)< numel(A(1,:)) && p_y{j}(jj,1)>0 && p_y{j}(jj,1)<numel(A(:,1))         
    prof_xy1(j,jj) = A(p_y{j}(jj,1),p_x{j}(jj,1));
    else
    prof_xy1(j,jj) = 0;
    end
    
    if p_x{j}(jj,1)>0 && p_x{j}(jj,1)< numel(A(1,:)) && p_yf{j}(jj,1)>0 && p_yf{j}(jj,1)<numel(A(:,1))         
    prof_xy2(j,jj) = A((p_yf{j}(jj,1)),p_x{j}(jj,1));
    else
    prof_xy2(j,jj) = 0;
    end
    
end
%x
[pks, locs, wds] = findpeaks(nonzeros(prof_x(j,:)),'MinPeakProminence',0.15);
pos = (locs < pleng/2+pleng*loctor).*(locs > pleng/2-pleng*loctor);
pos=pos>0;

if sum(pos)>0
[Fwidth(1,j),ID] = max(wds(pos)); 
Fwidth_h(1,j) = pks(ID); 
else
Fwidth(1,j) = 0; 
Fwidth_h(1,j) = 0; 
end    
%y     
[pks, locs, wds] = findpeaks(nonzeros(prof_y(j,:)),'MinPeakProminence',0.15);
pos = (locs < pleng/2+pleng*loctor).*(locs > pleng/2-pleng*loctor);
pos=pos>0;

if sum(pos)>0
[Fwidth(2,j),ID] = max(wds(pos)); 
Fwidth_h(2,j) = pks(ID); 
else
Fwidth(2,j) = 0; 
Fwidth_h(2,j) = 0; 
end               

          [pks, locs, wds] = findpeaks(nonzeros(prof_xy1(j,:)),'MinPeakProminence',0.15);

pos = (locs < pleng/2+pleng*loctor).*(locs > pleng/2-pleng*loctor);
pos=pos>0;

if sum(pos)>0
[Fwidth(3,j),ID] = max(wds(pos)); 
Fwidth(3,j) =Fwidth(3,j)*2^0.5;
Fwidth_h(3,j) = pks(ID); 
else
Fwidth(3,j) = 0; 
Fwidth_h(3,j) = 0; 
end           
                                
                    
[pks, locs, wds] = findpeaks(nonzeros(prof_xy2(j,:)),'MinPeakProminence',0.15);

pos = (locs < pleng/2+pleng*loctor).*(locs > pleng/2-pleng*loctor);
pos=pos>0;

if sum(pos)>0
[Fwidth(4,j),ID] = max(wds(pos)); 
Fwidth_h(4,j) = pks(ID); 
Fwidth(4,j) =Fwidth(4,j)*2^0.5;
else
Fwidth(4,j) = 0; 
Fwidth_h(4,j) = 0; 
end                                
end
Fwidth(Fwidth==0) = NaN;
Rmax = max(Fwidth,[],1);
Rmin = min((Fwidth),[],1);
Rmean = mean(Fwidth,1,'omitnan');

imagesc(A)
colormap(jet)
hold on
plot(xy_clean{i}(:,1),xy_clean{i}(:,2),'ko')
for j =1:numel(heights)
plot(p_x{j}(:,1),p_y{j}(:,1))
plot(p_x{j}(:,1),p_yf{j}(:,1))
str = num2str(j);
text(xy_clean{i}(j,1),xy_clean{i}(j,2),str,'FontSize',14,'Color','w')
end

figure(2)
plot(Rmax,heights,'.')
hold on
plot(Rmean,heights,'o')
ylim([0 12])
stats = [ heights, Rmean'];