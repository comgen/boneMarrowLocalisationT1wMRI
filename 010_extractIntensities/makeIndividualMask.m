function makeIndividualMask(subject)

m=MRIread('aseg.nii.gz');
m2=squeeze(m.vol(:,120,:));
subplot(1,2,1)
imagesc(m2); 
hold on
axis square
% find upper part of cerebellum
[y1,x1]=find(m2==47);
p1=find(y1==min(y1));
x1=x1(p1(1)); y1=y1(p1(1));
plot(x1, y1,'r*')
% find upper part accumbens
[y2,x2]=find(m2==58);
p2=find(y2==min(y2));
x2=x2(p2(1)); y2=y2(p2(1));
plot(x2, y2,'r*')

slope=(y2-y1) / (x2-x1);
yL=slope * (1 -x1) + y1

mask=ones(size(m2));
for i=1:size(m2,2)
 y=round(slope * i + yL);
mask(y:end,i)=0;
end 

subplot(1,2,2)
imagesc(mask)
axis square

set(gcf,'PaperPositionMode','auto','Position',[10 10 600 300]);                                                                                                                               
print(gcf, '-r300', '-dpng', [subject '_mask.png'])


m3=m.vol;
for i=1:size(m3,2)
 m3(:,i,:)=mask;
end     
m.vol=m3;                                                                                                                                                                                                                                                                                                                                                                                                      
MRIwrite(m, [subject '_individual_mask.mgz'],'mgz')    
