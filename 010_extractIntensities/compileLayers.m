function [mask,bms] = compileLayers(subject)

  all=[];
for i=0.5:0.5:25
	f=MRIread(sprintf('steps_%1.1fmm.nu.mgh',i));
f=f.vol;
all=cat(1,all, f);
end
allRaw=all;

dlmwrite([subject '_layers.csv'], all,';')
