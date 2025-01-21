%close all
%clear all

% This m-file is for making movies of wave propagation.

%%%%%%%%%%%%%%%%%%%%%%%% INPUT PARAMETER %%%%%%%%%%%%%%%%%%%%%%%%%
% Here the basic name of the binary snapshot files must be given:
% (The default extension is *.bin. The increasing number of the snapshot
% is added to the basic name, e.g. if ten snapshots were computed
% and the basic name is snap, then the filenames for the snapshots
% are snap1.bin, snap2.bin, ..., snap10.bin.
file1='../par/snap/hh_ve_1.bin.vx';
file2='../par/snap/hh_ve_1.bin.vy';


% gridsize and grid spacing (as specified in parameter-file) 
NX1=1; NX2=400;
NY1=1; NY2=400; 
IDX=1; IDY=1;
dh=2.0;

% time increment for snapshots:
TSNAPINC=0.01; TSNAP1=0.01;
FW=40.0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% grid size
nx=NX2-NX1+1; ny=NY2-NY1+1;

% plot range and increment
xp1=NX1*dh; xp2=NX2*dh; yp1=NY1*dh; yp2=NY2*dh; 

% Computing range for axis and subscript range for the movie
x=xp1:IDX*dh:xp2;
y=yp1:IDY*dh:yp2;


firstframe=13;
lastframe=13;
snapno=0;

caxis_value=1.0e-10;

%load /home/tbohlen/util/seismic.map
colormap(gray)

% load model
%  file='model/ktb10.pi';
%  disp([' loading file ' file]);
%  fid=fopen(file,'r','ieee-le');
%  model=fread(fid,[ny,nx],'float');
%  fclose(fid);
 %model=model(1:2:size(model,1),1:2:size(model,2));
   
 disp(['opening file ' file1]);
 fid1=fopen(file1,'r','ieee-le');
 
 disp(['opening file ' filediv]);
 fid2=fopen(file2,'r','ieee-le');

 for i=firstframe:1:lastframe,
 
   disp(['loading snapshot no ',int2str(i)]);
   % loading data:
    tsnap=(i-1)*TSNAPINC+TSNAP1;
   offset=4*nx*ny*(i-1);
   fseek(fid1,offset,-1);
   fseek(fid2,offset,-1);
   vx=fread(fid_div,[ny,nx],'float');
   vy=fread(fid_rot,[ny,nx],'float');
   
   subplot(121), 
      imagesc(x,y,vx);  
		 hold on
    	caxis([-caxis_value caxis_value]);
				
		set(text(1.4,-0.0,['T=',sprintf('%1.2f',tsnap),' s']),...
		 'FontSize',12,'FontWeight','bold','color','k');
       title('Vx');
       xlabel('X (m)')
       ylabel('Y (m)')
       set(gca,'DataAspectRatio',[1 1 1]);
       set(get(gca,'title'),'FontSize',12);
       set(get(gca,'title'),'FontWeight','bold');
       set(get(gca,'Ylabel'),'FontSize',12);
       set(get(gca,'Ylabel'),'FontWeight','bold');
       set(get(gca,'Xlabel'),'FontSize',12);
       set(get(gca,'Xlabel'),'FontWeight','bold');
       set(gca,'FontSize',12);
       set(gca,'FontWeight','bold');
       set(gca,'Linewidth',1.0);
       set(gca,'Box','on');
    
		 
	
   subplot(122), 
	%set(gcf,'Color',[1 1 1])
       imagesc(x,y,vy);  
	   caxis([-caxis_value caxis_value]);


       title('Vy')
         xlabel('X (m)')
       %ylabel('Y (m)')
       set(gca,'DataAspectRatio',[1 1 1]);
       set(get(gca,'title'),'FontSize',12);
       set(get(gca,'title'),'FontWeight','bold');
       set(get(gca,'Ylabel'),'FontSize',12);
       set(get(gca,'Ylabel'),'FontWeight','bold');
       set(get(gca,'Xlabel'),'FontSize',12);
       set(get(gca,'Xlabel'),'FontWeight','bold');
       set(gca,'FontSize',12);
       set(gca,'FontWeight','bold');
       set(gca,'Box','on');
       set(gca,'Linewidth',1.0);
     
 
	%brighten(0.5)
	

	
   % Saving the snapshot:
    
    %jpgfile=['eps/ktb11_',int2str(i),'.jpg'];
    %eval(['print -djpeg100 ' jpgfile]);

 end
fclose(fid1);
fclose(fid2);











