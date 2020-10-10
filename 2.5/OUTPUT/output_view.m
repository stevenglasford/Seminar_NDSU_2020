function output_view()

fprintf('ROCHJ v2.5.3 - matlab/octave postprocessing. ');
clear %all
clear xmini xmaxi ymini ymaxi
clear GRAPH_CONTOUR GRAPH_CONTOUR_EX
clear x
%close all

global TEST  %- default is TEST=0;  use TEST=1 for specific tests.

global level_set 
global MOVIE_PAUSE MOVIE_nbSaves MOVIE_INPUTNAME MOVIE_INPUTNAMEex
global MOVIE_PLOT_REACHABLE
global ITERMAX
global ALPHA BETA

global nn VALUE_PB
global OCTAVE
global dim filePrefix
global T DT
global xmin xmax cdd dx MESH
global Zmin Zmax	%- to be adjusted depending of [min,max] values of data
global dimcoupe

global VFile0		%- used in output_view3d()
global VexFile0
global V_EXACT
global VIDEO videoFile VIDEOfig
global xmin3d xmax3d nn3d
global nbSaves nbSavesCoupe

% 2018 for plot_reachable()
global PLOT_REACHABLE PLOT_CONTOUR PLOT_TOPT
global color_contour color_contour_ex lw_contour %- linewidth for contours
global legend_pos legend_tit_ex legend_tit_scheme
global xmini xmaxi ymini ymaxi empty
global AXIS_EQUAL
global xmesh ymesh n1 n2
global x val val_ex
global GRAPH_CONTOUR GRAPH_CONTOUR_EX GRAPH_CONTOUR_VF0
global tit1 tit2 tit3 %- xlabel and ylabel for graph
global TITLE_FILE
global PLOT_3d_ex
global BINARY
global FORMAT_FULLDATA
global dx1 dx2 
global xmini xmaxi ymini ymaxi
global xmeshi ymeshi


%- TEST variable
TEST=0;

%- OCTAVE detect	%- Default is OCTAVE=0. Special commands for octave : (graphs, legends...)
OCTAVE = (exist ('OCTAVE_VERSION', 'builtin') > 0); if OCTAVE fprintf('(OCTAVE DETECTED)\n'); else fprintf('(MATLAB DETECTED)\n'); end;

%-------------
%- fig 1:
%-------------
PLOT_CONTOUR   =1;	%- plot contours
PLOT_REACHABLE =1;	%- plot reachable region (points where scheme value <= level_set)
level_set      =0.00;	%- 
TARGET_PLOT    =1;	%- target plot = green dots where initial value VF0.dat is <= 0. For instance use in the data file:
			% const int     SAVE_VF_ALL        = 1; 
			% const int     SAVE_VF_ALL_STEP   = 10000;  //- or smaller
			%- in order to get a file VF0.dat
			%- 2019: for High Definition screens : 
factorMarkerSize=5*1;   %- 2019: might be necessary to ajust (From 5*1 to 5*2, or more) - this is for 'MarkerSize' in reachable set plots, etc.
factorLineWidth =  1;   %- 2019: might be necessary to ajust (From 1   to 5, or more) - this is for linewidth


%-------------
%- fig 2:
%-------------
PLOT_3d    =1;		%- Set figure 2 for plots: plot of VF or topt 
PLOT_TOPT  =0;		%- PLOT_3d=1 and PLOT_TOPT=1 ==> use topt.dat (optimal time function) instead of VF.dat  [toptFile0 / VFile0: defined below]


%-------------
%- fig 3:
%-------------
PLOT_3d_ex =1;		%- plot 3d of Vex (if furthermore compute_vex=1..)

%-------------
%- for figs 1 and 2 (plus, possibly fig 4 and 10) :
%-------------
AXIS_EQUAL=1;		%- will make axis equal for fig 1 (or fig 1 and 2); default is axis square
%ALPHA=65; BETA=35;	%- angles for the 3d plot (figs 2 and 3)
%ALPHA=25; BETA=50;	%- angles for the 3d plot (figs 2 and 3)
%ALPHA=45; BETA=0;	%- angles for motion dyn=(1,1) and u0=u0(x+y);
%ALPHA=-38; BETA=33;	%- angles standards
%ALPHA=0;  BETA=0;
%ALPHA=-70; BETA=30;
ALPHA=50; BETA=30;

%ALPHA=atan(2/1)*360/(2*pi); BETA=0;	%- direction (1,2) , angle is artg(2/1)*360/(2.pi)=63.43
%ALPHA=atan(1/2)*360/(2*pi); BETA=0;	%- direction (1,2) , angle is artg(2/1)*360/(2.pi)=63.43
%ALPHA=0; BETA=0;

%-------------
%- for figs 4 and 10 
%-------------
%TRAJECTORY=1;		%- plots of traj-xx.dat (first 2 components)
SUPPLEMENT_TRAJECTORY=0;%- each state (and possibly each control) is plotted vs time  (see also PRINT_TRAJ and PRINT_CONTROL parameters inside)

%-------------
%- Special parameters I : Movie, Video
%-------------
MOVIE=0; 		%- 0/1 (default 0): 1d/2d movie (can be used in particular if pb is 1d (if SAVE_VFALL=1), or if pb is 2d
MOVIE_PAUSE=0; 		%- 0/1 1:for graphic pause (0:no need to type "return" after each plot).
MOVIE_PLOT_REACHABLE=1;	%- 0/1 (default 0): surface graph(0) or level set (1)-in fig 1- (if level set then may plot also exact level set if available).

VIDEO=0;		%- Default is VIDEO=0. If VIDEO=1, creates a video .avi file of the animation (->videoFile)
VIDEOfig=1;             %- Choose figure for the movie (default is 2 for VF;  use 1 for reachable sets...)
videoFile='video2.avi'; %- avi filename

%- note: for MOVIE mode, FORCE TO SWITCH to VF files if dim==2, nbSavesCoupe=0 and nbSaves>0 
%-       however nbSavesCoupe should be > 0, we should rather plot from  "coupexx" files and not from "VFxx.dat" files.
MOVIE_nbSaves	  ='nbSavesCoupe'; %- The total number of iterations: should be one of {'nbSaves','nbSavescoupe'}
MOVIE_INPUTNAME   ='coupe';	%- Choice of input data  : from "coupe{it}.dat" 
%MOVIE_INPUTNAME   ='coupetopt';	%- Choice of input data  : from "coupe{it}.dat" 
MOVIE_INPUTNAMEex ='coupeex';	%- Choice of input data  : from "coupeex{it}.dat" 

%- uncomment here to change the filenames for movie plots
%MOVIE_nbSaves     ='nbSaves';	%- The total number of iterations: should be one of {'nbSaves','nbSavescoupe'}
%MOVIE_INPUTNAME   ='VF';	%- Choice of input data  : 'VF' = from "VF{it}.dat"  files
%MOVIE_INPUTNAMEex ='Vex';	%- Choice of input data  : 'Vex'= from "Vex{it}.dat"  files

%-------------
%- Special parameters II : 3d isoview
%-------------
PLOT_ISOVIEW=1;		%- Default is 0. Set 1 to obtain a 3d isoview instead of surface plot 
			%- Default is on figure 2.

%VALUE_PB=0; 		%- (Default : 0) This parameters is now read in data.data. 
			%- Indicates plot from VF.dat (if VALUE_PB=0) or from Value.dat (if VALUE_PB=1);

%-------------
%- Default parameters : filenames 
%-------------
%-  To adapt to the names given in stdafx.h
VFile0     ='coupe';		%- TODO 2018 : adapter aux noms de stdafx.h (sans .dat)
VexFile0   ='coupeex';
toptFile0  ='coupetopt';
ValueFile0 ='Value';
%VexFile0   ='Vex';

%-------------
%- Default parameters : labels pour les axes des coupes
%-------------
coupelabels={'x_1','x_2','x_3','x_4','x_5','x_6'};
%Front_ex=@Front_ex_xpp;
UPDATE_AXIS=1; %- parameter to tell if the axis have to be updated

%-----------------------------------------
%- Notes relative to some of the examples:
%-----------------------------------------
%- data_basicmodel.h:	set PLOT_3d=1; PLOT_TOPT=0 (or 1); TRAJECTOIRE=1; 
%- data_FD_zermelo.h:	set PLOT_3d=1; PLOT_TOPT=0 or 1;  AXIS_EQUAL=0 or 1;
%- data_FD_stat.h:	set PLOT_3d=1; PLOT_TOPT=0;  AXIS_EQUAL=0;
%- data_FD_dubinscar.h: 
%- data_FD_value.h:	set VALUE_PB=1; PLOT_TOPT=0 (or 1); TRAJECTORY=1;
%-     assumes also that following parameters are used at the begining of main.cpp file:
%-     VALUE_PB = 1; 
%-     SAVE_VALUE_FINAL = 1; 

%----------------------------------------------------
%- loading data.dat and initialisation of parameters
%----------------------------------------------------
%global nn VALUE_PB


filePrefix='';
[zz1,zz2]=unix('ls filePrefix.dat');
if zz1==0   % there is a filePrefix.dat file and we need to get the prefix name
  FILE=fopen('filePrefix.dat'); 
  filePrefix=fscanf(FILE,'%s');
  fclose(FILE);
end
dataFile=[filePrefix, 'data.dat'];

FILE=dataFile;
fprintf('.. Reading file: "%s" ...',FILE); 
[status, result]=unix(strcat('ls ./',FILE));
if status==2; 
  fprintf(strcat('.. No file ./',FILE,'. Abort.\n'));
  return; 
else
  param=load(dataFile,'-ascii'); %- problem & mesh parameters
  fprintf('\n');
end

i=1;
dim =param(i); i=i+1;
MESH=param(i); i=i+1;
dx=zeros(1,dim); xmin=zeros(1,dim); xmax=zeros(1,dim); nn=zeros(1,dim);
for d=1:dim
  dx(d)=param(i); i=i+1;
end
for d=1:dim
  xmin(d)=param(i); i=i+1;
  xmax(d)=param(i); i=i+1;
  fprintf('d=%i, xmin(d)=%6.3f, xmax(d)=%6.3f\n',d,xmin(d),xmax(d));
end
for d=1:dim
  nn(d)=param(i); i=i+1;
end

T=param(i); i=i+1;
V_EXACT = param(i); i=i+1; %- put 1 if we know the exact front (otherwise 0);
                           %- then the corresponding Front_exXX.m should be programmed.
compute_vex=V_EXACT;
PLOT_3d_ex=min(PLOT_3d_ex,compute_vex); 

%- SET_COUPE should be 0 (no coupe) or 1 (use of coupe)
dimcoupe = param(i); i=i+1; %-here param(i) contains the dimension of the cut

SET_COUPE=0; 
if dimcoupe>=1
  SET_COUPE=1; %- 2017: should be always SET_COUPE=1; excepted if dimcoupe==0
end

if (dimcoupe==3 || dim==3)
  xmin3d=xmin(1:3);
  xmax3d=xmax(1:3);
  nn3d=nn(1:3);
  %- This may be corrected afterwards if other coupes_dims(.) are choosen
end

if (1)
  fprintf('T=%6.3f\n', T);
  fprintf('V_EXACT=%2i \n',V_EXACT);
  fprintf('SET_COUPE=%2i [default should be 1: plots using "coupexx.dat files (excepted for dim=1)]\n',SET_COUPE);
  for d=1:dim; fprintf('nn(%i)=%3i; ',d,nn(d)); end; 
  fprintf('dim=%i; dimcoupe=%i; ',dim,dimcoupe);
  fprintf('\n');
end 



if (1)  % SET_COUPE == 1
  coupe_dims=zeros(dim,1); coupe_vals=zeros(dim,1);
  for d=1:dim,
    coupe_dims(d) = param(i); i=i+1;
    fprintf('coupe_dims(%i)=%2i; ',d,coupe_dims(d));
  end
  fprintf('\n');

  for d=1:dim 
    coupe_vals(d) = param(i); i=i+1; %Nb : only values for which coupe_dims(d)=1 are meaningful
  end
else
  %- default is 2d / MAR 2014
  cdd=[1 2];
end

if SET_COUPE == 1
  cdd = find(coupe_dims==1); %- dimensions kept fater cutting.
  %cv = find(coupe_dims==0); %- dimensions cutted
end


if dimcoupe==2
  fprintf('Directions of the cuts: cdd=[%i %i]\n',cdd(1),cdd(2));
  d1=cdd(1);     d2=cdd(2);
  n1=nn(cdd(1)); n2=nn(cdd(2)); 
  tit1=coupelabels(d1); tit2=coupelabels(d2);
  %- in case isoview 3d is used while coupe are also used:
  if (dim==3); tit3='x_d'; end;

  %xmini=xmin(d1); xmaxi=xmax(d1); 
  %ymini=xmin(d2); ymaxi=xmax(d2); 

  %fprintf('(n1=%2i;  n2=%2i;)\n',n1,n2); 
  %dx=dx(cdd); 
  %xmin=xmin(cdd); 
  %xmax=xmax(cdd);

elseif (dimcoupe==3) %- TODO 2018 
  fprintf('Directions of the cuts: cdd=[%i %i %i]\n',cdd(1),cdd(2),cdd(3));
  d1=cdd(1);     d2=cdd(2);     d3=cdd(3);
  %n1=nn(cdd(1)); n2=nn(cdd(2)); n3=nn(cdd(3)); 
  xmin3d=xmin([d1,d2,d3]);
  xmax3d=xmax([d1,d2,d3]);
  nn3d  =nn  ([d1,d2,d3]);
  tit1=coupelabels(d1); tit2=coupelabels(d2); tit3=coupelabels(d3);

elseif (dimcoupe==1) 
  %d1=1; d2=0;
  %n1=1; n2=0;
  d1=cdd(1); d2=0;
  n1=nn(cdd(1)); n2=0;
  tit1=coupelabels(d1); 
  %dx=dx(cdd); 
  %xmin=xmin(cdd); 
  %xmax=xmax(cdd);
else
  fprintf('Pb: dimcoupe(=%i)<>1 and <>2 not working. Abort.\n',dimcoupe);
  if dimcoupe==3; fprintf('Note: dimcoupe(=%i) could work with isoview.\n',dimcoupe); end
  return
end
if length(cdd)>2
  fprintf('.. the cut is more that 2d: will retain only first two directions for plots\n');
  cdd=cdd(1:2);
end   
if (dim~=1 && length(cdd)~= dimcoupe)
  fprintf('.. Error: dim<>1 and length(cdd)<>dimcoupe. Abort.\n'); 
  return; 
end

%%------------------- 
if 0
  dx=dx(cdd); 
  xmin=xmin(cdd); 
  xmax=xmax(cdd);
end
%%------------------- 

save_vfall=param(i); i=i+1;
nbSaves=param(i); i=i+1; 


nbTraj=param(i); i=i+1;
fprintf('nbTraj (nombre de trajectoires)=%2i\n',nbTraj);
TRAJECTORY=(nbTraj>=1);

VALUE_PB=param(i); i=i+1;
%VALUE_PB=1; fprintf('Warning: FORCING VALUE_PB=1 HERE FOR PLOTS'); %- ICI
if(VALUE_PB==1)
  VFile0=ValueFile0;
  fprintf('VALUE_PB=%2i : special settings for this case\n',VALUE_PB);
  fprintf(' .. setting VFile = %s ..\n',VFile0);
end

%- 2018
SAVE_COUPE_ALL=param(i); i=i+1;
nbSavesCoupe=param(i); i=i+1; 
if (nbSavesCoupe==0 && nbSaves>0 && dim<=2)
  MOVIE_nbSaves     ='nbSavesCoupe'; %- The total number of iterations: should be one of {'nbSaves','nbSavescoupe'}
  MOVIE_INPUTNAME   ='coupe';	%- Choice of input data  : from "coupe{it}.dat" 
  MOVIE_INPUTNAMEex ='coupeex';	%- Choice of input data  : from "coupeex{it}.dat" 
  %---------------------------
  %- force switch to VF files
  %---------------------------
  %fprintf(' ** WARNING : no coupexx.dat savings ==> force to "VF" names for MOVIE mode ..\n')
  %MOVIE_nbSaves     ='nbSaves';	%- The total number of iterations: should be one of {'nbSaves','nbSavescoupe'}
  %MOVIE_INPUTNAME   ='VF';	%- Choice of input data  : 'VF' = from "VF{it}.dat"  files
  %MOVIE_INPUTNAMEex ='Vex';	%- Choice of input data  : 'Vex'= from "Vex{it}.dat"  files
end

%- 2018
%- Supplementary variables (if exists...)
if i<=size(param,1);
  BINARY=param(i); i=i+1;
  fprintf('BINARY=%i\n',BINARY);
else
  BINARY=0;
  fprintf('Warning: param(i), i=%i does not exists. Continuing with BINARY=0 (default)\n',i);
end

if i<=size(param,1);
  FORMAT_FULLDATA=param(i); i=i+1;
  fprintf('FORMAT_FULLDATA=%i\n',FORMAT_FULLDATA);
else
  FORMAT_FULLDATA=1;
  fprintf('Warning: param(i), i=%i does not exists. Continuing with FORMAT_FULLDATA=1 (default)\n',i);
end
if (BINARY==1 && FORMAT_FULLDATA==1)
  FORMAT_FULLDATA=0;
  fprintf('Warning: BINARY MODE => forcing FORMAT_FULLDATA=0\n');
end

%----------------------------
%- END READING of param()
%----------------------------

if strcmp(MOVIE_nbSaves,'nbSaves')
  ITERMAX=nbSaves;
  if save_vfall==0
    ITERMAX=0;
  end
  fprintf('ITERMAX(=nbSaves)=%2i\n', ITERMAX);
elseif strcmp(MOVIE_nbSaves,'nbSavesCoupe')
  ITERMAX=nbSavesCoupe; %- 
  if SAVE_COUPE_ALL==0
    ITERMAX=0;
  end
  fprintf('ITERMAX(=nbSavesCoupe)=%2i\n', ITERMAX);
else 
  fprintf('MOVIE_nbSaves = %s not defined. Setting ITERMAX=0.\n',MOVIE_nbSaves);
  ITERMAX=0;
end

%-------------------
%- COMPLEMENT 2016
%-------------------
dtFile=[filePrefix, 'Dt.dat'];
Dt=load(dtFile,'-ascii'); %- problem & mesh parameters

% DT step for movie mode
%global DT
DT=Dt;

if strcmp(MOVIE_nbSaves,'nbSaves')
  if save_vfall>0
    DT=T/ITERMAX;
    %- HERE WE COULD ALSO ADAPT with Dtplot
  end
elseif strcmp(MOVIE_nbSaves,'nbSavesCoupe')
  if SAVE_COUPE_ALL>0
    DT=T/ITERMAX;
    %-------------------
    %- COMPLEMENT 2018 : the used Dtplot is adapted to SAVE_COUPE_ALL_STEP
    %- (ie is such that Dtplot = Dt*SAVE_COUPE_ALL_STEP)
    %-------------------
    dtplotFile=[filePrefix, 'Dtplot.dat'];
    [zz1,zz2]=unix(strcat('ls \ ',dtplotFile));
    if zz1==0  
      Dtplot=load(dtplotFile,'-ascii'); %- problem & mesh parameters
      fprintf('Loading Dtplot from %s: ..\n',dtplotFile);
      DT=Dtplot;
    else
      fprintf('File %s not found. Skip. Will use DT=%5.3f instead\n',dtplotFile,DT);
    end
  end
else 
   fprintf('MOVIE_nbSaves = %s not defined. DT set to 0.\n',MOVIE_nbSaves);
end

%-------------------------
%- begining graphics
%-------------------------

%-------------------------
%- LEGEND PARAMETERS & MISCELLANEOUS
%-------------------------
if dimcoupe==1
  %- LEGEND PARAMETERS (FOR 1d)
  legend_pos='northeast'; 
  %legend_pos='southeast';
  %legend_tit_ex='Trajectory';
  legend_tit_ex='Exact';
  legend_tit_scheme='Scheme';

elseif (dimcoupe==2)
  %- LEGEND PARAMETERS (FOR 2d)
  %legend_pos='northeast';
  legend_pos='southeast';
  %legend_tit_ex='Trajectory';
  legend_tit_ex='Exact';
  legend_tit_scheme='Reachable set';

  color_contour   ='r--';
  color_contour_ex='k--';
  lw_contour=1*factorLineWidth;	%- linewidth for contours

elseif (dimcoupe==0)
  fprintf('Legends: the case dimcoupe=%i is not specified ==> proceed as for dim=1;\n',dimcoupe); 
  legend_pos='northeast'; 
  legend_tit_ex='Exact';
  legend_tit_scheme='Scheme';
else
  fprintf('Legends: the case dimcoupe>=3 is not specified. Abort.\n');
  return; 
end



%------------------
%- MOVIE 1d or 2d:
%------------------
if (VIDEO && OCTAVE)
  fprintf('VIDEO(%i) mode not programmed in OCTAVE. Switching to VIDEO=0.\n',VIDEO);
  VIDEO=0;
end

%-----------
%- MOVIE 1d 
%-----------
if (dim==1 || dimcoupe==1) && MOVIE
  %%- ajusting for 1D only
  d1=cdd(1);
  dx1=dx(d1);
  xmin=xmin(d1); 
  xmax=xmax(d1);
  %if MESH==1;
  %  x=xmin + (0:n1)'*dx;
  %end
  %if MESH==0;
  %  x=xmin + ((1:n1)-0.5)'*dx;
  %end
  xmesh=xmin(d1) + (0:n1-1)'*dx1 + (1-MESH)*dx1/2;
  x=xmesh;

  Zmin=-1.5; Zmax=3.5;
  figure(1);
  drawperso();
  output_view_movie1d;
  input('exiting movie1d mode [waiting for "Enter"]');
  %return;
end

%-----------
%- MOVIE 2d 
%-----------
if ((dim>=2 && dimcoupe==2) || (dim==2 && dimcoupe==0)) && MOVIE
  %%- ajusting for 2D only
  %------------------------------
  %- mesh: following assumes dimcoupe=2
  %------------------------------
  xmini=xmin(d1); xmaxi=xmax(d1); 
  ymini=xmin(d2); ymaxi=xmax(d2); 
  dx1=dx(d1);     dx2=dx(d2);
  n1 =nn(d1);     n2 =nn(d2);

  %- TEST - ENCOURS 2019 - ISO-MAILLAGE 2d avec Dx(j)=1:
  %TEST=1;
  if (TEST)
    dx1=1.0;
    dx2=dx1;
    xmini=0; xmaxi=xmini+(n1-MESH)*dx1;
    ymini=0; ymaxi=ymini+(n2-MESH)*dx2;
  end

  xmesh=xmini + (0:n1-1)'*dx1 + (1-MESH)*dx1/2;
  ymesh=ymini + (0:n2-1)'*dx2 + (1-MESH)*dx2/2;
  xmeshi=xmesh;   ymeshi=ymesh;
  %- constructs x s.t. x(i,1) x(i,2) contain the list of x1 and x2 mesh points of the full 2d grid.
  ox=ones(size(xmesh));
  oy=ones(size(ymesh));
  xtemp=xmesh*oy'; xtemp=xtemp(:);
  ytemp=ox*ymesh'; ytemp=ytemp(:);
  x=[xtemp, ytemp]; %- liste 2d de toutes les coordonnees x et de toutes les coordonnees y


  figure(1);
  drawperso();
  %Zmin=-0.2; Zmax=0.2;
  %Zmin=-0.5; Zmax=0.5;
  output_view_movie2d();
  UPDATE_AXIS=0;
  input('exiting movie2d mode [waiting for "Enter"]');
  %return;
end


if (dim==1 || dimcoupe==1)
  %-------------------------------------
  %- Special part for case dim=1
  %-------------------------------------

  %%- ajusting for 1D only
  d1=cdd(1);
  dx=dx(d1);
  xmin=xmin(d1); 
  xmax=xmax(d1);
  %xmesh=xmin(d1) + (0:n1-1)'*dx1 + (1-MESH)*dx1/2;
  x=xmin + (0:n1-1)'*dx  + (1-MESH)*dx/2;
  %if MESH==1;
  %  x=xmin + (0:n1)'*dx;
  %end
  %if MESH==0;
  %  x=xmin + ((1:n1)-0.5)'*dx; % FALSE
  %end

  %- fig 1:
  if dim==1; fprintf('dim=1 !\n'); end;
  if dim>=2; fprintf('CAREFUL: COUPE is 1D !\n'); end;


  %- default files are "coupe" files, except if dim==1
  FILE  =[filePrefix 'coupe.dat'];
  FILEEX=[filePrefix 'coupeex.dat'];
  %FILE  =[filePrefix VFile0 '.dat'];
  %FILEEX=[filePrefix VexFile0 '.dat'];
  TITLE_FILE=FILE;

  %if dim==1;
  %  %FILE  =[filePrefix 'VF.dat'];
  %  %FILEEX=[filePrefix 'Vex.dat'];
  %  %FILE  =[filePrefix VFile0 '.dat'];
  %  %FILEEX=[filePrefix VexFile0 '.dat'];
  %end
  %fprintf('loading %s ...\n',FILE);

  %- New:
  %- ICI
  message='fig1';
  [val,err]=loadFile(FILE,message,n1,dimcoupe);
  if (err==1); return; end;
  %- calcul de x,nn (=n1),
  
  % if BINARY
  %   dataF=fopen(FILE,'rb');
  %   data=fread(dataF, [n1,1] ,'double');
  % else
  %   data=load(FILE); 
  % end

  %----------------------------------------------------------
  %- MAR 2014: adapted for case dim>=2 and dimcoupe=1
  %----------------------------------------------------------
  % if MESH==1;
  %   nn=size(data,1)-1;
  %   x=xmin + (0:nn)'*dx;
  % end
  % if MESH==0;
  %   nn=size(data,1);
  %   x=xmin + ((1:nn)-0.5)'*dx;
  % end
  % if FORMAT_FULLDATA==1
  %   val=data(:,2);		%- here if FORMAT_FULLDATA=1
  % else
  %   val=data(:,1); 		%- here if FORMAT_FULLDATA=0
  % end


  figure(1);
  clf;
  drawperso();

  %color_app='';
  %if OCTAVE; color_app='b*-'; color_exa='k*-';
  %else;      color_app='b.-'; color_exa='k.-';
  %end
  color_app='b.-'; 
  color_exa='k.-';
  global GRAPH GRAPH_EX 
  GRAPH=plot(x,val,color_app); hold on;

  if  V_EXACT==1

    message='(valex)';
    [valex,err]=loadFile(FILEEX,message,n1,dimcoupe);
    if (err==1); return; end;

    % fprintf('loading %s ...\n',FILEEX);
    % %data=load(FILEEX);
    % if BINARY
    %   dataF=fopen(FILEEX,'rb');
    %   data=fread(dataF, [n1,1] ,'double');
    % else
    %   data=load(FILEEX); 
    % end
    %     
    % %valex=data(:,2);
    % if FORMAT_FULLDATA==1
    %   valex=data(:,2);		%- FORMAT_FULLDATA=1
    % else
    %   valex=data(:,1); 		%- FORMAT_FULLDATA=0
    % end

    GRAPH_EX=plot(x,valex,color_exa); 
  end

  %- 2020 complement error
  COMPLEMENT_ERROR_1d=1;
  if COMPLEMENT_ERROR_1d
    color_err='r.-';
    scale=100;
    legend_tit_err=strcat('Error * ',num2str(scale));
    GRAPH_ERR=plot(x,scale*abs(val-valex),color_err); 
  end

  Zmin=min(val)-0.1*(max(val)-min(val)); 
  Zmax=max(val)+0.1*(max(val)-min(val));
  if Zmax<=Zmin; Zmax=abs(Zmax)*1.1; Zmin=-Zmax; end; %- to prevent bug and to have 0 in [Zmin,Zmax]
  axis([xmin,xmax,Zmin,Zmax]);
  grid();


  if (OCTAVE)
    if V_EXACT==1
      legend(legend_tit_scheme,legend_tit_ex,'Location',legend_pos);
      if COMPLEMENT_ERROR_1d
        legend(legend_tit_scheme,legend_tit_ex,legend_tit_err,'Location',legend_pos);
      end
    else
      legend(legend_tit_scheme,'Location',legend_pos);
    end
  else
    if V_EXACT==1
      legend([GRAPH(1), GRAPH_EX(1)],legend_tit_scheme,legend_tit_ex,'Location',legend_pos);
      if COMPLEMENT_ERROR_1d
        legend([GRAPH(1), GRAPH_EX(1), GRAPH_ERR(1)],legend_tit_scheme,legend_tit_ex,legend_tit_err,'Location',legend_pos);
      end
    else
      legend([GRAPH(1)],legend_tit_scheme,'Location',legend_pos);
    end
  end

  xlabel(tit1); 
  %TITLE=strcat('t=',num2str(T),'; level set=', num2str(level_set));
  ti=T;
  tiFile=['  (',TITLE_FILE,')'];
  title(strcat('t=',num2str(min(T,ti)),tiFile));

  fprintf('ALL DONE (1d)\n');
  return 
end

%-----------------------------
%-------- DIM >= 2 -----------
%-----------------------------
if (dim>=2 && SET_COUPE==1)

    if VALUE_PB
      FILE=[filePrefix 'Value.dat'];
    else
      %FILE=[filePrefix 'coupe.dat'];
      FILE  =[filePrefix VFile0   '.dat'];
      FILEex=[filePrefix VexFile0 '.dat'];
      %- can be changed : coupe.dat -->  topt.dat, VF.dat (if 2d problems) ...
    end
    TITLE=FILE;
    TITLEex=FILEex;


    %----------------------------------------------------------
    %- LOADING : TEX/BIN 
    %----------------------------------------------------------
    message='fig1';

else
  fprintf('Error: this case not treated in output_view\n'); 
  return; 
end


if VALUE_PB
  fprintf(' ** WARNING : VALUE PB ==> reformating for (Value.dat), to remove "INF" values \n');
  INF=1.e5;
  %i=find(data(:,3)<INF);
  %ValueMax=max(data(i,3));
  i=find(val<INF); %- val should be only 1 column
  ValueMax=max(val(i));
  level_set=ValueMax*(1-1e-2);
  fprintf(' ** WARNING :          ==> setting level_set =  %5.2f (max(|Value|) for |Value|<INF)\n',level_set);
end


%- OK for dim=2: cdd=[d1,d2].
%- <OBSOLETE>
% if FORMAT_FULLDATA==1
%   xi=data(:,1:2);
%   val=data(:,3); 
% else
%   val=data(:,1);
% end

if 1
  %- <B> <SHOULD REMAIN>
  %------------------------------
  %- following assumes dimcoupe=2
  %------------------------------
  xmini=xmin(d1); xmaxi=xmax(d1); 
  ymini=xmin(d2); ymaxi=xmax(d2); 
  dx1=dx(d1);     dx2=dx(d2);
  n1 =nn(d1);     n2 =nn(d2);

  %- TEST - ENCOURS 2019 - ISO-MAILLAGE 2d avec Dx(j)=1:
  %TEST=1;
  if (TEST)
    dx1=1.0;
    dx2=dx1;
    xmini=0; xmaxi=xmini+n1*dx1;
    ymini=0; ymaxi=ymini+n2*dx2;
  end

  xmesh=xmini + (0:n1-1)'*dx1 + (1-MESH)*dx1/2;
  ymesh=ymini + (0:n2-1)'*dx2 + (1-MESH)*dx2/2;
  xmeshi=xmesh;   ymeshi=ymesh;
  %- constructs x s.t. x(i,1) x(i,2) contain the list of x1 and x2 mesh points of the full 2d grid.
  ox=ones(size(xmesh));
  oy=ones(size(ymesh));
  xtemp=xmesh*oy'; xtemp=xtemp(:);
  ytemp=ox*ymesh'; ytemp=ytemp(:);
  x=[xtemp, ytemp];

end


empty=0; %- security variable: empty=1 means there is no level set, then no legends are drawns (to prevent error)

%----------------------------------------
%- REACHABLE PLOT (2019)
%----------------------------------------
if (1 && (PLOT_REACHABLE==1 || PLOT_CONTOUR==1))

  %- NEW 2019
  ifig=1;
  figure(1);
  clf;
  drawperso();

  hold on;

  message='fig1';
  [val,err]=loadFile(FILE,message,n1*n2,dimcoupe);
  if (err==1); return; end;


  if (V_EXACT==1)
    %------------------------------------------------------------------
    %- 2018 : Attempt to furthermore read exact file corresponding to current data
    %-        preparation of "val_ex" before calling plot_reachable()
    %------------------------------------------------------------------
    %FILEex=strcat(filePrefix,VexFile0,int2str(i),'.dat');
    %FILEex=[filePrefix VexFile0 int2str(i) '.dat'];
    ii=-1;
    if (ii==-1) 
      str='';
    else     
      srt=int2str(ii);
    end
    FILEex=[filePrefix VexFile0 str '.dat'];

    if (1)
      fprintf('loading %s .. ',FILEex); 
    end

    [status, result]=unix(strcat('ls ./',FILEex));
    if status==2; 
      fprintf(strcat('no file ./',FILEex,'. Skip'));
      %- trying with no {it} index 
      %- This may be used for the last iteration
      FILEex=[filePrefix VexFile0 '.dat'];
      [status, result]=unix(strcat('ls ./',FILEex));
      if status==2 
        fprintf(strcat('; ..no file ./',FILEex,'. Skip; Aborting\n'));
        return
      else
        fprintf('; using %s instead.',FILEex);
        %if ~PRINTF; fprintf('\n'); end;
        %i=ITERMAX;  %- in order to break for the for loop
      end
    end

    %----------------
    %- LOADING FILEex
    %----------------
    %data=load(FILEex);
    %val_ex=data(:,3);
    message='';
    [val_ex,err]=loadFile(FILEex,message,n1*n2,dimcoupe); % -dimcoupe should be 2!-
    if err==1; return; end


    TITLE_FILEex=FILEex;


    %------------------------------------------------------------------
    %- 2018 : end 
    %------------------------------------------------------------------
  end

  plot_reachable(-1);

  %------------------
  %- auxiliary plots
  %------------------
  FILE=['complement.m'];
  fprintf('.. Reading file: "%s" ... ',FILE); 
  [status, result]=unix(strcat('ls ./',FILE));
  if status==2; 
    fprintf(strcat('no file ./',FILE,'. Skipping.\n'));
  else
    fprintf('DONE.\n');
    complement;
  end
  %------------------
  %- end of auxiliary plots
  %------------------

end
%-----------------------
%- ENDS REACHABLE PLOT (2019)
%-----------------------

if (0 && (PLOT_REACHABLE==1 || PLOT_CONTOUR==1))

  ifig=1;
  figure(1);
  clf;
  drawperso();

  hold on;

  message='fig1';
  [val,err]=loadFile(FILE,message,n1*n2,dimcoupe);
  if (err==1); return; end;

  if PLOT_REACHABLE
    i=find(val<=level_set);
    xg=x(i,1);
    yg=x(i,2);
    color_app='b.';
    if length(i)>0
      GRAPH_REACHABLE=plot(xg,yg,color_app,'MarkerSize',factorMarkerSize);
    else
      fprintf('no level set ==> no level set plot !\n');
    end

    %------------------
    %- auxiliary plots
    %------------------
    COMPLEMENT=0;
    if COMPLEMENT
      %complement1;
      %complement_xpp; %- err. estimates, case of "data_xpp.h"
      %complement2;
      complement_test;
      %axis equal; %- to have normal scalings 
    end 

    if 0
    FILE=['complement.m'];
    fprintf('.. Reading file: "%s" ... ',FILE); 
    [status, result]=unix(strcat('ls ./',FILE));
    if status==2; 
      fprintf(strcat('no file ./',FILE,'. Skipping.\n'));
    else
      fprintf('DONE.\n');
      complement;
    end
    end
    %------------------
    %- end of auxiliary plots
    %------------------
  end

  grid on;

  %--------------------------------------------
  %- THE FOLLOWING WILL PLOT THE CONTOUR
  %--------------------------------------------

  if PLOT_CONTOUR
    V=zeros(n1,n2); V(:)=val;
    %level=min(min(abs(val))); %- this better for value problem
    level=[level_set, level_set];
    [zz,GRAPH_CONTOUR]=contour(xmesh,ymesh,V',level,color_contour,'LineWidth',lw_contour);
    axis([xmini,xmaxi,ymini,ymaxi]);
    hold on;
  end

  if V_EXACT
    %----------------------------------------------------------
    %- LOADING EXACT VALUE IN "vexa"
    %----------------------------------------------------------
    %FILE=[filePrefix 'coupeex.dat']; 
    %FILE=[filePrefix 'Vex.dat']; 
    FILEex=[filePrefix VexFile0 '.dat'];
    TITLEex=FILEex;
    message='fig1';
    [vexa,err]=loadFile(FILEex,message,n1*n2,dimcoupe);
    if err==1; 
      fprintf('.. V_EXACT=%i, but found no file %s. Going on with V_EXACT=0 (figure %i)\n',V_EXACT,FILEex,ifig);
      V_EXACT=0;
    end
  end

  %%- 2018 : plot of "val" and "val_ex" 
  if V_EXACT
    hold on        
    % figure(1);  %- (remettre graphique sur la fig 1 sans detruire / recreer fenetre)
    level=[level_set, level_set];
    V=zeros(n1,n2);

    if PLOT_CONTOUR
      V(:)=val;
      %level=min(min(abs(val))); %- this better for value problem
      [zz,GRAPH_CONTOUR]=contour(xmesh,ymesh,V',level,color_contour,'LineWidth',lw_contour);  % 2019 (3->)
      empty=(length(zz)==0);

      %FILE=[filePrefix 'coupeex.dat']; 
      %fprintf(strcat('loading file for fig1: (',FILE,', for exact contour) ...'));
      %vexa=load(FILE);
      %fprintf('DONE\n');

      V(:)=vexa;
      %level=min(min(abs(vexa))); %- this better for value problem
      [zz,GRAPH_CONTOUR_EX]=contour(xmesh,ymesh,V',level,color_contour_ex,'LineWidth',lw_contour); % 2019 (2->)
      xlabel(tit1);
      ylabel(tit2);
      emptyloc=(length(zz)==0);
      empty=min(empty,emptyloc);

      axis([xmini,xmaxi,ymini,ymaxi]);
    end

    if exist('Front_ex');
      t=T;
      X=Front_ex(t);
      xe=X(:,1); ye=X(:,2);
      [zz,GRAPH_CONTOUR_EX]=plot(xe,ye,'k--','LineWidth',2);
    end

  end


  if AXIS_EQUAL
    axis equal;
    axis([xmini,xmaxi,ymini,ymaxi]);
  end
  grid on;

  xlabel(tit1); ylabel(tit2);
  TITLE=strcat('t=',num2str(T),'; level set=', num2str(level_set),' (',FILE,')');
  title(TITLE);

end %- end of if (PLOT_REACHABLE==1 | PLOT_CONTOUR==1)
%-----------------------
%- ENDS REACHABLE PLOT
%-----------------------

%----------------------------------------
%- VF SURFACE PLOT (OR TOPT SURFACE PLOT)
%----------------------------------------
if PLOT_3d==1

  figure(2);
  clf;
  drawperso();

  %FILE='VF.dat';
  %VFile0='coupe.dat';
  FILE=[filePrefix VFile0 '.dat'];

  %if dim>=3; FILE='coupe.dat'; end;
  if (dim>=3) || SET_COUPE; FILE=[filePrefix 'coupe.dat']; end;
  
  if VALUE_PB; 
    if dim==3
      %FILE=[filePrefix 'Value.dat']; 
    else
      fprintf('ENCOURS... dim=%i(<>3) and VALUE_PB=1 not permitted for PLOT_3d=1\n',dim); 
      return;
    end
  end

  if PLOT_TOPT
    if dim>2; 
      fprintf('PLOT_TOPT=1 and dim>2 not fixed - abort this part for fig 2\n');
      return;
    end
    %toptFile0='toptcoupe.dat';
    %toptFile0='coupetopt.dat'; %- TODO 2018 : adapter aux noms de stdafx.h
    FILE=[filePrefix toptFile0 '.dat'];
    [status,result]=unix(strcat('ls ./',FILE));
    if status==2; 
      fprintf('PLOT_TOPT=1, but no file %s . Abort this plot.\n',FILE);
      return; 
    end
  end

  %----------------------------------------------------------
  %- LOADING TEX or BIN
  %----------------------------------------------------------
  message='fig2';
  [val,err]=loadFile(FILE,message,n1*n2,dimcoupe);
  if err==1; return; end



  if VALUE_PB; val=min(val,ValueMax); end;

  if PLOT_TOPT==0
    %xmesh=x(1:n1,1);
    %ymesh=x(n1*(1:n2),2);
    V=zeros(n1,n2);
    V(:)=val;
    colormap hsv;
    surf(xmesh,ymesh,V');

    title(TITLE);
    xlabel(tit1);
    ylabel(tit2);
    view(ALPHA,BETA);
    if UPDATE_AXIS
      vmin=min(V(:)); vmax=max(V(:));
      if (vmax-vmin)/(0.001+max(abs(V(:))))<=0.02e-1; %-special case when V=const !
         Zmin=vmin;
         Zmax=vmin+1;
         fprintf('Careful : SPECIAL CASE Zmin close to Zmax ==> setting Zmax=Zmin+1.\n');
         %axis([xmin(1),xmax(1),xmin(2),xmax(2),Zmin,Zmax]);
      else
        vmin=min(V(:)); vmax=max(V(:));
        Zmin=vmin-1.2*(vmax-vmin);
        Zmax=vmax+1.2*(vmax-vmin);
      end
      %Zmin=-0.5; Zmax=0.5;
    end
    axis([xmini,xmaxi,ymini,ymaxi,Zmin,Zmax]);

  end 


  if PLOT_TOPT==1
    %- typical for topt plot when topt is a minimal time 
    %- Will bound max(topt) and set to 0 all values s.t T=TMAX (that is, when T=+INF)
    TMAX=1.5*T;
    val=min(TMAX,val);	      

    i_val0=find(val==0); 

    i=find(val==TMAX); val(i)=0.0;

    %level=linspace(0,T,30); %- obsolete
    %xmesh=x(1:n1,1);
    %ymesh=x(n1*(1:n2),2);
    V=zeros(n1,n2);
    V(:)=val;
    %contour(xmesh,ymesh,V',level);
    %V=min(V,3);

    xmeshi=xmesh;   ymeshi=ymesh;   Vi=V;
    xmini=xmin(d1); xmaxi=xmax(d1); 
    ymini=xmin(d2); ymaxi=xmax(d2); 

    %- modif: plot in a subdomain:
    SUBDOMAIN=0;
    if SUBDOMAIN
      xmini=-1.0;  xmaxi=1.0;
      ymini=-1.0;  ymaxi=1.0;
      i=find(xmesh>=xmini  & xmesh<=xmaxi); xmeshi=xmesh(i);
      j=find(ymesh>=ymini  & ymesh<=ymaxi); ymeshi=ymesh(j);
      Vi=V(i,j);
    end
    colormap hsv;
    surf(xmeshi,ymeshi,Vi'); hold on;
    %- end of modif

    if length(i_val0)>0
      xx=x(i_val0,1); yy=x(i_val0,2); 
      plot3(xx,yy,zeros(size(xx)),'g.');
    end

    title(TITLE);
    %if OCTAVE==0; axis equal; end;
    zmin=0.0; zmax=TMAX;
    axis([xmini,xmaxi,ymini,ymaxi,zmin,zmax]);
    view(ALPHA,BETA);
    clear Vi
  end

end %- end if PLOT_3d==1



if (PLOT_3d_ex==1) %- plot of Vex.dat

  if V_EXACT

    figure(3);
    clf;
    %if OCTAVE; set(gcf,'visible','off'); set(gcf,'visible','on'); end
    drawperso();

    %FILE=[filePrefix 'Vex.dat'];
    FILE=[filePrefix VexFile0 '.dat'];
    %if dim>=3; FILE=[filePrefix 'coupeex.dat']; end;
    TITLE=FILE;

    message='fig3';
    [val,err]=loadFile(FILE,message,n1*n2,dimcoupe);
    if (err==1); return; end


    V=zeros(n1,n2);
    V(:)=val;
    colormap hsv;
    surf(xmesh,ymesh,V');

    title(TITLE);
    xlabel(tit1)
    ylabel(tit2)
    view(ALPHA,BETA);
    if UPDATE_AXIS 
      vmin=min(V(:)); vmax=max(V(:));
      if (vmax-vmin)/(0.001+max(abs(V(:))))<=0.02e-1; %-special case when V=const !
        Zmin=vmin;
        Zmax=vmin+1;
        fprintf('Careful : SPECIAL CASE Zmin close to Zmax ==> setting Zmax=Zmin+1.\n');
      else
        Zmin=vmin-1.2*(vmax-vmin);
        Zmax=vmax+1.2*(vmax-vmin);
      end
    end
    axis([xmini,xmaxi,ymini,ymaxi,Zmin,Zmax]);

  else %- V_EXACT==0
    %fprintf('Vex doesnot seem not exist - no Vex surface plot\n');
  end
else
  %- removing the 3d_ex plot if it exists 
  if (ishandle(3)); delete(3); end;
end


%----------------------------------
%- ADDING TARGET PLOT
%----------------------------------
while TARGET_PLOT && PLOT_ISOVIEW==0

  if (VALUE_PB==0 && dim~=2) || (VALUE_PB==1 && dim~=3)
    fprintf(' ** WARNING : TARGET_PLOT=1 may not work - only programmed for dim=2 or dimcoupe==2 (here dim=%i,dimcoupe=%i)\n',dim,dimcoupe);
    fprintf(' ** WARNING :                            - or dim=3 if VALUE_PB=1;\n'); 
    if ~(dimcoupe==2)
      fprintf('Skip targ. plot;\n');
      abort; 
    end;
  end

  %- Careful this may not work if VF0.dat not present or if FORMAT_FULLDATA=0
  for ifig=[1] 
  %for ifig=[1 2] 
  %if PLOT_REACHABLE==1;

  figure(ifig);
  %if OCTAVE; set(gcf,'visible','off'); set(gcf,'visible','on'); end
  drawperso();
  hold on;

  FILE=[filePrefix VFile0 '0' '.dat'];
  if VALUE_PB; 
    FILE=[filePrefix 'Value0.dat'];
  end
  level_set0=0.0;

  % if ifig==1; fprintf('TARGET_PLOT: .. '); end;
  % [status, result]=unix(strcat('ls ./',FILE));
  % if status==2; 
  %   fprintf(strcat('no file ./',FILE,'. Skipping.\n'));
  %   break;
  % else
  %   if ifig==1; 
  %     fprintf(strcat('loading VF0  (',FILE,') for green plotting of {x,VF0<=0}...'));
  %   end
  % end

  % data=load(FILE);           %- 2d
  % if size(data,2)==1
  %   fprintf('no standard format (may be FORMAT_FULL_DATA=0?). Aborting Target Plot.\n'); 
  %   return; 
  % else
  %   if ifig==1; fprintf('DONE\n'); end;
  % end
  % %- this works only for 2d problems 
  % val=data(:,3); 

  message='';
  if ifig==1; 
    message='fig1-target-plot-from-VF0';
  end
  [val,err]=loadFile(FILE,message,n1*n2,dimcoupe);
  if (err==1); 
    fprintf('ERR: file %s not present. Skip TARGET_PLOT.\n',FILE)
    break %- leaving "while TARGET_PLOT"
  end

  i=find(val<=level_set0);
  %level_temp=min(min(val)); i=find(val<=level_temp); %- this better for value problem
  xg=x(i,1);
  yg=x(i,2);
  
  if OCTAVE; color_target='g*'; color_contour='g--';
  else     ; color_target='g.'; color_contour='g--';
  end 
  %color_target='g.'; color_contour='g--';

  %- green region {VF0<=0}
  %factor_points=2;
  %GRAPH_VF0=plot(xg,yg,color_target,'MarkerSize',5*factor_points);
  %[zz,GRAPH_VF0]=plot(xg,yg,color_target,'MarkerSize',5*factor_points);
  legend_tit_VF0='Target set';
  plot(xg,yg,color_target,'MarkerSize',2*factorMarkerSize);

  %xmesh=x(1:n1,1);
  %ymesh=x(n1*(1:n2),2);
  V=zeros(n1,n2);
  V(:)=val;
  level=[level_set, level_set];
  color_contour_target='g-';
  [zz,GRAPH_CONTOUR_VF0]=contour(xmesh,ymesh,V',level,color_contour_target,'LineWidth',lw_contour); % 2019 (3->)

  %- adapt legend
  if (~empty)
  if V_EXACT==1
   %legend([GRAPH_REACH(1), GRAPH_CONTOUR_EX(1), GRAPH_VF0(1)],legend_tit_scheme,legend_tit_ex,legend_tit_VF0,'Location',legend_pos);
   legend([GRAPH_CONTOUR(1), GRAPH_CONTOUR_EX(1), GRAPH_CONTOUR_VF0(1)],legend_tit_scheme,legend_tit_ex,legend_tit_VF0,'Location',legend_pos);
  else
   %legend([GRAPH_REACHABLE, GRAPH_VF0],legend_tit_scheme,legend_tit_VF0,'Location',legend_pos);
   legend([GRAPH_CONTOUR(1), GRAPH_CONTOUR_VF0(1)],legend_tit_scheme,legend_tit_VF0,'Location',legend_pos);
   %legend([GRAPH_CONTOUR],legend_tit_scheme,'Location',legend_pos);
   %legend(GRAPH_CONTOUR_VF0(1),legend_tit_VF0,'Location',legend_pos);
  end
  end
  % ax=gca
  % ax.FontSize=14


  %- adding green contour {VF0==0} (tested for cdd=[1,2])
  if 0
    V=zeros(n1,n2); V(:)=val;
    level=[level_set, level_set];
    [zz,GRAPH_CONTOUR]=contour(xmesh,ymesh,V',level,color_contour,'LineWidth',3);
  end
  hold on;

  end %- end of [for ifig=[1]]

  break

end


%------------------------
%- DEC 2012 : trajectory
%------------------------
if TRAJECTORY;

  legend off;

  fprintf('TRAJECTORY PLOT (first 2 components on 2d coupe):');
  fprintf('\n');

  SPEC=0; %- set SPEC=1 for special graphs (2d) / zoom change

  %- reading file "successTrajectories.dat"

  FILE=[filePrefix 'successTrajectories.dat'];
  fprintf('.. Reading file: "%s" ... ',FILE); 
  [status, result]=unix(strcat('ls ./',FILE));
  if status==2; 
    fprintf(strcat('no file ./',FILE,'. Skipping.\n'));
    return; 
  end
  fprintf('DONE\n');
  success=load(FILE,'-ascii'); 
  success=success(:,end);
  %- end of reading file "successTrajectories.dat"

  fprintf('.. ploting traj : '); 

  for ifig=[1] 
  %for ifig=[1 2]  %- BUG 2015
  figure(ifig);
  %drawperso();
  hold on;

  GRAPH_TRAJ=zeros(nbTraj,1);


  for i=1:nbTraj
    
  nametraj=strcat(filePrefix,'traj-',num2str(i),'.dat');
  %fprintf('reading file nametraj'' \n');
  [status,result]=unix(strcat('ls ./',nametraj));
  if status==2; 
    fprintf('TRAJECTORY=1, but no file %s..\n',nametraj);
    return; 
  end
  traj=load(nametraj);
  xx=traj(:,1);
  yy=traj(:,2);


  %linemarkersize=5*2;
  %markersize    =5*4;
  factorLineMarkerSize=5;
  factorMarkerSizeTraj=5;
  linemarkersize=factorLineMarkerSize*1.5;
  markersize    =factorMarkerSizeTraj*4;

  lw_success=2;		%- linewidth for successful traj (fig 1)
  lw_unsuccess=2;	%- linewidth for unsuccessful traj (fig 1)

  couli='k.'; 		%- color of trajectory 
  coul_init='r.';	%- color of initial point 
  coul_fina='b.';	%- color of terminal point
  lw=lw_success;
  if success(i)==0;     
    couli='k-.';        %- color of trajectory     (unsuccessful trajectories)
    coul_init='r.';	%- color of initial point  (unsuccessful trajectories)
    coul_fina='b.';	%- color of terminal point (unsuccessful trajectories)
    lw=lw_unsuccess;
  end

  %if (SPEC && i==2); couli='-ok'; end
  %[GRAPH_TRAJ(i)]=plot3(xx,yy,zz,couli,'LineWidth',2);
  %plot3(xx(1),yy(1),zz(1)      ,'g.','MarkerSize',20);
  %plot3(xx(end),yy(end),zz(end),'r.','MarkerSize',20);
  [GRAPH_TRAJ(i)]=plot(xx,yy,couli,'LineWidth',lw,'MarkerSize',linemarkersize);
  plot(xx(end),yy(end),coul_fina,'MarkerSize',markersize);
  plot(xx(1),yy(1)    ,coul_init,'MarkerSize',markersize);
  %- SPECIAL TWO PLAYERS, dim=2+2=4
  %xx2=traj(:,3);
  %yy2=traj(:,4);
  %zz2=0*ones(size(xx2)); 
  %GRAPH_TRAJ2=plot3(xx2,yy2,zz2,'b.','LineWidth',2);
  %plot3(xx2(1),yy2(1),zz2(1)      ,'g.','MarkerSize',20);
  %plot3(xx2(end),yy2(end),zz2(end),'r.','MarkerSize',20);

  %if VALUE_PB && ifig==1
    %zz=traj(:,3); 
    %zz=zz-zz(end);
    %plot3(xx,yy,zz,'k.','LineWidth',2);% 'MarkerSize',15);
  %end	    
  grid on;
  hold on;

  end % end of nbTraj


  if SPEC
    legend_pos='northeast';
    if OCTAVE
      legend(['Traj2'],legend_pos);
    else
      legend([GRAPH_TRAJ(1)],['Traj2'],'Location',legend_pos);
    end
    %legend([GRAPH_TRAJ(1), GRAPH_TRAJ(2)],['Traj1'; 'Traj2'],1);
    %legend([GRAPH_TRAJ(2), GRAPH_TRAJ(1)],['Traj1'; 'Traj2'],'Location',1);
    %legend('Traj1', 'Traj2');%GRAPH_TRAJ(1),GRAPH_TRAJ(2),1);%,['Traj1','Traj2'],1);
    axis([-2, 0.5, -0.5 , 2.0]);
  end

  end % end ifig

  fprintf('DONE\n'); 

  if SUPPLEMENT_TRAJECTORY
    %- SUPPLEMENTARY GRAPH WITH TRAJECTORIES(time) and CONTROL(time) 
    PRINT_TRAJ=1;    	%- to show trajectories
    PRINT_CONTROL=1; 	%- to show corresponding control values
    trajmax=6;	 	%- maximum number of trajectories plotted
    nxmax=6;  		%- max number of traj. components shown
    ncmax=3;  		%- max number of control components shown
    ifig0=10;  		%- number of figure for for trajectories
    ifig0_control=20;	%- number of figure for corresponding controls
    lw_success=2;	%- linewidth for successful traj
    lw_unsuccess=1;	%- linewidth for unsuccessful traj
    %linestyle   ={'-','-','-','-','-','-'}; 	%- line style : with lines - also one can use one of '.','.-','-.','-', 'o'
    linestyle   ={'.','.','.','.','.','.'}; 	%- line style : with dots - advised style for periodic b.c.
    linestyle_unsuccessful={'-.'};		%- unsuccessful trajectory line style
    charstyle   ={'k','b', 'r', 'm', 'g', 'c'}; %- color line style: black, blue, red, magenta, green, cyan


    nametitle   ={'x1','x2','x3','x4','x5','x6'};
    controltitle={'control u1','control u2','control u3','control u4','control u5','control u6'};

    %fprintf('TRAJECTORY & CONTROL PLOTS (figs %i and %i) ...',ifig0,ifig0_control);
    fprintf('..Trajectory and control plots [Unsuccessful trajectories : with dotted lines]');

    if (dim>6)
      fprintf('Warning, only showing dimemsion <=6 of trajectories !\n');
    end
    nxmax=min(6,min(nxmax,dim)); %- max 6 trajectory plots here



    modeplot=1; %- set to 1 for 1 figure with subplots;  
    %modeplot=2; %- obsolete: for 1 figure per coordinate (for traj.only, not set for controls).

    %- clear previous drawings / set figures
    %figure('Name','Simulation Plot Window','NumberTitle','off')
    if modeplot==1 & ~OCTAVE
      if PRINT_TRAJ==1; 
        fig=figure(ifig0); clf;
        fig.Name='Trajectory components';
      end
      if PRINT_CONTROL==1; 
        fig=figure(ifig0_control); clf;
        fig.Name='Control components';
      end

    else %- obsolete
      % for j=1:nxmax;
      %   figure(ifig0+j-1); clf; 
      % end
      % figure(ifig0_control); clf;
    end


    ntraj=min(nbTraj,trajmax); %- max traj shown

    for i=1:ntraj

    %- loading traj
    nametraj=strcat(filePrefix,'traj-',num2str(i),'.dat');
    %fprintf('reading file nametraj'' \n');
    [status,result]=unix(strcat('ls ./',nametraj));
    if status==2; 
      fprintf('TRAJECTORY=1, but no file %s..\n',nametraj);
      return; 
    end
    traj=load(nametraj);


    %- traj plot.
    if PRINT_TRAJ

    figure(ifig0); 
    tt=traj(:,end); %- time
 

    cool=char(charstyle(i));
    ls=char(linestyle(i));
    cool=strcat(cool,ls);
    lw=lw_success;
    if success(i)==0;  %- special colors or style for unsuccessful trajectories:
      %cool='c';
      %linestyle(i)={'--'};
      cool=char(charstyle(i));
      ls=char(linestyle_unsuccessful); %- unsuccessful trajectory line style
      cool=strcat(cool,ls);
      lw=lw_unsuccess;
    end

    for j=1:nxmax;
      xx=traj(:,j); 
      %subplot(1,nxmax,j); plot(tt,xx,coul(i)); 
      %subplot(1,nxmax,j); plot(tt,xx,'LineStyle',coul(i),'LineWidth',2); 
      if modeplot==1
        figure(ifig0); 
        subplot(1,nxmax,j);
        %plot(tt,xx,cool,'LineStyle',char(linestyle(i)),'LineWidth',lw); hold on;
        %if length(ls)==0;
        %  plot(tt,xx,cool,'LineWidth',lw); hold on;
        %else
        %  plot(tt,xx,cool,'LineStyle',ls,'LineWidth',lw); hold on;
        %end
        plot(tt,xx,cool,'LineWidth',lw); hold on;
      else %- obsolete
        figure(ifig0+j-1); 
        cool=char(charstyle(i));
        plot(tt,xx,cool,'LineStyle',char(linestyle(i)),'LineWidth',lw); hold on;
      end
      %title(strcat(nametitle(j),num2str(j)));
      title(nametitle(j));
    end
    end %- end if PRINT_TRAJ

    %- control plot: traj contains number of colums dim (trajectory) + nc (controls) + 1 (time)
    if PRINT_CONTROL
      figure(ifig0_control); 
      %clf;
      ntot=size(traj,2); nc=ntot-1-dim; 
      tt=traj(:,end); %- time
      ncmax=min(nc,ncmax);
      for j=1:ncmax
        xx=traj(:,dim+j); 
        %subplot(1,ncmax,j); plot(tt,xx); title(strcat('control u',num2str(j)));
        subplot(1,ncmax,j); 
        %if length(ls)==0
        %  plot(tt,xx,cool,'LineWidth',lw); hold on;
        %else
        %  plot(tt,xx,cool,'LineStyle',ls,'LineWidth',lw); hold on;
        %end
        plot(tt,xx,cool,'LineWidth',lw); hold on;
        %cool='b.';
        %plot(tt,xnc,cool,'LineWidth',2); hold on;
        title(controltitle(j));
      end
    end %- end of SUPPLEMENT_TRAJECOTRY

    end %- end for i=1:nbTraj

  end

  %fprintf('\n');
end


if AXIS_EQUAL
  %for ifig=[2 1]
  for ifig=[1]
    figure(ifig);
    axis equal;
  end
end

if PLOT_ISOVIEW
  if ~(dim==3 || dimcoupe==3)
    fprintf('3d-ISOVIEW : dim(%i) <>3 or dimcoupe(%i)<>3 ==> Skip 3d isoview plot.\n',dim,dimcoupe);
  elseif OCTAVE
    fprintf('3d-ISOVIEW : 3d isoview plot not programmed with OCTAVE. Skip.\n');
  else 
    if ~(dimcoupe==3) %- then dim==3: force the use of the VF.dat file
      fprintf(' ** WARNING : Forcing the use of ''VF.dat'' file\n');
      VFile0='VF';
    end
    figure(4);
    drawperso();
    %FILE=[filePrefix VFile0 '.dat'];
    output_view3d(); 
    %input('exiting movie3d mode [waiting for "Enter"]');
  end
end


%-----------------------------------------
%- 2018: ADDING LEGEND ON FIG 1 
%-----------------------------------------
if (PLOT_REACHABLE==1 || PLOT_CONTOUR==1)
  %----------
  %- LEGEND:
  %----------
  figure(1);

  if PLOT_CONTOUR && ~empty
    if (OCTAVE)
      if exist('GRAPH_CONTOUR_EX');
        legend(legend_tit_scheme,legend_tit_ex,'Location',legend_pos);
      else
        legend(legend_tit_scheme,'Location',legend_pos);
      end
    else
      if exist('GRAPH_CONTOUR_EX(1)');
        if (length(GRAPH_CONTOUR_EX(1))>0) 
          legend([GRAPH_CONTOUR(1), GRAPH_CONTOUR_EX(1)],legend_tit_scheme,legend_tit_ex,'Location',legend_pos);
	end
      else
        legend([GRAPH_CONTOUR(1)],legend_tit_scheme,'Location',legend_pos);
        %legend(legend_tit_scheme,'Location',legend_pos);
	%return
      end
    end
  end

end

%------------------
%- COMPLEMENTS
%------------------
if 0
  figure(2);
  view(0,0);
end

if 0
  FILE=['complement.m'];
  fprintf('.. Reading file: "%s" ... ',FILE); 
  [status, result]=unix(strcat('ls ./',FILE));
  if status==2; 
    fprintf(strcat('no file ./',FILE,'. Skipping.\n'));
  else
    fprintf('DONE.\n');
    complement;
  end
end
%------------------
%- FIN COMPLEMENTS
%------------------

fprintf('ALL DONE\n');
return

end %- end of output_view()

%------------------------------------------
%- complement axes for example data_xppx.h
%------------------------------------------
%if (PLOT_CONTOUR || PLOT_REACHABLE)
%  figure(1);
%  plot([xmin(1);xmax(1)],[0;0],'k');
%  plot([0;0],[xmin(2);xmax(2)],'k');
%end

function drawperso()
  global OCTAVE
  if OCTAVE; set(gcf,'visible','off'); set(gcf,'visible','on'); end;
end

function [val,err]=loadFile(FILE,message,n,dim)
  %-------------------------------------------------------------------
  %- load tex or bin file to get value (val = 1st or 3rd column)
  %-  input : FILE, message, n = number of lines (should be n1 for 1d, or n1*n2 for 2d
  %-  input : dim : dimension (if TEXT file, will go to column dim+1)
  %-  output: val
  %-------------------------------------------------------------------
  global BINARY
  global FORMAT_FULLDATA

  PRINTloc=0;

  if length(message)==0;
    PRINTloc=0;
  else
    PRINTloc=1;
  end

  err=0;

  [status, result]=unix(strcat('ls ./',FILE));

  if status==2; 

    fprintf(strcat('.. no file ./',FILE,'. Skip; Aborting load.\n'));
    %break;
    val=0;
    err=1;
    return

  else

    if PRINTloc; fprintf('loading %12s (%s): .. ',FILE,message); end
    %data=load(FILE);

    if BINARY
      %----------------------------------
      %- check that FORMAT_FULLDATA=0 !?
      %----------------------------------
      if PRINTloc; fprintf('(BINARY)..'); end
      %FILE=[filePrefix 'VF.dat'];
      dataF=fopen(FILE,'rb');
      data=fread(dataF, [n,1] ,'double');
      FORMAT_FULLDATA=0; %- forcing - for the moment
    else
      data=load(FILE);	%- assumes file is 2d : list of i,j,val
    end
  
    if FORMAT_FULLDATA==1
      val=data(:,dim+1); 	%- assumes list of (i,val) or (i,j,val) or (i,j,k,val)
    else
      val=data(:,1);		%- assumes list of (val)
    end
    if PRINTloc; fprintf('DONE\n'); end

  end

end

