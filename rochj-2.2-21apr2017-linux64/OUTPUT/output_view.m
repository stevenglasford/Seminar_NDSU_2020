fprintf('ROCHJ v2.2 - matlab/octave postprocessing ');
clear
%close all

global level_set XPAUSE
global Zmin Zmax  %- to be adjusted depending of [min,max] values of data

%- OCTAVE test          %- Default is OCTAVE=0. Special commands for octave : (graphs, legends...)
OCTAVE = (exist ('OCTAVE_VERSION', 'builtin') > 0); if OCTAVE fprintf('(OCTAVE DETECTED)\n'); else fprintf('(MATLAB DETECTED)\n'); end;
%if OCTAVE; close all; end;

%-------------
%- fig 1:
%-------------
PLOT_CONTOUR   =1;	%- plot contours
PLOT_REACHABLE =1;	%- plot reachable region (points where scheme value <= level_set)
level_set      =0.00;
TARGET_PLOT    =1;	%- target plot = green dots where initial value VF0.dat is <= 0. For instance use in the data file:
			% const int     SAVE_VF_ALL        = 1; 
			% const int     SAVE_VF_ALL_STEP   = 10000;  //- or smaller
			%- in order to get a file VF0.dat


%-------------
%- fig 2:
%-------------
PLOT_3d    =1;		%- Set figure 2 for plots: plot of VF or topt 
PLOT_TMIN  =0; 		%- if PLOT_TMIN=1, will use topt.dat (optimal time function) instead of VF.dat 


%-------------
%- fig 3:
%-------------
PLOT_3d_ex =0;		%- plot 3d of VEX (if furthermore compute_vex=1..)

%-------------
%- for figs 1 and 2 (plus, possibly fig 4 and 10) :
%-------------
AXIS_EQUAL=0;		%- will make axis equal for fig 1 (or fig 1 and 2) 
%ALPHA=65; BETA=35;	%- angles for the 3d plot (figs 2 and 3)
ALPHA=25; BETA=50;	%- angles for the 3d plot (figs 2 and 3)
%ALPHA=0; BETA=0;	%- angles for the 3d plot (figs 2 and 3)
			%- (plotting tmin is only ok for dim=2)

%-------------
%- for figs 4 and 10 
%-------------
%TRAJECTORY=1;		%- plots of traj-xx.dat (first 2 components)
SUPPLEMENT_TRAJECTORY=1;%- each state (and possibly each control) is plotted vs time  (see also PRINT_TRAJ and PRINT_CONTROL parameters inside)

%-------------
%- Special parameters : Movie, value problem, etc.
%-------------
MOVIE=0; 		%- Default is 0. 1d/2d movie (can be used in particular if pb is 1d (if SAVE_VFALL=1), or if pb is 2d
XPAUSE=0; 		%- for graphic pause (no need to type "return" after each plot).
%VALUE_PB=0; 		%- Default is 0. Indicates plot from VF.dat (if VALUE_PB=0) or from Value.dat (if VALUE_PB=1);


%Front_ex=@Front_ex_xpp;

%-----------------------------------------
%- Notes relative to some of the examples:
%-----------------------------------------
%- data_basicmodel.h:	set PLOT_3d=1; PLOT_TMIN=0 (or 1); TRAJECTOIRE=1; 
%- data_FD_zermelo.h:	set PLOT_3d=1; PLOT_TMIN=0 or 1;  AXIS_EQUAL=0 or 1;
%- data_FD_stat.h:	set PLOT_3d=1; PLOT_TMIN=0;  AXIS_EQUAL=0;
%- data_FD_dubinscar.h: 
%- data_FD_value.h:	set VALUE_PB=1; PLOT_TMIN=0 (or 1); TRAJECTORY=1;
%-     assumes also that following parameters are used at the begining of main.cpp file:
%-     VALUE_PB = 1; 
%-     SAVE_VALUE_FINAL = 1; 

%----------------------------------------------------
%- loading data.dat and initialisation of parameters
%----------------------------------------------------
global nn VALUE_PB


filePrefix='';
[zz1,zz2]=unix('ls filePrefix.dat');
if zz1==0   % there is a filePrefix.dat file and we need to get the prefix name
  FILE=fopen('filePrefix.dat'); 
  filePrefix=fscanf(FILE,'%s');
  fclose(FILE);
end
dataFile=[filePrefix, 'data.dat'];

param=load(dataFile,'-ascii'); %- problem & mesh parameters

i=1;
dim =param(i); i=i+1;
MESH=param(i); i=i+1;
dx=zeros(1,dim); xmin=zeros(1,dim); xmax=zeros(1,dim); nn=zeros(1,dim);
for d=1:dim,
  dx(d)=param(i); i=i+1;
end
for d=1:dim,
  xmin(d)=param(i); i=i+1;
  xmax(d)=param(i); i=i+1;
  fprintf('d=%i, xmin(d)=%6.3f, xmax(d)=%6.3f\n',d,xmin(d),xmax(d));
end
for d=1:dim,
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
  SET_COUPE=1; %- 2017: should be always SET_COUPE=1 
end


if (1)
  fprintf('T=%6.3f\n', T);
  fprintf('V_EXACT=%2i \n',V_EXACT);
  fprintf('SET_COUPE=%2i [default should be 1: plots using "coupexx.dat files (excepted for dim=1)]\n',SET_COUPE);
  for d=1:dim; fprintf('nn(%i)=%3i; ',d,nn(d)); end; 
  fprintf('dim=%i; dimcoupe=%i; ',dim,dimcoupe);
  fprintf('\n');
end 



%- MAR 2014
global cdd

if SET_COUPE == 1
  coupe_dims=zeros(dim,1); coupe_vals=zeros(dim,1);
  for d=1:dim,
    coupe_dims(d) = param(i); i=i+1;
    fprintf('coupe_dims(%i)=%2i; ',d,coupe_dims(d));
  end
  fprintf('\n');

  cdd = find(coupe_dims==1); %- dimensions kept fater cutting.
  %cv = find(coupe_dims==0); %- dimensions cutted
  for d=1:dim 
    coupe_vals(d) = param(i); i=i+1; %Nb : only values for which coupe_dims(d)=1 are meaningful
  end
else
  %- default is 2d / MAR 2014
  cdd=[1 2];
end
if dim>1
  fprintf('Directions of the cuts: cdd=[%i %i]\n',cdd(1),cdd(2));
  d1=cdd(1);     d2=cdd(2);
  n1=nn(cdd(1)); n2=nn(cdd(2)); 
  %fprintf('(n1=%2i;  n2=%2i;)\n',n1,n2); 
else
  d1=1; d2=0;
  n1=1; n2=0;
end
if length(cdd)>2
  fprintf('.. the cut is more that 2d: will retain only first two directions for plots\n');
  cdd=cdd(1:2);
end   
if length(cdd)~= dimcoupe
  fprintf('.. Error: length(cdd)<>dimcoupe. Abort.\n'); 
  return; 
end


save_vfall=param(i); i=i+1;
nbSaves=param(i); i=i+1; 

ITERMAX=nbSaves;
if save_vfall==0
  ITERMAX=0;
end
fprintf('ITERMAX(=nbSaves)=%2i\n', ITERMAX);

nbTraj=param(i); i=i+1;
fprintf('nbTraj (nombre de trajectoires)=%2i\n',nbTraj);
TRAJECTORY=(nbTraj>=1);

VALUE_PB=param(i); i=i+1;
%VALUE_PB=1; fprintf('Warning: FORCING VALUE_PB=1 HERE FOR PLOTS'); %- ICI
if(VALUE_PB==1)
  fprintf('VALUE_PB=%2i : special settings for this case',VALUE_PB);
end

%-------------------
%- COMPLEMENT 2016
%-------------------
dtFile=[filePrefix, 'Dt.dat'];
Dt=load(dtFile,'-ascii'); %- problem & mesh parameters

% DT step for movie mode
global DT
DT=0;
if save_vfall>0
  DT=T/nbSaves;
end

%-------------------------
%- begining graphics
%-------------------------

%-------------------------
%- LEGEND PARAMETERS:
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
  lw_contour=1; 		%- linewidth for contours


else 
  fprintf('This dimcoupe=%i not programmed',dimcoupe); return; 
end



%-----------
%- MOVIE 1d 
%-----------
if (dim==1 || dimcoupe==1) && MOVIE
  Zmin=-1.5; Zmax=3.5;
  output_view_movie1d;
  input('exiting movie1d mode [waiting for "Enter"]');
  %return;
end

%-----------
%- MOVIE 2d 
%-----------
if (dim==2|| dimcoupe==2) && MOVIE

  figure(1);
  Zmin=-0.2; Zmax=0.2;
  output_view_movie2d;
  input('exiting movie2d mode [waiting for "Enter"]');
  %return;
end


if (dim==1 || dimcoupe==1)
%if (dim==1  |  (dim==2 & SET_COUPE==1))
 %-------------------------------------
 %- Special part for case dim=1
 %-------------------------------------
 %- fig 1:
 if dim==1; fprintf('dim=1 !\n'); end;
 if dim>=2; fprintf('CAREFUL: COUPE is 1D !\n'); end;

 if dim==1;
   FILE  =[filePrefix 'VF.dat'];
   FILEEX=[filePrefix 'VEX.dat'];
 else
   FILE  =[filePrefix 'coupe.dat'];
   FILEEX=[filePrefix 'coupeex.dat'];
   dx=dx(cdd);
   xmin=xmin(cdd);
 end
 fprintf('loading file %s ...\n',FILE);
 data=load(FILE);

 %----------------------------------------------------------
 %- MAR 2014: adapted for case dim>=2 and dimcoupe=1
 %----------------------------------------------------------
 if MESH==1;
   nn=size(data,1)-1;
   x=xmin + (0:nn)'*dx;
 end
 if MESH==0;
   nn=size(data,1);
   x=xmin + ((1:nn)-0.5)'*dx;
 end
 %val=data(:,1);	%- here if FORMAT_FULLDATA=0
 val=data(:,2);		%- here if FORMAT_FULLDATA=1

 figure(1);
 clf;

 color_app='';
 if 0&&OCTAVE; color_app='b*-'; color_exa='k*-';
 else;      color_app='b.-'; color_exa='k.-';
 end
 GRAPH=plot(x,val,color_app); hold on;
 %plot(x,zeros(size(x)),'k.');

 if  V_EXACT==1
   fprintf('loading file %s ...\n',FILEEX);
   data=load(FILEEX);
   valex=data(:,2);
   GRAPH_EX=plot(x,valex,color_exa); 
 end

 grid();


 if (OCTAVE)
   if V_EXACT==1
     legend(legend_tit_scheme,legend_tit_ex,'Location',legend_pos);
   else
     legend(legend_tit_scheme,'Location',legend_pos);
   end
 else
   if V_EXACT==1
     legend([GRAPH(1), GRAPH_EX(1)],legend_tit_scheme,legend_tit_ex,'Location',legend_pos);
   else
     legend([GRAPH(1)],legend_tit_scheme,'Location',legend_pos);
   end
 end

 fprintf('ALL DONE (1d)\n');
 return 
end


if (dim>=2 && SET_COUPE==1)
    if VALUE_PB; 
      FILE=[filePrefix 'Value.dat'];
    else
      FILE=[filePrefix 'coupe.dat'];
      %- can be changed : coupe.dat -->  topt.dat, VF.dat (if 2d problems) ...
    end
    TITLE=FILE;
    fprintf(strcat('loading file for fig-1: (',FILE,')... '));
    data=load(FILE);	%- assumes file is 2d : list of i,j,val
    fprintf('DONE\n');
else
  fprintf('Error: this case not treated in output_view\n'); 
  return; 
end


if VALUE_PB
  fprintf('WARNING : VALUE PB ==> reformating for (Value.dat), to remove "INF" values \n');
  INF=1.e5;
  i=find(data(:,3)<INF);
  ValueMax=max(data(i,3));
  level_set=ValueMax*(1-1e-2);
  fprintf('WARNING :          ==> setting level_set =  %5.2f (max(|Value|) for |Value|<INF)\n',level_set);
end


%- OK for dim=2: cdd=[d1,d2].
if (dim>=2)
  xi=data(:,1:2);
  val=data(:,3); 

  if (SET_COUPE==1)
    XMIN = repmat(xmin(cdd),length(xi),1);
    DX   = repmat(dx(cdd),length(xi),1);
  else
    XMIN = repmat(xmin(1:2),length(xi),1);
    DX   = repmat(dx(1:2),length(xi),1);
  end
  x = XMIN + xi.*DX + (1-MESH)*DX/2;
end

%----------------------------------
%- PREPARATION OF MESH FOR GRAPHICS
%----------------------------------
xmesh=x(1:n1,1);
ymesh=x(n1*(1:n2),2);
xmeshi=xmesh;   ymeshi=ymesh;
xmini=xmin(d1); xmaxi=xmax(d1); 
ymini=xmin(d2); ymaxi=xmax(d2); 


empty=0; %- security variable: empty=1 means there is no level set, then no legends are drawns (to prevent error)

%----------------------------------------
%- STARTS REACHABLE PLOT &/OR TMIN PLOT
%----------------------------------------
if (PLOT_REACHABLE==1 || PLOT_CONTOUR==1)
  figure(1);
  clf;

  if PLOT_REACHABLE
    i=find(val<=level_set);
    %level_temp=min(min(val)); i=find(val<=level_temp); %- this better for value problem
    xg=x(i,1);
    yg=x(i,2);
    if OCTAVE; color_app='b*';
    else     ; color_app='b.';
    end 
    if length(i)>0
      GRAPH_REACHABLE=plot(xg,yg,color_app,'MarkerSize',5);
    else
      fprintf('no level set ==> no level set plot !\n');
    end

    %title(strcat('level set=',num2str(level_set)));
    axis square;
    hold on;

    COMPLEMENT=0;
    %------------------
    %- auxiliary plots
    %------------------
    if COMPLEMENT
      %complement1;
      %complement_xpp; %- this only for err. estimates for the case of "data_xpp.h"
      complement2;
      axis equal; %- to have normal scalings 
    end 
    %------------------
    %- end of auxiliary plots
    %------------------
  end


  %--------------------------------------------
  %- THE FOLLOWING WILL ONLY PLOT THE CONTOURS
  %--------------------------------------------

  grid on;

  if V_EXACT==0
    if PLOT_CONTOUR
      xmesh=x(1:n1,1);
      ymesh=x(n1*(1:n2),2);
      V=zeros(n1,n2); V(:)=val;
      %level=min(min(abs(val))); %- this better for value problem
      level=[level_set, level_set];
      [zz,GRAPH_CONTOUR]=contour(xmesh,ymesh,V',level,color_contour,'LineWidth',lw_contour);
      axis square;
      axis([xmini,xmaxi,ymini,ymaxi]);
      hold on;
    end
  end


  if V_EXACT==1
        axis([xmin(cdd(1)),xmax(cdd(1)),xmin(cdd(2)),xmax(cdd(2))]);
        hold on        
        % figure(1);  %- (remettre graphique sur la fig 1 sans detruire / recreer fenetre)
        if PLOT_CONTOUR
            xmesh=x(1:n1,1);
            ymesh=x(n1*(1:n2),2);
            V=zeros(n1,n2);
	    
	    V(:)=val;
            %level=min(min(abs(val))); %- this better for value problem
            level=[level_set, level_set];
            [zz,GRAPH_CONTOUR]=contour(xmesh,ymesh,V',level,color_contour,'LineWidth',3);
            empty=(length(zz)==0);

            FILE=[filePrefix 'coupeex.dat']; 
	    fprintf(strcat('loading file for fig-1: (',FILE,')... '));
            vexa=load(FILE);
            fprintf('DONE\n');

            V(:)=vexa(:,3);
            %level=min(min(abs(vexa(:,3)))); %- this better for value problem
            level=[level_set, level_set];
            [zz,GRAPH_CONTOUR_EX]=contour(xmesh,ymesh,V',level,color_contour_ex,'LineWidth',2);
            emptyloc=(length(zz)==0);
            empty=min(empty,emptyloc);
             
        end
        if exist('Front_ex');
	    %disp('COUCOU');
            t=T;
            X=Front_ex(t);
            %xe=X(:,cdd(1)); ye=X(:,cdd(2));
            xe=X(:,1); ye=X(:,2);
	    %plot(xe,ye,'g.');
            [zz,GRAPH_CONTOUR_EX]=plot(xe,ye,'k--','LineWidth',2);
        end

  end


  if PLOT_CONTOUR && ~empty
    if (OCTAVE)
      if exist('GRAPH_CONTOUR_EX');
        legend(legend_tit_scheme,legend_tit_ex,'Location',legend_pos);
      else
        legend(legend_tit_scheme,'Location',legend_pos);
      end
    else
      if exist('GRAPH_CONTOUR_EX');
        legend([GRAPH_CONTOUR(1), GRAPH_CONTOUR_EX(1)],legend_tit_scheme,legend_tit_ex,'Location',legend_pos);
      else
        legend([GRAPH_CONTOUR(1)],legend_tit_scheme,'Location',legend_pos);
      end
    end
  end


  if AXIS_EQUAL
    axis equal; %pause(0.5);
    axis([xmini,xmaxi,ymini,ymaxi]);
  end
  grid on;

  %if (level_set==0) 
  %  xlabel('x'); ylabel('y'); title(strcat('t=',num2str(T)))
  %else
  xlabel('x'); ylabel('y');
  TITLE=strcat('t=',num2str(T),'; level set=', num2str(level_set));
  title(TITLE);
  %end

end %- end of if (PLOT_REACHABLE==1 | PLOT_CONTOUR==1)
%-----------------------
%- ENDS REACHABLE PLOT
%-----------------------


if PLOT_3d==1

    figure(2);
    clf;

    FILE='VF.dat';

    %if dim>=3; FILE='coupe.dat'; end;
    if (dim>=3) || SET_COUPE; FILE=[filePrefix 'coupe.dat']; end;
    
    if VALUE_PB; 
      if dim==3
        FILE=[filePrefix 'Value.dat']; 
      else
        fprintf('ENCOURS... dim=%i(<>3) and VALUE_PB=1 not permitted for PLOT_3d=1\n',dim); 
	return;
      end
    end

    if PLOT_TMIN
       if dim>2; 
	  fprintf('PLOT_TMIN=1 and dim>2 not fixed - abort this part for fig 2\n');
	  return;
       end
       tminFile0='toptcoupe.dat';
       tminFile=[filePrefix tminFile0];
       FILE=tminFile;
       [status,result]=unix(strcat('ls ./',tminFile));
       if status==2; 
         fprintf('TMIN=1, but no file topt.dat..\n');
         return; 
       end
    end

    TITLE=FILE;
    fprintf(strcat('loading file for fig-2: (',FILE,')... '));
    fprintf('DONE\n');
    data=load(FILE);
    val=data(:,3);
    if VALUE_PB; val=min(val,ValueMax); end;

    if PLOT_TMIN==0
      xmesh=x(1:n1,1);
      ymesh=x(n1*(1:n2),2);
      V=zeros(n1,n2);
      V(:)=val;
      colormap hsv;
      surf(xmesh,ymesh,V');


      if (max(V)-min(V))/(0.001+max(abs(V)))<=0.02e-1; %-special case when V=const !
  	 Zmin=floor(min(min(V))); 
	 Zmax=Zmin+1; 
	 fprintf('Careful : SPECIAL CASE Zmin close to Zmax ==> setting Zmax=Zmin+1.\n');
         %axis([xmin(1),xmax(1),xmin(2),xmax(2),Zmin,Zmax]);
      else
        Zmin=min(V(:));
        Zmax=max(V(:));
      end
      title(TITLE);
      xlabel('x')
      ylabel('y')
      view(ALPHA,BETA);
      %Zmin=-0.5; Zmax=0.5;
      %axis([xmin(cdd(1)),xmax(cdd(1)),xmin(cdd(2)),xmax(cdd(2)),Zmin,Zmax]);
      axis([xmini,xmaxi,ymini,ymaxi,Zmin,Zmax]);
    end 


    if PLOT_TMIN==1
      %- Will bound max(tmin) and set to 0 all values s.t T=TMAX (that is, when T=+INF)
      TMAX=1.5*T;
      val=min(TMAX,val);	      

      j=find(val==0); 

      i=find(val==TMAX); val(i)=0.0;

      level=linspace(0,T,30);
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
        xmini=-2.0;  xmaxi=2.0;
        ymini=-2.0;  ymaxi=2.0;
        i=find(xmesh>=xmini  & xmesh<=xmaxi); xmeshi=xmesh(i);
        j=find(ymesh>=ymini  & ymesh<=ymaxi); ymeshi=ymesh(j);
        Vi=V(i,j);
      end
      colormap hsv;
      surf(xmeshi,ymeshi,Vi'); hold on;
      %- end of modif

      if length(j)>0
        xx=x(j,1); yy=x(j,2); 
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



if (PLOT_3d_ex==1) %- plot of VEX.dat
  if V_EXACT==1

    figure(3);
    clf;

    FILE=[filePrefix 'VEX.dat'];
    if dim>=2; FILE=[filePrefix 'coupeex.dat']; end;

    TITLE=FILE;
    fprintf(strcat('loading file for fig-3: (',FILE,')... '));
    fprintf('DONE\n');
    data=load(FILE);
    val=data(:,3);

    V=zeros(n1,n2);
    V(:)=val;
    colormap hsv;
    surf(xmesh,ymesh,V');
    title(TITLE);
    xlabel('x')
    ylabel('y')
    view(ALPHA,BETA);
    Zmin=-0.5; Zmax=2.0;
    Zmin=min(V(:));
    Zmax=max(V(:));
    axis([xmini,xmaxi,ymini,ymaxi,Zmin,Zmax]);

  else %- V_EXACT==0
    fprintf('VEX doesnot seem not exist - abort 3d VEX plot\n');
  end
else
  %- removing the 3d_ex plot if it exists 
  if (ishandle(3)); delete(3); end;
end



if TARGET_PLOT

if (VALUE_PB==0 && dim~=2) || (VALUE_PB==1 && dim~=3)
    fprintf('TARGET_PLOT=1 only programmed for dim=2 (here dim=%i) (or dim=3 if VALUE_PB=1); Skipping target plot.\n',dim);
else

  %- Careful this may not work if VF0.dat not present or if FORMAT_FULLDATA=0
  for ifig=[1] 
  %for ifig=[1 2] 
  %if PLOT_REACHABLE==1;

  figure(ifig);
  hold on;

  FILE=[filePrefix 'VF0.dat'];
  if VALUE_PB; 
    FILE=[filePrefix 'Value0.dat'];
  end
  level_set0=0.00;
  if ifig==1; fprintf('TARGET_PLOT: .. '); end;
  [status, result]=unix(strcat('ls ./',FILE));
  if status==2; 
    fprintf(strcat('no file ./',FILE,'. Skipping.\n'));
    break;
  else
    if ifig==1; 
      fprintf(strcat('loading VF0  (',FILE,') for green plotting of {x,VF0<=0}...'));
    end
  end

  data=load(FILE);           %- 2d
  if size(data,2)==1
    fprintf('no standard format (may be FORMAT_FULL_DATA=0?). Aborting Target Plot.\n'); 
    return; 
  else
    if ifig==1; fprintf('DONE\n'); end;
  end

  %- this works only for 2d problems 
  val=data(:,3); 

  if 0
  xi(:,1:2)=data(:,1:2);
  if (SET_COUPE==1)
    XMIN = repmat(xmin(cdd),length(xi),1);
    DX   = repmat(dx(cdd),length(xi),1);
  else
    XMIN = repmat(xmin(1:2),length(xi),1);
    DX   = repmat(dx(1:2),length(xi),1);
  end
  x = XMIN + xi.*DX + (1-MESH)*DX/2;
  end

  i=find(val<=level_set0);
  %level_temp=min(min(val)); i=find(val<=level_temp); %- this better for value problem
  xg=x(i,1);
  yg=x(i,2);
  
  if 0 && OCTAVE; color_target='g*'; color_contour='g--';
  else     ; color_target='g.'; color_contour='g--';
  end 

  %- green region {VF0<=0}
  factor_points=2;
  legend_tit_VF0='Target set';
  GRAPH_VF0=plot(xg,yg,color_target,'MarkerSize',5*factor_points);

  %- adapt legend

  if (~empty)
  if V_EXACT==1
   %legend([GRAPH_REACH(1), GRAPH_CONTOUR_EX(1), GRAPH_VF0(1)],legend_tit_scheme,legend_tit_ex,legend_tit_VF0,'Location',legend_pos);
   legend([GRAPH_CONTOUR(1), GRAPH_CONTOUR_EX(1), GRAPH_VF0(1)],legend_tit_scheme,legend_tit_ex,legend_tit_VF0,'Location',legend_pos);
  else
   %legend([GRAPH_REACH, GRAPH_VF0],legend_tit_scheme,legend_tit_VF0,'Location',legend_pos);
   legend([GRAPH_CONTOUR, GRAPH_VF0],legend_tit_scheme,legend_tit_VF0,'Location',legend_pos);
  end
  end


  %- adding green contour {VF0==0} (tested for cdd=[1,2])
  if 0
    xmesh=x(1:n1,1);
    ymesh=x(n1*(1:n2),2);
    V=zeros(n1,n2); V(:)=val;
    level=[level_set, level_set];
    [zz,GRAPH_CONTOUR]=contour(xmesh,ymesh,V',level,color_contour,'LineWidth',3);
  end
  hold on;

  %end %- end of [if PLOT_REACHABLE==1]
  end %- end of [for ifig=[1]]

end
end


%------------------------
%- DEC 2012 : trajectory
%------------------------
if TRAJECTORY;

  fprintf('TRAJECTORY PLOT (first 2 components on 2d coupe) ...');
  SPEC=0; %- set SPEC=1 for special graphs (2d) / zoom change


  %- reading file "successTrajectories.dat"

  fprintf('\n');
  FILE=[filePrefix 'successTrajectories.dat'];
  fprintf('..Reading file: "%s" ...',FILE); 
  [status, result]=unix(strcat('ls ./',FILE));
  if status==2; 
    fprintf(strcat('no file ./',FILE,'. Skipping.\n'));
    return; 
  else
    %fprintf('\n');
  end
  success=load(FILE,'-ascii'); 
  success=success(:,end);
  %- end of reading file "successTrajectories.dat"


  for ifig=[1] 
  %for ifig=[1 2]  %- BUG 2015
  figure(ifig);
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


  linemarkersize=10;
  markersize=25;

  lw_success=2;		%- linewidth for successful traj (fig 1)
  lw_unsuccess=2;	%- linewidth for unsuccessful traj (fig 1)

  couli='k.'; 		%- color of trajectory 
  coul_init='r.';	%- color of initial point 
  coul_fina='b.';	%- color of terminal point
  lw=lw_success;
  if success(i)==0;     
    couli='k-.';         %- color of trajectory (unsuccessful trajectories)
    coul_init='r.';	%- color of initial point (unsuccessful trajectories)
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

  fprintf('\n'); 

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
      printf('Warning, only showing dimemsion <=6 of trajectories !\n');
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
    end

    end %- end for i=1:nbTraj

  end

  fprintf('\n');
end


if AXIS_EQUAL
  %for ifig=[2 1]
  for ifig=[1]
    figure(ifig);
    axis equal;
  end
end


fprintf('ALL DONE\n');

return
% grid(gca,'minor')

% complement axes pour example data_xppx.h
if (PLOT_CONTOUR || PLOT_REACHABLE)
  figure(1);
  plot([xmin(1);xmax(1)],[0;0],'k');
  plot([0;0],[xmin(2);xmax(2)],'k');
end


