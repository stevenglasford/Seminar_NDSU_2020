function output_view_movie2d()
%- movie 2d

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
global vidObj

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


PRINTF    =0;   	%- for more printings
TIME_PAUSE=0.01; 	%- time pause between 2 graphics

VFile0    = MOVIE_INPUTNAME;
VexFile0  = MOVIE_INPUTNAMEex;

%SPECIAL_TOPT=1;	%- special case : set 1 for plotting topt. Will truncate Zmax to some maximal finite value as in VFile0.
SPECIAL_TOPT=PLOT_TOPT;	%- special case : set 1 for plotting topt. Will truncate Zmax to some maximal finite value as in VFile0.
if SPECIAL_TOPT
  VFile0='coupetopt'; % 2018
  ValueMax=4;
end

if VALUE_PB==1;
  VFile0='Value';
end

if PRINTF; fprintf('dim=2\n'); end
if (dim~=2 && dimcoupe~=2);  fprintf('dim<>2 & dimcoupe<>2 !. Aborting.!\n'); return; end

if (1)
  fprintf('Entering 2d movie plot:\n');
  %- this is the loading of a file (like VF0.dat)
  %- in order to get first components  (obsolete)
  %-   <or> info on target
  %-   <or> info on min/max of value
  %- not used afterwards.
  FILE=[filePrefix VFile0 '0.dat'];
  fprintf('Loading %s for initial settings (Zmin/Zmax).. \n',FILE);

  % [status, result]=unix(strcat('ls ./',FILE));
  % if status==2; 
  %   fprintf(strcat('no file ./',FILE,'. Skip; Aborting\n'));
  %   %break;
  %   return;
  % else
  %   fprintf(strcat('loading (',FILE,'):..'));
  %   data=load(FILE);
  % end
  % fprintf('DONE\n');

  %----------------------------------------------------------
  message='0.dat';
  [val,err]=loadFile(FILE,message,n1*n2,dimcoupe); % -dimcoupe should be 2!-
  if err==1; return; end
  %----------------------------------------------------------



  if VALUE_PB
    fprintf(' ** WARNING : VALUE PB ==> reformating for (Value.dat), to remove "INF" values \n');
    INF=1.e5;
    %i=find(data(:,3)<INF);
    i=find(val<INF);
    %ValueMax=max(data(i,3));
    %level_set=ValueMax*(1-1e-2);
    ValueMax=3.0; fprintf(' ** WARNING : Forcing ValueMax=%10.5f;\n',ValueMax);
    level_set=ValueMax;
    fprintf(' ** WARNING :  Forcing level_set =  %5.2f (max(|Value|) for |Value|<INF)\n',level_set);
  end

  % OBSOLETE: xi,val,x,n1,n2
  %-special setting of Zmin and Zmax according to 'VF.dat'
  if 0
    xi =data(:,1:2);
    val=data(:,3); 
    % 2018 - some recomputed mesh parameters; 
    XMIN = repmat(xmin(cdd),length(xi),1);
    DX   = repmat(dx(cdd),length(xi),1);
    x = XMIN + xi.*DX + (1-MESH)*DX/2;
    n1=nn(cdd(1)); n2=nn(cdd(2));
    %- initialisation of xmesh,ymesh, V in a 2d setting
    %xmesh=x(1:n1,1);
    %ymesh=x(n1*(1:n2),2);
  end
  %- NEW jun 2018

  if 0
  d1 =cdd(1); d2 =cdd(2);
  dx1=dx(d1); dx2=dx(d2);
  n1 =nn(d1); n2 =nn(d2);
  %m1 =MESH(d1); m2=MESH(d2); 
  xmesh=xmin(d1) + (0:n1-1)*dx1 + (1-MESH)*dx1/2;
  ymesh=xmin(d2) + (0:n2-1)*dx2 + (1-MESH)*dx2/2;
  xmeshi=xmesh;   ymeshi=ymesh;
  end

  empty=0; %- security variable: empty=1 means there is no level set, then no legends are drawns (to prevent error)
  V=zeros(n1,n2);

  if VALUE_PB; val=min(val,ValueMax); end;

  if SPECIAL_TOPT 
    INF=1.e5;
    jj=find(val>=INF);
    val(jj)=0;
    Zmin=floor(min(val));
    Zmax=floor(max(val))+1;
  end

  %-special setting of Zmin and Zmax
  %maxmax=max(max(val),-min(val)); Zmax=maxmax; Zmin=-maxmax;
  %Zmin=min(val); 
  %Zmax=max(val);
  vmin=min(val);
  vmax=max(val);
  Zmin=vmin-0.2*(vmax-vmin);
  Zmax=vmax+0.2*(vmax-vmin);
  %if Zmax==Zmin; Zmax=Zmin+1; end; %- to prevent bug
  %if Zmax==Zmin; Zmax=abs(Zmax)*1.1; Zmin=-Zmax; end; %- to prevent bug and to have 0 in [Zmin,Zmax]
  if Zmax==Zmin; Zmax=Zmax+1.0; Zmin=Zmin-1.0; end; %- to prevent bug and to have 0 in [Zmin,Zmax]


end

%figure();
%drawperso();

for i=0:ITERMAX  % (ITERMAX, but may put ITERMAX-1 if not computing all steps to avoid a matlab error)

  ti=i*DT;
  FILE=strcat(filePrefix,VFile0,int2str(i),'.dat');
  TITLE_FILE=FILE;

  if (PRINTF||MOVIE_PAUSE) 
    fprintf('iteration i=%3i (t=%8.5f); ',i,ti);
    fprintf('loading %s .. ',FILE); 
  end

  PRINTloc=1;
  %if i>=2; PRINTloc=0; end

  [status, result]=unix(strcat('ls ./',FILE));
  if status==2; 
    if (PRINTloc); fprintf(strcat('no file ./',FILE,'. Skip')); end

    %- trying with no {it} index 
    %- This may be used for the last iteration
    FILE=[filePrefix VFile0 '.dat'];
    [status, result]=unix(strcat('ls ./',FILE));
    if status==2 
      if (PRINTloc); fprintf(strcat('; no file ./',FILE,'. Skip; Aborting\n')); end
      return
    else
      if (PRINTloc); fprintf('; using %s instead.',FILE); end;
      %if ~PRINTF; fprintf('\n'); end;
      i=ITERMAX;  %- in order to break from the for loop
    end

    %fprintf('; Aborting\n'));
    %return
  end


  %----------------
  %- LOADING FILE
  %----------------
  %data=load(FILE);
  %val=data(:,3);
  message='';
  [val,err]=loadFile(FILE,message,n1*n2,dimcoupe); % -dimcoupe should be 2!-
  if err==1; return; end

  if VALUE_PB; val=min(val,ValueMax); end;

  if SPECIAL_TOPT
    if i==0
       fprintf(' ** WARNING : SPECIAL FOR TOPT plot ==> reformating topt.dat, to remove "INF" values\n');
    end
    INF=1.e5;
    %----------------------------------------------------
    %- MISE A ZERO DES VALEURS >= INF (POUR LE GRAPHIQUE)
    %----------------------------------------------------
    %jj=find(data(:,3)>=INF);
    %val=data(:,3);
    jj=find(val>=INF);
    val(jj)=0; 
    %Zmax=floor(max(val))+1;
    Zmax=ValueMax;
  end

  clf;
  hold on;

  if length(val) ~= length(V(:));
    disp( [length(val) , length(V(:))])
    fprintf('\n!!! File lengths not corresponding !.. Abort movie mode; hint: may be use ./cleandat and restart ?\n');
    return;
  end

  %- PLOT_LEVEL_SET: 0/1 to decide between plotting level set (=call to plot reachable) or surface plot 
  %MOVIE_PLOT_LEVEL_SET=0; 

  if (V_EXACT==1)
    %------------------------------------------------------------------
    %- 2018 : Attempt to furthermore read exact file corresponding to current data
    %-        preparation of "val_ex" before calling plot_reachable()
    %------------------------------------------------------------------
    %FILEex=strcat(filePrefix,VexFile0,int2str(i),'.dat');
    FILEex=[filePrefix VexFile0 int2str(i) '.dat'];

    if (PRINTF||MOVIE_PAUSE) 
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

  %fprintf('\n'); disp([xmini,xmaxi]); 
  if MOVIE_PLOT_REACHABLE
    %- level set plot
    plot_reachable(i);
  end

  if 1 
    %- classical surf plot

    %- surface plot of "val"
    V(:)=val;
    %colormap  autumn; %parula; %winter;

    figure(2); %- normal way
    clf;
    %figure(1); %- test all videos in fig 1
    %subplot(1,2,2); %- if subplot then both fig2 and fig1 will be in fig 1
    if (i==0)
      drawperso();
    end

    %hold on;
    surf(xmesh,ymesh,V');
    xlabel(tit1)
    ylabel(tit2)

    %view(ALPHA,BETA);
    %axis([xmin(1),xmax(1),xmin(2),xmax(2),Zmin,Zmax]);

    axis([xmini,xmaxi,ymini,ymaxi,Zmin,Zmax]);
    view(ALPHA,BETA);
    grid on;

    ti=i*DT; %ti=floor(i*DT*1000)/1000;
    tiFile=['  (',TITLE_FILE,')'];
    title(strcat('t=',num2str(min(T,ti)),tiFile));


    %------------------
    %- auxiliary plots
    %------------------
    if 0
      FILE=['complement.m'];
      %fprintf('.. Reading file: "%s" ... ',FILE); 
      [status, result]=unix(strcat('ls ./',FILE));
      if status==2; 
        fprintf(strcat('no file ./',FILE,'. Skipping.\n'));
      else
        %fprintf('DONE.\n');
        complement;
      end
    end

    if TEST 
      figure(2); %- normal 
      %figure(1); %- test video
      view(0,0);
      %axis equal;
      axis square;
    end
    %------------------
    %- end of auxiliary plots
    %------------------
  end

  %--------------
  %- 2018 ENCOURS
  %--------------
  if PLOT_3d_ex 

    %- plot of Vex.dat/coupeex_xx.dat or related

    if V_EXACT && exist('val_ex');
    Vex=zeros(n1,n2);
    Vex(:)=val_ex;

    figure(3);
    clf;
    if (i==0)
      drawperso();
    end

    hold on;
    surf(xmesh,ymesh,Vex');
    xlabel(tit1)
    ylabel(tit2)

    axis([xmini,xmaxi,ymini,ymaxi,Zmin,Zmax]);
    view(ALPHA,BETA);
    %axis([xmin(1),xmax(1),xmin(2),xmax(2),Zmin,Zmax]);
    %axis([xmini,xmaxi,ymini,ymaxi,Zmin,Zmax]);
    grid on;

    ti=i*DT; %ti=floor(i*DT*1000)/1000;
    tiFile=['  (',TITLE_FILEex,')'];
    title(strcat('t=',num2str(min(T,ti)),tiFile));
    %if OCTAVE; set(gcf,'visible','on'); end;
    end
  end

  %--------------
  %- 2018 FIN ENCOURS
  %--------------

  pause(TIME_PAUSE);

  if VIDEO && i==0
    fprintf("VIDEO : using figure no %i.\n",VIDEOfig);
    figure(VIDEOfig)
    video_init(videoFile);
    for j=1:10; video_capture(); end
  end

  if VIDEO
    figure(VIDEOfig)
    video_capture();
  end

  if MOVIE_PAUSE; 
    aa=input('','s'); 
    if aa=='q'; return; end;
  else
    if PRINTF; fprintf('\n'); end;
  end

end


if VIDEO
  figure(VIDEOfig)
  video_end();
end

end %- end-of-function

%------------------
%- video functions
%------------------

function video_init(filename)
  global vidObj
  %global axvideo %TEST

  fprintf('VIDEO: video filename is: %s;\n',filename);

  % Prepare the new file.
  %vidObj = VideoWriter('network.avi');
  vidObj = VideoWriter(filename);
  vidObj.FrameRate=10;
  open(vidObj);
  % Create an animation.
  %Z = peaks; surf(Z); 
  %axis tight
  %plootvv(VV);

  set(gca,'nextplot','replacechildren');
  %axvideo=gca() %TEST

  %- 10 more of the first view
  for i=1:10; video_capture(); end

end

function video_capture()
  global vidObj
  %global axvideo %TEST

  currFrame = getframe;
  %currFrame = getframe(axvideo); %- TEST
  %myFrame = getframe(gca);
  %size(myFrame.cdata)
  %writeVideo(v,myFrame)
  %pause
  %a=size(currFrame.cdata);
  %fprintf("size(currFrame)=%3i %3i %3i\n",a(1),a(2),a(3))
  writeVideo(vidObj,currFrame);
end

function video_end()
  global vidObj
  global videoFile


  % 20 more of the last view
  for i=1:30; video_capture(); end

  % Close the file
  close(vidObj);

  fprintf('Writing video file in %s\n',videoFile);
end

%-----------------
%- OCTAVE figure() equivalent
%------------------
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

    fprintf(strcat('no file ./',FILE,'. Skip; Aborting load.\n'));
    %break;
    val=0;
    err=1;
    return

  else

    if PRINTloc; fprintf('loading %12s (%s):..',FILE,message); end
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

