function output_view_movie1d()
%- movie 1d

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

PLOT_DERIV=0; %- to show derivative of numerical solution
PRINTF    =0;   	%- for more printings
TIME_PAUSE=0.02; 	%- time pause between 2 graphics

VFile0    = MOVIE_INPUTNAME;
VexFile0  = MOVIE_INPUTNAMEex;

%- GENERAL PARAMETERS:
lw=1; %- linewidth

if VALUE_PB==1;
  VFile0='Value';
end

if PRINTF; fprintf('dim=1\n'); end
if (dimcoupe>1); fprintf('ordimcoupe(=%i)>1. Aborting.!\n',dimcoupe); return; end

if (1)
  fprintf('Entering 1d movie plot:\n');
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
  message='';
  [val,err]=loadFile(FILE,message,n1,dimcoupe); % -dimcoupe should be 2!-
  if err==1; return; end
  %----------------------------------------------------------



  if VALUE_PB
    fprintf(' ** WARNING : VALUE PB ==> reformating for (Value.dat), to remove "INF" values \n');
    INF=1.e5;
    %i=find(data(:,2)<INF);
    i=find(val<INF);
    %ValueMax=max(data(i,3));
    %level_set=ValueMax*(1-1e-2);
    ValueMax=3.0; fprintf(' ** WARNING : Forcing ValueMax=%10.5f;\n',ValueMax);
    level_set=ValueMax;
    fprintf(' ** WARNING :  Forcing level_set =  %5.2f (max(|Value|) for |Value|<INF)\n',level_set);
  end

  %-special setting of Zmin and Zmax according to 'VF.dat'
  %xi =data(:,1);
  %val=data(:,2); 
  if VALUE_PB; val=min(val,ValueMax); end;

  %-special setting of Zmin and Zmax
  %maxmax=max(max(val),-min(val)); Zmax=maxmax; Zmin=-maxmax;
  %Zmin=min(val); 
  %Zmax=max(val);
  vmin=min(val);
  vmax=max(val);
  Zmin=vmin-0.1*(vmax-vmin);
  Zmax=vmax+0.1*(vmax-vmin);
  %if Zmax==Zmin; Zmax=Zmin+1; end; %- to prevent bug
  %if Zmax==Zmin; Zmax=abs(Zmax)*1.1; Zmin=-Zmax; end; %- to prevent bug and to have 0 in [Zmin,Zmax]
  if Zmax==Zmin; Zmax=Zmax+1.0; Zmin=Zmin-1.0; end; %- to prevent bug and to have 0 in [Zmin,Zmax]

  %if MESH==1;
  %  nn=size(data,1)-1;
  %  x=xmin + (0:nn)'*dx;
  %  n1=nn+1;
  %else %- MESH==0;
  %  nn=size(data,1);
  %  x=xmin + ((1:nn)-0.5)'*dx;
  %  n1=nn;
  %end

  V=zeros(n1,1);

end

%figure();
drawperso();

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
  %val=data(:,2);
  message='';
  [val,err]=loadFile(FILE,message,n1,dimcoupe); % -dimcoupe should be 2!-
  if err==1; return; end


  if VALUE_PB; val=min(val,ValueMax); end;

  clf;
  hold on;

  if length(val) ~= length(V(:));
    disp( [length(val) , length(V(:))])
    fprintf('\n!!! File lengths not corresponding !.. Abort movie mode; restart using ./cleandat \n');
    return;
  end

  if OCTAVE
    %plot(x,val,'b-');  
    %plot(x,val,'b*','MarkerSize',5); 
    GRAPH=plot(x,val,'b.-','LineWidth',lw); 
  else
    GRAPH=plot(x,val,'b.-','LineWidth',lw); 
  end 

  %-------------------------------
  %- plot of a derivative estimate
  %-------------------------------
  if PLOT_DERIV;
    hold on;
    dx=x(2)-x(1);
    scaling=0.5;
    plot(x(2:end),scaling*(val(2:end)-val(1:end-1))/dx,'r.','LineWidth',1); 
  end

  axis([xmin,xmax,Zmin,Zmax]);
  grid on;
 

  if V_EXACT==1
    %FILE=[filePrefix VexFile0 '0.dat'];
    FILEex=strcat(filePrefix,VexFile0,int2str(i),'.dat');
    [status, result]=unix(strcat('ls ./',FILEex));
    if status==2
      fprintf(strcat('no file ./',FILEex,'. Skip'));
      %- trying with no {it} index:
      FILEex=[filePrefix VexFile0 '.dat'];
      [status, result]=unix(strcat('ls ./',FILEex));
      if status==2 
        fprintf(strcat('; ..no file ./',FILEex,'. Skip; Aborting\n'));
        return
      else
        fprintf('; using %s instead.',FILEex);
	PRINTF=1;
      end
    end

    %----------------
    %- LOADING FILEex
    %----------------
    %if PRINTF; fprintf('loading %s .. ',FILE); end;
    %data=load(FILE);
    %valex=data(:,2);
    message='';
    [valex,err]=loadFile(FILEex,message,n1,dimcoupe); % -dimcoupe should be 2!-
    if err==1; return; end


    if OCTAVE
      %plot(x,valex,'k*');  
      %GRAPH_EX=plot(x,valex,'k*','MarkerSize',5);
      GRAPH_EX=plot(x,valex,'k.-','MarkerSize',5,'LineWidth',lw);
    else
      GRAPH_EX=plot(x,valex,'k.-','LineWidth',lw); 
    end 

    %- 2020 complement error
    COMPLEMENT_ERROR_1d=1;
    if COMPLEMENT_ERROR_1d;
      color_err='r.-';
      scale=10;
      legend_tit_err=strcat('Error * ',num2str(scale));
      GRAPH_ERR=plot(x,scale*abs(val-valex),color_err); 
    end

  end

  %ti=floor(i*DT*1000)/1000;
  ti=i*DT;
  tiFile=['  (',TITLE_FILE,')'];
  title(strcat('t=',num2str(min(T,ti)),tiFile));

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

  if AXIS_EQUAL
    axis equal;
    axis([xmini,xmaxi,ymini,ymaxi]);
  end

  pause(TIME_PAUSE);

  if MOVIE_PAUSE; 
    aa=input(''); 
  else
    if PRINTF; fprintf('\n'); end;
  end

end

end %- end-of-function


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

    fprintf(strcat('no file ./',FILE,'. Skip; Aborting\n'));
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

