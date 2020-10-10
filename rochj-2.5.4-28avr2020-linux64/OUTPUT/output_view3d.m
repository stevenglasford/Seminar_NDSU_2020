function complement_isoview()
%- 3d ISOVIEW PLOT

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
global VIDEO videoFile
global vidObj


global tit1 tit2 tit3  %- xlabel and ylabel for graph

global xmin3d xmax3d nn3d

  n1=nn3d(1); n2=nn3d(2); n3=nn3d(3);

  fprintf('3d: isosurface ..\n');
  if ~(dim==3 || dimcoupe==3)
    fprintf('** WARNING : dim(%i) <>3 or dimcoupe(%i)<>3 ==> Skip 3d plot.\n',dim,dimcoupe);
    return
  end

  FILE=[filePrefix VFile0 '.dat'];
  TITLE=FILE;
  %----------------
  %- LOADING FILE
  %----------------
  %data=load(FILE);
  %val=data(:,4);
  message='fig-3d';
  [val,err]=loadFile(FILE,message,n1*n2*n3,3); % -dim should be 3!-
  if err==1; return; end


  clf;
  %drawperso();

  %[x,y,z,v] = flow;
  %p = patch(isosurface(x,y,z,v,-3));
  xmesh=linspace(xmin3d(1),xmax3d(1),nn3d(1))';
  ymesh=linspace(xmin3d(2),xmax3d(2),nn3d(2))';
  zmesh=linspace(xmin3d(3),xmax3d(3),nn3d(3))';
  V=zeros(nn3d(1),nn3d(2),nn3d(3));
  V(:)=val;
  %clear data; clear val;


  if OCTAVE

    fprintf('3d isoview not programmed yet in octave..');

  else %- MATLAB

  %for ll=[-0.02:0.01:0.02]
    %figure(2)
    %level_set=ll;
    %level_set
    %size(V)
    %size(xmesh)
    %size(ymesh)
    %size(zmesh)
    GRAPH3d=isosurface(ymesh,xmesh,zmesh,V,level_set);
    %GRAPH3d=isosurface(xmesh,ymesh,zmesh,V,level_set);
    %[X,Y,Z]=meshgrid(xmesh,ymesh,zmesh);
    %size(X)
    %size(Y)
    %size(Z)
    %GRAPH3d=isosurface(X,Y,Z,V,level_set);

    %GRAPH3d=isosurface(xmesh,ymesh,zmesh,V-ll,0.0);
    p=patch(GRAPH3d);
    axis([xmin3d(1),xmax3d(1),xmin3d(2),xmax3d(2),xmin3d(3),xmax3d(3)]);
    grid on; grid minor;
    %isonormals(x,y,z,v,p)
    p.FaceColor = 'red';
    p.EdgeColor = 'none';
    daspect([1,1,1]);
    %view(3); % axis tight
    camlight;
    lighting gouraud;

    axis square;

    xlabel(tit2)
    ylabel(tit1)
    zlabel(tit3)

    title('scheme')
    view(ALPHA,BETA);
    %title('scheme ENO2')
    %title('scheme ENO3')
    %title('scheme NBEE-OCOPT')
    %legend('NBEE scheme')
    %out=input('appuyer sur une touche:');
  %end

  end

  %fprintf('; DONE\n');

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

