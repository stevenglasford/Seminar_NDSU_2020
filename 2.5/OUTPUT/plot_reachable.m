%----------------------------------------
%- STARTS REACHABLE PLOT &/OR TMIN PLOT
%----------------------------------------
%- needed: x,val

function plot_reachable(ii)
%% if ii<= 0 then will call drawperso()
%% if ii==-1 then will load "coupe.dat" otherwise load coupe{ii}.dat for graphics

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

global vidObj
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




%fprintf('PLOT_REACHABLE = %i\n', PLOT_REACHABLE);
%return

empty=0; %- security variable: empty=1 means there is no level set, then no legends are drawns (to prevent error)

%----------------------------------------
%- STARTS REACHABLE PLOT &/OR TMIN PLOT
%----------------------------------------

if 0
if (OCTAVE &&  ~exist('xmini(1)')) %|| (~OCTAVE && ~exist('xmini'))
  fprintf('Re-Definition of xmini/xmaxi/ymini/ymax as xmin/xmax etc...:\n');
  xmini=xmin(cdd(1)); xmaxi=xmax(cdd(1));
  ymini=xmin(cdd(2)); ymaxi=xmax(cdd(2));
end
end


%----------------------------------------
%- REACHABLE PLOT
%----------------------------------------
if (PLOT_REACHABLE==1 || PLOT_CONTOUR==1)

  figure(1);
  clf;
  %subplot(1,2,1); %- movie option
  if (ii<=0)
    drawperso();
  end

  hold on;

  if PLOT_REACHABLE
    i=find(val<=level_set);
    xg=x(i,1);
    yg=x(i,2);
    color_app='b.';
    if length(i)>0
      axis([xmini,xmaxi,ymini,ymaxi]);
      GRAPH_REACHABLE=plot(xg,yg,color_app,'MarkerSize',5);
    else
      fprintf('no level set ==> no level set plot !\n');
    end

  end

  grid on;

  %--------------------------------------------
  %- THE FOLLOWING WILL PLOT THE CONTOUR
  %--------------------------------------------
  %ICI
  if PLOT_CONTOUR
    V=zeros(n1,n2); V(:)=val;
    %level=min(min(abs(val))); %- this better for value problem
    level=[level_set, level_set];
    [zz,GRAPH_CONTOUR]=contour(xmesh,ymesh,V',level,color_contour,'LineWidth',lw_contour);
    axis([xmini,xmaxi,ymini,ymaxi]);
    hold on;
  end

  if V_EXACT
    if PLOT_CONTOUR
      V=zeros(n1,n2); V(:)=val;
      %level=min(min(abs(val))); %- this better for value problem
      level=[level_set, level_set];
      [zz,GRAPH_CONTOUR]=contour(xmesh,ymesh,V',level,color_contour,'LineWidth',lw_contour);
      %axis square;
      axis([xmini,xmaxi,ymini,ymaxi]);
      hold on;
    end
  end

  %- 2018 : plot of "val" and "val_ex" 
  if V_EXACT
    hold on        
    figure(1);  %- (remettre graphique sur la fig 1 sans detruire / recreer fenetre)
    
    level=[level_set, level_set];
    V=zeros(n1,n2);

    if PLOT_CONTOUR
      V(:)=val;
      %level=min(min(abs(val))); %- this better for value problem
      [zz,GRAPH_CONTOUR]=contour(xmesh,ymesh,V',level,color_contour,'LineWidth',3);
      empty=(length(zz)==0);

      %FILEEX=[filePrefix 'coupeex.dat']; 
      %fprintf(strcat('loading file for fig-1: (',FILEEX,', for exact contour) ...'));
      %val_ex=load(FILEEX);
      %fprintf('DONE\n');

      V(:)=val_ex;
      %level=min(min(abs(valex(:,3)))); %- this better for value problem
      [zz,GRAPH_CONTOUR_EX]=contour(xmesh,ymesh,V',level,color_contour_ex,'LineWidth',2);
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


  if PLOT_CONTOUR && ~empty
    if (OCTAVE)
      if exist('GRAPH_CONTOUR_EX');
        legend(legend_tit_scheme,legend_tit_ex,'Location',legend_pos);
      else
        legend(legend_tit_scheme,'Location',legend_pos);
      end
    else
      if exist('GRAPH_CONTOUR_EX(1)');
        %if (length(GRAPH_CONTOUR_EX(1))>0) 
        legend([GRAPH_CONTOUR(1), GRAPH_CONTOUR_EX(1)],legend_tit_scheme,legend_tit_ex,'Location',legend_pos);
	%pause
	%end
      else
        legend([GRAPH_CONTOUR(1)],legend_tit_scheme,'Location',legend_pos);
      end
    end
  end


  if AXIS_EQUAL
    axis equal; 
    axis([xmini,xmaxi,ymini,ymaxi]);
  end
  grid on;

  %if (level_set==0) 
  %  xlabel('x'); ylabel('y'); title(strcat('t=',num2str(T)))
  %else
  %xlabel('x'); ylabel('y');
  %TITLE=strcat('t=',num2str(T),'; level set=', num2str(level_set));
  %title(TITLE);
  %end
  
  ti=ii*DT;
  tiFile=['  (',TITLE_FILE,')'];
  title(strcat('t=',num2str(min(T,ti)),tiFile,'; level set=',num2str(level_set)));
  xlabel(tit1)
  ylabel(tit2)

end %- end of if (PLOT_REACHABLE==1 | PLOT_CONTOUR==1)
return
end
%-----------------------
%- ENDS REACHABLE PLOT
%-----------------------

%-----------------
%- OCTAVE figure() equivalent
%------------------
function drawperso()
  global OCTAVE
  if OCTAVE; set(gcf,'visible','off'); set(gcf,'visible','on'); end;
end
