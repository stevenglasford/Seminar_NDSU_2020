%- movie 1d

global DT
global VALUE_PB level_set
global Zmin Zmax
global nn cdd 
global XPAUSE

PRINTF=0; %- for more printings

PLOT_DERIV=0; %- to show derivative of numerical solution

VFile='VF';
VFile0='VF0';
if VALUE_PB==1;
  VFile='Value';
  VFile0='Value0';
end

if PRINTF; fprintf('dim=1 !\n'); end
if dim>1;  fprintf('dim>1. Aboring.!\n'); return; end

if (1)
  fprintf('Entering 1d movie plot:\n');
  %- this is the loading of a file (VF.dat, or VF0.dat)
  %- in order to get the first too components 
  %- not used afterwards.
  FILE=[filePrefix VFile0 '.dat'];
  fprintf('loading file %s for initial settings (Zmin/Zmax)...\n',FILE);

  [status, result]=unix(strcat('ls ./',FILE));
  if status==2; 
    fprintf(strcat('no file ./',FILE,'. Skipping, Aborting\n'));
    %break;
    return;
  else
    fprintf(strcat('loading (',FILE,'): '));
    data=load(FILE);
  end
  fprintf('DONE\n');


  if VALUE_PB
    fprintf('WARNING : VALUE PB ==> reformating for (Value.dat), to remove "INF" values \n');
    INF=1.e5;
    i=find(data(:,2)<INF);
    %ValueMax=max(data(i,3));
    %level_set=ValueMax*(1-1e-2);
    ValueMax=3.0; fprintf('WARNING : Forcing ValueMax=%10.5f;\n',ValueMax);
    level_set=ValueMax;
    fprintf('WARNING :  Forcing level_set =  %5.2f (max(|Value|) for |Value|<INF)\n',level_set);
  end


  %-special setting of Zmin and Zmax according to 'VF.dat'
  xi=data(:,1);
  val=data(:,2); 
  if VALUE_PB; val=min(val,ValueMax); end;

  %-special setting of Zmin and Zmax
  %maxmax=max(max(val),-min(val)); Zmax=maxmax; Zmin=-maxmax;
  %Zmin=min(val); 
  %Zmax=max(val);
  Zmin=min(val)-0.1*(max(val)-min(val)); 
  Zmax=max(val)+0.1*(max(val)-min(val));
  if Zmax==Zmin; Zmax=Zmin+1; end; %- to prevent bug

  if MESH==1;
    nn=size(data,1)-1;
   x=xmin + (0:nn)'*dx;
    n1=nn+1;
  else %- MESH==0;
    nn=size(data,1);
    x=xmin + ((1:nn)-0.5)'*dx;
    n1=nn;
  end
  V=zeros(n1,1);

end


for i=0:ITERMAX  % (ITERMAX, but may put ITERMAX-1 if not computing all steps to avoid a matlab error)

  ti=i*DT;
  FILE=strcat(filePrefix,VFile,int2str(i),'.dat');
  if (PRINTF||XPAUSE) 
    fprintf('iteration i=%3i (t=%8.5f); ',i,ti);
    fprintf('loading file %10s ...',FILE); 
  end

  [status, result]=unix(strcat('ls ./',FILE));
  if status==2; 
    fprintf(strcat('no file ./',FILE,'. Skipping, Aborting\n'));
    return;
  else
    data=load(FILE);
  end
  val=data(:,2); 

  if VALUE_PB; val=min(val,ValueMax); end;

  figure(1);
  clf;
  hold on;

  if length(val) ~= length(V(:));
    fprintf('\n!!! File lengths not corresponding !.. Abort movie mode; restart using ./cleandat \n');
    return;
  end

  if OCTAVE
    plot(x,val,'b-');  
    plot(x,val,'b*','MarkerSize',5); 
  else
    plot(x,val,'b.-','LineWidth',1); 
    %- plot of a derivative estimate
    if PLOT_DERIV;
      hold on;
      dx=x(2)-x(1);
      scaling=0.5;
      plot(x(2:end),scaling*(val(2:end)-val(1:end-1))/dx,'r.','LineWidth',1); 
    end
  end 
  axis([xmin,xmax,Zmin,Zmax]);

  if V_EXACT==1
    VFile0='VEX';
    FILE=[filePrefix VFile0 '.dat'];
    if PRINTF; fprintf('loading file %s ...\n',FILE); end;
    data=load(FILE);
    valex=data(:,2);
    if OCTAVE
      plot(x,valex,'k-');  
      plot(x,valex,'k*','MarkerSize',5); 
    else
      plot(x,valex,'k.-','LineWidth',1); 
    end 
    grid();
  end

  %ti=floor(i*DT*1000)/1000;
  ti=i*DT;
  title(strcat('t=',num2str(min(T,ti))));
  pause(0.02);
  if XPAUSE; 
    aa=input(''); 
  end

end

