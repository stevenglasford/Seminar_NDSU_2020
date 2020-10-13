%- movie 2d

global DT
global VALUE_PB level_set
global Zmin Zmax
global nn cdd 
global XPAUSE

PRINTF=0; %- for more printings

VFile='VF';
VFile0='VF0';
if VALUE_PB==1;
  VFile='Value';
  VFile0='Value0';
end

if PRINTF; fprintf('dim=2\n'); end
if (dim~=2 && dimcoupe~=2);  fprintf('dim<>2! & dimcoupe<>2 !. Aborting.!\n'); return; end

if (1)
  fprintf('Entering 2d movie plot:\n');
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
    i=find(data(:,3)<INF);
    %ValueMax=max(data(i,3));
    %level_set=ValueMax*(1-1e-2);
    ValueMax=3.0; fprintf('WARNING : Forcing ValueMax=%10.5f;\n',ValueMax);
    level_set=ValueMax;
    fprintf('WARNING :  Forcing level_set =  %5.2f (max(|Value|) for |Value|<INF)\n',level_set);
  end


  %-special setting of Zmin and Zmax according to 'VF.dat'
  xi=data(:,1:2);
  val=data(:,3); 
  if VALUE_PB; val=min(val,ValueMax); end;

  %-special setting of Zmin and Zmax
  %maxmax=max(max(val),-min(val)); Zmax=maxmax; Zmin=-maxmax;
  %Zmin=min(val); 
  %Zmax=max(val);
  Zmin=min(val)-0.1*(max(val)-min(val)); 
  Zmax=max(val)+0.1*(max(val)-min(val));
  if Zmax==Zmin; Zmax=Zmin+1; end; %- to prevent bug

  %if (SET_COUPE==1)
  XMIN = repmat(xmin(cdd),length(xi),1);
  DX   = repmat(dx(cdd),length(xi),1);
  %else
  %XMIN = repmat(xmin(1:2),length(xi),1);
  %DX   = repmat(dx(1:2),length(xi),1);
  %end
  x = XMIN + xi.*DX + (1-MESH)*DX/2;

  n1=nn(cdd(1)); n2=nn(cdd(2));
  xmesh=x(1:n1,1);
  ymesh=x(n1*(1:n2),2);
  V=zeros(n1,n2);

end


for i=0:ITERMAX  % (ITERMAX, but may put ITERMAX-1 if not computing all steps to avoid a matlab error)

  ti=i*DT;
  FILE=strcat(filePrefix,VFile,int2str(i),'.dat');
  if (PRINTF||XPAUSE) 
    fprintf('iteration i=%3i (t=%8.5f); ',i,ti);
    fprintf('loading file %s ...',FILE); 
  end

  [status, result]=unix(strcat('ls ./',FILE));
  if status==2; 
    fprintf(strcat('no file ./',FILE,'. Skipping, Aborting\n'));
    return;
  else
    data=load(FILE);
  end
  val=data(:,3); 

  if VALUE_PB; val=min(val,ValueMax); end;

  figure(1);
  clf;
  hold on;

  if length(data) ~= length(V(:));
    fprintf('\n!!! File lengths not corresponding !.. Abort movie mode; restart using ./cleandat \n');
    return;
  end

  V(:)=val;
  %colormap  autumn; %parula; %winter;
  surf(xmesh,ymesh,V');
  xlabel('x')
  ylabel('y')
  view(ALPHA,BETA);
  axis([xmin(1),xmax(1),xmin(2),xmax(2),Zmin,Zmax]);
  grid on;

  %ti=floor(i*DT*1000)/1000;
  ti=i*DT;
  title(strcat('t=',num2str(min(T,ti))));
  pause(0.02);
  if XPAUSE; 
    aa=input(''); 
  end

end

