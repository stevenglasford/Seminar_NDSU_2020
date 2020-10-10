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

