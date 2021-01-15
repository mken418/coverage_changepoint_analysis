%This script was first written by Dr. Liumdmila Mainzer and edited to allow for more general use by Matt Kendzior 
%Callable on an input directory with per-base individual per base coverage files. At the moment this script only works with output from Sentieon's coverage metricss algorithm
%These files look like:
%Locus   Total_Depth     sample_Average_Depth    s_3595-ST-0709_Depth
%chr1:866431     47      47.00   47
%chr1:866432     48      48.00   48
%...


function [] = SignalDetection_Mar23_CL_args(putatVarStart, putatVarStop, filesDir, output)
% set the coordinates for the putative variant we are trying to detect
PutativeVariantStart=putatVarStart;
PutativeVariantStop=putatVarStop;

% find the coverage files
coveragefiles = dir(filesDir);
NumFiles=size(coveragefiles);


% create the output filer and print the headers
OutputFile = fopen(output,'w');
fprintf(OutputFile,'%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n','Sample','Chr','FirstChangePt','SecondChangePt','Residual','BeforeVarCoverageMean','BeforeVarCoverageStd',...
           'CoverageInVariantMean','CoverageInVariantStd','AfterVarCoverageMean','AfterVarCoverageStd','Zygocity');



% total count of individuals with the variant will be printed to the Matlab console
CountIndividualsWithPutativeVariant=0;

for file = 1:NumFiles(1)
   fid = fopen(strcat(coveragefiles(file).folder,'/',coveragefiles(file).name),'rt');
   InputFile = textscan(fid, '%s%u%f%u', 'MultipleDelimsAsOne',false, 'Delimiter','\t', 'HeaderLines',1);
   fclose(fid);

   % column 1 has coordinate
   % parse out chromosome and coord of the first element in the file
   FirstCoord=strsplit(InputFile{1}{1},':');

   % we are using column 3 for coverages
   % identify the length of region in the current file, in bases
   RegionLength=size(InputFile{3}(:,1));

   % plot the region coverage
   %pid=plot(1:RegionLength(1),InputFile{3}(:,1));
   %saveas(pid,strcat(coveragefiles(file).folder,'/',coveragefiles(file).name,'.png'));

   %find change points
   [TheChange,Residual]=findchangepts(InputFile{3}(:,1),'MaxNumChanges',2,'Statistic','mean');
   NumberOfChanges=size(TheChange);
   
   % if did not find change points, print a message; otherwise proceed with metrics
   if NumberOfChanges(1)~=2 
      fprintf(OutputFile,'%s\t%s\n',coveragefiles(file).name,'inappropriate number of change points found');

   % if both change points are inside the putative area, proceed with metrics; otherwise, print a message
     elseif ( ((str2num(FirstCoord{2})+TheChange(1)-1>PutativeVariantStart) & (str2num(FirstCoord{2})+TheChange(2)-1<PutativeVariantStop)) )

      %find mean and std coverage before, during and after the variant
      BeforeVarCoverage(1)=mean(InputFile{3}(1:TheChange(1),1));
      CoverageInVariant(1)=mean(InputFile{3}(TheChange(1):TheChange(2),1));
      AfterVarCoverage(1)=mean(InputFile{3}(TheChange(2):RegionLength(1),1));

      BeforeVarCoverage(2)=std(InputFile{3}(1:TheChange(1),1));
      CoverageInVariant(2)=std(InputFile{3}(TheChange(1):TheChange(2),1));
      AfterVarCoverage(2)=std(InputFile{3}(TheChange(2):RegionLength(1),1));

      if (CoverageInVariant(1) > BeforeVarCoverage(1))
	fprintf(OutputFile,'%s\t%s\n',coveragefiles(file).name,'higher coverage in var indicated');
	continue;
      end	


      %determine if this is a heterozygous variant
      if ((CoverageInVariant(1)/BeforeVarCoverage(1))>(0.5*(1-BeforeVarCoverage(2)/BeforeVarCoverage(1)))) & ...
         ((CoverageInVariant(1)/BeforeVarCoverage(1))<(0.5*(1+BeforeVarCoverage(2)/BeforeVarCoverage(1))))
         Zygocity='hetero';
      elseif (CoverageInVariant(1)/BeforeVarCoverage(1))<(0.5*(1-BeforeVarCoverage(2)/BeforeVarCoverage(1)))
         Zygocity='homo';
      else
         Zygocity='indeterminate';
      end 

      fprintf(OutputFile,'%s\t%s\t%d\t%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%s\n',...
           coveragefiles(file).name,FirstCoord{1},...
           str2num(FirstCoord{2})+TheChange(1)-1,str2num(FirstCoord{2})+TheChange(2)-1,Residual,...
           BeforeVarCoverage(1),BeforeVarCoverage(2),...
           CoverageInVariant(1),CoverageInVariant(2),AfterVarCoverage(1),AfterVarCoverage(2),Zygocity);
      
      CountIndividualsWithPutativeVariant=CountIndividualsWithPutativeVariant+1;
   else
      fprintf(OutputFile,'%s\t%s\n',coveragefiles(file).name,'change points are outside the putative variant region');      
   end

   clear BeforeVarCoverage CoverageInVariant AfterVarCoverage Zygocity
end
fclose(OutputFile);

% print the total count of individuals with the variant to the Matlab console
end
