%% load and format data
addpath('/Volumes/rds/RDS27534-NICAP/code/BCT/2019_03_03_BCT')

parc='HCP'
module='Yeo'

folders = dir(strcat('/Volumes/rds/RDS27534-NICAP/bids/derivatives/nandi/sMRI/age_structcov/age_windows/*',parc,'*'));
dirFlags = [folders.isdir] & ~strcmp({folders.name},'.') & ~strcmp({folders.name},'..');
folders = folders(dirFlags);

if parc == 'HCP'
    load(strcat('/Volumes/rds/RDS27534-NICAP/code/sMRI/structcov/map_modules/',module,'_HCP_modules.mat'));
    nreg=360;
elseif parc == 'DKT'
    load(strcat('/Volumes/rds/RDS27534-NICAP/code/sMRI/structcov/map_modules/',module,'_DKT_modules.mat'));
    nreg=62;
end

if module == "Yeo"
    nmod=7;
    mod=MYeo;
end

mod=str2double(mod(:,2));
    
%for each window specification, loop through bins and calculate graph metrics, and save output within each window folder
for f = 1:length(folders)  
    window = folders(f);
    cd([window.folder,'/',window.name])
    
    fileList = dir('*_bthresh.mat');
    
    [~, reindex] = sort( str2double( regexp( {fileList.name}, '\d+', 'match', 'once' )))
    fileList = fileList(reindex) ;
    
    %create sc array
    SC = zeros(nreg,nreg,length(fileList));

    for i = 1:length(fileList)
        file = fileList(i);
        sc = load([file.folder,'/',file.name]);
        SC(:,:,i) = sc.SC;
    end

    %set up arrays for output
    den = zeros(1,length(fileList));
    D = zeros(nreg,length(fileList));
    intramoduleD = zeros(nmod,3,length(fileList));
    intermoduleD = zeros(nmod,3,nmod,length(fileList));

    % calculate graph metrics
    for i = 1:length(fileList)
        W = SC(:,:,i);
        A = W > 0;

        %global density
        den(:,i) = density_und(A);
        
        %nodal degree & strength
        D(:,i) = degrees_und(A)';

        %intramodule (within-module) density
        for m = 1:nmod
            Mi = mod==m;
            intramod = zeros(nreg,nreg);
            for n = 1:nreg
                for o = 1:nreg
                    intramod(n,o) = (Mi(n) == 1 & Mi(o) == 1);                     
                end
            end
            intramod = intramod - diag(diag(intramod)); %removing diagonal connections in matrix
            intramoduleD(m,1,i) = nnz(intramod)/2; % number of all possible within module connections
            intramod_A=A.*intramod;
            intramoduleD(m,2,i) = sum(nonzeros(intramod_A))/2;
            intramoduleD(m,3,i) = intramoduleD(m,2,i) / intramoduleD(m,1,i);
        end

        %intermodule (between each pair of module) density 
        for m = 1:nmod
            for n = 1:nmod
                intermod = zeros(nreg,nreg);
                for o = 1:nreg
                    for p = 1:nreg
                        intermod(o,p) = (mod(o) == m & mod(p) == n) | (mod(o) == n & mod(p) == m);                     
                    end
                end
                intermoduleD(m,1,n,i) = nnz(intermod)/2; % number of all possible within module connections
                intermod_A=A.*intermod;
                intermoduleD(m,2,n,i) = sum(nonzeros(intermod_A))/2;
                intermoduleD(m,3,n,i) = intermoduleD(m,2,n,i) / intermoduleD(m,1,n,i);
            end
        end

        for m = 1:nmod
            x=intermoduleD(:,3,m,i);
            x(m,:)=[];
            intermoduleD_mean(m,i)=mean(x);
            intermoduleD(m,:,m,i) = NaN; 
        end
        
    end

    % save graph metrics for GAM modelling in R
    csvwrite([file.folder,'/global_density.csv'],den)
    csvwrite([file.folder,'/nodal_degree.csv'],D)
    csvwrite([file.folder,'/intramoduleD_',module,'.csv'],intramoduleD)
    csvwrite([file.folder,'/intermoduleD_',module,'.csv'],intermoduleD)
end
