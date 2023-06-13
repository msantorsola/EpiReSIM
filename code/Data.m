classdef Data


    properties
        name          % Name of the Data.
        class         % The label for real data
        pts           % List of data.
    end
    methods
        function D = Data(path,snpnum)
            [~, D.name, ~] = fileparts(path);
            fid = load(path);
            rng(0);
	    NUM  = round(rand(1,1)*(size(fid.pts,2)-snpnum-1));
            D.pts = fid.pts(:,NUM:NUM+snpnum-1);
            D.class =zeros(1,size(D.pts,1));
            for i = 1:size(D.pts,1)
                D.class(1,i)= fid.SampleInfo{i,5};%;%class(1,i);%{i,5};%class(1,i)
                i = i+1;
%                 D.class = fid.class;  
            end
            filename = 'results/input_snp_indices.txt'; % Specify the filename with the appropriate extension
	    fileID = fopen(filename, 'w');
            fprintf(fileID, 'start\tend\n');
            fprintf(fileID, '%d\t%d\n', NUM, NUM + snpnum - 1);
            fclose(fileID);

	end
    end
end


