function Main(input_file_name,output_file_name,WITH_MISMATCH)

if nargin<3
   WITH_MISMATCH = 0;
end

NUM_BASES_ADDED = 19;

T = readtable(input_file_name);

all_fasta_seq = T.fasta;
[unique_fasta_seq, i_all_fasta, i_unique_fasta] = unique(all_fasta_seq);

[homologies,mismatch_amount] = find_homologies(unique_fasta_seq,NUM_BASES_ADDED,WITH_MISMATCH);

%generate a list corresponding to the original list
all_homplogies = homologies(i_unique_fasta);
all_mismatch_amount = mismatch_amount(i_unique_fasta);

%go over the homologies and calculate the length 
del_start = T.start;
del_end = T.xEnd;
del_len = del_end-del_start+1;
MH_len = [];
status = [];
for i=1:length(all_homplogies)
   MH_len(i) = length(all_homplogies{i});
end
MH_len = MH_len';



%change names
names = T.Gene;
for i=1:length(names)
    cur_name = names{i};
    ind = strfind(cur_name,'ENST');
    if ~isempty(ind)
        names{i} = cur_name(1:ind-2);
    end
end

all_chr = T.chr;
all_start = T.start;
all_end = T.xEnd;
all_start_minus_19 = T.startMinus19;
all_end_plus_20 = T.endPlus20;
all_fasta_coordinates = T.fasta_coordinates;
all_fasta = T.fasta;


%generate the output file
fid = fopen(output_file_name,'w');

%insert header
fprintf(fid,'%s\t%s\t%s\t%s\t%s\t%s\t%s\t','Gene','chr','start','end','start_minus_19','end_plus_20','fasta_coordinates');
fprintf(fid,'%s\t%s\t%s\t%s\t%s\t%s\t%s\n','fasta','MH','MISMATCH_AMOUNT','MH_len','Del_len');
%insert vals
for i=1:length(all_fasta_seq)
    fprintf(fid,'%s\t%s\t%d\t%d\t%d\t%d\t%s\t',...
        names{i},all_chr{i},all_start(i),all_end(i),all_start_minus_19(i),all_end_plus_20(i),all_fasta_coordinates{i});
    fprintf(fid,'%s\t%s\t%d\t%d\t%d\t%f\t%d\n',all_fasta_seq{i},all_homplogies{i},all_mismatch_amount(i),MH_len(i),del_len(i));
end
fclose all;






