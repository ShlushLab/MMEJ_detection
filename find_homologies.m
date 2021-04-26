function [homologies,mismatch_amount] = find_homologies(unique_fasta_seq,NUM_BASES_ADDED,WITH_MISMATCH)

homologies = [];
mismatch_amount = [];
seq_amount = length(unique_fasta_seq);


if WITH_MISMATCH==0    
    for i=1:seq_amount

        current_seq = unique_fasta_seq{i};        
        %suppose the homolog is on the right side of the deletion (forword)
        cutted_seq = current_seq(NUM_BASES_ADDED+1 : end-NUM_BASES_ADDED-1);
        added_seq = current_seq(end-NUM_BASES_ADDED:end);
        MISMATCH = 0;
        base = 1;
        while MISMATCH<=1 & base<=length(cutted_seq) & base<=NUM_BASES_ADDED+1
            if added_seq(base) == cutted_seq(base)
                base = base+1;
            else
                MISMATCH = MISMATCH+1;                
            end
        end
        if base==1
            forward_homolog = [];
            mismatch_amount_for = -1000;
        else
            forward_homolog = cutted_seq(1:base-1);           
        end
        
        %suppose the homolog is on the left side of the deletion (reverse)
        cutted_seq = current_seq(NUM_BASES_ADDED+1 : end-NUM_BASES_ADDED-1);
        added_seq = current_seq(1:NUM_BASES_ADDED);
        MISMATCH = 0;
        base = 1;
        while MISMATCH<=1 & base<=length(cutted_seq) & base<=NUM_BASES_ADDED
            if added_seq(end-base+1) == cutted_seq(end-base+1)
                base = base+1;
            else
               MISMATCH = MISMATCH+1;               
            end
        end
        if base==1 
            reverse_homolog = [];
            mismatch_amount_rev = -1000;
        else
            reverse_homolog = cutted_seq(end-base+2:end);           
        end
        
        
        if isempty(reverse_homolog)
            homologies{i} = forward_homolog;
        elseif isempty(forward_homolog)
            homologies{i} = reverse_homolog;
        elseif length(forward_homolog) >= length(reverse_homolog)
            homologies{i} = forward_homolog;
        else
            homologies{i} = reverse_homolog;
        end
        mismatch_amount(i)  = 0;
    end
end






if WITH_MISMATCH==1  
    for i=1:seq_amount
        current_seq = unique_fasta_seq{i};        
        %suppose the homolog is on the right side of the deletion (forword)
        cutted_seq = current_seq(NUM_BASES_ADDED+1 : end-NUM_BASES_ADDED-1);
        added_seq = current_seq(end-NUM_BASES_ADDED:end);
        MISMATCH = 0;
        base = 1;
        while MISMATCH<=1 & base<=length(cutted_seq) & base<=NUM_BASES_ADDED+1
            if added_seq(base) == cutted_seq(base)
                base = base+1;
            else
                if MISMATCH<1
                    base = base+1;
                    MISMATCH = MISMATCH+1;
                elseif MISMATCH==1
                    MISMATCH = MISMATCH+1;
                end
            end
        end
        if base==1 | (base==2 & MISMATCH==1)
            forward_homolog = [];
            mismatch_amount_for = -1000;
        else
            forward_homolog = cutted_seq(1:base-1);
            if MISMATCH==0
                mismatch_amount_for = 0;
            else
                mismatch_amount_for = MISMATCH-1;
            end
        end
        
       
        %suppose the homolog is on the left side of the deletion (reverse)
        cutted_seq = current_seq(NUM_BASES_ADDED+1 : end-NUM_BASES_ADDED-1);
        added_seq = current_seq(1:NUM_BASES_ADDED);
        MISMATCH = 0;
        base = 1;
        while MISMATCH<=1 & base<=length(cutted_seq) & base<=NUM_BASES_ADDED
            if added_seq(end-base+1) == cutted_seq(end-base+1)
                base = base+1;
            else
                if MISMATCH<1
                    base = base+1;
                    MISMATCH = MISMATCH+1;
                elseif MISMATCH==1
                    MISMATCH = MISMATCH+1;
                end
            end
        end
        if base==1 | (base==2 & MISMATCH==1)
            reverse_homolog = [];
            mismatch_amount_rev = -1000;
        else
            reverse_homolog = cutted_seq(end-base+2:end);
            if MISMATCH==0
                mismatch_amount_rev = 0;
            else
                mismatch_amount_rev = MISMATCH-1;
            end
        end
        
        
        %end
        if isempty(reverse_homolog)
            homologies{i} = forward_homolog;
            mismatch_amount(i) = mismatch_amount_for;
        elseif isempty(forward_homolog)
            homologies{i} = reverse_homolog;
            mismatch_amount(i) = mismatch_amount_rev;
        elseif length(forward_homolog) >= length(reverse_homolog)
            homologies{i} = forward_homolog;
            mismatch_amount(i) = mismatch_amount_for;
        else
            homologies{i} = reverse_homolog;
            mismatch_amount(i) = mismatch_amount_rev;
        end
    end
end
