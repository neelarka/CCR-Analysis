clear
reads = fastaread('C:\Users\eajoseph\Desktop\UNCG\Projects\Cas9\emulsion\seq_list.fasta')
check = 'CTCCCATCACATCAATGGGCTTTGGAAAGG'
fullcheck = 'GACATCACCTCCCACAACGACGAAAAGAGGAGGAAGGGCCTGAGTCCGAGCAGAAGAAGAAGGG'
primer = 'GACATCACCTCCCACAACGA'
c = 0;
r1=0;
r2=0;
for x = 1:length(reads)
    s1=reads(x).Sequence;
    s2=seqrcomplement(reads(x).Sequence);
    L = length(s1);
    if L>250
      
      sc=swalign(s1(1:100),check);
       
      sr=swalign(s2(1:100),check);
      if (sr>70)&&(sr>sc)
          s1=s2;
          sc=sr;
      end
          
      if sc>70
             c=c+1;
          sc2=swalign(s1(1:100),fullcheck);
          scores(c)=sc2;
        if sc2>100
            r1=r1+1;
            L1(r1)=L;
        else
            if contains(s1(1:100),primer)==0
               r2=r2+1;
            L2(r2)=L;
            reads(x).Sequence=s1;
            bigcut(r2)=reads(x);
            end
        end
      end
      
    end
end
fastawrite('C:\Users\eajoseph\Desktop\UNCG\Projects\Cas9\emulsion\bigcut.fasta',bigcut)
