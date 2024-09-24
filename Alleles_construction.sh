for i in *.fastq ; do NanoFilt -q 15 ${i} >${i%.fastq}.q15.fastq.tmp; done
for i in *.fastq.tmp;
do
prefix=`echo ${i}|sed -e 's/20240820.//'|sed -e 's/.q15.fastq.tmp//'`
vsearch --orient ${i} --db OG0009774.fasta --fastqout orient/${prefix}.fastq
done
for i in *.fastq.tmp;
do
prefix=`echo ${i}|sed -e 's/20240820.//'|sed -e 's/.q15.fastq.tmp//'`
vsearch --cluster_size orient/${prefix}.fastq --id 0.99 --iddef 4 --consout orient/${prefix}_consensus.fna
grep '>' orient/${prefix}_consensus.fna | sed -e 's/=/ /g' |sort -k3,3n|sed -e 's/ /=/g'  |awk '{if($2>=50)print $0}'|sed -e 's/>//' >orient/${prefix}_orient.list
seqtk subseq orient/${prefix}_consensus.fna orient/${prefix}_orient.list > orient/${prefix}_consensus_size50.fasta
done
for i in *.fastq.tmp;
do
prefix=`echo ${i}|sed -e 's/20240820.//'|sed -e 's/.q15.fastq.tmp//'`
sed -i "s/>.*.;/>${prefix}_/" orient/${prefix}_consensus_size50.fasta
done
