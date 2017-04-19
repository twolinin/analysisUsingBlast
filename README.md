# analysisUsingBlast

makeblastdb -in reference.fasta -dbtype nucl -hash_index -parse_seqids -out reference.db

blastn -task blastn -query analysisData -db reference.db -evalue 0.0001 -outfmt 6 > out.blast

g++ analysisUsingBlast.cpp -o analysisUsingBlast

./analysisUsingBlast out.blast
