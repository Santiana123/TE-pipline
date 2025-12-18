#!/bin/bash
#PBS -j oe 
#PBS -q FAFU3
#PBS -V 
#PBS -l nodes=1:ppn=40

cd $PBS_O_WORKDIR

# Input parameters
fa=/public/home/yuejingjing/tyh/1.project/1.emb506/10.15-4YYh/7.busco/15-4YYh.hapA.Yh.fasta
genome_size=334768680
name=15-4YYh.hapA
rp_species=Carica_papaya
cpu=40

# Activate EDTA environment
source activate EDTA

# Create EDTA output directory if not exists
if [ -d 00.EDTA ]
then
      echo '00.EDTA existed'
else
      mkdir ./00.EDTA
      echo '00.EDTA new create'
fi

cd ./00.EDTA
# Create symbolic link to genome file
ln -s ${fa} genome.fasta
#########EDTA annotation
EDTA.pl --genome genome.fasta --species others --sensitive 0 --anno 1 --overwrite 0 \
        --threads ${cpu} --curatedlib /public/home/yuejingjing/tyh/db/MIPS-database/mipsREdat_9.3p_ALL_library2.fix.fa
cd ../

#########RepeatMasker annotation
# Create RepeatMasker output directory if not exists
if  [ -d 01.rpanno ]
then
      echo '01.rpanno existed'
else
      mkdir ./01.rpanno
      echo '01.rpanno new create'
fi

cd ./01.rpanno

# Copy genome file
cp ${fa} ./${name}.fa

### Build database for RepeatModeler
BuildDatabase -name ${name} ${name}.fa

### De novo repeat annotation using RepeatModeler
RepeatModeler -database ./${name} -pa ${cpu} -engine ncbi 1>repeatmodeler.log 2>&1

### Annotation using existing library with RepeatMasker
RepeatMasker -e rmblast -pa ${cpu} -gff -a -species ${rp_species} -dir ./rmout_${name} ${name}.fa 1>repeatmasker_A.log 2>&1

### Combine TE libraries (EDTA and RepeatModeler) for genome annotation
cat ../00.EDTA/genome.fasta.mod.EDTA.raw/genome.fasta.mod.EDTA.intact.fa ./${name}-families.fa > holistic.TE-lib.fa
RepeatMasker -e rmblast -pa ${cpu} -gff -a -html -poly -lib holistic.TE-lib.fa -dir ./rmout_denovo_${name} rmout_${name}/${name}.fa.masked 1>repeatmasker_denovo.log 2>&1

# Combine annotation results
if [ ! -d rmout_combine_${name} ]
then
      mkdir  rmout_combine_${name}
fi

perl /public/software/RepeatMasker/util/combineRMFiles.pl rmout_${name}/${name}.fa  ./rmout_denovo_${name}/${name}.fa.masked rmout_combine_${name}/${name}
# Convert to GFF format
/public/software/RepeatMasker/util/rmOutToGFF3.pl rmout_combine_${name}/${name}.out > rmout_combine_${name}/${name}.out.gff
# Softmask genome (repeat sequences in lowercase)
/public/software/RepeatMasker/util/maskFile.pl -fasta ${name}.fa -annotations rmout_combine_${name}/${name}.out -softmask
mv ${name}.fa.masked  rmout_combine_${name}/${name}.sm.fa
# Mask genome (repeat sequences replaced with N)
/public/software/RepeatMasker/util/maskFile.pl -fasta ${name}.fa -annotations rmout_combine_${name}/${name}.out
mv ${name}.fa.masked  rmout_combine_${name}/${name}.rm.fa

# For subsequent protein sequence annotation: remove annotations of simple repeats and generate GFF and softmasked fasta
egrep -v "Simple|Satellite|Low_complexity" rmout_combine_${name}/${name}.out > rmout_combine_${name}/${name}.filtered.out
# Convert to GFF
/public/software/RepeatMasker/util/rmOutToGFF3.pl rmout_combine_${name}/${name}.filtered.out > rmout_combine_${name}/${name}.filtered.out.gff
# Re-perform softmasking
/public/software/RepeatMasker/util/maskFile.pl -fasta ${name}.fa -annotations rmout_combine_${name}/${name}.filtered.out -softmask
mv ${name}.fa.masked rmout_combine_${name}/${name}.filtered.sm.fa

########Re-annotate unknown TEs
cd ./rmout_combine_${name}

# Create directory for unknown TEs if not exists
if [ ! -d ./Unknown_TE ]
then
      mkdir ./Unknown_TE
fi
cd ./Unknown_TE
# Process format of filtered out file (remove simple repeat annotations)
grep -v '*' ../${name}.filtered.out > ../${name}.filtered.f2.out
sed 's/^[ \t]*//' ../${name}.filtered.f2.out | sed 's/ \+/\t/g' > test.out

# Convert out file to BED format
perl /public/home/yuejingjing/tyh/protocol/rpanno-stat/rmout_to_bed.pl -input test.out -bed test.bed
# Extract unknown TE sequences
bedtools getfasta -fi ../../${name}.fa -bed test.bed -s > unknown.fa

# Generate tandem.all.bed to check tandem repeats
perl /public/home/yuejingjing/tyh/protocol/rpanno-stat/tandemRepeatFinder.pl -i ../../${name}.fa

# Activate DeepTE environment
source ~/.bashrc
source activate DeepTE

# Annotate unknown TEs using DeepTE
python /public/home/yuejingjing/biosofts/DeepTE-master/DeepTE.py -d working_dir \
-i unknown.fa -sp P -m_dir /public/home/yuejingjing/tyh/db/Plants_model

# Modify output format
sed  's/(+)__/#/g' opt_DeepTE.fasta > opt_DeepTE.fix.fa
# Rename transposons for TEclassify.pl recognition
python3  /public/home/yuejingjing/tyh/protocol/rpanno-stat/merge_deepTE.py -i opt_DeepTE.fix.fa -o DeepTE.fa
python /public/home/yuejingjing/tyh/protocol/rpanno-stat/update_annotation.py -inf DeepTE.fa -i test.out -o update.out

# Generate required file for TEclassify-2.0.pl
perl Scripts/all-types-stat.pl -input update.out -output tmp.TE.out

# Extract remaining unknowns after DeepTE for further statistics
perl /public/home/yuejingjing/tyh/protocol/rpanno-stat/rmout_to_bed.pl -input update.out -bed final-unknown.bed
bedtools getfasta -fi ../../${name}.fa  -bed final-unknown.bed -s > final.unknown.fa

# Run TEclassify-2.0.pl
perl /public/home/yuejingjing/tyh/protocol/rpanno-stat/TEclassify-2.0.pl -i update.out -u final.unknown.fa -g ${genome_size} > results.txt
# Convert update.out to BED format for subsequent methylation analysis
source ~/.bashrc
source activate EDTA
~/mambaforge/envs/EDTA/share/RepeatMasker/util/RM2Bed.py update.out
sort -k1,1V -s update_rm.bed > update_rm.sort.bed
awk 'BEGIN{OFS="\t"} {tmp=$10; $10=$11; $11=tmp; print}' update.out > update.out2
cd ../../../
ln -s 01.rpanno/rmout_combine_${name}/Unknown_TE/results.txt
/public/software/RepeatMasker/util/rmOutToGFF3.pl 01.rpanno/rmout_combine_${name}/Unknown_TE/update.out2 > update.out.gff

# Create TEsorter directory and run TEsorter
mkdir TEsorter
cd TEsorter
source ~/.bashrc
conda activate TEsorter
TEsorter ${fa} -genome -p ${cpu} -db rexdb-plant -prob 0.9 -cov 30 -eval 1e-5 -score 1 1>2 2>TEsorter.log
ln -s ../update.out.gff

# Process GFF files
sed -E 's/(([^\t]*\t){8}).*Target=([^ ;]+).*/\1\3/' update.out.gff > update.out.gff3
sed -E 's/(([^\t]*\t){8}).*\|([^:]+):.*/\1\3/' *.rexdb-plant.dom.gff3 > TEsorter.dom.gff3
bedtools intersect -a update.out.gff3 -b TEsorter.dom.gff3 -wa -wb | awk -F'\t' 'BEGIN{OFS="\t"} {map[$1","$4","$5","$7]=$18} END{while((getline < "update.out.gff3")>0){split($0, old, "\t"); key=old[1]","old[4]","old[5]","old[7]; if(key in map){$9=map[key]} print}}' > updated_combined.gff

# Standardize annotation names
awk 'BEGIN { FS=OFS="\t" }; $9 == "LTR/Gypsy" { $9 = "Class_I/LTR/Ty3_gypsy" }; $9 == "DNA/DTM" { $9 = "Class_II/subclass_1/TIR/MuDR/Mutator" }; $9 == "MITE/DTM" { $9 = "Class_II/subclass_1/TIR/MuDR/Mutator" }; $9 == "DNA/PIF-Harbinger" { $9 = "Class_II/subclass_1/TIR/PIF/Harbinger" }; $9 == "DNA/DTH" { $9 = "Class_II/subclass_1/TIR/PIF/Harbinger" }; $9 == "MITE/DTH" { $9 = "Class_II/subclass_1/TIR/PIF/Harbinger" }; $9 == "DNA/DTC" { $9 = "Class_II/subclass_1/TIR/EnSpm/CACTA" }; $9 == "MITE/DTC" { $9 = "Class_II/subclass_1/TIR/EnSpm/CACTA" }; $9 == "LTR/Copia" { $9 = "Class_I/LTR/Ty1_copia" }; $9 == "DNA" { $9 = "Class_II/subclass_1/TIR" }; $9 == "hAT" { $9 = "Class_II/subclass_1/TIR/hAT" }; $9 == "LINE/L1" { $9 = "Class_I/LINE" }; $9 == "ClassI_nLTR_LINE_I" { $9 = "Class_I/LINE" }; $9 == "ClassI_nLTR_LINE" { $9 = "Class_I/LINE" }; $9 == "DNA/Helitron" { $9 = "Class_II/subclass_2/Helitron" }; $9 == "RC/Helitron" { $9 = "Class_II/subclass_2/Helitron" }; $9 == "LTR" { $9 = "Class_I/LTR" }; $9 == "DNA/DTA" { $9 = "Class_II/subclass_1/TIR/hAT" }; $9 == "MITE/DTA" { $9 = "Class_II/subclass_1/TIR/hAT" }; $9 == "DNA/DTT" { $9 = "Class_II/subclass_1/TIR/Tc1/Mariner" }; $9 == "MITE/DTT" { $9 = "Class_II/subclass_1/TIR/Tc1/Mariner" }; $9 == "LTR/unknown" { $9 = "Class_I/LTR" }; $9 == "SINE" { $9 = "Class_I/SINE" }; $9 == "DNA/MULE-MuDR" { $9 = "Class_II/subclass_1/TIR/MuDR/Mutator" }; $9 == "tRNA" { next }; $9 == "Retroposon" { next }; $9 == "rRNA" { next }; $9 == "snRNA" { next }; { print }' updated_combined.gff > formatted.gff

echo "TE annotation completed"
