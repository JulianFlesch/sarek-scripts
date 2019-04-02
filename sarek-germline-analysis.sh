#PBS -N sarek_mrkh
#PBS -l nodes=1:ppn=2,walltime=100:00:00
#PBS -o LOG

cd #DUMMY-STRING#
nextflow run SciLifeLab/Sarek/main.nf --genome_base /beegfs/work/zxmbm36/reference/hg38 --genome GRCh38 -profile binac --sample samples.tsv --containerPath /beegfs/work/zxmbm36/containers/
touch .main.completed
nextflow run SciLifeLab/Sarek/germlineVC.nf --genome_base /beegfs/work/zxmbm36/reference/hg38 --genome GRCh38 -profile binac --sample Preprocessing/Recalibrated/recalibrated.tsv --containerPath /beegfs/work/zxmbm36/containers/ --tools HaplotypeCaller
touch .vc.completed
nextflow run SciLifeLab/Sarek/annotate.nf -profile binac --containerPath /beegfs/work/zxmbm36/containers/ --annotateTools HaplotypeCaller --genome_base /beegfs/work/zxmbm36/reference/hg38 --tools VEP
touch .annotate.completed
nextflow run SciLifeLab/Sarek/runMultiQC.nf -profile binac --containerPath /beegfs/work/zxmbm36/containers 
touch .multiqc.completed