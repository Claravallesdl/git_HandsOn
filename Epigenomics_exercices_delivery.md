# 4. EN‐TEx ATAC‐seq data: downstream analyses

## Task 4.1: Move to folder ATAC-seq, and create folders to store bigBed data files and peaks analyses files. Make sure the files are organized in a consistent way as done for ChIP-seq.
Initial clarification: The commands used will be found within code boxes, and the obtained results to be shown will be below the corresponding command highlighted in bold.

We are currently working in the directory epigenomics_uvic, which contains the following directories:

    ~/epigenomics_uvic ls
**`ATAC-seq  ChIP-seq  bin  docker  handsOn_images  nextflow  test`**
    
To replicate the same directory structure as the one in the ChIP-seq directory, we navigate to the ATAC-seq directory, create the necessary directories inside it, and check the contents afterward.
 
    cd ATAc-seq

    mkdir -p data annotation analyses
    ls
**`analyses  annotation  data`**

## Task 4.2: Retrieve from a newly generated metadata file ATAC-seq peaks (bigBed narrow, pseudoreplicated peaks, assembly GRCh38) for stomach and sigmoid_colon for the same donor used in the previous sections. Hint: have a look at what we did [here](https://github.com/bborsari/epigenomics_uvic/wiki/3.1.-EN%E2%80%90TEx-ChIP%E2%80%90seq-data:-how-to-navigate-the-portal-and-run-the-chipnf-pipeline#2-how-to-find-experiments-and-programmatically-download-files-from-the-encode-portal-). Make sure your md5sum values coincide with the ones provided by ENCODE.

From the [ENCODE](https://www.encodeproject.org/) website, we select the following:

Data > Epigenomes from four individuals (ENTEx) > Individual number 2 (ENCDO451RUA) > Functional genomics experiments from this donor (view all)

Select the parameters:
* Assay type: DNA accessibility
* Assay title: ATAC-seq
* Biosamples: stomach and sigmoid colon

After selecting the desired data, we obtain the following download link:

    https://www.encodeproject.org/metadata/?replicates.library.biosample.donor.uuid=d370683e-81e7-473f-8475-7716d027849b&status=released&status=submitted&status=in+progress&assay_title=ATAC-seq&biosample_ontology.term_name=sigmoid+colon&biosample_ontology.term_name=stomach&assay_slims=DNA+accessibility&type=Experiment

and we save download it:

../bin/download.metadata.sh https://www.encodeproject.org/metadata/?replicates.library.biosample.donor.uuid=d370683e-81e7-473f-8475-7716d027849b&status=released&status=submitted&status=in+progress&assay_title=ATAC-seq&biosample_ontology.term_name=sigmoid+colon&biosample_ontology.term_name=stomach&assay_slims=DNA+accessibility&type=Experiment

Now, we explore the row names of metadata.tsv to locate the information we want to know (IDs):

    head -1 metadata.tsv 
**`File_accession  File_format     File_type       File_format_type        Output_type     File_assembly   Experiment_accession    Assay   Donor(s)Biosample_term_id       Biosample_term_name     Biosample_type ...`**

After identifying the pertinent columns (ID and file type) within the metadata.tsv file, we proceed to filter and extract the IDs from the bigBed peak files. The provided command performs this task effectively, utilizing various filtering criteria to isolate ATAC-seq experiments with bigBed narrowPeak files generated from pseudoreplicated peaks in the GRCh38 assembly. Subsequently, the command extracts the essential columns (ID, biosample term ID, and file type), sorts them accordingly, and stores the resulting IDs in the analyses/bigBed.peaks.ids.txt file for further analysis.

    grep -F ATAC-seq metadata.tsv |\
    grep -F "bigBed_narrowPeak" |\
    grep -F "pseudoreplicated_peaks" |\
    grep -F "GRCh38" |\
    awk 'BEGIN{FS=OFS="\t"}{print $1, $11, $23}' |\
    sort -k2,2 -k1,1r |\
    sort -k2,2 -u > analyses/bigBed.peaks.ids.txt

Si accedim al directori analyses podem comprovar que el fitxer s'ha creat correctament:

    cd analyses
    ls
**`bigBed.peaks.ids.txt`**

    cat bigBed.peaks.ids.txt

**`ENCFF287UHP     sigmoid_colon`**

**`ENCFF762IFP     stomach`**

We navigate to the "data" directory and create a new subdirectory named "bigBed.files" within it to store the bigBed files that we are going to download.

    cd data
    mkdir bigBed.files

Now we procede to download the corresponding bigBedfiles:

    cut -f1 analyses/bigBed.peaks.ids.txt |\
    while read filename; do
     wget -P data/bigBed.files "https://www.encodeproject.org/files/$filename/@@download/$filename.bigBed"
    done

And we check them:

    cd data
    cd bigBed.files
    ls

**`ENCFF287UHP.bigBed`**

**`ENCFF762IFP.bigBed`**

Check the integrity of the downloaded files:
    for file_type in bigBed; do

      # retrieve original MD5 hash from the metadata
      ../bin/selectRows.sh <(cut -f1 analyses/"$file_type".*.ids.txt) metadata.tsv | cut -f1,46 > data/"$file_type".files/md5sum.txt

      # compute MD5 hash on the downloaded files 
      cat data/"$file_type".files/md5sum.txt |\
      while read filename original_md5sum; do 
        md5sum data/"$file_type".files/"$filename"."$file_type" |\
        awk -v filename="$filename" -v original_md5sum="$original_md5sum" 'BEGIN{FS=" "; OFS="\t"}{print filename, original_md5sum, $1}' 
      done > tmp 
      mv tmp data/"$file_type".files/md5sum.txt

      # make sure there are no files for which original and computed MD5 hashes differ
      awk '$2!=$3' data/"$file_type".files/md5sum.txt

    done

There is no output and it indicates that there are no discrepancies found between the original and computed MD5 hashes for the downloaded files. This suggests that the downloaded files are intact and have not been corrupted during the download process.

## Task 4.3:For each tissue, run an intersection analysis using BEDTools report:

### 1) the number of peaks that intersect promoter regions

    mkdir data/bed.files

Convert bigBed files to BED files with the `bigBedToBed` command:

    cut -f1 analyses/bigBed.peaks.ids.txt |\
    while read filename; do
      bigBedToBed data/bigBed.files/"$filename".bigBed data/bed.files/"$filename".bed
    done

We check the new files:

    cd data/bed.files
    ls

**`ENCFF287UHP.bed  ENCFF762IFP.bed`**

To carry out the retrieval of genes with peaks at the promoter region in each tissue, we will utilize the BED files generated in the ChIP-seq exercises. Specifically, we'll be working with two BED files: "gencode.v24.protein.coding.non.redundant.TSS.bed," containing non-redundant transcription start site (TSS) coordinates of protein-coding genes, and "gencode.v24.protein.coding.gene.body.bed," delineating the gene body coordinates of protein-coding genes. 

    cut -f-2 analyses/bigBed.peaks.ids.txt |\
    while read filename tissue; do 
      bedtools intersect -b annotation/gencode.v24.protein.coding.non.redundant.TSS.bed -a data/bed.files/"$filename".bed -u > 
    analyses/peaks.analysis/genes.with.peaks."$tissue".txt
    done

The the number of peaks that intersect promoter regions is:
    
    cd analyses/peaks.analysis
    ls

**`genes.with.peaks.sigmoid_colon.txt  genes.with.peaks.stomach.txt`**
    
    cd ..
    wc -l peaks.analysis/*

**`47871 genes.with.peaks.sigmoid_colon.txt`**

**`44749 genes.with.peaks.stomach.txt`**

  **`92620 total`**

### 2) the number of peaks that fall outside gene coordinates (whole gene body, not just the promoter regions). Hint: have a look at what we did [here](https://github.com/bborsari/epigenomics_uvic/wiki/3.2.-EN%E2%80%90TEx-ChIP%E2%80%90seq-data:-downstream-analyses#51-genes-with-peaks-of-h3k4me3-in-each-tissue) and [here](https://github.com/bborsari/epigenomics_uvic/wiki/3.2.-EN%E2%80%90TEx-ChIP%E2%80%90seq-data:-downstream-analyses#22-prepare-a-bed-file-with-gene-body-coordinates-of-protein-coding-genes).

    cut -f-2 analyses/bigBed.peaks.ids.txt |\
    while read filename tissue; do 
      bedtools intersect -b annotation/gencode.v24.protein.coding.gene.body.bed -a data/bed.files/"$filename".bed -v > 
    analyses/peaks.analysis/genes.with.peaks."$tissue".outside.txt
    done

    cd analyses/peaks.analysis
    wc -l outside.txt

  **`37035 genes.with.peaks.sigmoid_colon.outside.txt`**

  **`34537 genes.with.peaks.stomach.outside.txt`**

 **` 71572 total`**


# 5. Distal regulatory activity

## Task 5.1: Create a folder regulatory_elements inside epigenomics_uvic. This will be the folder where you store all your subsequent results.
To replicate the same directory structure as the one in the ChIP-seq and ATAC-seq directory, we create a new folder names `regulatory_elements`  , create the necessary directories inside it. We also copy the metadata that we obtain to ChIP-seq exercices.

    mkdir regulatory_elements
    cd regulatory_elements
    mkdir analyses
    mkdir data
    cd data
    mkdir bigBed.files 
    mkdir bed.files

    cp ChIP-seq/metadata.tsv regulatory_elements


## Task 5.2: Distal regulatory regions are usually found to be flanked by both H3K27ac and H3K4me1. From your starting catalogue of open regions in each tissue, select those that overlap peaks of H3K27ac AND H3K4me1 in the corresponding tissue. You will get a list of candidate distal regulatory elements for each tissue. How many are they?

We proceed by first creating the peaks ID files for the H3K27ac and H3K4me1 datasets, filtering the metadata.tsv file and then downloading the corresponding BigBed files.

**H3K27ac:**

    grep -F H3K27ac metadata.tsv |\
    grep -F "bigBed_narrowPeak" |\
    grep -F "pseudoreplicated_peaks" |\
    grep -F "GRCh38" |\
    awk 'BEGIN{FS=OFS="\t"}{print $1, $11, $23}' |\
    sort -k2,2 -k1,1r |\
    sort -k2,2 -u > analyses/bigBed.H3K27ac_peaks.ids.txt

    cut -f1 analyses/bigBed.H3K27ac_peaks.ids.txt |\
    while read filename; do
      wget -P data/bigBed.files "https://www.encodeproject.org/files/$filename/@@download/$filename.bigBed"
    done

**H3K4me1:**

    grep -F H3K4me1 metadata.tsv |\
    grep -F "bigBed_narrowPeak" |\
    grep -F "pseudoreplicated_peaks" |\
    grep -F "GRCh38" |\
    awk 'BEGIN{FS=OFS="\t"}{print $1, $11, $23}' |\
    sort -k2,2 -k1,1r |\
    sort -k2,2 -u > analyses/bigBed.H3K4me1_peaks.ids.txt

    cut -f1 analyses/bigBed.H3K4me1_peaks.ids.txt |\
    while read filename; do
      wget -P data/bigBed.files "https://www.encodeproject.org/files/$filename/@@download/$filename.bigBed"
    done

Now we can check which IDs corresponds to H3K27ac and H3K4me1 peaks:

    cd analyses
    ls
    cat *


**`ENCFF872UHN     sigmoid_colon   H3K27ac-human`**

**`ENCFF977LBD     stomach H3K27ac-human`**

**`ENCFF724ZOF     sigmoid_colon   H3K4me1-human`**

**`ENCFF844XRN     stomach H3K4me1-human`**

Check the integrity of the downloaded files:

    for file_type in bigBed; do

      # retrieve original MD5 hash from the metadata
      ../bin/selectRows.sh <(cut -f1 analyses/"$file_type".*.ids.txt) metadata.tsv | cut -f1,46 > data/"$file_type".files/md5sum.txt

      # compute MD5 hash on the downloaded files 
      cat data/"$file_type".files/md5sum.txt |\
      while read filename original_md5sum; do 
        md5sum data/"$file_type".files/"$filename"."$file_type" |\
        awk -v filename="$filename" -v original_md5sum="$original_md5sum" 'BEGIN{FS=" "; OFS="\t"}{print filename, original_md5sum, $1}' 
      done > tmp 
      mv tmp data/"$file_type".files/md5sum.txt

      # make sure there are no files for which original and computed MD5 hashes differ
      awk '$2!=$3' data/"$file_type".files/md5sum.txt

    done

Convert bigBed files to BED files with the `bigBedToBed` command:

**H3K27ac:**

    cut -f1 analyses/bigBed.H3K27ac_peaks.ids.txt |\
    while read filename; do
      bigBedToBed data/bigBed.files/"$filename".bigBed data/bed.files/"$filename".bed
    done

**H3K4me1:**

    cut -f1 analyses/bigBed.H3K4me1_peaks.ids.txt |\
    while read filename; do
      bigBedToBed data/bigBed.files/"$filename".bigBed data/bed.files/"$filename".bed
    done

And we check the files:

    cd data/bed.files
    ls

**`ENCFF724ZOF.bigBed  ENCFF844XRN.bigBed  ENCFF872UHN.bigBed  ENCFF977LBD.bigBed  md5sum.txt`**

To select those regions that overlap with peaks of both H3K27ac and H3K4me1 in the corresponding tissue from the initial catalog of open regions in each tissue (genes.with.peaks.sigmoid_colon.outside.txt and genes.with.peaks.stomach.outside.txt), so that the intersection is carried out simultaneously for H3K27ac and H3K4me1, we proceed as follows:

**List of files for sigmoid_colon and stomach:**

    sigmoid_colon_files=(regulatory_elements/data/bed.files/ENCFF872UHN.bed regulatory_elements/data/bed.files/ENCFF724ZOF.bed)
    stomach_files=(regulatory_elements/data/bed.files/ENCFF977LBD.bed regulatory_elements/data/bed.files/ENCFF844XRN.bed)

**Intersection for sigmoid_colon:**

    bedtools intersect -a "ATAC-seq/analyses/peaks.analysis/genes.with.peaks.sigmoid_colon.outside.txt" -b "${sigmoid_colon_files[@]}" -u > 
    "regulatory_elements/intersection/common_sigmoid_colon.bed"

**Perform the intersection for stomach:**

    bedtools intersect -a "ATAC-seq/analyses/peaks.analysis/genes.with.peaks.stomach.outside.txt" -b "${stomach_files[@]}" -u > 
    "regulatory_elements/intersection/common_stomach.bed"

**Count:**

    wc -l regulatory_elements/intersection/*


**`23130 regulatory_elements/intersection/common_sigmoid_colon.bed`**

**`18028 regulatory_elements/intersection/common_stomach.bed`**

**`41158 total`**

## Task 5.3: Focus on regulatory elements that are located on chromosome 1 (hint: to parse a file based on the value of a specific column, have a look at what we did [here](https://github.com/bborsari/epigenomics_uvic/wiki/3.1.-EN%E2%80%90TEx-ChIP%E2%80%90seq-data:-how-to-navigate-the-portal-and-run-the-chipnf-pipeline#22-how-to-parse-the-metadata-file-to-retrieve-data-of-a-specific-experiment)), and generate a file regulatory.elements.starts.tsv that contains the name of the regulatory region (i.e. the name of the original ATAC-seq peak) and the start (5') coordinate of the region.

    cd regulatory_elements
    cd intersection
    ls

**`common_sigmoid_colon.bed  common_stomach.bed`**

    head -1 common_sigmoid_colon.bed

**`chr1    778339  779193  Peak_175893     36      .       2.20588 3.66321 1.70835 781`**

Filter the regulatory regions located on chromosome 1 and extract the name of the regulatory region and its start coordinate for sigmoid_colon

    awk '$1 == "chr1" {print $4 "\t" $2}' intersection/common_sigmoid_colon.bed > starts.tsv/sigmoid_colon_starts.tsv

Filter the regulatory regions located on chromosome 1 and extract the name of the regulatory region and its start coordinate for stomach

    awk '$1 == "chr1" {print $4 "\t" $2}' intersection/common_stomach.bed > starts.tsv/stomach_starts.tsv

Combine the results into a single file using cat

    cat starts.tsv/sigmoid_colon_starts.tsv starts.tsv/stomach_starts.tsv > starts.tsv/regulatory.elements.starts.tsv

Finally we check the regulatory.elements.starts.tsv file

    head -3 starts.tsv/regulatory.elements.starts.tsv

**`Peak_175893     778339`**

**`Peak_203412     778339`**

**`Peak_2485       778339`**


## Task 5.4: Focus on protein-coding genes located on chromosome 1. From the BED file of gene body coordinates that you generated [here](https://github.com/bborsari/epigenomics_uvic/wiki/3.2.-EN%E2%80%90TEx-ChIP%E2%80%90seq-data:-downstream-analyses#22-prepare-a-bed-file-with-gene-body-coordinates-of-protein-coding-genes), prepare a tab-separated file called gene.starts.tsv which will store the name of the gene in the first column, and the start coordinate of the gene on the second column (REMEMBER: for genes located on the minus strand, the start coordinate will be at the 3'). Use the command below as a starting point:

    awk 'BEGIN{FS=OFS="\t"}{if ($6=="+"){start=$2} else {start=$3}; print $4, start}'

In our case, the correct commnd to focus on chr1 would be:

    awk 'BEGIN{FS=OFS="\t"} $1=="chr1" {if ($6=="+") {start=$2} else {start=$3}; print $4, start}' \
    ChIP-seq/annotation/gencode.v24.protein.coding.gene.body.bed > regulatory_elements/starts.tsv/gene.starts.tsv

And we obtain the number of genes in the list:

    cd regulatory_elements/starts.tsv
    wc -l gene.starts.tsv

**`2047 gene.starts.tsv`**



## Task 5.5: Download or copy [this python script](https://public-docs.crg.es/rguigo/Data/bborsari/UVIC/epigenomics_course/get.distance.py) inside the epigenomics_uvic/bin folder.Complete the python script so that for a given coordinate --start the script returns the closest gene, the start of the gene and the distance of the regulatory element. To make sure your script is working fine, run the following command: `python ../bin/get.distance.py --input gene.starts.tsv --start 980000` You should be getting this result: `ENSG00000187642.9	982093 2093`

    nano get.distance.py

    for line in open_input.readlines():  # For each line in the input file
        gene, start = line.strip().split('\t')  # Split the line into two columns based on a tab
        start = int(start)  # Convert start coordinate to integer

        # Compute the absolute value of the difference between position and enhancer_start
        distance = abs(start - enhancer_start)

        # If this absolute value is lower than the current minimum distance x
        if distance < x:
            # Update the current minimum distance x
            x = distance

            # Save gene as selectedGene
            selected_gene = gene

            # Save start coordinate as selectedGeneStart
            selected_gene_start = start

        # Print the closest gene, its start coordinate, and the distance from the regulatory element
        print("\t".join([selected_gene, str(selected_gene_start), str(x)]))

We move the get.distance.py to bin directory:

    mv get.distance.py bin/
    cd ..

And we run the command:

    python bin/get.distance.py --input regulatory_elements/starts.tsv/gene.starts.tsv --start 980000

**`ENSG00000187642.9       982093  2093`**

## Task 5.6. For each regulatory element contained in the file regulatory.elements.starts.tsv, retrieve the closest gene and the distance to the closest gene using the python script you created above. Use the command below as a starting point:
    
    cat regulatory.elements.starts.tsv | while read element start; do 
       python ../bin/get.distance.py ... # to be completed by you; 
    done > regulatoryElements.genes.distances.tsv

In our case:

    cat regulatory_elements/starts.tsv/regulatory.elements.starts.tsv | while read element start; do
       python bin/get.distance.py -i regulatory_elements/starts.tsv/gene.starts.tsv --start "$start";
    done > regulatory_elements/regulatoryElements.genes.distances.tsv

    cd regulatory_elements
    ls

**`analyses  data  intersection  metadata.tsv  regulatoryElements.genes.distances.tsv  starts.tsv`**

    wc -l  regulatoryElements.genes.distances.tsv

**`4488 regulatoryElements.genes.distances.tsv`**

    
## Task 5.7: Use R to compute the mean and the median of the distances stored in regulatoryElements.genes.distances.tsv.
    
    R
    
    # Read the distances file
    data <- read.csv("regulatoryElements.genes.distances.tsv", header = FALSE, sep = "\t")

    # Calculate the mean and median of the third column (distances)
    mean_distance <- mean(data$V3)
    median_distance <- median(data$V3)

    # Create a dataframe with the results
    result <- data.frame(
      Metric = c("Mean distance", "Median distance"),
      Value = c(mean_distance, median_distance)
    )

    # Write the results to a text file
    write.table(result, "mean_median_results.txt", sep = "\t", row.names = FALSE, col.names = FALSE)

    q()

Out of R:

    cat mean_median_results.txt

**`"Mean distance" 65868.7130124777`**

**`"Median distance"       31878.5`**
