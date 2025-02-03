# Rscripts
store Rscripts I made
# perform_edgeR.r
This function takes count data, group information, and a list of comparisons as input,
and performs edgeR analysis to identify differentially expressed genes for each comparison.

Usage:
1. Prepare your count data as a CSV file (e.g., counts.csv) where rows are genes and columns are samples.
   The first column should be gene names and the first row should be sample names.
   Example:
      gene,sample1,sample2,sample3,sample4,sample5,sample6
      gene1,100,120,50,60,200,220
      gene2,50,60,20,30,100,110
      ...

2. Prepare your group information as a CSV file (e.g., group.csv) with one column indicating the group for each sample.
   The order of groups in this file must correspond to the column order in your count data.
   Example:
      control
      control
      treat1
      treat1
      treat2
      treat2

3. Prepare a comparisons list as a CSV file (e.g., comparisons.csv) with two columns: 'target' and 'reference' groups for each comparison.
   Each row defines a comparison you want to perform (target vs reference).
   Example:
      treat1,control
      treat2,control
      treat2,treat1

4. Run this script from the command line using Docker (or Rscript directly if you have dependencies installed):
   docker run -v $(pwd):/data your-dockerhub-username/edger-generalized \
       Rscript script.R \
       -c /data/counts.csv \
       -g /data/group.csv \
       -v /data/comparisons.csv \
       -o /data/output

   Replace 'your-dockerhub-username/edger-generalized' with your Docker Hub image name.
   Make sure 'counts.csv', 'group.csv', and 'comparisons.csv' are in the same directory where you run the command.
   Results will be saved in the 'output' directory within the same directory.

Arguments:
  countdata: Count data frame (rows: gene names, columns: sample names).
  group:     Factor vector of sample groups (corresponding to the column order of countdata).
  comparisons: A list of comparisons. Each element is a vector of two group names to compare (e.g., list(c("groupA", "groupB"), c("groupA", "groupC"))).


