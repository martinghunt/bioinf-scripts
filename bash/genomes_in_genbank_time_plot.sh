#!/usr/bin/env bash
set -e

if [ $# -ne 1 ]
then
    echo "usage: $0 output_directory

Makes a plot of number of assemblies in genbank over time."
    exit
fi


mkdir $1
cd $1

wget ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_genbank.txt

awk -F"\t" '$15~/\// {print $15}' assembly_summary_genbank.txt | awk -F"/" '{print $1"-"$2"-01"}' | sort | uniq -c | awk 'BEGIN{print"Date\tNumber"} {print $2"\t"$1}' > count_by_month.tsv

echo 'library(ggplot2)
a = read.csv(file="count_by_month.tsv", sep="\t", header=TRUE, colClasses=c("Date", "numeric"))
ggplot(a, aes(x=Date, y=cumsum(Number))) +
  geom_line() +
  ylab("Cumulative number of genomes") +
  ggtitle("Assemblies in GenBank")
ggsave("plot.png", scale=0.7)' > script.R

R CMD BATCH script.R

