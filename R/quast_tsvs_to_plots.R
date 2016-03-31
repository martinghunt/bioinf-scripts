#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

if (length(args) < 2){
    write(paste("Usage:", commandArgs()[1], "<transposed_report.tsv> <outprefix> [gage_transposed_report.tsv]"), stderr())
    write("\nMakes some plots using the output of QUAST", stderr())
    quit(save="no", status=1)
}



quast_tsv = args[1]
out = args[2]

library(ggplot2)


write_bar_plot <- function(d, column, y_label, filename) {
    t = theme(text=element_text(size=15), axis.text.x = element_text(angle=45, hjust=1,size=10), legend.position="none")
    ggplot(data=d) +
        geom_bar(aes_string(x="Assembly",y=column), stat="identity", width=0.75, position=position_dodge()) +
        t +
        scale_x_discrete(labels=d$Assembly, limits=d$Assembly) +
        ylab(y_label)
    ggsave(filename=filename, scale=0.7)
}



quast = read.csv(file=quast_tsv, sep="\t", header=TRUE, na.strings="-")
write_bar_plot(quast, "X..contigs.....0.bp.", "Number of contigs", paste(out, "number_of_contigs.pdf", sep="."))
write_bar_plot(quast, "X..local.misassemblies", "Local missasemblies", paste(out, "local_misassmeblies.pdf", sep="."))

if (length(args) == 3) {
    gage = read.csv(file=args[3], sep="\t", header=TRUE, na.strings="-")
    write_bar_plot(gage, "Avg.idy", "Average identity", paste(out, "average_identity.pdf", sep="."))
    write_bar_plot(gage, "Compressed.reference.bases", "Compressed refrence bases", paste(out, "compressed_ref_bases.pdf", sep="."))
    write_bar_plot(gage, "Duplicated.reference.bases", "Duplicated refrence bases", paste(out, "duplicated_ref_bases.pdf", sep="."))
    write_bar_plot(gage, "Indels....5", "Indels at least 5bp long", paste(out, "indels_at_least_5bp.pdf", sep="."))
    write_bar_plot(gage, "Indels...5bp", "Indels shorter than 5bp", paste(out, "indels_shorter_than_5bp.pdf", sep="."))
    write_bar_plot(gage, "Inversions", "Inversions", paste(out, "inversions.pdf", sep="."))
    write_bar_plot(gage, "Relocation", "Relocations", paste(out, "relocations.pdf", sep="."))
    write_bar_plot(gage, "Translocation", "Translocations", paste(out, "relocations.pdf", sep="."))
    write_bar_plot(gage, "SNPs", "SNPs", paste(out, "snps.pdf", sep="."))
}

unlink("Rplots.pdf") # Aaaargh!
