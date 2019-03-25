#!/usr/bin/env Rscript

library(ggplot2)
library(RColorBrewer)

setwd(".")

args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)<4) {
  stop("Four arguments must be supplied: infile, outfile, region bed file, qual cutoff", call.=FALSE)
}

infile<-args[1]
out_plot<-args[2]
region_file<-args[3]
qual_set <- as.numeric(args[4])


region_df<-read.table(region_file, sep="\t", header=T)
region_df

chr = region_df[1, "chr"]
start = as.numeric(region_df[1, "start"])
end = as.numeric(region_df[1, "end"])

chr
start
end

dat<-read.table(infile, sep="\t", header=T)
dat_sub<-subset(dat, chr_x==paste0("chr",chr) & pos_x>start & pos_x<end & qual_x>=qual_set & qual_y>=qual_set)
dat_sub

theme_set(theme_bw(20))
p1 <- ggplot(dat_sub, aes(x=pos_y/1000,y=pos_x/1000000)) + 
  geom_point(size=0.5, color = rgb(255, 0, 0, 30, maxColorValue = 255)) + 
  scale_colour_manual(name = "chr_y") + #values=col_v_rndm
  theme(strip.background=element_blank(), 
        legend.title=element_blank(),
	legend.text=element_text(size=15),
	legend.background = element_rect(color = "black"),
        axis.title.x = element_text(size=18),
        axis.title.y = element_text(size=18)) +
  xlab("scaffold coordinate (kb)") +
  ylab(paste0("chr", chr," coordinate (Mb)")) +
  ylim(start/1000000, end/1000000)

#p1

ggsave(out_plot, plot=p1, width=4, height=5.5, units="in", dpi=100)

