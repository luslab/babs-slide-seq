#!/usr/bin/Rscript

library(dplyr)
library(readr)

annot <- readr::read_csv("annotations.csv")

struct <-
	readr::read_csv("read_structures.csv") %>%
	dplyr::mutate(run=as.character(run))


df <-
	readr::read_csv("samples.csv") %>%
	dplyr::select(sample_name, species, fastq_1, fastq_2) %>%
	dplyr::distinct() %>%
	dplyr::rename("name"=sample_name) %>%
	dplyr::mutate(puck=name) %>%
	dplyr::mutate(run=gsub("_.*", "", name)) %>%
	dplyr::left_join(struct, by="run") %>%
	dplyr::left_join(annot, by="species")

readr::write_csv(df, "design.csv")

