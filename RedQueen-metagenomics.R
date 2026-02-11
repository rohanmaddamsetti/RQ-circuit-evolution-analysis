## RedQueen-metagenomics.R by Rohan Maddamsetti.
## for now, just analyze Day 100 metagenomes (assembled in consensus mode, so mutations at 100% allele frequency).

library(tidyverse)
library(cowplot)
library(patchwork)
library(ggthemes)
library(viridis)
library(ggrepel)


## colorblind-friendly palette.
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

## assert that we are in the src directory, such that
## projdir is the parent of the current directory.
stopifnot(endsWith(getwd(), file.path("Hye-in_RedQueen","src")))
projdir <- file.path("..")

## This file *just* has the evolved populations (Day 9).
pop.clone.labels <- read.csv(
  file.path(projdir,
            "data/sample-metadata.csv"),
  stringsAsFactors=FALSE) %>%
    ## remove the Day 1 samples from the analysis.
    filter(Day != "D1")

## This is the key data file for the analysis.
evolved.mutations <- read.csv(
    file.path(projdir,
              "results/evolved_mutations.csv"),
    stringsAsFactors=FALSE) %>%
    mutate(Mbp.coordinate=Position/1000000)


###############################################
## Figure 1: make a stacked bar plot of the kinds of mutations in each treatment.

## This function sums mutations per replicate population.
make.mutation.class.df <- function(evolved.mutations.df) {
    evolved.mutations.df %>%
        ## give nicer names for mutation classes.
        mutate(Mutation=recode(Mutation,
                               MOB = "Mobile element transposition",
                               DEL = "Indel",
                               INS = "Indel",
                               SUB = "Multiple-base substitution",
                               nonsynonymous = "Nonsynonymous",
                               synonymous = "Synonymous",
                               nonsense = "Nonsense",
                               pseudogene = "Pseudogene",
                               intergenic = "Intergenic",
                               )) %>%
        group_by(Sample, Day, aTc, Spacer_type, Colony, Mutation) %>%
        summarize(Count=n()) %>%
        ungroup() %>%
        data.frame() %>%
        mutate(Mutation=as.factor(as.character(Mutation)))
}


plot.mutation.summary.stackbar <- function(mutation.class.df, leg=FALSE) {

    fig <- ggplot(mutation.class.df, aes(x=Sample, y=Count, fill=Mutation)) +
        ylab("Count") +
        facet_grid(.~Spacer_type) +
        geom_bar(stat='identity') +
        scale_fill_brewer(palette = "RdYlBu", direction=-1,drop=FALSE) +        
        theme_classic(base_family='Helvetica') +
        theme(axis.text.x=element_text(size=12,angle=45,hjust=1),
              axis.text.y=element_text(size=12),
              panel.border=element_blank(),
              strip.background = element_blank(),
              panel.spacing.x=unit(1, "cm"),
              panel.spacing.y=unit(0.5, "cm"))

    if (leg == TRUE) {
        fig <- fig +
            theme(legend.title=element_text(size=8, face="bold"),
                  legend.title.align=0.5,
                  legend.text=element_text(size=8),
                  legend.position="bottom")
    } else {
        fig <- fig + guides(fill = "none")
    }
    
    return(fig)
}

## Now make Figure 1.
mutation.class.df <- make.mutation.class.df(evolved.mutations)

Fig1 <- plot.mutation.summary.stackbar(mutation.class.df)
fig1.output <- "../results/mutation_stackbar.pdf"
ggsave(Fig1, file=fig1.output)


#################################################################################
## analysis of parallel evolution at the same nucleotide.
## discuss numbers and finding in the text (no figure.).
## This could be a Supplementary Table.

bp.parallel.mutations <- evolved.mutations %>% group_by(Position) %>%
summarise(count = n()) %>% filter(count>1) %>% inner_join(evolved.mutations)


parallel.DEL <- filter(bp.parallel.mutations,Mutation=='DEL')
parallel.INS <- filter(bp.parallel.mutations,Mutation=='INS')
parallel.dN <- filter(bp.parallel.mutations,Mutation=='nonsynonymous')
parallel.dS <- filter(bp.parallel.mutations,Mutation=='synonymous')

## examine parallel evolution at amino acid level (only one case, in robA).
parallel.AA.dN <- evolved.mutations %>% filter(Mutation=='nonsynonymous') %>% group_by(Position) %>% summarize(count=n()) %>% filter(count > 1)
parallel.dN.Table <- filter(evolved.mutations, Position %in% parallel.AA.dN$Position) %>% arrange(Position)

## check parallel evolution for synonymous mutations too.
parallel.dS.Table <- filter(evolved.mutations, Position %in% parallel.dS$Position) %>% arrange(Position)


#################################################################################
## analysis of parallel evolution at the gene level (including intergenic regions).

gene.level.parallel.mutations <- evolved.mutations %>% group_by(Gene) %>%
summarise(count = n()) %>% filter(count>1) %>% inner_join(evolved.mutations)

parallel.genes <- gene.level.parallel.mutations %>%
    select(Gene, count, Plasmid, Transposon, Tet) %>%
    distinct() %>%
    arrange(desc(count))

#################################################################################
### Figure 2: make a matrix plot of genes with mutations in two or more clones.
################################################################################
MakeMutCountMatrixFigure <- function(evolved.muts, show.all=FALSE) {

    ## First, make a mutation matrix for plotting.
    matrix.data <- evolved.muts %>%
        group_by(Gene, Sample, aTc, Spacer_type, Colony) %>%
        summarize(mutation.count = n()) 

    total.muts <- matrix.data %>%
        group_by(Gene) %>%
        summarize(total.mutation.count = sum(mutation.count))
    
    matrix.data <- left_join(matrix.data, total.muts)
    
    if (!show.all) { ## then filter out genes that are only hit in one sample.
        matrix.data <- matrix.data %>%
            filter(total.mutation.count > 1)
    }
    
    ## sort genes by number of mutations in each row, but put all the transposon mutations together.
    ## also check out the alternate sorting method that follows.
    gene.hit.sort <- matrix.data %>%
        group_by(Gene, .drop = FALSE) %>%
        summarize(hits=sum(mutation.count)) %>%
        arrange(desc(hits))
    
    ## now sort genes.
    matrix.data$Gene <- factor(matrix.data$Gene,levels=rev(gene.hit.sort$Gene))

    ## cast mutation.count into a factor for plotting.
    matrix.data$mutation.count <- factor(matrix.data$mutation.count)

    make.matrix.panel <- function(mdata, spacer.type, leg=FALSE) {
        panel.data <- filter(mdata, Spacer_type==spacer.type)
        fig <- ggplot(panel.data,
                      aes(x=Sample,
                          y=Gene,
                          fill=mutation.count,
                          frame=Spacer_type)
                      ) +
            geom_tile(color="black",size=0.1) +
            ggtitle(spacer.type) +
            theme_tufte(base_family='Helvetica') +
            theme(axis.ticks = element_blank(),
                  axis.text.x = element_text(size=10,angle=45,hjust=1),
                  axis.text.y = element_text(size=10,hjust=1,face="italic"),
                  axis.title.x = element_blank(),
                  axis.title.y = element_blank()) +
            scale_y_discrete(drop=FALSE) + ## don't drop missing genes.
            scale_fill_manual(name="Mutations",
                              values = c("#ffdf00", "#bebada", "#fb8072", "#80b1d3", "#fdb462"))
        
        if (leg == FALSE) {
            fig <- fig + guides(fill = "none")
        }
        return(fig)
    }


    ## get the legend from the first panel.
    sp1.matrix.panel <- make.matrix.panel(matrix.data, "sp1", leg=TRUE)

    Fig.legend <- get_legend(sp1.matrix.panel)

    ## now remove the legend from the first panel.
    sp1.matrix.panel <- sp1.matrix.panel + guides(fill = "none")
    
    ## Remove the gene labels from the second panel to save space.
    nt.matrix.panel <- make.matrix.panel(matrix.data, "nt") +
        theme(axis.text.y=element_blank())
    
    ## Using the patchwork library for layout.
    matrix.panels <-
        sp1.matrix.panel +
        nt.matrix.panel +
        Fig.legend +
        plot_layout(nrow = 1)

    ## hack to label x-axis from comments at: https://github.com/thomasp85/patchwork/issues/150
    matrix.panels.grob <- patchwork::patchworkGrob(matrix.panels)
    matrix.figure <- gridExtra::grid.arrange(matrix.panels.grob, left = "", bottom = "Evolved populations")
    
    return(matrix.figure)
}


## Figure 2
Fig2 <- MakeMutCountMatrixFigure(evolved.mutations, show.all=TRUE)
Fig2.outf <- "../results/mutation_matrix.pdf"
ggsave(Fig2.outf, Fig2, height=30, width=10)

Fig3 <- MakeMutCountMatrixFigure(evolved.mutations, show.all=FALSE)
Fig3.outf <- "../results/parallel_mutation_matrix.pdf"
ggsave(Fig3.outf, Fig3, height=30, width=10)



 
