# 16S amplicon analyses in R

### Lets start with some basics in statistics and touch up on the "Microbial Ecology" we are investigating

Microbial ecology is the study of microorganisms in their natural environments.

* What have we done up to here?
We sequenced the amplified product of the 16S gene of an entire community.
* What does this mean?
We have taken a gene possessed by all members of a community, and PCR amplified it to allow sequencing. 
Sequencing allowed us to get all the members in the communities' genes. 
* What now?

We covered rarefaction. Why do that?
This allows us to get an 'equal' sampling depth per community, not counting more in one sample than the other. 

## Our data

#### In R, we need to specify the packages we are going to use
This is a simple script which load all the packages we are going to use in one go:
```
# load pacakges into R 
libs <- c("ggplot2", "tidyverse", "vegan", "phyloseq", "gplots", "venneuler", "reshape")
lapply(libs, require, character.only = TRUE)
```

#### We first need to read our data into R 
```
#Set working directory
setwd("~/Software_carpentry_course/otu_table/")

#Load OTU table into R and check dimensions
otu_table <- read.delim("./ninja_otutable_format.txt", header = T, row.names = 1)
str(otu_table)

#Read last col, and remove if taxonomy
otu_table[,length(otu_table)]
otu_table <- as.data.frame(t(otu_table[,1:length(otu_table)-1]))

#Load mapping file containing groups for samples and make it available for R in the global environment
mapping_file <- read.delim("./mapping_file_format.txt", header = T, row.names = 1)
attach(mapping_file)
```
We will also be using a package called phyloseq, that requires the data to be loaded into an object specific to the package:
```
#Create a phyloseq object which we need for some of the analyses
#Read the otu table, mapping file (already done) and reference taxonomy sets - the assigned taxonomy by QIIME
taxonomy_mapping <- read.delim("./taxonomy_phyloseq.txt", row.names = 1)

#Make the reference taxa a matrix
taxonomy_mapping <- as.matrix(taxonomy_mapping)

#Remove X's in col names and correctly transpose dataframe
otu_table <- t(otu_table)
colnames(otu_table) <- row.names(SAMPL)

#Create otu table and tax for phyloseq objects and make the combined object
OTU <- otu_table(otu_table, taxa_are_rows = T)
TAX <- tax_table(taxonomy_mapping)
SAMPL <- sample_data(mapping_file)
otus_physeq = phyloseq(OTU, TAX, SAMPL)
otus_physeq 


```

## What questions can we ask with the data we have?


### Diversity based


#### Alpha diversity
Well first of all, we want to understand how many 'species' are present in a given environment.
	This allows us to compare it to other environments, determine if it is more or less diverse, and whether it is as diverse as expected.

Some of the measures used are as follow:

* Species richness: This is a measure taking only richness into consideration
```
alpha <-data.frame(specnumber(otu_table))
alpha.mean <-tapply(specnumber(otu_table), type, mean)
alpha.sd <-tapply(specnumber(otu_table), type, sd)
names(alpha) <-c("alpha")
```

##### Test the differences in diversities between variables
It is now possible to test if these diversities are significantly different between sampling sites using an analysis of variance.
```
attach(alpha)
ANOVA1 <-aov(alpha~type)
summary(ANOVA1)
TukeyHSD(ANOVA1)
```
#### **Challenge**
We get species number, which is what you tested. How would you change it to be for Shannon diversity, Simpson and Inverse Simpson?

These measures all give a good estimate of the number and distribution of species within a single sample. 
But what if we want to have a look at more than one sample?

#### Gamma diversty
Before assessing beta diversity, which I will explain next, it is important that we look at gamma diversity.
Gamma diversity is simply the total diversity of an environment, or sampling group:
```
gamma <-specnumber(otu_table, type)
gamma
```

#### Beta diversity
To assess how different samples are within sampling groups, we use beta diversity.
This measure gives an overall impression of how similar samples are, within a defined group. 

To make it simple, beta diversity is simply gamma/alpha. 
This means two things, each sample will have its own beta diversity, giving its similarity/dissimilarity to its defined group.

The second is better explained by an example:
If you have four samples with each 200 species, and the gamma diversity is also 200, it means that all four samples are identical, and beta = 1.
In a different scenario, you have four samples with each 200 species, and the gamma diversity is 800, it means that you have no sample with a single species in common, 
and beta will be four. 

```
#Calculate the average beta diversity
beta.whittaker <-gamma/alpha.mean

#First create separate dataframes with each of the sites
#Then calculate the beta diversity of each of the sites individually

het <- alpha[1:4,]
ko <- alpha[5:62,]
wt <- alpha[63:116,]

beta_het <- gamma[1]/het
beta_ko <- gamma[2]/ko
beta_wt <- gamma[3]/wt 
```

##### Test the diversities
It is again possible to test if these diversities are significantly different between sampling sites using an analysis of variance.
```
beta_all <-data.frame(c(beta_het, beta_ko, beta_wt), mapping_file)
names(beta_all)<-c("beta","type")
anova.beta <-aov(beta ~type, data=beta_all)
summary(anova.beta)
TukeyHSD(anova.beta)
```

With this you will note that beta is dependent on sample number and the group size. 
This makes the practice of comparing beta diversities between studies undesirable, and I would reccommend using it with caution. 


### Community structure based
These methods are still diversity based, assessing beta-diversity, but the underlying statistics uses whole community data.
#####Transformation
We need to transform community data, this allows us to compare 'apples with apples'
```
otu_transf <- decostand(otu_table, "hellinger")
```

### Other ways of visualizing beta diversity
While the values of beta diversity provides a great way to "look" at how dissimilar or similar your samples are, 
it is not intuative, and you will not be able to explain it to someone outside of microbial ecology. 

Within statistics there are numerous ways to visualize multidimentional data, which simplifies the description of beta diversity, or the differences between samples.

Generally PCA and NMDS is most widely used:
#### PCA
Principal Component Analysis is a commonly used dimentionality reduction technique. With the use of multiple regressions between each of the "dimensions" which is actually only your OTUs/species, the technique tries to explain as much of the variation as possible on as few new axes, called (principal components).

With this statistical test, the amount of variation explained by each principal component is given.

We typically only use the first two components, as this allows visualization on two axes, but three is also possible, and even more can be used in statistical tests. 
The only assumption of this test, as with many others, is normality.
Biological and specifically microbial community data are not assumed to be normal, and we use a non-parametric approach instead. 

#### NMDS
A nonparametric alternative to PCA is Non-metric Multidimensional Scaling.
This method models the data in multidimensional space, and tries to find a plane which best intersects the data, minimizing the variance around the axes. 

This variance around the axes is then explained as a stress value, which shows how well the data suits the model. 

This value has to be below 0.2 to remain statistically viable. 

```
#Bray curtis distance matrix in vegan
vegdist(otu_transf, "bray") -> d

#Perform the multidimensional scaling using metaMDS
fit <- metaMDS(d, "bray", k = 2, trymax = 1)

#Take a look at the results
fit

#Extract the scores used for the plot
data.scores <- as.data.frame(scores(fit))
data.scores

attach(mapping_file)
#Plot using ggplot
#Set the colors and shapes
cols = rainbow(3)
#Change shapes using pch numbers
shps = c(22, 21)

plot_1 <- ggplot() +
  geom_point(data = data.scores, aes(x = NMDS1, y = NMDS2, fill = type, shape = location), size = 6) + 
  scale_fill_manual(values = cols)+
  scale_shape_manual(values = shps) +
  coord_fixed(ratio = 1.05) + 
  theme_bw()  +
  annotate("text", x = 0.15, y = -0.55, label = "stress = 0.14", size = 7) +
  theme(axis.text.x = element_blank(), 
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_text(size = 18), 
        axis.title.y = element_text(size = 18), 
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank(),
        legend.text = element_text(size = 18), 
        legend.title = element_text(size = 18))
```
In both these methods the distance between two points indicate the relationship between samples. 

#### RDA
A third way of visualizing the data will be through Redundancy analysis. 
This method is also parametric, but allows the incorporation of environmental variables. 
This shows how much variation each of the variables explain.
However, for this we need continuous variables for explanation of biological data, and we do not have any for this dataset.

#### Venn
We can also assess the comonalities of microbial communities through venn diagrams.
This shows how members are shared between sampling sites/groups.
```
#Venn diagram

otu_table_t <- t(otu_table)
comb <-data.frame(cbind(rowSums(otu_table_t[,1:4]), rowSums(otu_table_t[,5:62]), rowSums(otu_table_t[,63:116]))>0)
names(comb)<-c("Het", "KO", "WT")
venn(comb)

#Or draw using venneuler

plot<- venneuler(comb > 0)
plot(plot)
```

#### Permanova
This is a permutational multi-variate analysis of variance. 
It allows us to assess if the entire community is different between sites, using the basics of analysis of variance, and applying it to multiple variables and non-parametric data. 

`adonis(otu_trans~type*location)`

We can fine tune the model, as we need a model:

`adonis(otu_trans~location)`

### Taxonomy based
#### Barplots
One of the main objectives of microbial ecology is to understand WHO is in an environment. 
This is allowed by using taxonomic databases, with assigned taxonomy. 
One of the most simple ways of looking at this is through barplots.
```
#Load summarized dataset
otu_table_tax_l2 <- read.delim("./ninja_otutable_L2.txt", header = T, row.names = 1)
mapping_file <- read.delim("./mapping_file_format.txt", header = T)

#Transpose table to be the correct orientation
otu_table_tax_l2 <- t(otu_table_tax_l2)

#Combine the OTU table and mapping file with the variables corresponding to the correct samples
otu_map <- cbind(otu_table_tax_l2, mapping_file)

#Simplify this into a format which makes plotting easier for R
otu_melted <- melt(otu_map)

#Ensure the correct order on the x axis
corr_order <- mapping_file$Sample
otu_melted$Sample <- factor(otu_melted$Sample, levels = corr_order)

attach(mapping_file)
filcols <- rainbow(9)
plot1 <- ggplot() + 
  geom_bar(aes(y = value, x = Sample, fill = variable), colour = "black", data = otu_melted, stat = "identity") +
  labs(x = "Sample Site", y = "Relative abundance (%)") +
  scale_fill_manual(values = filcols) +
  theme(axis.line.x = element_line(size=.5, colour = "black"),
        axis.line.y = element_line(size=.5, colour = "black"),
        axis.text.x=element_text(colour="black", size = 10, angle = 90, hjust = 1),
        axis.text.y=element_text(colour="black", size = 10),
        legend.key=element_rect(fill="white", colour="white"),
        legend.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(),
        panel.background = element_blank())
```

#### Multiple regressions
If we have environmental parameters, we could compare them to each of the taxa found, at various taxonomic levels. 
I.e., does age correlate with Proteobacteria?
Again we do not have continuous variables for the data, and can't do this.

### Microbial community interactions based
This section actually needs an entire course on its own, but I will try to make it as easy as possible. Please read the papers if you want to know more. 

We can infer putative microbial community interactions through co-occurrance. This has mostly been done through multiple regressions and correlations. Simply put, if a species have a signficant correlation, they cooccur.

This is best visualized using networks. With species representing nodes, and interactions being edges, which can be positive or negative depending on the method used. 

```
#First subset the otu table
otu_table_ko <- otu_table[,5:62]
otu_table_wt <- otu_table[,63:116]

#Create otu table and tax for phyloseq objects and make the combined object
TAX <- tax_table(taxonomy_mapping)
SAMPL <- sample_data(mapping_file)
OTU_KO <- otu_table(otu_table_ko, taxa_are_rows = T)
OTU_WT <- otu_table(otu_table_wt, taxa_are_rows = T)
otus_ko_physeq = phyloseq(otu_table_ko, TAX, SAMPL)
otus_wt_physeq = phyloseq(otu_table_ko, TAX, SAMPL)

#Make network for KO
spiec.out=spiec.easi(otus_ko_physeq , method="mb", lambda.min.ratio = 1e-2, nlambda = 9,icov.select.params=list(rep.num=20))
spiec.graph=adj2igraph(spiec.out$refit, vertex.attr=list(name=taxa_names(otus_ko_physeq )))
plot_network(spiec.graph, otus_ko_physeq , type='taxa', color="Phylum", label=NULL)

#Make network for WT
spiec.out=spiec.easi(otus_wt_physeq , method="mb", lambda.min.ratio = 1e-2, nlambda = 9,icov.select.params=list(rep.num=20))
spiec.graph=adj2igraph(spiec.out$refit, vertex.attr=list(name=taxa_names(otus_wt_physeq )))
plot_network(spiec.graph, otus_wt_physeq , type='taxa', color="Phylum", label=NULL)
```
And this output is best visualized in Cytoscape.
Remember to do this for both networks.

```
write.graph(spiec.graph,file="spieceasi.ncol.hmp.txt",format="ncol") 
write.table(TAX,file="taxonomy_network.txt",sep="\t", quote=FALSE)
```

