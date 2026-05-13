#Day 2 of Proteomics workshop
#Boxes and Jitters: alternatives to Dynamite Plots

#Load the needed library
library(ggplot2)

#If the library is not install use the install command on terminal or from packages -> type install and install
#the needed library

### Read df_transcripts.csv into R
Transcripts <- read.table(file="df_transcripts.csv",sep=",",header=T)

### Pull out just the transcripts with qvalues below 0.05; require a gene name, too.
NamedOnly <- subset(Transcripts, geneNames != ".")
Diffs = subset(NamedOnly, qval < 0.05)

### Remove the gene that has intensities SO MUCH BIGGER than the others.
Diffs = subset(Diffs, geneIDs != "MSTRG.506")

### How many genes were called differential? #with using [1] I am taking the first element -> num of rows 
NumDiffs = dim(Diffs)[1]

### Create a vector to specify males and females from "geuvadis_phenodata.csv"
FASTQ = c("ERR188044","ERR188104","ERR188234","ERR188245","ERR188257","ERR188273","ERR188337","ERR188383","ERR188401","ERR188428","ERR188454","ERR204916")
Sex = c("male","male","female","female","male","female","female","male","male","female","male","female")
Sex = rep(Sex,each=NumDiffs)
Population = c("YRI","YRI","YRI","GBR","GBR","YRI","GBR","GBR","GBR","GBR","YRI","YRI")
Population = rep(Population,each=NumDiffs)
Population

### Create a long differences table from the wide one, including just geneIDs and FPKM measurements
LongDiffs = reshape(Diffs[c(1,3:14)], direction="long", varying=list(2:13), times=FASTQ)
### Handle some renaming stuff that results from reshape
colnames(LongDiffs) = c("geneNames","FASTQ","FPKM","id")
rownames(LongDiffs) = NULL

### Append the metadata columns and rename
ExtDiffs = cbind(Sex,Population, LongDiffs)
ExtDiffs

### About the "Grammar of Graphics"
# ggplot2 embodies a _layered_ approach to describe and construct visualizations in a structured manner.
# (https://towardsdatascience.com/a-comprehensive-guide-to-the-grammar-of-graphics-for-effective-visualization-of-multi-dimensional-1f92b4ed4149)
# For each image, expect to define each of the following:
# Data: what object contains the information for plotting?
# Aesthetics: what field governs positioning and axes?
# Scale: will we need to mathematically transform values or add breaks in an axis?
# Geometrics Objects: what visual items depict the data?
# Statistics: are we plotting the data or values derived from the data?
# Facets: will we subdivide the visual space to multiple subplots?
# Coordinate system: are we sticking with Cartesian coordinates or getting crazy with polar?

### "Dynamite plot" via ggplot2
# This site has excellent information on plotting summaries: https://ggplot2tutor.com/tutorials/summary_statistics
# In this case, we are not plotting our data directly but rather a statistical summary based on our data
ggplot(ExtDiffs, aes(x=geneNames, y=FPKM, fill=Sex)) +
  stat_summary(fun="mean", geom="bar", position="dodge") +
  stat_summary(fun.data="mean_se", geom="errorbar", position="dodge")


### Boxplot via ggplot2
# http://www.sthda.com/english/wiki/ggplot2-box-plot-quick-start-guide-r-software-and-data-visualization
# Again, we are plotting just a statistical summary based on our data
ggplot(ExtDiffs, aes(x=geneNames, y=FPKM, fill=Sex)) +
  geom_boxplot()



### QUESTIONS:
# 1) What do the pieces of each box represent:
# -A Top of upper whisker
# -B Top edge of box
# -C Midline of box
# -D Lower edge of box
# -E Bottom of lower whisker
# 2) What do the dots outside the boxes mean?
# 3) Does this plot visualize the mean values?

### "Jittered dot plot" via ggplot2
# https://hbiostat.org/bbr/descript.html
# In this case, we are plot our data directly plus a statistical summary based on our data
ggplot(ExtDiffs, aes(x=geneNames, y=FPKM, fill=Sex)) +
  geom_boxplot(position="dodge") +
  geom_dotplot(binaxis="y", stackdir="center",binwidth=5)

# We can also omit the statistical summary, showing the data directly.
ggplot(ExtDiffs, aes(x=geneNames, y=FPKM, fill=Sex)) +
  geom_dotplot(binaxis="y", stackdir="center",binwidth=5)


### QUESTIONS:
# 1) Do the data points look like we would have expected on the basis of the box positions? of the means and SEs?
# 2) Why are so many blue points appearing in a horizonal line at the bottom?
# 3) How many points could you practically visualize this way?

### Save image in a specific format, with user given dimensions. 
ggsave("jitter.pdf",width=14,height=12,units="cm")
ggsave("jitter.png",width=1120,height=960,units="px")

#units = pixels




