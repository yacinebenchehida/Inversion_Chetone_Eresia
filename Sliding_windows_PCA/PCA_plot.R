##############################
# Install required libraries #
##############################
library(plotly)
library(ggplot2)
library(Cairo)
options(bitmapType='cairo') # Necessary to make ggplotly works in the cluster
library(htmlwidgets)
library(rmarkdown)

###############
# Import data #
###############
# Eigenvectors
dat = commandArgs(trailingOnly=TRUE)
subset = read.table(dat[1])
subset = subset[,c(1,3,4)]
colnames(subset)=c("Ind","PC1","PC2")

# Eigenvalues
contribution = read.table(dat[2])
contribution = as.numeric(unlist(contribution))
  PC1_cont = (contribution[1]/sum(contribution))*100
  PC2_cont = (contribution[2]/sum(contribution))*100

# Phenotypes
pheno = read.table(dat[3])

subset = cbind(subset,pheno)
colnames(subset) = c("Ind","PC1","PC2","phenotype")
numb_color = length(unique(subset$phenotype))

###################
# Get output name #
###################
out_name = "Marsaeus"
print(paste(dat[4],"/",out_name,".html",sep=""))

########
# plot #
########
p1 <- ggplot(subset, aes(x=PC1, y=PC2,fill=phenotype)) + 
  geom_point(size = 3,color="black",shape=21,aes(text = Ind)) + 
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values=rainbow(numb_color)) +
  theme(legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-10,20,-5,0)) +
    xlab(paste("PC1 (",round(x=PC1_cont,2),"%)",sep="")) +
  ylab(paste("PC2 (",round(x=PC2_cont,2),"%)",sep=""))  

pdf(paste(dat[4],"/",out_name,".pdf",sep=""),12,6)
p1
dev.off()

#########################################################
# Create an html widget to have an interactive document #
#########################################################
widg = ggplotly(p1,tooltip = c("Ind","phenotype"))
htmlwidgets::saveWidget(as_widget(widg),paste(dat[4],"/",out_name,".html",sep=""))

