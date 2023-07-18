\name{GSEAfact}
\alias{GSEAfact}
\title{
Gene Set Enrichment Analysis of Gene Ontology / KEGG pathway
}
\description{
Gene Set Enrichment Analysis of Gene Ontology / KEGG pathway
}
\usage{
gsearesult <- GSEAfact(geneList, org = "human", ont = "KEGG", ver = "latest", dir = system.file(package = "GOfact"),
                       exponent = 1, nPerm = 1000, minGSSize = 1, maxGSSize = 20000,
                       pvalueCutoff = 0.05, pAdjustMethod = "BH", verbose = TRUE, seed=FALSE)
}
\arguments{
  \item{geneList}{
order ranked geneList
}
  \item{org}{
organism information, one of "human","mouse" and "rat"
}
  \item{ont}{
ontology information, one of "MSigDB","GO" and "KEGG"
}
  \item{ver}{
ontology files version, update ontology information by GOupdate & KEGGupdate
}
  \item{dir}{
ontology files directory
}
  \item{exponent}{
weight of each step
}
  \item{nPerm}{
permutation numbers
}
  \item{minGSSize}{
minimal size of each geneSet for analyzing
}
  \item{maxGSSize}{
maximal size of genes annotated for testing
}
  \item{pvalueCutoff}{
pvalue Cutoff
}
  \item{pAdjustMethod}{
pvalue adjustment method
}
  \item{verbose}{
print message or not
}
  \item{seed}{
logical
}
}

\value{
gseaResult object
}
\references{
Yu G, Wang L, Han Y, et al. clusterProfiler: an R Package for Comparing Biological Themes Among Gene Clusters[J]. Omics A Journal of Integrative Biology, 2012, 16(5): 284-287.
}
\author{
{\url{lidong.bprc@foxmail.com}}

{\url{jiangp@geneworks.cn}}
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
gsearesult <- GSEAfact(geneList,
                       ont = "KEGG",
                       ver = "latest",
                       dir = system.file(package = "GOfact"),
                       exponent = 1,
                       nPerm = 1000,
                       minGSSize = 1,
                       maxGSSize = 20000,
                       pvalueCutoff = 0.05,
                       pAdjustMethod = "BH",
                       verbose = TRUE,
                       seed=FALSE)

gseatable <- gsearesult@result
}
