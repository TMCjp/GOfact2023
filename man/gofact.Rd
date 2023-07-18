\name{gofact}
\alias{gofact}
\title{
Gene Ontology Enrichment Analysis
}
\description{
Ontology Enrichment Analysis of a gene set. Given a vector of genes or a data.frame which with multicol genes, this function will return the enrichment Ontology categories after FDR control.
}
\usage{
gofact(gene, slim = "kegg", org = "human",
       ver = "latest", dir = system.file(package = "GOfact"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gene}{
a vector of genes or a data.frame which with multicol genes
}
  \item{slim}{
ontology slim information
}
  \item{org}{
organism information, one of "human","mouse" and "rat"
}
  \item{ver}{
ontology files version, update ontology information by GOupdate & KEGGupdate
}
  \item{dir}{
ontology files directory
}
}
\value{
A enrichResult instance.
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
genelist <- c("4312","8318","10874","55143","55388","991","6280","2305","9493","1062","3868")
goresult <- gofact(genelist, slim = "kegg", org = "human", ver = "latest")

exportHtml(goresult, "test.html")
exportExcel(goresult, "test.xlsx", sheetname = "test")
}

