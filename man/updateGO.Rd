\name{updateGO}
\alias{updateGO}

\title{
update gene ontology hierarchy
}
\description{
update gene ontology hierarchy
}
\usage{
updateGO(ver = format(Sys.time(), "\%Y\%m\%d"), dir = system.file(package = "GOfact"))
}

\arguments{
  \item{ver}{
Gene ontology hierarchy version information, default is the date of the day
}
  \item{dir}{
Gene ontology hierarchy directory information
}
}

\value{
output Rdata to directory
}
\references{
Carlson M (2019). GO.db: A set of annotation maps describing the entire Gene Ontology. R package version 3.8.2.
}

\author{
{\url{lidong.bprc@foxmail.com}}

{\url{jiangp@geneworks.cn}}
}

\examples{
updateGO(ver = "20190618", dir = './')
}
