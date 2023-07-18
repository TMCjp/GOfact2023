\name{updateGOA}
\alias{updateGOA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
update gene ontology annotation
}
\description{
update gene ontology annotation
}
\usage{
updateGOA(org = "human", ver = format(Sys.time(), "\%Y\%m\%d"), dir = system.file(package = "GOfact"))
}

\arguments{
  \item{org}{
organism information, one of "human","mouse" and "rat"
}
  \item{ver}{
Gene ontology annotation version information, default is the date of the day
}
  \item{dir}{
Gene ontology annotation directory information
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
updateGOA(org = "human", ver = "20190618", dir = './')
}
