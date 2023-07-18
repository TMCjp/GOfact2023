\name{updateKOA}
\alias{updateKOA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
update kegg pathway annotation
}
\description{
update kegg pathway annotation
}
\usage{
updateKOA(org = "human", ver = format(Sys.time(), "\%Y\%m\%d"), dir = system.file(package = "GOfact"))
}

\arguments{
  \item{org}{
organism information, one of "human","mouse" and "rat"
}
  \item{ver}{
KEGG pathway annotation version information, default is the date of the day
}
  \item{dir}{
KEGG pathway annotation directory information
}
}

\value{
output Rdata to directory
}
\references{
{\url{https://www.kegg.jp/kegg/rest/keggapi.html}}
}

\author{
{\url{lidong.bprc@foxmail.com}}

{\url{jiangp@geneworks.cn}}
}

\examples{
updateKOA(org = "human", ver = "20190618", dir = './')
}
