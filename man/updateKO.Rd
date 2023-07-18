\name{updateKO}
\alias{updateKO}

\title{
update KEGG pathway hierarchy
}
\description{
update KEGG pathway hierarchy
}
\usage{
updateKO(ver = format(Sys.time(), "\%Y\%m\%d"), dir = system.file(package = "GOfact"))
}

\arguments{
  \item{ver}{
KEGG pathway hierarchy version information, default is the date of the day
}
  \item{dir}{
KEGG pathway hierarchy directory information
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
updateKO(ver = "20190618", dir = './')
}
