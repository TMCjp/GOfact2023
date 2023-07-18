\name{OUTtype}
\alias{OUTtype}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
gofact analysis result output type
}
\description{
output gofact analysis result
}
\usage{
exportHtml(result, htmlfile = "test.html")

exportExcel(result, xls.file = "test.xlsx", sheetname = "test", multisheet = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{result}{
gofact analysis result
}
  \item{xls.file}{
output xlsx file
}
  \item{sheetname}{
sheet name, if multisheet = TRUE, sheetname is invalid.
}
  \item{multisheet}{
Logical, whether output multiple sheet xlsx file
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
html file;
xlsx file;
}
\references{
%% ~put references to the literature/web site here ~
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
