\name{PSMnorm}
\alias{appendTabellona}

\title{ Append single patients results in Tabellona
}
\description{
This function updates Tabellona.xlsx file from two .txt patient's files
(TargetPeptideSpectrumMatch and TargetProtein).
}
\usage{
appendTabellona(targetPeptide_name, targetProtein_name, Tabellona_name)
}
\arguments{
  \item{targetPeptide_name}{ TargetPeptideSpectrumMatch.txt filename}
  \item{targetProtein_name}{ targetProtein_name.txt filename}
  \item{Tabellona_name}{ Tabellona_name.xlsx filename}
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\item{TabellonaUPDATED }{.xlsx file with patient's data appended in the last row}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Eleonora Fresi
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
targetPeptide <- "1234 - AB-11-111 TargetPeptideSpectrumMatch.txt"
targetProtein <- "1234 - AB-11-111 TargetProtein.txt"
Tabellona_name <- "Tabellona.xlsx"
appendTabellona(targetPeptide, targetProtein, Tabellona_name)
# the function updates "Tabellona.xlsx" file in the defined working directory
}




