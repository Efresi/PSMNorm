\name{PSMnorm}
\alias{appendTabellona}

\title{ Append single patients results in Tabellona
}
\description{
This function updates Tabellona.xlsx file from two .txt patient's files
(TargetPeptideSpectrumMatch and TargetProtein).
}
\usage{
appendTabellona(targetPeptide_name, targetProtein_name, Tabellona_PATH, id_pz_escl = NA)
}
\arguments{
  \item{targetPeptide_name}{ TargetPeptideSpectrumMatch.txt filename}
  \item{targetProtein_name}{ targetProtein_name.txt filename}
  \item{Tabellona_PATH}{ a file path for Tabellona.xlsx}
  \item{id_pz_escl}{ vector of patient ID to be excluded from Tabellona (Default=NA inserts every patient)}
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\item{Tabellona }{.xlsx file with patient's data appended in the last row}
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
Tabellona_PATH <- "C:/Users/utente/prova/Tabellona.xlsx"
appendTabellona(targetPeptide, targetProtein, Tabellona_PATH)
# the function updates "Tabellona.xlsx" file
}
