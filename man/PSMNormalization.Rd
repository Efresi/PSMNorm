\name{PSMnorm}
\alias{PSMnormalization}

\title{ Peptide-Spectrum Match (PSM) normalization
}
\description{
This function takes two .txt files (TargetPeptideSpectrumMatch and TargetProtein)
of a patient and gives a .xlsx formatted file with normalized PSM counts.
It provides also the sum of PSM amyloidogenic proteins:
(Apolipoprotein E, Apolipoprotein A-IV, Serum amyloid P-component)
}
\usage{
PSMnormalization(targetPeptide_name, targetProtein_name)
}
\arguments{
  \item{targetPeptide_name}{ TargetPeptideSpectrumMatch.txt filename}
  \item{targetProtein_name}{ targetProtein_name.txt filename}

}
}
\value{
\item{normalizedOUTPUT }{.xlsx formatted file with normalized PSMs and amyloidogenic proteins}
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

\examples{
targetPeptide <- "1234 - AB-11-111 TargetPeptideSpectrumMatch.txt"
targetProtein <- "1234 - AB-11-111 TargetProtein.txt"
PSMnormalization(targetPeptide, targetProtein)
# the function creates "AB-11-111.xlsx" file in the defined working directory
}

