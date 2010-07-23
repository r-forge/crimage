\name{processAperio}
\alias{processAperio}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Cellularity Calculation of Aperio TX Scanner}
\description{
  Procession of Aperio TX Slides.
}
\usage{
processAperio(classifier = classifier, inputFolder = inputFolder, outputFolder = outputFolder, identifier = identifier, numSlides = numSlides, cancerIdentifier = cancerIdentifier)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{classifier}{ The classifier.}
  \item{inputFolder}{ The path to the image folder. }
  \item{outputFolder}{ The path to the output folder. }
  \item{identifier}{ The identifier of the files ("Ss" or "Da") }
  \item{numSlides}{ The number of sections in the image. }
  \item{cancerIdentifier}{ The identifier of the cancer class}
}
\details{
  The function processes images of Aperio TX scanners. The images have to be saved in the CWS format. 

}
\value{
 Four folders are created in the output folder.
  \item{Files}{Cellularity values and cell numbers are saved in the file}
  \item{classifiedImage}{Subimages wit marked tumour and non tumour cells}
  \item{tumourDensity}{Cancer heatmaps for every subimage}
  \item{cellCoordinates}{Coordinates and cell class for every cell in the subimage}
}
\references{ ~put references to the literature/web site here ~ }
\author{ Henrik Failmezger }

}
\examples{
#classify aperio
f = system.file("data",  package="CRImage")
load("/Users/failmezger/package/classifier")
classifier=model
dir.create("AperiOutput")
processAperio(classifier,f,"AperiOutput","Da",2,"c")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{misc}
