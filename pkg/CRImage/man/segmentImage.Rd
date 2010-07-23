\name{segmentImage}
\alias{segmentImage}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
 The function segments cells or cell nuclei in the image.
}
\usage{
segmentImage(filename = "", image = NA, maxShape = NA, minShape =NA, failureRegion = NA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{filename}{ A path to an image}
  \item{image}{ An 'image' object, if no filename is specified.}
  \item{maxShape}{Maximum size of cell nuclei}
  \item{minShape}{ Minimum size of cell nuclei }
  \item{failureRegion}{ minimum size of failure regions}
}
\details{
The image is converted to greyscale and thresholded.
 Clutter is deleted using morphological operations. Clustered objects are sperated using watershed algorithm.
 Segmented Cell nuclei, which exceed the maximum size are thresholded and segmented again. Cell nuclei which fall below the minimum size are
 deleted. Dark regions which exceed the parameter failureRegion are considered as artefacts and deleted.
 If the parameters are not defined, the operations will not be executed.
 Features are generated for every segmented object.
}
\value{
  A list is returned containing
  \item{image}{The original image}
  \item{segmented image}{Them segmented image}
  ...
}
\references{EBImage, 'http://www.bioconductor.org/packages/release/bioc/html/EBImage.html' }
\author{ Henrik Failmezger }

}
\examples{
#segment image
f = system.file("data", "exImg.jpg", package="CRImage")
segmentationValues=segmentImage(f)
image=segmentationValues[[1]]
segmentedImage=segmentationValues[[2]]
imageFeatures=segmentationValues[[3]]
}
\keyword{misc}
