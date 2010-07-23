\name{CRImage-package}
\alias{CRImage-package}
\alias{CRImage}
\docType{package}
\title{
CRImage is a package to analyze images and classify cells.
}
\description{
CRImage allows classification of cells in biological images. It offers methods to segment cells or cell nuclei
in biological images for example HE stained images. It offers methods to create classifiers and to classify cells in these images.
Furthermore it allows the calculation of tumour cellularities for large microscopical slides.

CRImage makes use of the image processing package EBImage, which uses the 'ImageMagick' library for image I/O operations and the 'GTK' library to display images.
}
\details{
\tabular{ll}{
Package: \tab CRImage\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2010-04-27\cr
License: \tab What license is it under?\cr
LazyLoad: \tab yes\cr
}
~~ An overview of how to use the package, including the most important ~~
~~ functions ~~
}
\author{
Henrik Failmezger, failmezger@cip.ifi.lmu.de
Yinyin Yuan
Florian Markowetz

CRI Cambridge
Ludwigs-Maximilians University of Munich
Technical University of Munich


}
\references{
~~ Literature or other references for background information ~~
}

\keyword{misc}

\examples{
example(segmentImage)
example(createClassifier)
example(classifyImage)
}
