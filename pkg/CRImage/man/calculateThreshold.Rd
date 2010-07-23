\name{calculateThreshold}
\alias{calculateThreshold}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  Calculates the grey value which separates the grey-level image best in a binary image.
  To calculate the threshold, the Otsu Method is used.
}
\usage{
calculateThreshold(allGreyValues)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{allGreyValues}{ a vector of grey values. }
}
\details{
  ~~ If necessary, more details than the description above ~~
}
\value{
  The calculated threshold.
}
\references{ Otsu, N. A threshold selection method from gray level histograms IEEE Trans. Systems, Man and Cybernetics, 1979, 9, 62-66}
\author{ Henrik Failmezger }
\note{ ~~further notes~~ 

 ~Make other sections like Warning with \section{Warning }{....} ~
}
\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

f= system.file("data", "exImg.jpg", package="CRImage")
img=readImage(f)
#convert to grayscale
imgB=channel(img,"gray")
#find white pixels and exclude them from thresholding(if white is background)
indexWhitePixel=which(img[,,1]>0.85 &img[,,2]>0.85 & img[,,3]>0.85)
#calculate threshold
t=calculateThreshold(as.vector(imgB[-indexWhitePixel]))
#create binary image
imgB[imgB>t]=-1
imgB[imgB != -1]=0
imgB[imgB == -1]=1
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{misc}
