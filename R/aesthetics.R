#Functions which create plots



#' Plot distribution
#'
#' Creates a plot the distribution of a performance measure across sample size
#'
#' @param data a dataframe with each column being the result of one experiment
#' @param n1 the smallest sample size to display on the graph
#' @param n the largest sample size to display on the graph
#' @param xlabel the x axis label
#' @param ylabel the y axis label
#'
#' @return None
#'
#' @examples
#' shadeplot(data, 10, 100, "Sample size", "D-optimality")
#'
#' @export
#'
shadeplot <- function(data, n1= 4, n=100, xlabel="", ylabel="", ... ){


  atkins.bal.quant<- apply(data, 1, FUN=quantile, probs=c(0.1, 0.4, 0.5, 0.6, 0.9), na.rm=TRUE)
  atkins.bal1 <- atkins.bal.quant[1,]
  atkins.bal25 <- atkins.bal.quant[2,]
  atkins.bal50  <- atkins.bal.quant[3,]
  atkins.bal75  <-atkins.bal.quant[4,]
  atkins.bal100 <- atkins.bal.quant[5,]


  plot(x=n1:n, y=atkins.bal50 ,type="l", xlab=xlabel,ylab=ylabel,  cex.lab=1.7, cex.lab=1.7, cex.axis=1.5, cex.main=1.5, ... )
  polygon(c(n1:n,n:n1),c(atkins.bal1,rev(atkins.bal100)),col="gray60",border=NA)
  polygon(c(n1:n,n:n1),c(atkins.bal25,rev(atkins.bal75)),col="gray30",border=NA)
  lines(n1:n, atkins.bal50,lty=1,lwd=2, col="black")

}
