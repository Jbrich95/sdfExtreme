#' Exponent functions
#' 
#' Exponent function for the Brown-Resnick process and its first- and second-order partial derivatives. For use in likelihoods.
#'
#' @param x \code{x} component.
 #' @param y \code{y} component.
 #' @param a \eqn{\sqrt{2\gamma(s_x-s_y)}}, where \eqn{s_x} and \eqn{s_y} are the sampling locations for the \eqn{x} and \eqn{y} components.
#' @return Value

#' @rdname Vterms
#' @export
Vterm=function(x,y,a){
  return( 1/x*pnorm(a/2-1/a*log(x/y))+1/y*pnorm(a/2-1/a*log(y/x)))
}
#' @rdname Vterms
#' @export
Vpart=function(x,y,a){
  return(-1/x^2*pnorm(a/2-1/a*log(x/y))-1/(a*x^2)*dnorm(a/2-1/a*log(x/y))+1/(a*y*x)*dnorm(a/2-1/a*log(y/x)))
}
#' @rdname Vterms
#' @export
Vpart2=function(x,y,a){
  return(-1/(x^2*a*y)*dnorm(a/2-1/a*log(x/y))-1/(a^2*x^2*y)*(-(a/2-1/a*log(x/y)))*dnorm(a/2-1/a*log(x/y))
         -1/(a*y^2*x)*dnorm(a/2-1/a*log(y/x))-1/(a^2*y^2*x)*(-(a/2-1/a*log(y/x)))*dnorm(a/2-1/a*log(y/x)))
}
