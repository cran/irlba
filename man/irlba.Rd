\name{irlba}
\alias{irlba}
\title{Find a few approximate singular values and singular vectors of a matrix. }
\description{
The implicitly restarted Lanczos bi-diagonalization (IRLBA) algorithm 
finds a few approximate singular values and corresponding singular vectors
of a matrix using the augmented implicitly-restarted Lanczos bidiagonalization
method of Baglama and Reichel.  It is a fast and memory-efficient way to 
compute a partial SVD.
}
\usage{
irlba(A, nu = 5, nv = 5, adjust = 3, aug = "ritz", sigma = "ls", maxit = 1000, m_b = 20, reorth = 2, tol = 1e-06, V = NULL, matmul = NULL)
}
\arguments{
  \item{A}{A double-precision real or complex matrix or real sparse matrix}
  \item{nu}{Number of desired left singular vectors}
  \item{nv}{Number of desired right singular vectors}
  \item{adjust}{Number of extra approximate singular values to compute to enhance convergence}
  \item{aug}{"ritz" for Ritz "harm" for harmonic Ritz vector augmentation}
  \item{sigma}{"ls" for largest few singular values, "ss" for smallest}
  \item{maxit}{Maximum number of iterations}
  \item{m_b}{Size of the projected bidiagonal matrix}
  \item{reorth}{Either 1 or 2: full (2) or one-sided (1) reorthogonalization}
  \item{tol}{Convergence is determined when ||A*V - U*S|| <= tol*||A||, where ||A|| is approximated by the largest singular value of all projection matrices. }
  \item{V}{Optional matrix of approximate right singular vectors}
  \item{matmul}{Option matrix multiply function--if specified, matmul must
                be a function that takes three arguments:
                \code{A, B, transpose} and computes \code{A \%*\% B} if
                \code{transpose=FALSE} or \code{t(A) \%*\% B} if
                \code{transpose=TRUE}.
   }
}
\value{
  \item{d}{min (nu, nv) approximate singular values}
  \item{u}{nu approximate left singular vectors}
  \item{v}{nv approximate right singular vectors}
  \item{iter}{The number of Lanczos iterations carried out}
}
\details{
The syntax of \code{irlba} largely conforms to \code{svd}, with an
important exception. The R \code{svd} function always returns a complete
set of singular values, even if the number of singular vectors \code{nu}
or \code{nv} is set less than the maximum. The \code{irlba} function
returns a number of singular vectors conforming to the number of
specified singular vectors \code{nu} and \code{nv}.
}
\references{"Augmented Implicitly Restarted Lanczos Bidiagonalization Methods",
     J. Baglama and L. Reichel, SIAM J. Sci. Comput. 2005.  }
\author{ Adapted for R by B. W. Lewis <Bryan@revolution-computing.com>}
\examples{
A <- matrix(runif(100*100),100,100)
S <- irlba (A)
S$d
}
\keyword{array}
\keyword{algebra}
