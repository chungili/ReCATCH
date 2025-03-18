\name{CRC}
\alias{CRC}
\title{Colorectal Cancer Dataset}
\description{
  The dataset consists of 71 colorectal cancer (CRC) patients who underwent surgery and adjuvant chemotherapy. The study aims to classify cancer recurrence based on adjusted computed tomography (CT) scan radiomics features and immune genomic expression.
}
\usage{
data(CRC)
}
\format{
  For each observation in the dataset, there are three types of data: radiomics features in a 33x32 matrix, gene expression levels in a 30-dimensional vector, and a class label indicating cancer recurrence status.
  \describe{
    \item{RF}{
      A list with 71 elements, each containing radiomics features in a 33x32 matrix.
    }
    \item{recurrence}{
      A vector with 71 elements, each containing a binary class label of 2 or 1, indicating whether cancer recurrence occurs or not, respectively.
    }
    \item{DE}{
      A matrix with dimensions 71x30, where each row contains the gene expression levels.
    }
  }
}
\details{
  These data come from a colorectal cancer research project and are used for constructing a classification model. The data have been preprocessed and can be directly used for statistical analysis.
}
\source{
  Huang, Y. C., Tsai, Y. S., Li, C. I., Chan, R. H., Yeh, Y. M., Chen, P. C., M.R. Shen and Lin, P. C. (2022). Adjusted CT image-based radiomic features combined with immune genomic expression achieve accurate prognostic classification and identification of therapeutic targets in Stage III Colorectal Cancer. Cancers(Basel), 14(8), 1895.
}
\examples{
  library(ReCatch)
  data(CRC)
  results <- ReCatch(x = RF, U = DE, y = recurrence, setseed = 123)
  results$pred
}
