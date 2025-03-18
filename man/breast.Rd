\name{breast}
\alias{breast}
\title{Breast Cancer Dataset}
\description{
  The study included 410 female breast cancer patients who were treated and followed up between 2011 and 2021. All participants underwent a post-dual-energy X-ray absorptiometry (DXA) examination to assess bone mineral density (BMD) at the lumbar spine and a post-thoracic CT scan that included the L1 vertebra. The aim of the study was to evaluate the correlation between CT imaging features and DXA results to develop a model for predicting osteoporosis after treatment in breast cancer patients.
}
\usage{
data(breast)
}
\format{
  The data breast contains training dataset and testing dataset, each dataset has its own radiomics features, covariate features and class label.
  \describe{
    \item{trainX}{
      A list of length 310, each element containing radiomics features in a 22x22 matrix.
    }
    \item{trainOP}{
      A vector with 310 elements, each containing a binary class label of 2 or 1, indicating the presence or absence of osteoporosis, respectively.
    }
    \item{trainU}{
      A matrix with dimensions 310x2, where each row contains the age and body weight (in kg) of a patient.
    }
    \item{testX}{
      A list of length 100, each element containing radiomics features in a 22x22 matrix.
    }
    \item{testOP}{
      A vector with 100 elements, each containing a binary class label of 2 or 1, indicating the presence or absence of osteoporosis, respectively.
    }
    \item{testU}{
      A matrix with dimensions 100x2, where each row contains the age and body weight (in kg) of a patient.
    }
  }
}
\details{
  These data come from a breast cancer research project and are used for constructing a classification model. The data have been preprocessed and can be directly used for statistical analysis.
}
\source{
  Lai, Y. H., Tsai, Y. S., Su, P. F., Li, C. I., & Chen, H. H. (2024). A computed tomography radiomics-based model for predicting osteoporosis after breast cancer treatment. Physical and Engineering Sciences in Medicine, 47(1), 239-248.
}
\examples{
  library(ReCatch)
  data(breast)
  results <- CatchPlus(x = trainX, U = trainU, y = trainOP,
           testx = testX, testU = testU, testy = testOP,
           setseed = 1234)
  results$pred
}
