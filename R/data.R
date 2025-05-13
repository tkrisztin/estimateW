#' Covid incidences data
#'
#' COVID-19 data set provided by Johns Hopkins University (Dong et al., 2020). The
#' database contains information on (official) daily infections for a large panel of
#' countries around the globe in the very beginning of the outbreak
#' from 17 February to 20 April 2020.
#'
#' Data is provided for countries: Australia (AUS), Bahrain (BHR), Belgium (BEL),
#' Canada (CAN), China (CHN), Finland (FIN), France (FRA), Germany (DEU), Iran (IRN), Iraq (IRQ),
#' Israel (ISR), Italy (ITA), Japan (JPN), Kuwait (KWT), Lebanon (LBN), Malaysia (MYS), Oman (OMN),
#' Republic of Korea (KOR), Russian Federation (RUS), Singapore (SGP), Spain (ESP), Sweden (SWE),
#' Thailand (THA), United Arab Emirates (ARE), United Kingdom (GBR), United States of America (USA),
#' and Viet Nam (VNM).
#'
#' The dataset includes daily data on the country specific maximum measured temperature (Temperature) and
#' precipitation levels (Precipitation) as additional covariates (source: Dark Sky API).
#' The stringency index (Stringency) put forward by Hale et al. (2020), which summarizes country-specific
#' governmental policy measures to contain the spread of the virus. We use the biweekly average of the
#' reported stringency index.
#'
#' @name covid
#' @keywords covid infections stringency
#'
#' @docType data
#'
#' @format A \code{data.frame} object.
#'
#' @references
#'   Dong, E., Du, H., and Gardner, L. (2020). An interactive web-based dashboard to track
#'   COVID-19 in real time. \emph{The Lancet Infectious Diseases}, \bold{20(5)}, 533–534.
#'   \doi{10.1016/S1473-3099(20)30120-1}.
#'
#'   Hale, T., Petherick, A., Phillips, T., and Webster, S. (2020). Variation in government
#'   responses to COVID-19. Blavatnik School of Government Working Paper, 31, 2020–2011.
#'   \doi{10.1038/s41562-021-01079-8}.
#'
#'   Krisztin, T., and Piribauer, P. (2022). A Bayesian approach for the estimation
#'   of weight matrices in spatial autoregressive models, \emph{Spatial Economic Analysis},
#'   1-20. \doi{10.1080/17421772.2022.2095426}.
#'
#'   Krisztin, T., Piribauer, P., and Wögerer, M. (2020). The spatial econometrics of the
#'   coronavirus pandemic. \emph{Letters in Spatial and Resource Sciences}, \bold{13 (3)}, 209-218.
#'   \doi{10.1007/s12076-020-00254-1}.
#'
#'   Dong, E., Du, H., and Gardner, L. (2020). An interactive web-based dashboard to track
#'   COVID-19 in real time. \emph{The Lancet Infectious Diseases}, \bold{20(5)}, 533–534.
#'   \doi{10.1016/S1473-3099(20)30120-1}.
"covid"

#' Regional growth data
#'
#' Regional growth data set contains information on annual growth rates of GVA per worker
#' (labor productivity) for 90 European NUTS-1 regions for the period 1999-2019.
#’ The data set moreover contains initial log-levels of labor productivity as well as information
#’ on the share of low- and tertiary education attainment on working age population.
#’ Data comes from the Annual Regional Database of the European Commission's Directorate General for Regional and Urban Policy (ARDECO), and the Eurostat regional database.
#’ The data set can be viewed as a reduced form of the application in Piribauer et al. (2023).
#'
#' Data is provided for NUTS-1 regions of the following countries: ...
#'
#' The dataset includes growth of GDP per worker
#'
#' @name nuts1growth
#' @keywords regional growth
#'
#' @docType data
#'
#' @format A \code{data.frame} object.
#'
#' @references
#'   Piribauer, P., Glocker, C., & Krisztin, T. (2023). Beyond distance: The spatial relationships of European regional economic growth. Journal of Economic Dynamics and Control, 155, 104735.
"nuts1growth"
