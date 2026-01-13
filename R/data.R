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

#' Regional growth data for European NUTS-1 regions
#'
#' Annual growth rates of GVA per worker (labor productivity) for 90 NUTS-1 regions, 2001-2019 (T=19).
#' Explanatory variables are lagged by one year.
#'
#' The dataset contains annual regional data for 26 European Union countries, disaggregated at the NUTS-1 level. The countries covered are: Austria, Belgium, Bulgaria, Cyprus, Czechia, Denmark, Estonia, Finland, France, Germany, Hungary, Ireland, Italy, Lithuania, Luxembourg, Latvia, Malta, Netherlands, Poland, Portugal, Romania, Sweden, Slovenia, and Slovakia. The NUTS-1 codes identify the first-level administrative regions within these countries.
#'
#' @format A data frame with 1710 observations (90 regions × 19 years) and 6 variables:
#' \describe{
#'   \item{\code{NUTS1}}{NUTS-1 region code}
#'   \item{\code{year}}{Year (2001-2019)}
#'   \item{\code{growth_gdp_pw}}{Annual growth rate of GVA per worker (log differences), for year \code{year} vs. \code{year-1}}
#'   \item{\code{init_gdp_pw}}{Log initial GVA per worker (measured at \code{year-1})}
#'   \item{\code{loweduc}}{Share of working-age population with low education (ISCED 0-2, at \code{year-1})}
#'   \item{\code{higheduc}}{Share of working-age population with tertiary education (ISCED 5-8, at \code{year-1})}
#' }
#'
#' @source ARDECO (European Commission DG REGIO) and Eurostat regional database.
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
