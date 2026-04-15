#' California Proposition 99 Tobacco Control Data (Common Timing)
#'
#' Panel data for evaluating the effect of California's Proposition 99
#' tobacco control program (1989) on cigarette sales. This is a classic
#' common timing difference-in-differences dataset where California is
#' the sole treated state and all other states serve as controls.
#'
#' Originally used in Abadie, Diamond, and Hainmueller (2010) for
#' synthetic control analysis. Used in Lee and Wooldridge (2025, 2026)
#' as the common timing DiD benchmark.
#'
#' ATT benchmark (demean, RA, homoskedastic SE):
#' \eqn{\hat{\tau} \approx -27.35} (log cigarette sales reduction).
#'
#' @format A data frame with 1209 rows and 12 columns:
#' \describe{
#'   \item{state}{State identifier (integer, 1-39). In Stata .dta the
#'     column carries value labels (e.g. 3 = "California"). In the
#'     Python CSV the column contains string state names. In this R
#'     package the column is a plain integer without labels.}
#'   \item{year}{Calendar year (integer, 1970-2000).}
#'   \item{cigsale}{Per-capita cigarette sales in packs (numeric).}
#'   \item{lcigsale}{Log per-capita cigarette sales, i.e.
#'     \code{log(cigsale)} (numeric). Primary outcome variable for
#'     DiD analysis.}
#'   \item{retprice}{Average retail price of cigarettes (numeric).}
#'   \item{lretprice}{Log retail price, i.e. \code{log(retprice)}
#'     (numeric).}
#'   \item{lnincome}{Log per-capita state personal income (numeric).
#'     Contains some NA values.}
#'   \item{beer}{Per-capita beer consumption (numeric). Contains some
#'     NA values.}
#'   \item{age15to24}{Proportion of population aged 15-24 (numeric).
#'     Contains some NA values.}
#'   \item{d}{Treatment indicator (integer, 0 or 1). Equals 1 for
#'     California (state = 3) in all years, 0 otherwise. This is the
#'     unit-level treatment assignment, not the interaction d*post.}
#'   \item{post}{Post-treatment indicator (integer, 0 or 1). Equals 1
#'     for years 1989-2000, 0 for years 1970-1988.}
#'   \item{treat}{Treatment-post interaction (integer, 0 or 1). Equals
#'     \code{d * post}. Equals 1 only for California in post-1989
#'     years.}
#' }
#'
#' @usage data(smoking)
#'
#' @source
#' Abadie, A., Diamond, A., and Hainmueller, J. (2010).
#' "Synthetic Control Methods for Comparative Case Studies:
#' Estimating the Effect of California's Tobacco Control Program."
#' \emph{Journal of the American Statistical Association}, 105(490),
#' 493-505.
#'
#' Lee, S. and Wooldridge, J. M. (2025).
#' "A Simple Transformation Approach to Difference-in-Differences
#' Estimation for Panel Data." SSRN No. 4516518.
#'
#' Lee, S. and Wooldridge, J. M. (2026).
#' "Simple Difference-in-Differences Estimation in Fixed Effects
#' Models." \emph{Journal of Econometrics} (forthcoming).
#'
#' @seealso \code{\link{castle}}, \code{\link{lwdid}},
#'   \code{\link{simulate_panel_data}}
#'
#' @examples
#' data(smoking)
#' str(smoking)
#' # Common timing DiD:
#' # result <- lwdid(data = smoking, y = "lcigsale", ivar = "state",
#' #                tvar = "year", d = "d", post = "post")
"smoking"


#' Castle Doctrine Law Data (Staggered Adoption)
#'
#' Panel data for evaluating the effect of Castle Doctrine / Stand Your
#' Ground laws on homicide rates across US states. This is a staggered
#' adoption dataset where 21 states adopted the law between 2005 and
#' 2009, while 29 states serve as never-treated controls.
#'
#' Originally used in Cheng and Hoekstra (2013). Used in Lee and
#' Wooldridge (2025, 2026) as the staggered DiD benchmark (Section 6).
#'
#' Cohort structure (number of states per adoption year):
#' 2005: 1 (Florida), 2006: 13, 2007: 4, 2008: 2, 2009: 1, NA: 29.
#'
#' @format A data frame with 550 rows and 18 columns:
#' \describe{
#'   \item{sid}{State numeric identifier (integer, 1-51, skipping 9).}
#'   \item{state_name}{State name (character). Renamed from \code{state}
#'     in the original Stata .dta to avoid confusion with numeric IDs.}
#'   \item{year}{Calendar year (integer, 2000-2010).}
#'   \item{lhomicide}{Log homicide rate, i.e. \code{log(homicide)}
#'     (numeric). Primary outcome variable for DiD analysis.}
#'   \item{homicide}{Homicide rate per 100,000 population (numeric).}
#'   \item{population}{State population (numeric).}
#'   \item{police}{Police employment per capita (numeric).}
#'   \item{unemployrt}{Unemployment rate (numeric).}
#'   \item{income}{Per-capita income (numeric).}
#'   \item{poverty}{Poverty rate (numeric).}
#'   \item{blackm_15_24}{Proportion of Black males aged 15-24
#'     (numeric).}
#'   \item{whitem_15_24}{Proportion of White males aged 15-24
#'     (numeric).}
#'   \item{prisoner}{Incarceration rate (numeric).}
#'   \item{cdl}{Castle Doctrine Law intensity (numeric, 0.0-1.0).
#'     Despite the Stata label "=1 if CDL is effective", this is a
#'     continuous proportion variable with intermediate values between
#'     0 and 1, not a binary indicator.}
#'   \item{post}{Post-treatment indicator (integer, 0 or 1).}
#'   \item{effyear}{Effective year of Castle Doctrine adoption
#'     (integer or NA). NA for never-treated states. Values in
#'     \{2005, 2006, 2007, 2008, 2009\}.}
#'   \item{gvar}{Treatment cohort variable (integer or NA). Derived
#'     from \code{effyear}: \code{gvar = as.integer(effyear)}. NA
#'     indicates never-treated units (R convention). In the Python
#'     implementation, never-treated units are marked with
#'     \code{gvar = 0} instead.}
#'   \item{dinf}{Never-treated indicator (integer, 0 or 1). Equals 1
#'     if the state never adopted Castle Doctrine (\code{is.na(gvar)}),
#'     0 otherwise.}
#' }
#'
#' @usage data(castle)
#'
#' @source
#' Cheng, C. and Hoekstra, M. (2013).
#' "Does Strengthening Self-Defense Law Deter Crime or Escalate
#' Violence? Evidence from Expansions to Castle Doctrine."
#' \emph{Journal of Human Resources}, 48(3), 821-854.
#'
#' Lee, S. and Wooldridge, J. M. (2025).
#' "A Simple Transformation Approach to Difference-in-Differences
#' Estimation for Panel Data." SSRN No. 4516518.
#'
#' Cunningham, S. (2021). \emph{Causal Inference: The Mixtape}.
#' Yale University Press.
#'
#' @seealso \code{\link{smoking}}, \code{\link{lwdid}},
#'   \code{\link{simulate_panel_data}}
#'
#' @examples
#' data(castle)
#' str(castle)
#' table(castle$gvar[!duplicated(castle$sid)], useNA = "always")
#' # Staggered DiD:
#' # result <- lwdid(data = castle, y = "lhomicide", ivar = "sid",
#' #                tvar = "year", gvar = "gvar")
"castle"


#' Walmart Store Openings and Local Labor Markets (Staggered Adoption)
#'
#' County-level panel data for evaluating the effect of Walmart store
#' openings on local labor markets. This is a staggered adoption dataset
#' where 886 counties experienced a first Walmart opening between 1986
#' and 1999, while 394 counties serve as never-treated controls.
#'
#' Originally constructed by Brown and Butts (2025), who follow Basker
#' (2005) in using County Business Patterns (CBP) data from 1964 and
#' 1977--1999. Used in Lee and Wooldridge (2025) Section 6 as the
#' large-scale staggered DiD application.
#'
#' The dataset is limited to counties with more than 1,500 total
#' employment in 1964 and non-negative employment growth during
#' 1964--1977. Counties whose first Walmart opened before 1985 are
#' excluded, yielding a balanced panel of 1,280 counties over 23 years
#' (1977--1999).
#'
#' @format A data frame with 29,440 rows and 21 columns:
#' \describe{
#'   \item{fips}{County FIPS code (integer). Unit identifier.}
#'   \item{state_fips}{State FIPS code (integer).}
#'   \item{county_fips}{County FIPS code within state (integer).}
#'   \item{year}{Calendar year (integer, 1977--1999).}
#'   \item{retail_emp}{County-level retail employment (numeric).}
#'   \item{nonretail_emp}{Non-retail employment (numeric).}
#'   \item{emp1964}{Total employment in 1964 (numeric).}
#'   \item{emp1977}{Total employment in 1977 (numeric).}
#'   \item{any_open}{Whether any Walmart store is open (logical).}
#'   \item{n_open}{Number of Walmart stores open (integer).}
#'   \item{g}{Treatment cohort variable (numeric). Year of first
#'     Walmart opening. \code{Inf} for never-treated counties.}
#'   \item{rel_year}{Relative year to first Walmart opening (numeric).
#'     \code{-Inf} for never-treated counties.}
#'   \item{log_retail_emp}{Log retail employment (numeric). Primary outcome.}
#'   \item{log_nonretail_emp}{Log non-retail employment (numeric).}
#'   \item{log_wholesale_emp}{Log wholesale employment (numeric).}
#'   \item{log_manufacturing_emp}{Log manufacturing employment (numeric).}
#'   \item{log_construction_emp}{Log construction employment (numeric).}
#'   \item{total_pop}{Total population, 1980 Census (numeric).}
#'   \item{retail_emp_share}{Retail employment share (numeric).}
#'   \item{nonretail_emp_share}{Non-retail employment share (numeric).}
#'   \item{balanced}{Whether the county is in the balanced panel (logical).}
#' }
#'
#' @usage data(walmart)
#'
#' @source
#' Brown, N. and Butts, K. (2025).
#' "Dynamic Treatment Effect Estimation with Interactive Fixed Effects
#' and Short Panels."
#' \emph{Journal of Econometrics} (forthcoming).
#'
#' Basker, E. (2005).
#' "Job Creation or Destruction? Labor Market Effects of Wal-Mart
#' Expansion."
#' \emph{Review of Economics and Statistics}, 87(1), 174--183.
#'
#' Lee, S. and Wooldridge, J. M. (2025).
#' "A Simple Transformation Approach to Difference-in-Differences
#' Estimation for Panel Data." SSRN No. 4516518.
#'
#' @seealso \code{\link{castle}}, \code{\link{smoking}},
#'   \code{\link{lwdid}}
#'
#' @examples
#' data(walmart)
#' dim(walmart)
#' # Cohort distribution:
#' table(walmart$g[!duplicated(walmart$fips)])
#' # Staggered DiD on balanced subsample:
#' # bal <- walmart[walmart$balanced, ]
#' # result <- lwdid(data = bal, y = "log_retail_emp", ivar = "fips",
#' #                tvar = "year", gvar = "g", rolling = "detrend")
"walmart"
